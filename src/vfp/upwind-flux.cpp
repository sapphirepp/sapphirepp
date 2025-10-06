// -----------------------------------------------------------------------------
//
// Copyright (C) 2023 by the Sapphire++ authors
//
// This file is part of Sapphire++.
//
// Sapphire++ is free software: you can redistribute it and/or modify it under
// the terms of the GNU Lesser General Public License as published by the Free
// Software Foundation, either version 3 of the License, or (at your option) any
// later version.
//
// Sapphire++ is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
// A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more
// details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with Sapphire++. If not, see <https://www.gnu.org/licenses/>.
//
// -----------------------------------------------------------------------------

/**
 * @file upwind-flux.cpp
 * @author Nils Schween (nils.schween@mpi-hd.mpg.de)
 * @author Florian Schulze (florian.schulze@mpi-hd.mpg.de)
 * @brief Implement @ref sapphirepp::VFP::UpwindFlux
 */

#include "upwind-flux.h"

#include <deal.II/base/exceptions.h>

#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/lapack_support.h>
// direct access to the Fortran driver routines of Lapack
#include <deal.II/lac/lapack_templates.h>
#include <deal.II/lac/vector.h>

#include <algorithm>
// for testing arbitrary coefficients in matrix sum and the corresponding
// eigenvalue/eigenvector computations
#include <iomanip>
#include <random>

#include "config.h"



template <unsigned int dim, bool has_momentum, bool logarithmic_p>
sapphirepp::VFP::UpwindFlux<dim, has_momentum, logarithmic_p>::UpwindFlux(
  const PDESystem          &system,
  const VFPParameters<dim> &solver_control,
  const PhysicalParameters &physical_parameters)
  : pde_system{system}
  , matrix_size{pde_system.system_size}
  , matrix_size_blas{static_cast<dealii::types::blas_int>(matrix_size)}
  // NOLINTBEGIN(google-readability-casting)
  , advection_matrices(3, std::vector<double>(matrix_size * matrix_size))
  , adv_mat_products(6, std::vector<double>(matrix_size * matrix_size))
  // NOLINTEND(google-readability-casting)
  , matrix_sum(matrix_size * matrix_size)
  , eigenvalues(matrix_size)
  , eigenvectors(matrix_size * matrix_size)
  , positive_eigenvalues(matrix_size)
  , negative_eigenvalues(matrix_size)
  , eigenvalues_advection_matrices(matrix_size)
  // NOLINTBEGIN(google-readability-casting)
  , eigenvectors_advection_matrices(3,
                                    std::vector<double>(matrix_size *
                                                        matrix_size))
  // NOLINTEND(google-readability-casting)
  , background_velocity_field(physical_parameters)
  , particle_velocity_func(solver_control.mass)
  , particle_gamma_func(solver_control.mass)
  , mass(solver_control.mass)
  , charge(solver_control.charge)
  , velocity(solver_control.velocity)
  , isuppz(2 * matrix_size)
  , jobz{&dealii::LAPACKSupport::V}
  , range{&dealii::LAPACKSupport::A}
  , uplo{&dealii::LAPACKSupport::U}
  , abstol{0.}
  , // see Documentation of xsyever
  int_dummy{&dealii::LAPACKSupport::one}
  , double_dummy{1.}
{
  // NOTE: Since we very often call compute_matrix_sum and the matrices classes
  // of dealii do not allow unchecked access to there raw data, we create copies
  // of the matrices in the hope that this additional memory consumption is made
  // up for by performance gains.
  const std::vector<dealii::LAPACKFullMatrix<double>> &adv_mat =
    pde_system.get_advection_matrices();
  for (unsigned int k = 0; k < 3; ++k)
    {
      for (unsigned int i = 0; i < matrix_size; ++i)
        for (unsigned int j = 0; j < matrix_size; ++j)
          // NOTE: When computing the eigenvalues of matrix_sum, I call a
          // Fortran routine which requires the matrix elements to be stored in
          // column-major order. But since the resulting matrix is symmetric,
          // there is no difference between column-major order and row-major
          // order. I stick with the C convention and use row-major order.
          advection_matrices[k][i * matrix_size + j] = adv_mat[k](i, j);
    }
  const std::vector<dealii::LAPACKFullMatrix<double>> &adv_mat_prod =
    pde_system.get_adv_mat_products();
  for (unsigned int k = 0; k < 6; ++k)
    {
      for (unsigned int i = 0; i < matrix_size; ++i)
        for (unsigned int j = 0; j < matrix_size; ++j)
          adv_mat_products[k][i * matrix_size + j] = adv_mat_prod[k](i, j);
    }

  prepare_work_arrays_for_lapack();
  prepare_upwind_fluxes();
}



template <unsigned int dim, bool has_momentum, bool logarithmic_p>
void
sapphirepp::VFP::UpwindFlux<dim, has_momentum, logarithmic_p>::set_time(
  double time)
{
  background_velocity_field.set_time(time);
}



template <unsigned int dim, bool has_momentum, bool logarithmic_p>
void
sapphirepp::VFP::UpwindFlux<dim, has_momentum, logarithmic_p>::
  compute_upwind_fluxes(
    const std::vector<dealii::Point<dim>>     &q_points,
    const std::vector<dealii::Tensor<1, dim>> &normals,
    std::vector<dealii::FullMatrix<double>>   &positive_flux_matrices,
    std::vector<dealii::FullMatrix<double>>   &negative_flux_matrices)
{
  // Determine if we are computing the flux of an interior face whose normal
  // points into the x,y,z or p direction. Since we are using a rectangular
  // grid, the normal is the same for all quadrature points and it points either
  // in the x, y, z or p direction. Hence it can determined outside the loop
  // over the quadrature points
  unsigned int          component    = 0;
  [[maybe_unused]] bool found_normal = false;
  for (unsigned int i = 0; i < dim; ++i)
    {
      // NOTE: Such a simple comparison of doubles is only possible, because I
      // know that I am comparing 0. with 1. or something very close to one with
      // one. If this ever fails a more sophisticated comparison is needed, e.g.
      // https://floating-point-gui.de/errors/comparison/
      // std::abs is nessecary because of the boundary faces, i.e. <-| n_k = -1
      // for some faces.
      if (std::abs(1. - std::abs(normals[0][i])) <
          VFPParameters<dim>::epsilon_d)
        {
          component    = i;
          found_normal = true;
          break;
        }
    }
  // Check in Debug mode, if the normal was found. The double comparison could
  // fail.
  Assert(
    found_normal,
    dealii::ExcMessage(
      "When computing the flux through a cell surface, it was not possible to determine the normal of the cell interface."));

  // Fluxes in the spatial directions
  if (component < dim_cs)
    {
      // Return zero fluxes for case with no spatial advection term
      // but dim_cs != 0
      if constexpr ((vfp_flags & VFPFlags::spatial_advection) == VFPFlags::none)
        return;

      // Background velocity field
      // Get the value of the velocity field at every quadrature point
      // TODO: Value list for a specific component would be faster.
      std::vector<dealii::Vector<double>> velocities(q_points.size(),
                                                     dealii::Vector<double>(3));
      background_velocity_field.vector_value_list(q_points, velocities);

      // Particle
      std::vector<double> particle_velocities(q_points.size());
      if constexpr (has_momentum)
        {
          // if distribution functions depends on p, the particle velocity
          // depends on the position in the phase space and needs to be computed
          particle_velocity_func.value_list(q_points, particle_velocities);
        }
      else
        {
          std::fill(particle_velocities.begin(),
                    particle_velocities.end(),
                    velocity);
        }
      // Compute flux at every quadrature point
      for (unsigned int q_index = 0; q_index < q_points.size(); ++q_index)
        compute_flux_in_space_directions(component,
                                         normals[q_index][component],
                                         velocities[q_index][component],
                                         particle_velocities[q_index],
                                         positive_flux_matrices[q_index],
                                         negative_flux_matrices[q_index]);
    }
  else
    { // Fluxes in the p direction
      // NOTE: If momentum terms are included, then the last component of
      // normals, points etc. are the ones corresponding to the momentum
      // direction
      std::vector<dealii::Vector<double>> material_derivative_vel(
        q_points.size(), dealii::Vector<double>(3));
      background_velocity_field.material_derivative_list(
        q_points, material_derivative_vel);

      std::vector<dealii::FullMatrix<double>> jacobians_vel(
        q_points.size(), dealii::FullMatrix<double>(3));
      background_velocity_field.jacobian_list(q_points, jacobians_vel);

      std::vector<double> particle_gammas(q_points.size());
      particle_gamma_func.value_list(q_points, particle_gammas);
      for (unsigned int q_index = 0; q_index < q_points.size(); ++q_index)
        compute_flux_in_p_direction(normals[q_index][component],
                                    q_points[q_index][component],
                                    particle_gammas[q_index],
                                    material_derivative_vel[q_index],
                                    jacobians_vel[q_index],
                                    positive_flux_matrices[q_index],
                                    negative_flux_matrices[q_index]);
    }
}



template <unsigned int dim, bool has_momentum, bool logarithmic_p>
void
sapphirepp::VFP::UpwindFlux<dim, has_momentum, logarithmic_p>::
  compute_local_lax_friedrichs_fluxes(
    const std::vector<dealii::Point<dim>>     &q_points,
    const std::vector<dealii::Tensor<1, dim>> &normals,
    std::vector<dealii::FullMatrix<double>>   &flux_matrices,
    std::vector<double>                       &max_eigenvalues)
{
  bool space_direction = true;
  if constexpr (has_momentum)
    {
      if (std::abs(1. - std::abs(normals[0][dim - 1])) <
          VFPParameters<dim>::epsilon_d)
        space_direction = false;
    }

  if (space_direction)
    {
      // Return zero fluxes for case with no spatial advection term
      // but dim_cs != 0
      if constexpr ((vfp_flags & VFPFlags::spatial_advection) == VFPFlags::none)
        return;

      // Background velocity field
      // Get the value of the velocity field at every quadrature point
      // TODO: Value list for a specific component would be faster.
      std::vector<dealii::Vector<double>> velocities(q_points.size(),
                                                     dealii::Vector<double>(3));
      background_velocity_field.vector_value_list(q_points, velocities);

      // Particle
      std::vector<double> particle_velocities(q_points.size());
      if constexpr (has_momentum)
        {
          // if distribution functions depends on p, the particle velocity
          // depends on the position in the phase space and needs to be computed
          particle_velocity_func.value_list(q_points, particle_velocities);
        }
      else
        {
          std::fill(particle_velocities.begin(),
                    particle_velocities.end(),
                    velocity);
        }

      for (unsigned int q_index = 0; q_index < q_points.size(); ++q_index)
        {
          for (unsigned int d = 0; d < dim_cs; ++d)
            for (unsigned int i = 0; i < matrix_size; ++i)
              for (unsigned int j = 0; j < matrix_size; ++j)
                flux_matrices[q_index](i, j) +=
                  (normals[q_index][d] * particle_velocities[q_index] *
                     advection_matrices[d][i * matrix_size + j] +
                   ((i == j) ? normals[q_index][d] * velocities[q_index][d] :
                               0.)); // NOTE: integer divide is on purpose )

          double normal_velocity = 0;
          for (unsigned int d = 0; d < dim_cs; ++d)
            normal_velocity += normals[q_index][d] * velocities[q_index][d];

          // eigenvalues_advection_matrices are in ascending order,
          // see LAPACK.dsyevr()
          // https://www.netlib.org/lapack/explore-html/d1/d56/group__heevr_gaa334ac0c11113576db0fc37b7565e8b5.html#gaa334ac0c11113576db0fc37b7565e8b5
          const double min_eigenvalue = eigenvalues_advection_matrices[0];
          const double max_eigenvalue =
            eigenvalues_advection_matrices[matrix_size - 1];
          const double max_eigenvalue_advection_matrices =
            std::max(std::abs(min_eigenvalue), std::abs(max_eigenvalue));

          /** @todo Precompute max eigenvalue */
          max_eigenvalues[q_index] =
            particle_velocities[q_index] * max_eigenvalue_advection_matrices +
            std::abs(normal_velocity);
        }
    }
  else
    { // Fluxes in the p direction
      // NOTE: If momentum terms are included, then the last component of
      // normals, points etc. are the ones corresponding to the momentum
      // direction
      std::vector<dealii::Vector<double>> material_derivative_vel(
        q_points.size(), dealii::Vector<double>(3));
      background_velocity_field.material_derivative_list(
        q_points, material_derivative_vel);

      std::vector<dealii::FullMatrix<double>> jacobians_vel(
        q_points.size(), dealii::FullMatrix<double>(3));
      background_velocity_field.jacobian_list(q_points, jacobians_vel);

      std::vector<double> particle_gammas(q_points.size());
      particle_gamma_func.value_list(q_points, particle_gammas);


      for (unsigned int q_index = 0; q_index < q_points.size(); ++q_index)
        {
          // compute the matrix sum at the point q at time t. Overwrites the
          // member variable matrix_sum
          compute_matrix_sum(normals[q_index][dim - 1],  // n_p
                             q_points[q_index][dim - 1], // momentum
                             particle_gammas[q_index],
                             material_derivative_vel[q_index],
                             jacobians_vel[q_index]);

          // copy the flux matrices
          for (unsigned int i = 0; i < matrix_size; ++i)
            {
              for (unsigned int j = 0; j < matrix_size; ++j)
                {
                  flux_matrices[q_index](i, j) =
                    matrix_sum[i * matrix_size + j];
                }
            }


          dealii::syevr(&dealii::LAPACKSupport::N, // compute eigenvalues only
                        &dealii::LAPACKSupport::A, // compute all eigenvalues
                        uplo,
                        &matrix_size_blas,
                        matrix_sum.data(),
                        &matrix_size_blas,
                        &double_dummy,
                        &double_dummy,
                        int_dummy,
                        int_dummy,
                        &abstol,
                        &num_eigenvalues,
                        eigenvalues.data(),
                        eigenvectors.data(),
                        &matrix_size_blas,
                        isuppz.data(),
                        work.data(),
                        &lwork,
                        iwork.data(),
                        &liwork,
                        &info);
          Assert(info >= 0, dealii::ExcInternalError());
          if (info != 0)
            std::cerr /** @todo Better error handling */
              << "The computation of the eigenvalues and eigenvectors for the "
                 "flux in p direction failed: LAPACK error in syevr"
              << std::endl;

          const double min_eigenvalue = eigenvalues[0];
          const double max_eigenvalue = eigenvalues[matrix_size - 1];
          max_eigenvalues[q_index] =
            std::max(std::abs(min_eigenvalue), std::abs(max_eigenvalue));
        }
    }
}



template <unsigned int dim, bool has_momentum, bool logarithmic_p>
void
sapphirepp::VFP::UpwindFlux<dim, has_momentum, logarithmic_p>::
  prepare_work_arrays_for_lapack()
{
  // Preparations for the eigenvalue and eigenvector computations
  //
  // The computation of the upwind flux in the p direction requires to the
  // compute eigenvalues and eigenvectors for every face with n_p = 1 (or -1).
  // Lapack needs correctly sized work arrays to do these computations. It can
  // determine the sizes it needs with a specific call to the used eigenvalue
  // routine. The work arrays will be allocated only once and then reused
  // every time eigenvalues and eigenvectors are computed.

  // https://stackoverflow.com/questions/46618391/what-is-the-use-of-the-work-parameters-in-lapack-routines
  // Quote: "If dynamic allocation is available, the most common use of LAPACK
  // functions is to first perform a workspace query using LWORK = -1, then
  // use the return value to allocate a WORK array of the correct size and
  // finally call the routine of LAPACK to get the expected result. High-end
  // wrappers of LAPACK such as LAPACKE features function doing just that:
  // take a look at the source of LAPACKE for function LAPACKE_dsyev()! It
  // calls twice the function LAPACKE_dsyev_work, which calls LAPACK_dsyev
  // (wrapping dsyev()).

  // Wrappers still feature functions such as LAPACKE_dsyev_work(), where the
  // arguments work and lwork are still required. The number of allocations
  // can therefore be reduced if the routine is called multiple times on
  // similar sizes by not deallocating WORK between calls, but the user must
  // do that himself (see this example). In addition, the source of ILAENV,
  // the function of LAPACK called to compute the optimzed size of WORK,
  // features the following text:

  //     This version provides a set of parameters which should give good, but
  //     not optimal, performance on many of the currently available
  //     computers. Users are encouraged to modify this subroutine to set the
  //     tuning parameters for their particular machine using the option and
  //     problem size information in the arguments.

  // As a result, testing sizes of WORK larger than the size returned by the
  // workspace query could improve performances."

  // Compute the size of the working arrays for the Lapack routine xsyevr

  // The optimal size of the work arrays will be stored in their first element
  work.resize(1);
  iwork.resize(1);
  lwork  = -1;
  liwork = -1;
  // create dummy variables for the matrix sum
  double                     n_p_dummy   = 1.;
  double                     gamma_dummy = 1.;
  double                     p_dummy     = 1.;
  dealii::Vector<double>     material_derivative_dummy{1., 1., 1.};
  dealii::FullMatrix<double> jacobian_dummy(3);
  for (unsigned int i = 0; i < jacobian_dummy.m(); ++i)
    for (unsigned int j = 0; j < jacobian_dummy.n(); ++j)
      jacobian_dummy[i][j] = 1.;

  compute_matrix_sum(
    n_p_dummy, p_dummy, gamma_dummy, material_derivative_dummy, jacobian_dummy);
  // call Lapack routine
  dealii::syevr(jobz,
                range,
                uplo,
                &matrix_size_blas,
                matrix_sum.data(),
                &matrix_size_blas,
                &double_dummy,
                &double_dummy,
                int_dummy,
                int_dummy,
                &abstol,
                &num_eigenvalues,
                eigenvalues.data(),
                eigenvectors.data(),
                &matrix_size_blas,
                isuppz.data(),
                work.data(),
                &lwork,
                iwork.data(),
                &liwork,
                &info);
  Assert(info == 0, dealii::ExcInternalError());

  // resize the working arrays
  lwork  = static_cast<dealii::types::blas_int>(work[0]);
  liwork = iwork[0];
  work.resize(static_cast<unsigned int>(lwork));
  iwork.resize(static_cast<unsigned int>(liwork));
}



template <unsigned int dim, bool has_momentum, bool logarithmic_p>
void
sapphirepp::VFP::UpwindFlux<dim, has_momentum, logarithmic_p>::
  prepare_upwind_fluxes()
{
  // The eigenvalues of A_x are also the eigenvalues of A_y and A_z. The
  // eigenvalues of A_x are the roots of the associated legendre polynomials.
  // Since the associated Legendre Polynomials are orthogonal polynomials, it
  // can be shown that all roots are contained in in the intervall [-1., 1.].
  // Physically this makes sense, because c = 1 and the eigenvalues encode the
  // speed of 'information' transport.
  //
  // Helper variables
  const char *const V  = &dealii::LAPACKSupport::V;
  const double      vl = -1.1;
  const double      vu = 1.1;
  // NOTE: We need another copy of the advection matrix, because it gets
  // destroyed when the xsyevr is called.
  std::vector<double> A_x;
  A_x = advection_matrices[0];
  dealii::syevr(jobz,
                V,
                uplo,
                &matrix_size_blas,
                A_x.data(),
                &matrix_size_blas,
                &vl,
                &vu,
                int_dummy,
                int_dummy,
                &abstol,
                &num_eigenvalues,
                eigenvalues_advection_matrices.data(),
                eigenvectors_advection_matrices[0].data(),
                &matrix_size_blas,
                isuppz.data(),
                work.data(),
                &lwork,
                iwork.data(),
                &liwork,
                &info);
  // NOTE: The eigenvalues of A_x can computed very efficiently because A_x is
  // (when ordered correctly) a tridiagonal symmetric matrix, i.e. it is in
  // Hessenberg-Form. At the moment it is not ordered correctly. The
  // computation has only to be done once.
  Assert(info >= 0, dealii::ExcInternalError());
  if (info != 0)
    std::cerr << "The computation of the eigenvalues and eigenvectors of A_x "
                 "failed: LAPACK error in syevr"
              << std::endl;
  std::vector<double> A_y;
  A_y = advection_matrices[1];
  dealii::syevr(jobz,
                V,
                uplo,
                &matrix_size_blas,
                A_y.data(),
                &matrix_size_blas,
                &vl,
                &vu,
                int_dummy,
                int_dummy,
                &abstol,
                &num_eigenvalues,
                eigenvalues_advection_matrices.data(),
                eigenvectors_advection_matrices[1].data(),
                &matrix_size_blas,
                isuppz.data(),
                work.data(),
                &lwork,
                iwork.data(),
                &liwork,
                &info);
  Assert(info >= 0, dealii::ExcInternalError());
  if (info != 0)
    std::cerr << "The computation of the eigenvalues and eigenvectors of A_y "
                 "failed: LAPACK error in syevr"
              << std::endl;

  std::vector<double> A_z;
  A_z = advection_matrices[2];
  dealii::syevr(jobz,
                V,
                uplo,
                &matrix_size_blas,
                A_z.data(),
                &matrix_size_blas,
                &vl,
                &vu,
                int_dummy,
                int_dummy,
                &abstol,
                &num_eigenvalues,
                eigenvalues_advection_matrices.data(),
                eigenvectors_advection_matrices[2].data(),
                &matrix_size_blas,
                isuppz.data(),
                work.data(),
                &lwork,
                iwork.data(),
                &liwork,
                &info);
  Assert(info >= 0, dealii::ExcInternalError());
  if (info != 0)
    std::cerr << "The computation of the eigenvalues and eigenvectors of A_z "
                 "failed: LAPACK error in syevr"
              << std::endl;

  // TODO: Rotate the eigenvectors of A_x to get the eigenvectors of A_y and
  // A_z, i.e. implement the rotation matrices e^{-i\Omega_x pi/2} and e^{-i
  // \Omega_z pi/2} For now we will three times compute the same eigenvalues
  // to get the eigenvectors of A_y and A_z
}



template <unsigned int dim, bool has_momentum, bool logarithmic_p>
void
sapphirepp::VFP::UpwindFlux<dim, has_momentum, logarithmic_p>::test()
{
  saplog << "Eigenvalues: \n";
  for (auto &lambda : eigenvalues_advection_matrices)
    saplog << lambda << " ";
  saplog << std::endl;

  dealii::FullMatrix<double> test_positive_flux_matrix(matrix_size);
  dealii::FullMatrix<double> test_negative_flux_matrix(matrix_size);

  // Test space directions
  compute_flux_in_space_directions(
    1, 1., -0.5, 0.9, test_positive_flux_matrix, test_negative_flux_matrix);

  // saplog << "positive_flux_matrix \n";
  // test_positive_flux_matrix.print_formatted(saplog);
  // saplog << "negative_flux_matrix \n";
  // test_negative_flux_matrix.print_formatted(saplog);

  // Test momentum direction
  // Random number generator
  class RandomNumber
  {
  public:
    RandomNumber(double low, double high)
      : m_uniform_dist{low, high}
    {}

    double
    operator()()
    {
      return m_uniform_dist(m_re);
    }
    void
    seed(unsigned int s)
    {
      m_re.seed(s);
    }

  private:
    std::default_random_engine             m_re;
    std::uniform_real_distribution<double> m_uniform_dist; // Default intervall
  };
  RandomNumber       rnd_number_generator{-1.e4, 1.e4};
  std::random_device seed; // Seeded with /dev/urandom
  rnd_number_generator.seed(seed());

  // Fill the arguments of compute_flux_in_p_directon with randon numbers
  double                     test_n_p   = 1.;
  double                     test_gamma = rnd_number_generator();
  double                     test_p     = rnd_number_generator();
  dealii::Vector<double>     test_material_derivative{rnd_number_generator(),
                                                  rnd_number_generator(),
                                                  rnd_number_generator()};
  dealii::FullMatrix<double> test_jacobian(3);
  for (unsigned int i = 0; i < test_jacobian.m(); ++i)
    for (unsigned int j = 0; j < test_jacobian.n(); ++j)
      test_jacobian[i][j] = rnd_number_generator();

  // Print out the random values:
  saplog << "n_p: " << test_n_p << "\n";
  saplog << "gamma: " << test_gamma << "\n";
  saplog << "momentum: " << test_p << "\n";
  saplog << "Material derivative: \n";
  test_material_derivative.print(saplog);
  saplog << "Jacobian: "
         << "\n";
  test_jacobian.print(saplog);

  // Python array ouput to ease testing
  saplog << std::fixed << std::setprecision(9) << "[" << test_n_p << ", "
         << test_gamma << ", " << test_material_derivative[0] << ", "
         << test_material_derivative[1] << ", " << test_material_derivative[2]
         << ", " << test_p << ", " << test_jacobian[0][0] << ", "
         << test_jacobian[1][1] << ", " << test_jacobian[2][2] << ", "
         << test_jacobian[0][1] << ", " << test_jacobian[1][0] << ", "
         << test_jacobian[0][2] << ", " << test_jacobian[2][0] << ", "
         << test_jacobian[1][2] << ", " << test_jacobian[2][1] << "]"
         << std::endl;

  // reset matrices
  test_positive_flux_matrix = 0;
  test_negative_flux_matrix = 0;
  // Compute flux matrices
  compute_flux_in_p_direction(test_n_p,
                              test_p,
                              test_gamma,
                              test_material_derivative,
                              test_jacobian,
                              test_positive_flux_matrix,
                              test_negative_flux_matrix);

  saplog << "positive_flux_matrix \n";
  test_positive_flux_matrix.print_formatted(saplog);
  saplog << "negative_flux_matrix \n";
  test_negative_flux_matrix.print_formatted(saplog);
}



template <unsigned int dim, bool has_momentum, bool logarithmic_p>
void
sapphirepp::VFP::UpwindFlux<dim, has_momentum, logarithmic_p>::
  compute_flux_in_space_directions(
    const unsigned int          component,
    const double                n_component,
    const double                background_velocity,
    const double                particle_velocity,
    dealii::FullMatrix<double> &positive_flux_matrix,
    dealii::FullMatrix<double> &negative_flux_matrix)
{
  // Create a copy of the eigenvalues of the advection matrices to compute the
  // positive and negative fluxes at point q at time t
  eigenvalues = eigenvalues_advection_matrices;
  // Update the eigenvalues
  std::for_each(eigenvalues.begin(), eigenvalues.end(), [=](double &lambda) {
    lambda = n_component * (background_velocity + particle_velocity * lambda);
  });

  // split eigenvalues in positive and negative ones
  std::replace_copy_if(
    eigenvalues.begin(),
    eigenvalues.end(),
    positive_eigenvalues.begin(),
    [](const double &val) { return val < 0.0; },
    0.0);
  std::replace_copy_if(
    eigenvalues.begin(),
    eigenvalues.end(),
    negative_eigenvalues.begin(),
    [](const double &val) { return val > 0.0; },
    0.0);

  // compute the flux matrices: triple product V * Lambda_{+/-} * V^{T})
  //
  // This is reimplementation of a part of Lapacks
  // dgemm routine to allow for a diagonal matrix in between the two matrix
  // factors.
  // https://netlib.org/lapack/explore-html/d1/d54/group__double__blas__level3_gaeda3cbd99c8fb834a60a6412878226e1.html
  double temp_positive = 0;
  double temp_negative = 0;
  for (unsigned int j = 0; j < matrix_size; ++j)
    {
      for (unsigned int l = 0; l < matrix_size; ++l)
        {
          temp_positive =
            positive_eigenvalues[l] *
            eigenvectors_advection_matrices[component][l * matrix_size + j];
          temp_negative =
            negative_eigenvalues[l] *
            eigenvectors_advection_matrices[component][l * matrix_size + j];
          for (unsigned int i = 0; i < matrix_size; ++i)
            {
              positive_flux_matrix(i, j) +=
                temp_positive *
                eigenvectors_advection_matrices[component][l * matrix_size + i];
              negative_flux_matrix(i, j) +=
                temp_negative *
                eigenvectors_advection_matrices[component][l * matrix_size + i];
            }
        }
    }
}



template <unsigned int dim, bool has_momentum, bool logarithmic_p>
void
sapphirepp::VFP::UpwindFlux<dim, has_momentum, logarithmic_p>::
  compute_matrix_sum(const double                      n_p,
                     const double                      momentum,
                     const double                      gamma,
                     const dealii::Vector<double>     &material_derivative,
                     const dealii::FullMatrix<double> &jacobian)
{
  if (logarithmic_p)
    {
      double p = std::exp(momentum);
      for (unsigned int i = 0; i < matrix_size * matrix_size; ++i)
        matrix_sum[i] =
          -n_p * (gamma / p *
                    (material_derivative[0] * advection_matrices[0][i] +
                     material_derivative[1] * advection_matrices[1][i] +
                     material_derivative[2] * advection_matrices[2][i]) +
                  jacobian[0][0] * adv_mat_products[0][i] +
                  jacobian[1][1] * adv_mat_products[3][i] +
                  jacobian[2][2] * adv_mat_products[5][i] +
                  (jacobian[0][1] + jacobian[1][0]) * adv_mat_products[1][i] +
                  (jacobian[0][2] + jacobian[2][0]) * adv_mat_products[2][i] +
                  (jacobian[1][2] + jacobian[2][1]) * adv_mat_products[4][i]);
    }
  else
    for (unsigned int i = 0; i < matrix_size * matrix_size; ++i)
      matrix_sum[i] =
        -n_p * (gamma * mass *
                  (material_derivative[0] * advection_matrices[0][i] +
                   material_derivative[1] * advection_matrices[1][i] +
                   material_derivative[2] * advection_matrices[2][i]) +
                momentum *
                  (jacobian[0][0] * adv_mat_products[0][i] +
                   jacobian[1][1] * adv_mat_products[3][i] +
                   jacobian[2][2] * adv_mat_products[5][i] +
                   (jacobian[0][1] + jacobian[1][0]) * adv_mat_products[1][i] +
                   (jacobian[0][2] + jacobian[2][0]) * adv_mat_products[2][i] +
                   (jacobian[1][2] + jacobian[2][1]) * adv_mat_products[4][i]));
}



template <unsigned int dim, bool has_momentum, bool logarithmic_p>
void
sapphirepp::VFP::UpwindFlux<dim, has_momentum, logarithmic_p>::
  compute_flux_in_p_direction(const double                  n_p,
                              const double                  momentum,
                              const double                  gamma,
                              const dealii::Vector<double> &material_derivative,
                              const dealii::FullMatrix<double> &jacobian,
                              dealii::FullMatrix<double> &positive_flux_matrix,
                              dealii::FullMatrix<double> &negative_flux_matrix)
{
  // compute the matrix sum at the point q at time t. Overwrites the member
  // variable matrix_sum
  compute_matrix_sum(n_p, momentum, gamma, material_derivative, jacobian);
  // compute eigenvalues and eigenvectors. Overwrites the member variables
  // eigenvalues and eigenvectors
  dealii::syevr(jobz,
                range,
                uplo,
                &matrix_size_blas,
                matrix_sum.data(),
                &matrix_size_blas,
                &double_dummy,
                &double_dummy,
                int_dummy,
                int_dummy,
                &abstol,
                &num_eigenvalues,
                eigenvalues.data(),
                eigenvectors.data(),
                &matrix_size_blas,
                isuppz.data(),
                work.data(),
                &lwork,
                iwork.data(),
                &liwork,
                &info);
  Assert(info >= 0, dealii::ExcInternalError());
  if (info != 0)
    std::cerr << "The computation of the eigenvalues and eigenvectors for the "
                 "flux in p direction failed: LAPACK error in syevr"
              << std::endl;

  // split eigenvalues in positive and negative ones
  std::replace_copy_if(
    eigenvalues.begin(),
    eigenvalues.end(),
    positive_eigenvalues.begin(),
    [](const double &val) { return val < 0.0; },
    0.);
  std::replace_copy_if(
    eigenvalues.begin(),
    eigenvalues.end(),
    negative_eigenvalues.begin(),
    [](const double &val) { return val > 0.0; },
    0.);

  // compute the flux matrices
  // see comment in compute_flux_in_space_directions
  double temp_positive = 0;
  double temp_negative = 0;
  for (unsigned int j = 0; j < matrix_size; ++j)
    {
      for (unsigned int l = 0; l < matrix_size; ++l)
        {
          temp_positive =
            positive_eigenvalues[l] * eigenvectors[l * matrix_size + j];
          temp_negative =
            negative_eigenvalues[l] * eigenvectors[l * matrix_size + j];
          for (unsigned int i = 0; i < matrix_size; ++i)
            {
              positive_flux_matrix(i, j) +=
                temp_positive * eigenvectors[l * matrix_size + i];
              negative_flux_matrix(i, j) +=
                temp_negative * eigenvectors[l * matrix_size + i];
            }
        }
    }
}



// explicit instantiation
template class sapphirepp::VFP::UpwindFlux<1, true, true>;
template class sapphirepp::VFP::UpwindFlux<1, true, false>;
template class sapphirepp::VFP::UpwindFlux<1, false, true>;
template class sapphirepp::VFP::UpwindFlux<1, false, false>;
template class sapphirepp::VFP::UpwindFlux<2, true, true>;
template class sapphirepp::VFP::UpwindFlux<2, true, false>;
template class sapphirepp::VFP::UpwindFlux<2, false, true>;
template class sapphirepp::VFP::UpwindFlux<2, false, false>;
template class sapphirepp::VFP::UpwindFlux<3, true, true>;
template class sapphirepp::VFP::UpwindFlux<3, true, false>;
template class sapphirepp::VFP::UpwindFlux<3, false, true>;
template class sapphirepp::VFP::UpwindFlux<3, false, false>;
