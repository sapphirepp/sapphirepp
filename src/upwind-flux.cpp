#include "upwind-flux.h"

#include <deal.II/base/exceptions.h>
#include <deal.II/lac/lapack_support.h>
#include <deal.II/lac/lapack_templates.h>  // direct access to the Fortran driver routines of Lapack
#include <deal.II/lac/vector.h>

#include <algorithm>

template <int dim>
VFPEquation::UpwindFlux<dim>::UpwindFlux(const PDESystem &system, bool momentum)
    : pde_system{system},
      matrix_size(pde_system.system_size()),
      advection_matrices(3),
      adv_mat_products(6),
      matrix_sum(matrix_size * matrix_size),
      eigenvalues(matrix_size),
      eigenvectors(matrix_size * matrix_size),
      positive_eigenvalues(matrix_size),
      negative_eigenvalues(matrix_size),
      eigenvectors_advection_matrices(
          3, std::vector<double>(matrix_size * matrix_size)),
      eigenvalues_advection_matrices(matrix_size),
      isuppz(2 * matrix_size),
      jobz{&dealii::LAPACKSupport::V},
      range{&dealii::LAPACKSupport::A},
      uplo{&dealii::LAPACKSupport::U},
      abstol{0.},  // see Documentation of xsyever
      int_dummy{&dealii::LAPACKSupport::one},
      double_dummy{1.},
      momentum{momentum} {
  // NOTE: Since we very often call compute_matrix_sum and the matrixes classes
  // of dealii do not allow unchecked access to there raw data, we create copies
  // of the matrices in the hope that this additional memory consumption is made
  // up for by perfomance gains.
  const std::vector<dealii::LAPACKFullMatrix<double>> &adv_mat =
      pde_system.get_advection_matrices();
  for (unsigned int k = 0; k < 3; ++k) {
    for (int i = 0; i < matrix_size; ++i)
      for (int j = 0; j < matrix_size; ++j)
        // NOTE: When computing the eigenvalues of matrix_sum, I call a Fortran
        // routine which requires the matrix elements to be stored in
        // column-major order. But since the resulting matrix is symmteric,
        // there is no difference between column-major order and row-major
        // order. I stick with the C convention and use row-major order.
        advection_matrices[k][i * matrix_size + j] = adv_mat[k](i, j);
  }
  const std::vector<dealii::LAPACKFullMatrix<double>> &adv_mat_prod =
      pde_system.get_adv_mat_products();
  for (unsigned int k = 0; k < 6; ++k) {
    for (int i = 0; i < matrix_size; ++i)
      for (int j = 0; j < matrix_size; ++j)
        adv_mat_products[k][i * matrix_size + j] = adv_mat_prod[k](i, j);
  }

  prepare_work_arrays_for_lapack();
  prepare_upwind_fluxes();
}

template <int dim>
void VFPEquation::UpwindFlux<dim>::set_time(double time) {
  background_velocity_field.set_time(time);
}

template <int dim>
void VFPEquation::UpwindFlux<dim>::compute_upwind_fluxes(
    const std::vector<dealii::Point<dim>> &q_points,
    const std::vector<dealii::Tensor<1, dim>> &normals,
    std::vector<dealii::FullMatrix<double>> &positive_flux_matrices,
    std::vector<dealii::FullMatrix<double>> &negative_flux_matrices) {
  // Determine if we are computing the flux of an interior face whose normal
  // points into the x,y,z or p direction. Since we are using a rectangular
  // grid, the normal is the same for all quadrature points and it points either
  // in the x, y, z or p direction. Hence it can determined outside the loop
  // over the quadrature points

  unsigned int dim_cs = dim - momentum;
  Coordinate component = none;  // Produces an error if not overwritten

  for (unsigned int i = 0; i < dim_cs; ++i)
    // NOTE: Such a simple comparison of doubles is only possible, because I
    // know that I am comparing 0. with 1. or something very close to one with
    // one. If this ever fails a more sophisticated comparison is needed, e.g.
    // https://floating-point-gui.de/errors/comparison/
    // std::abs is nessecary because of the boundary faces, i.e. <-| n_k = -1
    // for some faces.
    if (std::abs(1. - std::abs(normals[0][i])) < 1e-5) component = i;

  // If the momentum terms are included, then the last component of normals is
  // the momentum direction
  if (momentum && std::abs(1. - std::abs(normals[0][dim])) < 1e-5)
    component = p;

  // Flux in the spatial directions
  if (component == x || component == y || component == z) {
    // Background velocity field
    // Get the value of the velocity field at every quadrature point
    // TODO: Value list for a specific component would be faster.
    std::vector<dealii::Vector<double>> velocities(q_points.size(),
                                                   dealii::Vector<double>(3));
    background_velocity_field.vector_value_list(q_points, velocities);

    // Particle
    std::vector<double> particle_velocities(q_points.size());
    if (momentum) {
      // if distribution functions depends on p, the particle velocity depends
      // on the position in the phase space and needs to be computed
      particle_velocity_func.value_list(q_points, particle_velocities);
    } else {
      std::fill(particle_velocities.begin(), particle_velocities.end(),
                particle_properties.velocity);
    }
    // Compute flux at every quadrature point
    for (unsigned int q_index; q_index < q_points.size(); ++q_index)
      compute_flux_in_space_directions(
          normals[component], velocities[q_index], particle_velocities[q_index],
          positive_flux_matrices[q_index], negative_flux_matrices[q_index]);
  } else if (component == p) {
    std::vector<dealii::Vector<double>> material_derivative_vel(
        q_points.size(), dealii::Vector<double>(3));
    background_velocity_field.material_derivative_list(q_points,
                                                       material_derivative_vel);

    std::vector<std::vector<dealii::Vector<double>>> jacobians_vel(
        q_points.size(),
        std::vector<dealii::Vector<double>>(3, dealii::Vector<double>(3)));
    background_velocity_field.jacobian_list(q_points, jacobians_vel);

    std::vector<double> particle_gammas(q_points.size());
    particle_gamma_func.value_list(q_points, particle_gammas);
    for (unsigned int q_index = 0; q_index < q_points.size(); ++q_index)
      compute_flux_in_p_direction(
          normals[p], q_points[p], particle_gammas[q_index],
          material_derivative_vel[q_index], jacobians_vel[q_index],
          positive_flux_matrices[q_index], negative_flux_matrices[q_index]);
  }
}

template <int dim>
void VFPEquation::UpwindFlux<dim>::prepare_work_arrays_for_lapack() {
  // Preparations for the eigenvalue and eigenvector computations
  //
  // The computation of the upwind flux in the p direction requires to the
  // compute eigenvalues and eigenvectors for every face with n_p = 1 (or -1).
  // Lapack needs correctly sized work arrays to do these computations. It can
  // determine the sizes it needs with a specific call to the used eigenvalue
  // routine. The work arrays will be allocated only once and then reused every
  // time eigenvalues and eigenvectors are computed.

  // https://stackoverflow.com/questions/46618391/what-is-the-use-of-the-work-parameters-in-lapack-routines
  // Quote: "If dynamic allocation is available, the most common use of LAPACK
  // functions is to first perform a workspace query using LWORK = -1, then use
  // the return value to allocate a WORK array of the correct size and finally
  // call the routine of LAPACK to get the expected result. High-end wrappers of
  // LAPACK such as LAPACKE features function doing just that: take a look at
  // the source of LAPACKE for function LAPACKE_dsyev()! It calls twice the
  // function LAPACKE_dsyev_work, which calls LAPACK_dsyev (wrapping dsyev()).

  // Wrappers still feature functions such as LAPACKE_dsyev_work(), where the
  // arguments work and lwork are still required. The number of allocations can
  // therefore be reduced if the routine is called multiple times on similar
  // sizes by not deallocating WORK between calls, but the user must do that
  // himself (see this example). In addition, the source of ILAENV, the function
  // of LAPACK called to compute the optimzed size of WORK, features the
  // following text:

  //     This version provides a set of parameters which should give good, but
  //     not optimal, performance on many of the currently available computers.
  //     Users are encouraged to modify this subroutine to set the tuning
  //     parameters for their particular machine using the option and problem
  //     size information in the arguments.

  // As a result, testing sizes of WORK larger than the size returned by the
  // workspace query could improve performances."

  // Compute the size of the working arrays for the Lapack routine xsyevr

  // The optimal size of the work arrays will be stored in their first element
  work.resize(1);
  iwork.resize(1);
  lwork = -1;
  liwork = -1;
  // create dummy variables for the matrix sum
  double n_p_dummy = 1.;
  double gamma_dummy = 1.;
  double p_dummy = 1.;
  dealii::Vector<double> material_derivative_dummy{1., 1., 1.};
  std::vector<dealii::Vector<double>> jacobian_dummy(
      3, dealii::Vector<double>{1., 1., 1.});
  compute_matrix_sum(n_p_dummy, p_dummy, gamma_dummy, material_derivative_dummy,
                     jacobian_dummy);
  // call Lapack routine
  dealii::syevr(jobz, range, uplo, &matrix_size, matrix_sum.data(),
                &matrix_size, &double_dummy, &double_dummy, int_dummy,
                int_dummy, &abstol, &num_eigenvalues, eigenvalues.data(),
                eigenvectors.data(), &matrix_size, isuppz.data(), work.data(),
                &lwork, iwork.data(), &liwork, &info);
  Assert(info == 0, dealii::ExcInternalError());

  // resize the working arrays
  lwork = static_cast<dealii::types::blas_int>(work[0]);
  liwork = iwork[0];
  work.resize(lwork);
  iwork.resize(liwork);
}

template <int dim>
void VFPEquation::UpwindFlux<dim>::prepare_upwind_fluxes() {
  // The eigenvalues of A_x are also the eigenvalues of A_y and A_z. The
  // eigenvalues of A_x are the roots of the associated legendre polynomials.
  // Since the associated Legendre Polynomials are orthogonal polynomials, it
  // can be shown that all roots are contained in in the intervall [-1., 1.].
  // Physically this makes sense, because c = 1 and the eigenvalues encode the
  // speed of 'information' transport.
  //
  // Helper variables
  const char *const V = &dealii::LAPACKSupport::V;
  const double vl = -1.1;
  const double vu = 1.1;
  // NOTE: We need another copy of the advection matrix, because it gets
  // destroyed when the xsyevr is called.
  std::vector<double> A_x;
  A_x = advection_matrices[0];
  dealii::syevr(jobz, V, uplo, &matrix_size, A_x.data(), &matrix_size, &vl, &vu,
                int_dummy, int_dummy, &abstol, &num_eigenvalues,
                eigenvalues_advection_matrices.data(),
                eigenvectors_advection_matrices[0].data(), &matrix_size,
                isuppz.data(), work.data(), &lwork, iwork.data(), &liwork,
                &info);
  // NOTE: The eigenvalues of A_x can computed very efficiently because A_x is
  // (when ordered correctly) a tridiagonal symmetric matrix, i.e. it is in
  // Hessenberg-Form. At the moment it is not ordered correctly. The computation
  // has only to be done once.
  Assert(info >= 0, dealii::ExcInternalError());
  if (info != 0)
    std::cerr << "The computation of the eigenvalues and eigenvectors of A_x "
                 "failed: LAPACK error in syevr"
              << std::endl;
  std::vector<double> A_y;
  A_y = advection_matrices[1];
  dealii::syevr(jobz, V, uplo, &matrix_size, A_y.data(), &matrix_size, &vl, &vu,
                int_dummy, int_dummy, &abstol, &num_eigenvalues,
                eigenvalues_advection_matrices.data(),
                eigenvectors_advection_matrices[1].data(), &matrix_size,
                isuppz.data(), work.data(), &lwork, iwork.data(), &liwork,
                &info);
  Assert(info >= 0, dealii::ExcInternalError());
  if (info != 0)
    std::cerr << "The computation of the eigenvalues and eigenvectors of A_y "
                 "failed: LAPACK error in syevr"
              << std::endl;

  std::vector<double> A_z;
  A_z = advection_matrices[2];
  dealii::syevr(jobz, V, uplo, &matrix_size, A_z.data(), &matrix_size, &vl, &vu,
                int_dummy, int_dummy, &abstol, &num_eigenvalues,
                eigenvalues_advection_matrices.data(),
                eigenvectors_advection_matrices[2].data(), &matrix_size,
                isuppz.data(), work.data(), &lwork, iwork.data(), &liwork,
                &info);
  Assert(info >= 0, dealii::ExcInternalError());
  if (info != 0)
    std::cerr << "The computation of the eigenvalues and eigenvectors of A_z "
                 "failed: LAPACK error in syevr"
              << std::endl;

  // TODO: Rotate the eigenvectors of A_x to get the eigenvectors of A_y and
  // A_z, i.e. implement the rotation matrices e^{-i\Omega_x pi/2} and e^{-i
  // \Omega_z pi/2} For now we will three times compute the same eigenvalues to
  // get the eigenvectors of A_y and A_z
}

template <int dim>
void VFPEquation::UpwindFlux<dim>::compute_flux_in_space_directions(
    const Coordinate component, const double n_component,
    const double background_velocity, const double particle_velocity,
    dealii::FullMatrix<double> &positive_flux_matrix,
    dealii::FullMatrix<double> &negative_flux_matrix) {
  // Create a copy of the eigenvalues of the advection matrices to compute the
  // positive and negative fluxes at point q at time t
  eigenvalues = eigenvalues_advection_matrices;
  // Update the eigenvalues
  std::for_each(eigenvalues.begin(), eigenvalues.end(), [=](double &lambda) {
    lambda = n_component * (background_velocity + particle_velocity * lambda);
  });

  // split eigenvalues in positive and negative ones
  std::replace_copy_if(
      eigenvalues.begin(), eigenvalues.end(), positive_eigenvalues.begin(),
      std::bind(std::less<double>(), std::placeholders::_1, 0.), 0.);
  std::replace_copy_if(
      eigenvalues.begin(), eigenvalues.end(), negative_eigenvalues.begin(),
      std::bind(std::greater<double>(), std::placeholders::_1, 0.), 0.);

  // compute the flux matrices
  for (unsigned int i = 0; i < matrix_size; ++i)
    for (unsigned int j = 0; j < matrix_size; ++j) {
      // NOTE: We compute the triple matrix product = V * Lambda_{+/-} V^T. The
      // eigenvectors V are computed with the Lapack routine xsyevr. Since this
      // routine is a Fortran routine it returns an array V in column-major
      // order instead of row-major order. This is why it looks like the roles
      // of V^T and V are exchanged.
      positive_flux_matrix(i, j) =
          eigenvectors_advection_matrices[component][j * matrix_size + i] *
          positive_eigenvalues[i] *
          eigenvectors_advection_matrices[component][i * matrix_size + j];
      negative_flux_matrix(i, j) =
          eigenvectors_advection_matrices[component][j * matrix_size + i] *
          negative_eigenvalues[i] *
          eigenvectors_advection_matrices[component][i * matrix_size + j];
    }
}

template <int dim>
void VFPEquation::UpwindFlux<dim>::compute_matrix_sum(
    const double n_p, const double momentum, const double gamma,
    const dealii::Vector<double> &material_derivative,
    const std::vector<dealii::Vector<double>> &jacobian) {
  for (int i = 0; i < matrix_size * matrix_size; ++i)
    matrix_sum[i] =
        -n_p *
        (gamma * (material_derivative[0] * advection_matrices[0][i] +
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

template <int dim>
void VFPEquation::UpwindFlux<dim>::compute_flux_in_p_direction(
    const double n_p, const double momentum, const double gamma,
    const dealii::Vector<double> &material_derivative,
    const std::vector<dealii::Vector<double>> &jacobian,
    dealii::FullMatrix<double> &positive_flux_matrix,
    dealii::FullMatrix<double> &negative_flux_matrix) {
  // compute the matrix sum at the point q at time t. Overwrites the member
  // variable matrix_sum
  compute_matrix_sum(n_p, momentum, gamma, material_derivative, jacobian);
  // compute eigenvalues and eigenvectors. Overwrites the member variables
  // eigenvalues and eigenvectors
  dealii::syevr(jobz, range, uplo, &matrix_size, matrix_sum.data(),
                &matrix_size, &double_dummy, &double_dummy, int_dummy,
                int_dummy, &abstol, &num_eigenvalues, eigenvalues.data(),
                eigenvectors.data(), &matrix_size, isuppz.data(), work.data(),
                &lwork, iwork.data(), &liwork, &info);
  Assert(info >= 0, dealii::ExcInternalError());
  if (info != 0)
    std::cerr << "The computation of the eigenvalues and eigenvectors for the "
                 "flux in p direction failed: LAPACK error in syevr"
              << std::endl;

  // split eigenvalues in positive and negative ones
  std::replace_copy_if(
      eigenvalues.begin(), eigenvalues.end(), positive_eigenvalues.begin(),
      std::bind(std::less<double>(), std::placeholders::_1, 0.), 0.);
  std::replace_copy_if(
      eigenvalues.begin(), eigenvalues.end(), negative_eigenvalues.begin(),
      std::bind(std::greater<double>(), std::placeholders::_1, 0.), 0.);

  // compute the flux matrices
  for (unsigned int i = 0; i < matrix_size; ++i)
    for (unsigned int j = 0; j < matrix_size; ++j) {
      // NOTE: see comment in compute_flux_in_space_directions
      positive_flux_matrix(i, j) = eigenvectors[j * matrix_size + i] *
                                   positive_eigenvalues[i] *
                                   eigenvectors[i * matrix_size + j];
      negative_flux_matrix(i, j) = eigenvectors[j * matrix_size + i] *
                                   negative_eigenvalues[i] *
                                   eigenvectors[i * matrix_size + j];
    }
}
