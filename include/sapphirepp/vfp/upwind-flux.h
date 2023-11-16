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

#ifndef VFP_UPWINDFLUX_H
#define VFP_UPWINDFLUX_H

#include <deal.II/base/point.h>
#include <deal.II/base/tensor.h>

#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/lapack_full_matrix.h>
#include <deal.II/lac/lapack_support.h>
#include <deal.II/lac/vector.h>

#include <vector>

#include "config.h"
#include "particle-functions.h"
#include "pde-system.h"
#include "vfp-parameters.h"

namespace Sapphire
{
  namespace VFP
  {
    template <unsigned int dim, bool has_momentum, bool logarithmic_p>
    class UpwindFlux
    {
    private:
      static constexpr unsigned int dim_cs = dim - has_momentum;

    public:
      UpwindFlux(const PDESystem          &system,
                 const VFPParameters<dim> &solver_control,
                 const PhysicalParameters &physical_parameters);
      void
      set_time(double time);
      void
      compute_upwind_fluxes(
        const std::vector<dealii::Point<dim>>     &q_points,
        const std::vector<dealii::Tensor<1, dim>> &normals,
        std::vector<dealii::FullMatrix<double>>   &positive_flux_matrices,
        std::vector<dealii::FullMatrix<double>>   &negative_flux_matrices);
      void
      test();

    private:
      void
      prepare_work_arrays_for_lapack();
      void
      prepare_upwind_fluxes();
      void
      compute_flux_in_space_directions(
        unsigned int                component,
        const double                n_component,
        const double                background_velocity,
        const double                particle_velocity,
        dealii::FullMatrix<double> &positive_flux_matrix,
        dealii::FullMatrix<double> &negative_flux_matrix);
      void
      compute_matrix_sum(const double                  n_p,
                         const double                  momentum,
                         const double                  gamma,
                         const dealii::Vector<double> &material_derivative,
                         const std::vector<dealii::Vector<double>> &jacobian);
      void
      compute_flux_in_p_direction(
        const double                               n_p,
        const double                               momentum,
        const double                               gamma,
        const dealii::Vector<double>              &material_derivative,
        const std::vector<dealii::Vector<double>> &jacobian,
        dealii::FullMatrix<double>                &positive_flux_matrix,
        dealii::FullMatrix<double>                &negative_flux_matrix);
      // PDE System matrices
      const PDESystem              &pde_system;
      const dealii::types::blas_int matrix_size;

      std::vector<std::vector<double>> advection_matrices;
      std::vector<std::vector<double>> adv_mat_products;

      // Will be used and overwritten in flux computations
      std::vector<double> matrix_sum;
      std::vector<double> eigenvalues;
      std::vector<double> eigenvectors;
      int                 num_eigenvalues;
      std::vector<double> positive_eigenvalues;
      std::vector<double> negative_eigenvalues;

      // Are used in the flux computations in space directions
      std::vector<double>              eigenvalues_advection_matrices;
      std::vector<std::vector<double>> eigenvectors_advection_matrices;

      // Physical input
      BackgroundVelocityField<dim>         background_velocity_field;
      ParticleVelocity<dim, logarithmic_p> particle_velocity_func;
      ParticleGamma<dim, logarithmic_p>    particle_gamma_func;

      // Particle properties
      const double mass;
      const double charge;
      // TransportOnly
      const double velocity;

      // Arguments for the Lapack routine xsyevr
      // Documentation:
      // https://netlib.org/lapack/explore-html/d3/d88/group__real_s_yeigen_ga24155d2da67fb4a896c5f8257589b19f.html#ga24155d2da67fb4a896c5f8257589b19f
      // https://www.intel.com/content/www/us/en/develop/documentation/onemkl-developer-reference-fortran/top/lapack-routines/lapack-least-squares-and-eigenvalue-problem/lapack-least-squares-eigenvalue-problem-driver/symmetric-eigenvalue-problems-lapack-driver/syevr.html#syevr

      // Lapack's error code
      dealii::types::blas_int info;
      // working array containting doubles
      std::vector<double> work;
      // length of the working array -> to be  determined by the first call of
      // syevr
      dealii::types::blas_int lwork;
      // working array containing integers
      std::vector<dealii::types::blas_int> iwork;
      dealii::types::blas_int              liwork;
      // SUPPZ is INTEGER array, dimension at least ( 2*max(1,matrix_size) )
      // The support of the eigenvectors in "eigenvectors" array, i.e., the
      // indices indicating the nonzero elements in "eigenvectors" array
      std::vector<dealii::types::blas_int> isuppz;
      // jobz = 'V' means to compute eigenvalues and eigenvectors
      const char *const jobz;
      // range = 'A' means to compute *all* eigenvalues and eigenvectors
      const char *const range;
      // the input array is an upper triangular matrix
      const char *const uplo;
      // The absolute error tolerance for the eigenvalues.
      //
      // An approximate eigenvalue is accepted as converged when it is
      // determined to lie in an interval [a,b] of width less than or equal to
      //                 ABSTOL + EPS *   max( |a|,|b| ) ,
      // where EPS is the machine precision. If ABSTOL is less than or equal to
      // zero, then EPS*|T| will be used in its place, where |T| is the 1-norm
      // of the tridiagonal matrix obtained by reducing A to tridiagonal form.
      const double abstol;
      // If range = 'A', no upper, nor lower bounds for the eigenvalues are
      // needed. It is still necessary to hand over dummies, which are not
      // referenced when the routine is called.
      const dealii::types::blas_int *const int_dummy;
      const double                         double_dummy;
    };
  } // namespace VFP
} // namespace Sapphire
#endif
