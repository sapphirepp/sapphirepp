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
 * @file vfp-solver.h
 * @author Nils Schween (nils.schween@mpi-hd.mpg.de)
 * @author Florian Schulze (florian.schulze@mpi-hd.mpg.de)
 * @brief Define @ref sapphirepp::VFP::VFPSolver
 */

#ifndef VFP_VFPSOLVER_H
#define VFP_VFPSOLVER_H

#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/function.h>
#include <deal.II/base/index_set.h>
#include <deal.II/base/mpi.h>
#include <deal.II/base/quadrature.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/tensor_function.h>
#include <deal.II/base/timer.h>
#include <deal.II/base/types.h>

#include <deal.II/distributed/shared_tria.h>
#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/mapping_q1.h>

#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/lapack_full_matrix.h>
#include <deal.II/lac/petsc_sparse_matrix.h>
#include <deal.II/lac/petsc_vector.h>
#include <deal.II/lac/sparsity_pattern.h>

#include <deal.II/numerics/vector_tools.h>

#include <mpi.h>

#include <type_traits>

#include "config.h"
#include "output-parameters.h"
#include "pde-system.h"
#include "upwind-flux.h"
#include "vfp-flags.h"
#include "vfp-parameters.h"



namespace sapphirepp
{
  namespace VFP
  {
    using namespace dealii;


    /**
     * @brief This class solves the Vlasov-Fokker-Planck equation
     *
     * Solve the Vlasov-Fokker-Planck equation specified with by the @ref
     * VFPFlags. We distinguish between two different dimensions of the problem:
     *
     *  - The dimension of the \f$ (\mathbf{x}) \f$-space, called configuration
     *    space (`dim_cs`). Depending on the setup, this can be 1, 2 or 3.
     *  - The total dimension of the problem in the reduced phase space \f$
     *    (\mathbf{x}, p) \f$, `dim_ps`. In the momentum depended case this
     *    equals `dim_ps = dim_cs + 1`, in the momentum independent
     *    (transport-only) case `dim_ps = dim_cs`.
     *
     * Since we decompose the \f$ \mathbf{p} \f$ momentum part into spherical
     * harmonics, we always solve for the full 3 dimensions in momentum space.
     *
     * Note that `dim_ps` equals the dimension of the numerical problem, i.e.
     * the grid and solution are of dimension `dim_ps`.
     *
     * @note Due to limitations of @dealii, the total dimension of the problem
     *       must be smaller or equal to 3, `dim_ps <= 3`.
     *
     *
     * @tparam dim Dimension of the reduced phase space \f$ (\mathbf{x}, p) \f$,
     *         `dim_ps`
     */
    template <unsigned int dim>
    class VFPSolver
    {
    public:
      /**
       * Type of triangulation.
       *
       * @note `parallel::distributed:Triangulation` does not allow 1D. To
       *       maintain the possibility to compute 1D scenarios, we replace it
       *       with `parallel::shared::Triangulation`.
       */
      using Triangulation =
        typename std::conditional<dim != 1,
                                  parallel::distributed::Triangulation<dim>,
                                  parallel::shared::Triangulation<dim>>::type;

      /** Is the momentum term activated? */
      static constexpr bool momentum =
        (vfp_flags & VFPFlags::momentum) != VFPFlags::none ? true : false;
      /** Dimension in reduced phase space */
      static constexpr int dim_ps = dim;
      /** Dimension of the configuration space */
      static constexpr int dim_cs = dim - momentum;
      /** Do we use a logarithmic momentum variable? */
      static constexpr bool logarithmic_p =
        (vfp_flags & VFPFlags::linear_p) == VFPFlags::none;



      /**
       * @brief Constructor
       *
       * Constructs the VFP solver with the given parameters. It does not setup
       * the system yet. This is done in the @ref run() method.
       *
       * @param vfp_parameters Parameters for the VFP equation
       * @param physical_parameters User defined parameters of the problem
       * @param output_parameters Parameters for the output
       */
      VFPSolver(const VFPParameters<dim_ps>   &vfp_parameters,
                const PhysicalParameters      &physical_parameters,
                const Utils::OutputParameters &output_parameters);



      /**
       * @brief Solve the VFP equation
       *
       * This method solves the VFP equation with the parameters given in the
       * constructor. It is the main method of this class.
       */
      void
      run();



      /** @{ */
      /**
       * @brief Compute the total global error of the solution with respect to
       *        the exact solution
       *
       * The total global error is computed as
       *
       * \f[
       *    d = \lVert \mathbf{d}_K \rVert_X
       * \f]
       *
       * with the global norm \f$ X \f$ and the cell-wise error
       *
       * \f[
       *    \mathbf{d}_K = \lVert \mathbf{u} - \mathbf{u_h} \rVert_Y
       * \f]
       *
       * where \f$ Y \f$ is the cell-wise norm, \f$ u \f$ denotes the exact
       * solution and \f$ u_h \f$ numerical approximation.
       *
       * Note that this function can only be called **after** @ref run().
       *
       * For more details see
       * @dealref{VectorTools::integrate_difference,namespaceVectorTools,aec4da3324bbce54d7c12dd54c59dd915}
       * and
       * @dealref{VectorTools::compute_global_error,namespaceVectorTools,a21eb62d70953182dcc2b15c4e14dd533}.
       *
       * @param exact_solution Exact/comparison solution of the VFP equation
       *        \f$ u \f$
       * @param cell_norm Cell-wise norm \f$ Y \f$
       * @param global_norm Global norm \f$ X \f$
       * @param weight The additional argument weight allows to evaluate
       *        weighted norms. For details see
       *        @dealref{VectorTools::integrate_difference,namespaceVectorTools,aec4da3324bbce54d7c12dd54c59dd915}.
       * @return The total global error \f$ d \f$
       */
      double
      compute_global_error(
        const Function<dim_ps>         &exact_solution,
        const VectorTools::NormType    &cell_norm,
        const VectorTools::NormType    &global_norm,
        const Function<dim_ps, double> *weight = nullptr) const;

      /**
       * @brief Compute the (weighted) norm of the solution
       *
       * Note that this function can only be called **after** @ref run().
       *
       * For more details see
       * @ref compute_global_error(),
       * @dealref{VectorTools::integrate_difference,namespaceVectorTools,aec4da3324bbce54d7c12dd54c59dd915}
       * and
       * @dealref{VectorTools::compute_global_error,namespaceVectorTools,a21eb62d70953182dcc2b15c4e14dd533}.
       *
       * @param cell_norm Cell-wise norm
       * @param global_norm Global norm
       * @param weight The additional argument weight allows to evaluate
       *        weighted norms. For details see
       *        @dealref{VectorTools::integrate_difference,namespaceVectorTools,aec4da3324bbce54d7c12dd54c59dd915}.
       * @return The weighted norm of the solution
       */
      double
      compute_weighted_norm(
        const VectorTools::NormType    &cell_norm,
        const VectorTools::NormType    &global_norm,
        const Function<dim_ps, double> *weight = nullptr) const;
      /** @} */



      /** @{ */
      /**
       * @brief Get the triangulation object
       *
       * @return const Triangulation&
       */
      const Triangulation &
      get_triangulation() const;

      /**
       * @brief Get the dof handler object
       *
       * @return const DoFHandler<dim>&
       */
      const DoFHandler<dim> &
      get_dof_handler() const;

      /**
       * @brief Get the current (locally relevant) solution
       *
       * @return const PETScWrappers::MPI::Vector&
       */
      const PETScWrappers::MPI::Vector &
      get_current_solution() const;

      /**
       * @brief Get the timer object
       *
       * @return const TimerOutput&
       */
      const TimerOutput &
      get_timer() const;
      /** @} */



    private:
      /** @{ */
      /** VFP parameter */
      const VFPParameters<dim_ps> vfp_parameters;
      /** User defined parameter */
      const PhysicalParameters physical_parameters;
      /**
       * Output parameter
       *
       * @note The object should be constant, but HDF5 output requires non-const
       */
      Utils::OutputParameters output_parameters;
      /** @} */

      /** @{ */
      /** @ref PDESystem */
      PDESystem pde_system;
      /** Maximum expansion order \f$ l_{\rm max} \f$ */
      const unsigned int expansion_order;
      /** Number of expansion coefficients, same as `pde_system.size`*/
      const unsigned int num_exp_coefficients;
      /** @} */

      /** @ref UpwindFlux */
      UpwindFlux<dim_ps, momentum, logarithmic_p> upwind_flux;

      /** MPI communicator */
      const MPI_Comm mpi_communicator;

      /** @dealref{Triangulation}, i.e. Grid for the problem */
      Triangulation triangulation;

      /** @{ */
      /** @dealref{DoFHandler} */
      DoFHandler<dim_ps> dof_handler;

      /** Set of locally owned dofs */
      IndexSet locally_owned_dofs;
      /** Set of locally relevant  dofs */
      IndexSet locally_relevant_dofs;
      /** @}  */

      /**
       * @dealref{Mapping}
       *
       * @note The explicit use of a mapping is most likely related to the usage
       *       of mesh_loop as well
       */
      const MappingQ1<dim_ps> mapping;

      /** @dealref{FESystem} */
      const FESystem<dim_ps> fe;

      /** @{ */
      /**
       * Cell quadrature
       *
       * @note Quadratures are members of this class (and not e.g. part of the
       *       assemble_system method), because I am using the mesh_loop
       * function
       */
      const QGauss<dim_ps> quadrature;
      /** Face quadrature */
      const QGauss<dim_ps - 1> quadrature_face;
      /** @} */

      /**
       * @brief @dealref{AffineConstraints}
       *
       * The constraints object is used for the @ref
       * distribute_local_to_global() function
       */
      const AffineConstraints<double> constraints;


      /** @{ */
      /** @dealref{SparsityPattern} */
      SparsityPattern sparsity_pattern;
      /** Mass matrix */
      PETScWrappers::MPI::SparseMatrix mass_matrix;
      /** DG matrix */
      PETScWrappers::MPI::SparseMatrix dg_matrix;
      /** System matrix, depends on time stepping method */
      PETScWrappers::MPI::SparseMatrix system_matrix;
      /** Source */
      PETScWrappers::MPI::Vector locally_owned_current_source;
      /** System right hand side, depends on time stepping method */
      PETScWrappers::MPI::Vector system_rhs;
      /** @} */

      /** @{ */
      /** Previous solution */
      PETScWrappers::MPI::Vector locally_owned_previous_solution;
      /** Current solution */
      PETScWrappers::MPI::Vector locally_relevant_current_solution;
      /** @} */

      /** @{ */
      /**
       * Output stream for timer output
       *
       * @note Should only be used for TimerOutput
       */
      ConditionalOStream pcout;
      /** @dealref{TimerOutput} */
      TimerOutput timer;
      /** @} */



      /** @{ */
      /**
       * @brief Creates the grid
       *
       * Setup @ref triangulation
       */
      void
      make_grid();

      /**
       * @brief Setup data structure for the linear system
       *
       * Setup the @ref dof_handler
       */
      void
      setup_system();

      /**
       * @brief Assemble the mass matrix
       *
       * Setup and compute the @ref mass_matrix
       */
      void
      assemble_mass_matrix();

      /**
       * @brief Project a function onto the finite element space
       *
       * @param f Function to project
       * @param projected_function Vector returning the projected functions
       */
      void
      project(const Function<dim>        &f,
              PETScWrappers::MPI::Vector &projected_function);
      /** @} */



      /** @{ */
      /**
       * @brief Assemble the DG matrix
       *
       * Compute the @ref dg_matrix
       *
       * @param time Time of the current time step
       */
      void
      assemble_dg_matrix(const double time);

      /**
       * @brief Compute the source term
       *
       * Compute the @ref locally_owned_current_source
       *
       * @param source_function Source function
       */
      void
      compute_source_term(const Function<dim> &source_function);
      /** @} */



      /** @{ */
      /**
       * @brief Calculate one time step with the theta method
       *
       * @param time Current time
       * @param time_step Time step size
       */
      void
      theta_method(const double time, const double time_step);

      /**
       * @brief Calculate one time step with the fourth order explicit
       *        Runge-Kutta method
       *
       * @param time Current time
       * @param time_step Time step size
       */
      void
      explicit_runge_kutta(const double time, const double time_step);

      /**
       * @brief Calculate one time step with the low-storage explicit
       *        Runge-Kutta method
       *
       * @param time Current time
       * @param time_step Time step size
       *
       * @todo The method is not forth order yet
       */
      void
      low_storage_explicit_runge_kutta(const double time,
                                       const double time_step);
      /** @} */


      /**
       * @brief Output the results
       *
       * @param time_step_number time step number
       *
       * @note This function should be a const member, but HDF5 output requires
       *       non-const
       */
      void
      output_results(const unsigned int time_step_number);
    };

  } // namespace VFP
} // namespace sapphirepp
#endif
