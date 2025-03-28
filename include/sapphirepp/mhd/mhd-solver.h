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
 * @file mhd-solver.h
 * @author Florian Schulze (florian.schulze@mpi-hd.mpg.de)
 * @brief Define @ref sapphirepp::MHD::MHDSolver
 */

#ifndef MHD_MHDSOLVER_H
#define MHD_MHDSOLVER_H

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
#include "mhd-equations.h"
#include "mhd-flags.h"
#include "mhd-parameters.h"
#include "numerical-flux.h"
#include "output-parameters.h"
#include "slope-limiter.h"



namespace sapphirepp
{
  namespace MHD
  {
    using namespace dealii;



    /**
     * @brief This class solves the MHD equations
     *
     * @tparam dim Dimension of the configuration space \f$ (\mathbf{x}) \f$,
     *         `dim`
     */
    template <unsigned int dim>
    class MHDSolver
    {
    public:
      /**
       * Use Lagrange multiplier \f$ \psi \f$
       * for hyperbolic divergence cleaning.
       */
      static constexpr bool divergence_cleaning =
        (mhd_flags & MHDFlags::hyperbolic_divergence_cleaning) !=
        MHDFlags::none;

      /** @ref MHDEquations::n_components */
      static constexpr unsigned int n_components =
        MHDEquations<dim, divergence_cleaning>::n_components;
      /** @ref MHDEquations::density_component */
      static constexpr unsigned int density_component =
        MHDEquations<dim, divergence_cleaning>::density_component;
      /** @ref MHDEquations::first_momentum_component */
      static constexpr unsigned int first_momentum_component =
        MHDEquations<dim, divergence_cleaning>::first_momentum_component;
      /** @ref MHDEquations::energy_component */
      static constexpr unsigned int energy_component =
        MHDEquations<dim, divergence_cleaning>::energy_component;
      /** @ref MHDEquations::first_magnetic_component */
      static constexpr unsigned int first_magnetic_component =
        MHDEquations<dim, divergence_cleaning>::first_magnetic_component;
      /** @ref MHDEquations::divergence_cleaning_component */
      static constexpr unsigned int divergence_cleaning_component =
        MHDEquations<dim, divergence_cleaning>::divergence_cleaning_component;

      /** @ref MHDEquations::state_type */
      using state_type =
        typename MHDEquations<dim, divergence_cleaning>::state_type;
      /** @ref MHDEquations::flux_type */
      using flux_type =
        typename MHDEquations<dim, divergence_cleaning>::flux_type;


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



      /**
       * @brief Constructor
       *
       * Constructs the MHD solver with the given parameters. It does not setup
       * the system yet. This is done in the @ref setup() method.
       *
       * @param mhd_parameters Parameters for the MHD equation
       * @param physical_parameters User defined parameters of the problem
       * @param output_parameters Parameters for the output
       */
      MHDSolver(const MHDParameters<dim>      &mhd_parameters,
                const PhysicalParameters      &physical_parameters,
                const Utils::OutputParameters &output_parameters);



      /** @{ */
      /**
       * @brief Setup the MHD solver
       *
       * This methods must be called before any other function. It sets up the
       * grid, the dof_handler, the linear system, and the inital conditions.
       */
      void
      setup();

      /**
       * @brief Solve the MHD equations
       *
       * This method solves the MHD equations with the parameters given in the
       * constructor. It is the main method of this class.
       */
      void
      run();
      /** @} */



      /**
       * @name Time stepping methods
       * @{
       */
      /**
       * @brief Calculate one time step with the explicit Euler method
       *
       * @param time Current time
       * @param time_step Time step size
       */
      void
      forward_euler_method(const double time, const double time_step);

      /**
       * @brief Calculate one time step with a strong stability-preserving
       *        Runge-Kutta method.
       *
       * @param time Current time
       * @param time_step Time step size
       */
      void
      explicit_runge_kutta(const double time, const double time_step);
      /** @} */



      /**
       * @brief Output the results
       *
       * @param time_step_number time step number
       * @param cur_time Simulation time
       *
       * @note This function should be a const member, but HDF5 output requires
       *       non-const
       */
      void
      output_results(const unsigned int time_step_number,
                     const double       cur_time);



      /**
       * @name Utility functions
       * @{
       */
      /**
       * @brief Project a function onto the finite element space
       *
       * Note that this function can only be called **after** the setup in
       * @ref setup().
       *
       * @param f Function to project
       * @param projected_function Vector returning the projected functions
       */
      void
      project(const Function<dim>        &f,
              PETScWrappers::MPI::Vector &projected_function) const;


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
       * Note that this function can only be called **after** @ref setup().
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
      compute_global_error(const Function<dim>         &exact_solution,
                           const VectorTools::NormType &cell_norm,
                           const VectorTools::NormType &global_norm,
                           const Function<dim, double> *weight = nullptr) const;

      /**
       * @brief Compute the (weighted) norm of the solution
       *
       * Note that this function can only be called **after** @ref setup().
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
        const VectorTools::NormType &cell_norm,
        const VectorTools::NormType &global_norm,
        const Function<dim, double> *weight = nullptr) const;
      /** @} */



      /**
       * @name Accessors
       * @{
       */
      /**
       * @brief Get the MHDEquations object
       *
       * @return const MHDEquations&
       */
      const MHDEquations<dim, divergence_cleaning> &
      get_mhd_equations() const;

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
       * @brief Get the cell average of the solution
       *
       * @return const std::vector<Vector<double>>&
       */
      const std::vector<Vector<double>> &
      get_cell_average() const;

      /**
       * @brief Get the cell average of one component the solution
       *
       * @return const Vector<double>
       */
      Vector<double>
      get_cell_average_component(unsigned int component) const;

      /**
       * @brief Get the shock_indicator
       *
       * @return const Vector<double>&
       */
      const Vector<double> &
      get_shock_indicator() const;

      /**
       * @brief Get the positivity_limiter_indicator
       *
       * @return const Vector<float>&
       */
      const Vector<float> &
      get_positivity_limiter_indicator() const;

      /**
       * @brief Get the magnetic_divergence
       *
       * @return const Vector<double>&
       */
      const Vector<double> &
      get_magnetic_divergence() const;

      /**
       * @brief Get the maximum CFL time step on each cell
       *
       * @return const Vector<double>&
       */
      const Vector<double> &
      get_cell_dt() const;

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
      /** MHD parameter */
      const MHDParameters<dim> mhd_parameters;
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
      /** @ref MHDEquations */
      MHDEquations<dim, divergence_cleaning> mhd_equations;
      /** @ref NumericalFlux */
      NumericalFlux<dim, divergence_cleaning> numerical_flux;
      /** @ref SlopeLimiter */
      SlopeLimiter<dim, divergence_cleaning> slope_limiter;
      /** @} */

      /** MPI communicator */
      const MPI_Comm mpi_communicator;

      /** @dealref{Triangulation}, i.e. Grid for the problem */
      Triangulation triangulation;

      /** @{ */
      /** @dealref{DoFHandler} */
      DoFHandler<dim> dof_handler;

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
      const MappingQ1<dim> mapping;

      /** @dealref{FESystem} */
      const FESystem<dim> fe;

      /** @{ */
      /**
       * Cell quadrature
       *
       * @note Quadratures are members of this class (and not e.g. part of the
       *       assemble_system method), because I am using the mesh_loop
       * function
       */
      const QGauss<dim> quadrature;
      /** Face quadrature */
      const QGauss<dim - 1> quadrature_face;
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
      /** DG right-hand-side */
      PETScWrappers::MPI::Vector dg_rhs;
      /** System right hand side, depends on time stepping method */
      PETScWrappers::MPI::Vector system_rhs;
      /** @} */

      /** @{ */
      /** Locally owned solution */
      PETScWrappers::MPI::Vector locally_owned_solution;
      /** Current solution */
      PETScWrappers::MPI::Vector locally_relevant_current_solution;
      /** Component wise cell average of the solution. */
      std::vector<Vector<double>> cell_average;
      /** Indicator for a shock in the solution. */
      Vector<double> shock_indicator;
      /** Indicator where the solution is positivity limited. */
      Vector<float> positivity_limiter_indicator;
      /** Magnetic divergence \f$ \nable \cdot \mathbf{B} \f$ per cell. */
      Vector<double> magnetic_divergence;
      /** @} */

      /** @{ */
      /** Maximum CFL time step on each cell. */
      Vector<double> cell_dt;
      /** Global maximum CLF time step. */
      double global_dt_cfl;
      /** Global minimum cell size. */
      double global_dx_min;
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
      /** @} */



      /** @{ */
      /**
       * @brief Assemble the right hand side of the DG system.
       *
       * Compute the @ref dg_rhs using the @ref
       * locally_relevant_current_solution.
       *
       * @param time Time of the current time step
       */
      void
      assemble_dg_rhs(const double time);
      /** @} */



      /** @{ */
      /**
       * @brief Compute CFL condition for the system at the current time step.
       *
       * The CFL condition is given by,
       * \f[
       *   \Delta t_{\mathrm{CFL}} =
       *   \mathcal{C} \, \frac{h}{(2k+1) \lambda_{\mathrm{max}}} \,,
       * \f]
       * where \f$ h \f$ is the cell size,
       * \f$ k \f$ the polynomial degree of the basis functions,
       * \f$ \lambda_{\mathrm{max}} \f$ the maximum eigenvalue
       * and \f$ \mathcal{C} \lesssim 1 \f$ the Courant/CFL number.
       *
       * This function calculates the CFL condition for each cell, @ref cell_dt,
       * and updates @ref global_dt_cfl
       * as well as @ref global_dx_min accordingly.
       */
      void
      compute_cfl_condition();
      /** @} */


      /**
       * @name Limiter
       * @{
       */
      /**
       * @brief Computes @ref cell_average of the current solution
       */
      void
      compute_cell_average();

      /**
       * @brief Computes the @ref shock_indicator using the KXRCF method
       *
       * @note This function can only be called after updating the
       *       @ref cell_average using  @ref compute_cell_average().
       */
      void
      compute_shock_indicator();


      /**
       * @brief Indicates if the solution needs positivity limiting.
       *
       * @param states Vector of states on the cell.
       * @return bool `true` if one of the state is in an non-admissible state,
       *              `false` otherwise.
       */
      bool
      indicate_positivity_limiting(const std::vector<state_type> &states) const;


      /**
       * @brief Retrieves the precomputed cell average for a given cell.
       *
       * @param cell The cell for which to compute the average.
       * @return Returns the cell average of the
       *         @ref state_type "MHD state".
       */
      inline state_type
      get_cell_average(
        const typename dealii::DoFHandler<dim>::cell_iterator &cell) const
      {
        if (cell->is_active())
          {
            Assert(!cell->is_artificial(),
                   dealii::ExcMessage("Cell average is only calculated for "
                                      "locally active or ghost cells."));

            return cell_average[cell->active_cell_index()];
          }
        else
          {
            state_type avg(n_components);
            double     measure = 0;
            for (const auto &child_cell : cell->child_iterators())
              {
                avg.add(child_cell->measure(), get_cell_average(cell));
                measure += child_cell->measure();
              }
            avg /= measure;
            return avg;
          }
      }


      /**
       * @brief Applys a slope limiter to the current solution.
       */
      void
      apply_limiter();
      /** @} */



      /**
       * @brief Computes @ref magentic_divergence of the current solution
       */
      void
      compute_magnetic_divergence();
    };

  } // namespace MHD
} // namespace sapphirepp
#endif
