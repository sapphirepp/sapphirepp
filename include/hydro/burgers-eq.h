/**
 * @file burgers-eq.h
 * @author Florian Schulze (florian.schulze@mpi-hd.mpg.de)
 * @brief Solve Burgers' equation
 * @version 0.1
 * @date 2023-07-12
 *
 *
 * Implement the class BurgersEq that solves the 1d Burgers' equation
 * \f$ \frac{\partial u}{\partial t} + \beta \, \frac{\partial
 * u^2}{\partial x} = 0 \f$.
 */

#ifndef HYDROSOLVER_BURGERSEQ_H
#define HYDROSOLVER_BURGERSEQ_H

#include <deal.II/base/function.h>
#include <deal.II/base/mpi.h>
#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/smartpointer.h>
#include <deal.II/base/tensor_function.h>
#include <deal.II/base/timer.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_q1.h>

#include <deal.II/grid/tria.h>

#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/vector.h>

#include <deal.II/numerics/data_out.h>

#include <mpi.h>

#include "flux.h"
#include "numerics.h"
#include "output-module.h"

namespace Sapphire
{
  namespace Hydro
  {
    using namespace dealii;
    using namespace Sapphire::Utils;

    /**
     * @brief Solve Burgers' equation.
     *
     * This class solves the (1d) Burgers' equation
     * \f$ \frac{\partial u}{\partial t} + \beta \, \frac{\partial
     * u^2}{\partial x} = 0 \f$
     * where \f$ u(\mathbf{x}, t) \f$ is the solution and the flux is
     * \f$ \mathbf{f}(u) = \beta u^2 \f$, with $\beta$ a constant. The
     * initial condition is given by \f$ u(\mathbf{x}, 0) = u_0(\mathbf{x}) \f$
     * and the boundary condition is given by \f$ u(\mathbf{x}, t) =
     * u_b(\mathbf{x}, t) \f$.
     *
     * \tparam dim dimension of the problem (must be 1)
     */
    template <int dim>
    class BurgersEq
    {
    public:
      /**
       * @brief Construct a new Burgers Eq object
       *
       *
       * @param initial_condition initial condition \f$ u_0(\mathbf{x}) \f$
       * @param boundary_values boundary values \f$ u_b(\mathbf{x}, t) \f$
       * @param exact_solution exact solution for comparison \f$ u(\mathbf{x}, t)
       * \f$
       */
      BurgersEq(Function<dim>           *initial_condition,
                Function<dim>           *boundary_values,
                Function<dim>           *exact_solution,
                const ParameterParser   &prm,
                const OutputModule<dim> &output_module,
                const double             beta = 1.0);

      void
      init();
      void
      do_timestep();
      void
      run();

    private:
      void
      make_grid();
      void
      setup_system();
      void
      assemble_mass_matrix();
      void
      assemble_dg_vector();
      void
      assemble_system();
      void
      slope_limiter();
      void
      perform_time_step();
      /**
       * @brief Solve the linear system
       *
       * This function solves the linear system \f$ M \mathbf{u} = \mathbf{b}
       * \f$ where \f$ M \f$ is the mass matrix, \f$ \mathbf{u} \f$ is the
       * solution vector and \f$ \mathbf{b} \f$ is the right hand side vector.
       */
      void
      solve_linear_system();
      void
      output_results();
      void
      process_results();

      const SmartPointer<Function<dim>> initial_condition;
      const SmartPointer<Function<dim>> boundary_values;
      const SmartPointer<Function<dim>> exact_solution;

      HDSolverControl   hd_solver_control;
      Flux<dim>         flux;
      OutputModule<dim> output_module;

      const double beta; //< factor in front of the flux

      MPI_Comm mpi_communicator;

      Triangulation<dim>   triangulation;
      const MappingQ1<dim> mapping;

      FE_DGQ<dim>           fe;
      DoFHandler<dim>       dof_handler;
      const QGauss<dim>     quadrature_formula;
      const QGauss<dim - 1> face_quadrature_formula;

      AffineConstraints<double> constraints;

      SparsityPattern      sparcity_pattern;
      SparseMatrix<double> mass_matrix;
      SparseMatrix<double> system_matrix;

      Vector<double> solution;

      // solution at the current intermediate time step, used to calculate the
      // flux
      Vector<double> current_solution;

      // solution of the last time step
      Vector<double> old_solution;

      // vector of the fluxes at the current intermediate time step
      Vector<double> dg_vector;

      // right hand side of the linear system \f$ \mathbf{b} \f$
      Vector<double> system_rhs;

      std::vector<float> error_with_time;
      Vector<float>      mark_cell_for_limiter;

      double       time;
      double       current_time;
      unsigned int timestep_number;

      ConditionalOStream pcout;
      TimerOutput        computing_timer;
    };

  } // namespace Hydro
} // namespace Sapphire
#endif