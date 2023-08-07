/**
 * @file hd-solver.h
 * @author Florian Schulze (florian.schulze@mpi-hd.mpg.de)
 * @brief Implement HDSolver class to solve the Euler equations
 * @version 0.1
 * @date 2023-08-07
 */

#ifndef HYDROSOLVER_HDSOLVER_H
#define HYDROSOLVER_HDSOLVER_H

#include <deal.II/base/function.h>
#include <deal.II/base/mpi.h>
#include <deal.II/base/parameter_handler.h>
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

#include "config.h"
#include "flux.h"
#include "hd-solver-control.h"
#include "output-module.h"

namespace Sapphire
{
  namespace Hydro
  {
    using namespace dealii;
    using namespace Sapphire::Utils;

    template <int dim>
    class HDSolver
    {
    public:
      HDSolver(const ParameterParser   &prm,
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
      perform_time_step();
      void
      output_results();
      void
      process_results();

      HDInitialCondition<dim> initial_condition;
      HDBoundaryValues<dim>   boundary_values;
      HDExactSolution<dim>    exact_solution;


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

      double       time;
      double       current_time;
      unsigned int timestep_number;

      TimerOutput computing_timer;
    };

  } // namespace Hydro
} // namespace Sapphire
#endif