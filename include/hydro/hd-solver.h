/**
 * @file hd-solver.h
 * @author Florian Schulze (florian.schulze@mpi-hd.mpg.de)
 * @brief Implement HDSolver class to solve the Euler equations
 * @version 0.1
 * @date 2023-08-07
 */

#ifndef HYDROSOLVER_HDSOLVER_H
#define HYDROSOLVER_HDSOLVER_H

#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/function.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/time_stepping.h>
#include <deal.II/base/timer.h>
#include <deal.II/base/utilities.h>
#include <deal.II/base/vectorization.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_system.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/la_parallel_vector.h>

#include <deal.II/matrix_free/fe_evaluation.h>
#include <deal.II/matrix_free/matrix_free.h>
#include <deal.II/matrix_free/operators.h>

#include <deal.II/numerics/data_out.h>

#include <fstream>
#include <iomanip>
#include <iostream>

#include "config.h"
#include "euler-operator.h"
#include "flux.h"
#include "hd-solver-control.h"
#include "output-module.h"
#include "parameter-parser.h"
#include "postprocessor.h"
#include "sapphire-logstream.h"
#include "time-stepping.h"

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
               const OutputModule<dim> &output_module);

      // void
      // do_timestep(); //TODO_HD: implement this function
      void
      run();

    private:
      void
      make_grid_and_dofs();

      void
      output_results(const unsigned int result_number);

      HDInitialCondition<dim> initial_condition; // unused
      HDExactSolution<dim>    exact_solution;    // unused

      HDSolverControl   hd_solver_control; // unused
      OutputModule<dim> output_module;

      LinearAlgebra::distributed::Vector<Number> solution;

      ConditionalOStream pcout; // TODO_HD: replace with saplog

#if dim > 1
      parallel::distributed::Triangulation<dim> triangulation;
#else
      Triangulation<dim> triangulation;
#endif

      FESystem<dim>   fe;
      MappingQ<dim>   mapping;
      DoFHandler<dim> dof_handler;

      TimerOutput timer;

      EulerOperator<dim, fe_degree, n_q_points_1d> euler_operator;

      double time, time_step;
    };

  } // namespace Hydro
} // namespace Sapphire
#endif