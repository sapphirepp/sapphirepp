#include "hd-solver.h"

#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_interface_values.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_q1.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/tria.h>

#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/precondition_block.h>
#include <deal.II/lac/solver_richardson.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/vector.h>

#include <deal.II/meshworker/mesh_loop.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/vector_tools.h>

#include <fstream>
#include <iostream>

#include "hd-solver-control.h"
#include "sapphire-logstream.h"

namespace Sapphire
{
  namespace Hydro
  {
    using namespace dealii;
    template <int dim>
    class ScratchDataDG
    {
    public:
      // Constructor
      ScratchDataDG(
        const Mapping<dim>        &mapping,
        const FiniteElement<dim>  &fe,
        const Quadrature<dim>     &quadrature,
        const Quadrature<dim - 1> &face_quadrature,
        const UpdateFlags update_flags = update_values | update_gradients |
                                         update_quadrature_points |
                                         update_JxW_values,
        const UpdateFlags face_update_flags = update_values |
                                              update_quadrature_points |
                                              update_JxW_values |
                                              update_normal_vectors,
        const UpdateFlags neighbor_face_update_flags = update_values)
        : fe_values(mapping, fe, quadrature, update_flags)
        , fe_face_values(mapping, fe, face_quadrature, face_update_flags)
        , fe_face_values_neighbor(mapping,
                                  fe,
                                  face_quadrature,
                                  neighbor_face_update_flags)
      {}

      // Copy constructor
      ScratchDataDG(const ScratchDataDG &scratch_data)
        : fe_values(scratch_data.fe_values.get_mapping(),
                    scratch_data.fe_values.get_fe(),
                    scratch_data.fe_values.get_quadrature(),
                    scratch_data.fe_values.get_update_flags())
        , fe_face_values(scratch_data.fe_face_values.get_mapping(),
                         scratch_data.fe_face_values.get_fe(),
                         scratch_data.fe_face_values.get_quadrature(),
                         scratch_data.fe_face_values.get_update_flags())
        , fe_face_values_neighbor(
            scratch_data.fe_face_values_neighbor.get_mapping(),
            scratch_data.fe_face_values_neighbor.get_fe(),
            scratch_data.fe_face_values_neighbor.get_quadrature(),
            scratch_data.fe_face_values_neighbor.get_update_flags())
      {}

      FEValues<dim>     fe_values;
      FEFaceValues<dim> fe_face_values;
      FEFaceValues<dim> fe_face_values_neighbor;
    };

    struct CopyDataFaceDG
    {
      Vector<double> cell_vector_1;
      Vector<double> cell_vector_2;

      std::vector<types::global_dof_index> local_dof_indices;
      std::vector<types::global_dof_index> local_dof_indices_neighbor;

      template <typename Iterator>
      void
      reinit(const Iterator &cell,
             const Iterator &neighbor_cell,
             unsigned int    dofs_per_cell)
      {
        cell_vector_1.reinit(dofs_per_cell);
        cell_vector_2.reinit(dofs_per_cell);

        local_dof_indices.resize(dofs_per_cell);
        cell->get_dof_indices(local_dof_indices);

        local_dof_indices_neighbor.resize(dofs_per_cell);
        neighbor_cell->get_dof_indices(local_dof_indices_neighbor);
      }
    };

    struct CopyDataDG
    {
      Vector<double>                       cell_vector;
      std::vector<types::global_dof_index> local_dof_indices;
      std::vector<CopyDataFaceDG>          face_data;

      template <typename Iterator>
      void
      reinit(const Iterator &cell, unsigned int dofs_per_cell)
      {
        cell_vector.reinit(dofs_per_cell);

        local_dof_indices.resize(dofs_per_cell);
        cell->get_dof_indices(local_dof_indices);
      }
    };
  } // namespace Hydro
} // namespace Sapphire

template <int dim>
Sapphire::Hydro::HDSolver<dim>::HDSolver(const ParameterParser   &prm,
                                         const OutputModule<dim> &output_module,
                                         const double             beta)
  : initial_condition(prm)
  , boundary_values(prm)
  , exact_solution(prm)
  , hd_solver_control(prm)
  , flux(prm, beta)
  , output_module(output_module)
  , beta(beta)
  , mpi_communicator(MPI_COMM_WORLD)
  , mapping()
  , fe(hd_solver_control.fe_degree)
  , dof_handler(triangulation)
  , quadrature_formula(fe.tensor_degree() + 1)
  , face_quadrature_formula(fe.tensor_degree() + 1)
  , computing_timer(mpi_communicator,
                    std::cout,
                    TimerOutput::never,
                    TimerOutput::wall_times)
{
  LogStream::Prefix p("HDSolver", saplog);
  AssertDimension(dim, 1);
  saplog << "Create HDSolver" << std::endl;
}

template <int dim>
void
Sapphire::Hydro::HDSolver<dim>::make_grid()
{
  TimerOutput::Scope t(computing_timer, "Make grid");
  LogStream::Prefix  p("HDSolver", saplog);
  saplog << "Make grid" << std::endl;

  GridGenerator::hyper_cube(triangulation, -1, 1);
  triangulation.refine_global(hd_solver_control.refinement_level);
  saplog << "Number of active cells:\t" << triangulation.n_active_cells()
         << std::endl;
}

template <int dim>
void
Sapphire::Hydro::HDSolver<dim>::setup_system()
{
  TimerOutput::Scope t(computing_timer, "Setup system");
  saplog << "Setup system" << std::endl;

  dof_handler.distribute_dofs(fe);

  DynamicSparsityPattern dsp(dof_handler.n_dofs());
  DoFTools::make_flux_sparsity_pattern(dof_handler, dsp);
  sparcity_pattern.copy_from(dsp);

  mass_matrix.reinit(sparcity_pattern);
  system_matrix.reinit(sparcity_pattern);

  solution.reinit(dof_handler.n_dofs());
  current_solution.reinit(dof_handler.n_dofs());
  old_solution.reinit(dof_handler.n_dofs());
  dg_vector.reinit(dof_handler.n_dofs());
  system_rhs.reinit(dof_handler.n_dofs());

  VectorTools::interpolate(mapping, dof_handler, initial_condition, solution);

  constraints.clear();
  constraints.close();
}

template <int dim>
void
Sapphire::Hydro::HDSolver<dim>::assemble_mass_matrix()
{
  TimerOutput::Scope t(computing_timer, "Assemble mass matrix");
  saplog << "Assemble mass matrix" << std::endl;

  MatrixCreator::create_mass_matrix(mapping,
                                    dof_handler,
                                    quadrature_formula,
                                    mass_matrix);
}

template <int dim>
void
Sapphire::Hydro::HDSolver<dim>::assemble_dg_vector()
{
  TimerOutput::Scope t(computing_timer, "Assemble DG vector");
  saplog << "Assemble DG vector" << std::endl;

  dg_vector = 0;
  boundary_values.set_time(current_time);

  using Iterator = typename DoFHandler<dim>::active_cell_iterator;

  const auto cell_worker = [&](const Iterator     &cell,
                               ScratchDataDG<dim> &scratch_data,
                               CopyDataDG         &copy_data) {
    FEValues<dim> &fe_values = scratch_data.fe_values;

    fe_values.reinit(cell);
    const unsigned int n_dofs = fe_values.get_fe().n_dofs_per_cell();

    Tensor<1, dim>      flux_value;
    std::vector<double> current_solution_values(n_dofs);
    fe_values.get_function_values(current_solution, current_solution_values);

    copy_data.reinit(cell, n_dofs);

    for (const unsigned int q_index : fe_values.quadrature_point_indices())
      {
        flux.flux(current_solution_values[q_index], flux_value);

        for (const unsigned int i : fe_values.dof_indices())
          {
            copy_data.cell_vector(i) -= flux_value *
                                        fe_values.shape_grad(i, q_index) *
                                        fe_values.JxW(q_index);
          }
      }
  };

  const auto boundary_worker = [&](const Iterator     &cell,
                                   const unsigned int &face_no,
                                   ScratchDataDG<dim> &scratch_data,
                                   CopyDataDG         &copy_data) {
    scratch_data.fe_face_values.reinit(cell, face_no);
    const FEFaceValuesBase<dim> &fe_face_values = scratch_data.fe_face_values;

    const unsigned int n_dofs = fe_face_values.get_fe().n_dofs_per_cell();

    double              boundary_value;
    std::vector<double> current_solution_values(n_dofs);
    fe_face_values.get_function_values(current_solution,
                                       current_solution_values);

    for (const unsigned int q_index : fe_face_values.quadrature_point_indices())
      {
        const double f_dot_n = beta * current_solution_values[q_index] *
                               current_solution_values[q_index] *
                               fe_face_values.normal_vector(q_index)[0];

        if (f_dot_n > 0.0)
          { // outflow boundary
            for (const unsigned int i : fe_face_values.dof_indices())
              {
                copy_data.cell_vector(i) +=
                  f_dot_n * fe_face_values.shape_value(i, q_index) *
                  fe_face_values.JxW(q_index);
              }
          }
        else
          { // inflow boundary
            for (const unsigned int i : fe_face_values.dof_indices())
              {
                boundary_value = boundary_values.value(
                  fe_face_values.quadrature_point(q_index));
                const double boundary_f_dot_n =
                  beta * boundary_value * boundary_value *
                  fe_face_values.normal_vector(q_index)[0];

                copy_data.cell_vector(i) +=
                  boundary_f_dot_n * fe_face_values.shape_value(i, q_index) *
                  fe_face_values.JxW(q_index);
              }
          }
      }
  };

  const auto face_worker = [&](const Iterator     &cell,
                               const unsigned int &face_no,
                               const unsigned int &subface_no,
                               const Iterator     &neighbor_cell,
                               const unsigned int &neighbor_face_no,
                               const unsigned int &neighbor_subface_no,
                               ScratchDataDG<dim> &scratch_data,
                               CopyDataDG         &copy_data) {
    // supress unused variable warning
    (void)subface_no;
    (void)neighbor_subface_no;

    FEFaceValues<dim> &fe_face_values = scratch_data.fe_face_values;
    fe_face_values.reinit(cell, face_no);

    FEFaceValues<dim> &fe_face_values_neighbor =
      scratch_data.fe_face_values_neighbor;
    fe_face_values_neighbor.reinit(neighbor_cell, neighbor_face_no);

    copy_data.face_data.emplace_back();
    CopyDataFaceDG &copy_data_face = copy_data.face_data.back();

    const unsigned int n_dofs = fe_face_values.get_fe().n_dofs_per_cell();
    copy_data_face.reinit(cell, neighbor_cell, n_dofs);

    std::vector<double> current_solution_values_1(n_dofs);
    fe_face_values.get_function_values(current_solution,
                                       current_solution_values_1);
    std::vector<double> current_solution_values_2(n_dofs);
    fe_face_values_neighbor.get_function_values(current_solution,
                                                current_solution_values_2);

    double flux_dot_n;

    for (const unsigned int q_index : fe_face_values.quadrature_point_indices())
      {
        flux_dot_n = flux.numerical_flux(current_solution_values_1[q_index],
                                         current_solution_values_2[q_index],
                                         fe_face_values.normal_vector(q_index));

        for (const unsigned int i : fe_face_values.dof_indices())
          {
            copy_data_face.cell_vector_1(i) +=
              flux_dot_n * fe_face_values.shape_value(i, q_index) *
              fe_face_values.JxW(q_index);

            copy_data_face.cell_vector_2(i) -=
              flux_dot_n * fe_face_values_neighbor.shape_value(i, q_index) *
              fe_face_values.JxW(q_index);
          }
      }
  };

  const auto copier = [&](const CopyDataDG &c) {
    constraints.distribute_local_to_global(c.cell_vector,
                                           c.local_dof_indices,
                                           dg_vector);

    for (auto &cdf : c.face_data)
      {
        constraints.distribute_local_to_global(cdf.cell_vector_1,
                                               cdf.local_dof_indices,
                                               dg_vector);
        constraints.distribute_local_to_global(cdf.cell_vector_2,
                                               cdf.local_dof_indices_neighbor,
                                               dg_vector);
      }
  };

  ScratchDataDG<dim> scratch_data(mapping,
                                  fe,
                                  quadrature_formula,
                                  face_quadrature_formula);
  CopyDataDG         copy_data;

  const auto filtered_iterator_range =
    dof_handler.active_cell_iterators() | IteratorFilters::LocallyOwnedCell();
  MeshWorker::mesh_loop(filtered_iterator_range,
                        cell_worker,
                        copier,
                        scratch_data,
                        copy_data,
                        MeshWorker::assemble_own_cells |
                          MeshWorker::assemble_boundary_faces |
                          MeshWorker::assemble_own_interior_faces_once,
                        boundary_worker,
                        face_worker);
}

template <int dim>
void
Sapphire::Hydro::HDSolver<dim>::assemble_system()
{
  TimerOutput::Scope t(computing_timer, "Assemble system");
  saplog << "Assemble system" << std::endl;

  system_matrix = 0;
  system_rhs    = 0;

  /** Nothing to do here, RHS of equation is zero */
}

template <int dim>
void
Sapphire::Hydro::HDSolver<dim>::perform_time_step()
{
  TimerOutput::Scope t(computing_timer, "Time step");
  saplog << "Time step" << std::endl;
  LogStream::Prefix p("TimeStep", saplog);

  old_solution     = solution;
  current_solution = old_solution;
  Vector<double> tmp(dof_handler.n_dofs());

  SolverControl solver_control(hd_solver_control.max_iterations,
                               hd_solver_control.tolerance);
  SolverRichardson<Vector<double>> solver(solver_control);

  PreconditionSSOR<SparseMatrix<double>> preconditioner;

  double time_step = hd_solver_control.time_step;

  switch (hd_solver_control.scheme)
    {
      case TimeSteppingScheme::ForwardEuler:
        {
          assemble_system();
          assemble_dg_vector();

          system_rhs.add(-1.0, dg_vector);
          system_rhs *= time_step;

          mass_matrix.vmult(tmp, current_solution);
          system_rhs.add(1.0, tmp);

          system_matrix.copy_from(mass_matrix);
          preconditioner.initialize(system_matrix);

          {
            TimerOutput::Scope t(computing_timer, "Solve linear system");
            saplog << "Solve linear system" << std::endl;
            LogStream::Prefix p("Solve", saplog);
            solver.solve(system_matrix, solution, system_rhs, preconditioner);
            saplog << "Solver converged in " << solver_control.last_step()
                   << " iterations." << std::endl;
          }

          break;
        }

      case TimeSteppingScheme::ExplicitRK:
        {
          //  Butcher's array
          Vector<double> a({0.5, 0.5, 1.});
          Vector<double> b({1. / 6, 1. / 3, 1. / 3, 1. / 6});
          Vector<double> c({0., 0.5, 0.5, 1.});
          int            i = 0;

          Vector<double> k1(dof_handler.n_dofs());
          i                = 0;
          tmp              = 0;
          current_time     = time + c[i] * time_step;
          current_solution = old_solution;
          assemble_system();
          assemble_dg_vector();
          system_matrix.copy_from(mass_matrix);
          system_rhs.add(-1.0, dg_vector);
          {
            TimerOutput::Scope t(computing_timer, "Solve linear system");
            saplog << "Solve linear system" << std::endl;
            LogStream::Prefix p("Solve", saplog);
            preconditioner.initialize(system_matrix);
            solver.solve(system_matrix, solution, system_rhs, preconditioner);
            saplog << "Solver converged in " << solver_control.last_step()
                   << " iterations." << std::endl;
          }
          k1 = solution;

          Vector<double> k2(dof_handler.n_dofs());
          i                = 1;
          tmp              = 0;
          current_time     = time + c[i] * time_step;
          current_solution = old_solution;
          current_solution.add(a[i - 1] * time_step, k1);
          assemble_system();
          assemble_dg_vector();
          system_matrix.copy_from(mass_matrix);
          system_rhs.add(-1.0, dg_vector);
          {
            TimerOutput::Scope t(computing_timer, "Solve linear system");
            saplog << "Solve linear system" << std::endl;
            LogStream::Prefix p("Solve", saplog);
            solver.solve(system_matrix, solution, system_rhs, preconditioner);
            saplog << "Solver converged in " << solver_control.last_step()
                   << " iterations." << std::endl;
          }
          k2 = solution;

          Vector<double> k3(dof_handler.n_dofs());
          i                = 2;
          tmp              = 0;
          current_time     = time + c[i] * time_step;
          current_solution = old_solution;
          current_solution.add(a[i - 1] * time_step, k2);
          assemble_system();
          assemble_dg_vector();
          system_matrix.copy_from(mass_matrix);
          system_rhs.add(-1.0, dg_vector);
          {
            TimerOutput::Scope t(computing_timer, "Solve linear system");
            saplog << "Solve linear system" << std::endl;
            LogStream::Prefix p("Solve", saplog);
            solver.solve(system_matrix, solution, system_rhs, preconditioner);
            saplog << "Solver converged in " << solver_control.last_step()
                   << " iterations." << std::endl;
          }
          k3 = solution;

          Vector<double> k4(dof_handler.n_dofs());
          i                = 3;
          tmp              = 0;
          current_time     = time + c[i] * time_step;
          current_solution = old_solution;
          current_solution.add(a[i - 1] * time_step, k3);
          assemble_system();
          assemble_dg_vector();
          system_matrix.copy_from(mass_matrix);
          system_rhs.add(-1.0, dg_vector);
          {
            TimerOutput::Scope t(computing_timer, "Solve linear system");
            saplog << "Solve linear system" << std::endl;
            LogStream::Prefix p("Solve", saplog);
            solver.solve(system_matrix, solution, system_rhs, preconditioner);
            saplog << "Solver converged in " << solver_control.last_step()
                   << " iterations." << std::endl;
          }
          k4 = solution;

          solution = old_solution;
          solution.add(b[0] * time_step, k1);
          solution.add(b[1] * time_step, k2);
          solution.add(b[2] * time_step, k3);
          solution.add(b[3] * time_step, k4);
          break;
        }

      default:
        Assert(false, ExcNotImplemented());
        break;
    }

  time += time_step;
}

template <int dim>
void
Sapphire::Hydro::HDSolver<dim>::output_results()
{
  saplog << "Output results" << std::endl;
  // LogStream::Prefix p("OutputResults", saplog);

  Vector<double> exact_solution_values(dof_handler.n_dofs());
  exact_solution.set_time(time);
  VectorTools::interpolate(mapping,
                           dof_handler,
                           exact_solution,
                           exact_solution_values);

  DataOut<dim> data_out;

  data_out.attach_dof_handler(dof_handler);
  data_out.add_data_vector(solution, "Solution");
  data_out.add_data_vector(exact_solution_values, "ExactSolution");

  data_out.build_patches(fe.tensor_degree());

  output_module.write_results(data_out, timestep_number);
}

template <int dim>
void
Sapphire::Hydro::HDSolver<dim>::process_results()
{
  TimerOutput::Scope t(computing_timer, "Process results");
  saplog << "Process results" << std::endl;
  LogStream::Prefix p("ProcessResults", saplog);

  Vector<float> difference_per_cell(triangulation.n_active_cells());
  exact_solution.set_time(time);

  // Use different quadrature for error computation
  const QTrapezoid<1>  q_trapez;
  const QIterated<dim> q_iterated(q_trapez, fe.degree * 2 + 1);

  VectorTools::integrate_difference(mapping,
                                    dof_handler,
                                    solution,
                                    exact_solution,
                                    difference_per_cell,
                                    q_iterated,
                                    VectorTools::L2_norm);
  float L2_error = VectorTools::compute_global_error(triangulation,
                                                     difference_per_cell,
                                                     VectorTools::L2_norm);
  saplog << "L2 error:\t\t" << L2_error << std::endl;

  VectorTools::integrate_difference(mapping,
                                    dof_handler,
                                    solution,
                                    exact_solution,
                                    difference_per_cell,
                                    q_iterated,
                                    VectorTools::Linfty_norm);
  float Linf_error =
    VectorTools::compute_global_error(triangulation,
                                      difference_per_cell,
                                      VectorTools::Linfty_norm);
  saplog << "L-infinity error:\t" << Linf_error << std::endl;

  error_with_time.push_back(L2_error);
}

template <int dim>
void
Sapphire::Hydro::HDSolver<dim>::init()
{
  LogStream::Prefix p("HDSolver", saplog);
  saplog << "Init HDSolver" << std::endl;
  LogStream::Prefix p2("Init", saplog);
  time            = 0.0;
  timestep_number = 0;
  error_with_time.clear();
  error_with_time.reserve(
    (unsigned int)(hd_solver_control.end_time / hd_solver_control.time_step) /
    output_module.output_frequency);

  make_grid();
  setup_system();
  assemble_mass_matrix();
  // assemble_system();
  {
    TimerOutput::Scope t(computing_timer, "Output results");
    output_results();
  }
}

template <int dim>
void
Sapphire::Hydro::HDSolver<dim>::do_timestep()
{
  LogStream::Prefix p("HDSolver", saplog);
  saplog << "Timestep " << timestep_number + 1
         << " (time = " << time + hd_solver_control.time_step << "/"
         << hd_solver_control.end_time << ")" << std::endl;
  LogStream::Prefix p2("TimeStep", saplog);

  perform_time_step();
  timestep_number++;
  if (timestep_number % output_module.output_frequency == 0)
    {
      {
        TimerOutput::Scope t(computing_timer, "Output results");
        output_results();
      }
      process_results();
    }
}

template <int dim>
void
Sapphire::Hydro::HDSolver<dim>::run()
{
  {
    LogStream::Prefix p("HDSolver", saplog);
    saplog << "Run HDSolver" << std::endl;
  }
  init();
  while (time < hd_solver_control.end_time)
    {
      do_timestep();
    }

  computing_timer.print_summary();
  computing_timer.reset();

  {
    LogStream::Prefix p("HDSolver", saplog);
    LogStream::Prefix p2("Run", saplog);

    saplog << "L2 error with time:" << std::endl;
    char buffer[100];
    for (const auto &e : error_with_time)
      {
        snprintf(buffer, 100, "%.2e ", e);
        saplog << buffer;
      }
    saplog << std::endl;
  }
}

// explicit instantiation
template class Sapphire::Hydro::HDSolver<1>;
