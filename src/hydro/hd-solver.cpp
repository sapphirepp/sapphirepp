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

template <int dim>
Sapphire::Hydro::HDSolver<dim>::HDSolver(const ParameterParser   &prm,
                                         const OutputModule<dim> &output_module)
  : initial_condition(prm)
  , exact_solution(prm)
  , hd_solver_control(prm)
  , output_module(output_module)
  , pcout(std::cout, Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
#if dim > 1
  , triangulation(MPI_COMM_WORLD)
#endif
  , fe(FE_DGQ<dim>(fe_degree), dim + 2)
  , mapping(fe_degree)
  , dof_handler(triangulation)
  , timer(pcout, TimerOutput::never, TimerOutput::wall_times)
  , euler_operator(timer)
  , time(0)
  , time_step(0)
{
  LogStream::Prefix p("HDSolver", saplog);
  saplog << "Create HDSolver" << std::endl;
}


template <int dim>
void
Sapphire::Hydro::HDSolver<dim>::make_grid_and_dofs()
{
  // TODO_HD: read grif from file
  switch (testcase)
    {
      case 0:
        {
          Point<dim> lower_left;
          for (unsigned int d = 1; d < dim; ++d)
            lower_left[d] = -5;

          Point<dim> upper_right;
          upper_right[0] = 10;
          for (unsigned int d = 1; d < dim; ++d)
            upper_right[d] = 5;

          GridGenerator::hyper_rectangle(triangulation,
                                         lower_left,
                                         upper_right);
          triangulation.refine_global(2);

          euler_operator.set_inflow_boundary(
            0, std::make_unique<HDExactSolution<dim>>(0));

          break;
        }

      case 1:
        {
          GridGenerator::channel_with_cylinder(triangulation, 0.03, 1, 0, true);

          euler_operator.set_inflow_boundary(
            0, std::make_unique<HDExactSolution<dim>>(0));
          euler_operator.set_subsonic_outflow_boundary(
            1, std::make_unique<HDExactSolution<dim>>(0));

          euler_operator.set_wall_boundary(2);
          euler_operator.set_wall_boundary(3);

          if (dim == 3)
            euler_operator.set_body_force(
              std::make_unique<Functions::ConstantFunction<dim>>(
                std::vector<double>({0., 0., -0.2})));

          break;
        }

      case 2:
        {
          GridGenerator::hyper_cube(triangulation, -1, 1);
          triangulation.refine_global(n_global_refinements);

          euler_operator.set_inflow_boundary(
            0, std::make_unique<HDExactSolution<dim>>(0));
          euler_operator.set_inflow_boundary(
            1, std::make_unique<HDExactSolution<dim>>(0));

          break;
        }

      default:
        Assert(false, ExcNotImplemented());
    }

  triangulation.refine_global(n_global_refinements);

  dof_handler.distribute_dofs(fe);

  euler_operator.reinit(mapping, dof_handler);
  euler_operator.initialize_vector(solution);

  std::locale s = pcout.get_stream().getloc();
  pcout.get_stream().imbue(std::locale(""));
  pcout << "Number of degrees of freedom: " << dof_handler.n_dofs()
        << " ( = " << (dim + 2) << " [vars] x "
        << triangulation.n_global_active_cells() << " [cells] x "
        << Utilities::pow(fe_degree + 1, dim) << " [dofs/cell/var] )"
        << std::endl;
  pcout.get_stream().imbue(s);
}

template <int dim>
void
Sapphire::Hydro::HDSolver<dim>::output_results(const unsigned int result_number)
{
  // TODO_HD: Use OutputModule
  const std::array<double, 3> errors =
    euler_operator.compute_errors(HDExactSolution<dim>(time), solution);
  const std::string quantity_name =
    (testcase == 0 or testcase == 2) ? "error" : "norm";

  pcout << "Time:" << std::setw(8) << std::setprecision(3) << time
        << ", dt: " << std::setw(8) << std::setprecision(2) << time_step << ", "
        << quantity_name << " rho: " << std::setprecision(4) << std::setw(10)
        << errors[0] << ", rho * u: " << std::setprecision(4) << std::setw(10)
        << errors[1] << ", energy:" << std::setprecision(4) << std::setw(10)
        << errors[2] << std::endl;

  {
    TimerOutput::Scope t(timer, "output");

    Postprocessor<dim> postprocessor;
    DataOut<dim>       data_out;

    DataOutBase::VtkFlags flags;
    if (dim > 1)
      flags.write_higher_order_cells = true;
    data_out.set_flags(flags);

    data_out.attach_dof_handler(dof_handler);
    {
      std::vector<std::string> names;
      names.emplace_back("density");
      for (unsigned int d = 0; d < dim; ++d)
        names.emplace_back("momentum");
      names.emplace_back("energy");

      std::vector<DataComponentInterpretation::DataComponentInterpretation>
        interpretation;
      interpretation.push_back(
        DataComponentInterpretation::component_is_scalar);
      for (unsigned int d = 0; d < dim; ++d)
        interpretation.push_back(
          DataComponentInterpretation::component_is_part_of_vector);
      interpretation.push_back(
        DataComponentInterpretation::component_is_scalar);

      data_out.add_data_vector(dof_handler, solution, names, interpretation);
    }
    data_out.add_data_vector(solution, postprocessor);

    LinearAlgebra::distributed::Vector<Number> reference;
    if ((testcase == 0 && dim == 2) or (testcase == 2))
      {
        reference.reinit(solution);
        euler_operator.project(HDExactSolution<dim>(time), reference);
        reference.sadd(-1., 1, solution);
        std::vector<std::string> names;
        names.emplace_back("error_density");
        for (unsigned int d = 0; d < dim; ++d)
          names.emplace_back("error_momentum");
        names.emplace_back("error_energy");

        std::vector<DataComponentInterpretation::DataComponentInterpretation>
          interpretation;
        interpretation.push_back(
          DataComponentInterpretation::component_is_scalar);
        for (unsigned int d = 0; d < dim; ++d)
          interpretation.push_back(
            DataComponentInterpretation::component_is_part_of_vector);
        interpretation.push_back(
          DataComponentInterpretation::component_is_scalar);

        data_out.add_data_vector(dof_handler, reference, names, interpretation);
      }

    Vector<double> mpi_owner(triangulation.n_active_cells());
    mpi_owner = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);
    data_out.add_data_vector(mpi_owner, "owner");

    data_out.build_patches(mapping,
                           fe.degree,
                           DataOut<dim>::curved_inner_cells);

    const std::string filename = "../results/solution_" +
                                 Utilities::int_to_string(result_number, 3) +
                                 ".vtu";
    data_out.write_vtu_in_parallel(filename, MPI_COMM_WORLD);
  }
}

template <int dim>
void
Sapphire::Hydro::HDSolver<dim>::run()
{
  // TODO_HD: split into individual time steps
  {
    const unsigned int n_vect_number = VectorizedArray<Number>::size();
    const unsigned int n_vect_bits   = 8 * sizeof(Number) * n_vect_number;

    pcout << "Running with " << Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD)
          << " MPI processes" << std::endl;
    pcout << "Vectorization over " << n_vect_number << ' '
          << (std::is_same<Number, double>::value ? "doubles" : "floats")
          << " = " << n_vect_bits << " bits ("
          << Utilities::System::get_current_vectorization_level() << ')'
          << std::endl;
  }

  make_grid_and_dofs();

  const LowStorageRungeKuttaIntegrator integrator(lsrk_scheme);

  LinearAlgebra::distributed::Vector<Number> rk_register_1;
  LinearAlgebra::distributed::Vector<Number> rk_register_2;
  rk_register_1.reinit(solution);
  rk_register_2.reinit(solution);

  euler_operator.project(HDExactSolution<dim>(time), solution);

  double min_vertex_distance = std::numeric_limits<double>::max();
  for (const auto &cell : triangulation.active_cell_iterators())
    if (cell->is_locally_owned())
      min_vertex_distance =
        std::min(min_vertex_distance, cell->minimum_vertex_distance());
  min_vertex_distance =
    Utilities::MPI::min(min_vertex_distance, MPI_COMM_WORLD);

  time_step = courant_number * integrator.n_stages() /
              euler_operator.compute_cell_transport_speed(solution);
  pcout << "Time step size: " << time_step
        << ", minimal h: " << min_vertex_distance
        << ", initial transport scaling: "
        << 1. / euler_operator.compute_cell_transport_speed(solution)
        << std::endl
        << std::endl;

  output_results(0);

  unsigned int timestep_number = 0;

  while (time < final_time - 1e-12)
    {
      ++timestep_number;
      if (timestep_number % 5 == 0)
        time_step = courant_number * integrator.n_stages() /
                    Utilities::truncate_to_n_digits(
                      euler_operator.compute_cell_transport_speed(solution), 3);

      {
        TimerOutput::Scope t(timer, "rk time stepping total");
        integrator.perform_time_step(euler_operator,
                                     time,
                                     time_step,
                                     solution,
                                     rk_register_1,
                                     rk_register_2);
      }

      time += time_step;

      if (static_cast<int>(time / output_tick) !=
            static_cast<int>((time - time_step) / output_tick) ||
          time >= final_time - 1e-12)
        output_results(
          static_cast<unsigned int>(std::round(time / output_tick)));
    }

  timer.print_wall_time_statistics(MPI_COMM_WORLD);
  pcout << std::endl;
}

// explicit instantiation
template class Sapphire::Hydro::HDSolver<Sapphire::Hydro::dimension>;
// TODO_HD: add all instantiations
//  template class Sapphire::Hydro::HDSolver<1>;
//  template class Sapphire::Hydro::HDSolver<2>;
//  template class Sapphire::Hydro::HDSolver<3>;
