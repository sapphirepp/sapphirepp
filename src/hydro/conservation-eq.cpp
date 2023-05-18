#include "conservation-eq.h"

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
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/vector.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/vector_tools.h>

#include <fstream>
#include <iostream>

void Sapphire::Hydro::ExactSolution::vector_value(
    const Point<1> &p, Vector<double> &values) const {
  AssertDimension(values.size(), 1);
  values(0) = std::sin(numbers::PI * (p(0) - a * this->get_time()));
}

void Sapphire::Hydro::InitialCondition::vector_value(
    const Point<1> &p, Vector<double> &values) const {
  ExactSolution(a, 0.0).vector_value(p, values);
}

Sapphire::Hydro::ConservationEq::ConservationEq()
    : mapping(), fe(1), dof_handler(triangulation),
      quadrature_formula(fe.tensor_degree() + 1),
      face_quadrature_formula(fe.tensor_degree() + 1) {
  std::cout << "Setup conservation equation" << std::endl;
  time = 0.0;
  time_step = 0.1;
  timestep_number = 0;
}

void Sapphire::Hydro::ConservationEq::make_grid() {
  std::cout << "Make grid" << std::endl;

  GridGenerator::hyper_cube(triangulation, -1, 1);
  triangulation.refine_global(5);
  std::cout << "  Number of active cells:       "
            << triangulation.n_active_cells() << std::endl;
}

void Sapphire::Hydro::ConservationEq::setup_system() {
  std::cout << "Setup system" << std::endl;

  dof_handler.distribute_dofs(fe);

  DynamicSparsityPattern dsp(dof_handler.n_dofs());
  DoFTools::make_flux_sparsity_pattern(dof_handler, dsp);
  sparcity_pattern.copy_from(dsp);

  mass_matrix.reinit(sparcity_pattern);
  system_matrix2.reinit(sparcity_pattern);
  solution.reinit(dof_handler.n_dofs());
  old_solution.reinit(dof_handler.n_dofs());
  system_rhs.reinit(dof_handler.n_dofs());

  VectorTools::interpolate(dof_handler, InitialCondition(a), solution);

  constraints.clear();
  constraints.close();
}

void Sapphire::Hydro::ConservationEq::assemble_system() {
  std::cout << "Assemble system" << std::endl;

  FEValues<dim> fe_values(mapping, fe, quadrature_formula,
                          update_values | update_gradients |
                              update_quadrature_points | update_JxW_values);

  FullMatrix<double> cell_mass_matrix(fe.dofs_per_cell, fe.dofs_per_cell);
  FullMatrix<double> cell_system_matrix(fe.dofs_per_cell, fe.dofs_per_cell);
  Vector<double> cell_rhs(fe.dofs_per_cell);

  std::vector<types::global_dof_index> local_dof_indices(fe.dofs_per_cell);

  for (const auto &cell : dof_handler.active_cell_iterators()) {
    cell_mass_matrix = 0;
    cell_system_matrix = 0;
    cell_rhs = 0;

    fe_values.reinit(cell);

    for (const unsigned int q_index : fe_values.quadrature_point_indices()) {
      for (const unsigned int i : fe_values.dof_indices()) {
        for (const unsigned int j : fe_values.dof_indices()) {
          cell_mass_matrix(i, j) +=
              (fe_values.shape_value(i, q_index) *
               fe_values.shape_value(j, q_index) * fe_values.JxW(q_index));
          cell_system_matrix(i, j) +=
              a * (fe_values.shape_value(i, q_index) *
                   fe_values.shape_value(j, q_index) * fe_values.JxW(q_index));
        }

        cell_rhs(i) += 0.0; // TODO: missing the Lax-Friedrichs flux
      }
    }

    cell->get_dof_indices(local_dof_indices);
    constraints.distribute_local_to_global(cell_mass_matrix, local_dof_indices,
                                           mass_matrix);
    constraints.distribute_local_to_global(cell_system_matrix,
                                           local_dof_indices, system_matrix2);
    constraints.distribute_local_to_global(cell_rhs, local_dof_indices,
                                           system_rhs);
  }

  Vector<double> tmp(dof_handler.n_dofs());

  system_matrix2.vmult(tmp, solution);
  system_rhs.add(1.0, tmp);
  system_rhs *= time_step;
  mass_matrix.vmult(tmp, old_solution);
  system_rhs.add(1.0, tmp);

  // system_rhs = tmp; // constant solution
}

void Sapphire::Hydro::ConservationEq::solve() {
  std::cout << "Solve" << std::endl;

  SolverControl solver_control(1000, 1e-12);
  SolverCG<Vector<double>> solver(solver_control);

  // PreconditionSSOR<SparseMatrix<double>> preconditioner;
  // preconditioner.initialize(system_matrix, 1.2);

  PreconditionIdentity preconditioner;

  solver.solve(mass_matrix, solution, system_rhs, preconditioner);

  // constraints.distribute(solution);

  std::cout << "   " << solver_control.last_step() << " CG iterations."
            << std::endl;
}

void Sapphire::Hydro::ConservationEq::output_results() const {
  std::cout << "Output results" << std::endl;

  Vector<double> exact_solution(dof_handler.n_dofs());
  VectorTools::interpolate(dof_handler, ExactSolution(a, time + time_step),
                           exact_solution);
  // VectorTools::interpolate(dof_handler, ExactSolution(a, time),
  // exact_solution);

  DataOut<dim> data_out;

  data_out.attach_dof_handler(dof_handler);
  data_out.add_data_vector(solution, "Solution");
  data_out.add_data_vector(exact_solution, "ExactSolution");

  data_out.build_patches();

  const std::string filename = "../results/solution-" +
                               Utilities::int_to_string(timestep_number, 3) +
                               ".vtu";
  DataOutBase::VtkFlags vtk_flags;
  vtk_flags.compression_level =
      DataOutBase::VtkFlags::ZlibCompressionLevel::best_speed;
  data_out.set_flags(vtk_flags);
  std::ofstream output(filename);
  data_out.write_vtu(output);
}

void Sapphire::Hydro::ConservationEq::run() {
  std::cout << "Run conservation equation" << std::endl;
  make_grid();
  setup_system();
  output_results();
  timestep_number++;

  for (; time < 1; time += time_step, ++timestep_number) {
    old_solution = solution;
    assemble_system();
    std::cout << "  Time: " << time << std::endl;
    solve();
    output_results();
  }
}
