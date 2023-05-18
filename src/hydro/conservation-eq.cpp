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
#include <deal.II/lac/precondition_block.h>
#include <deal.II/lac/solver_richardson.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/vector.h>

#include <deal.II/meshworker/mesh_loop.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/derivative_approximation.h>
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
  ExactSolution(0.0).vector_value(p, values);
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

  system_matrix.reinit(sparcity_pattern);
  solution.reinit(dof_handler.n_dofs());
  old_solution.reinit(dof_handler.n_dofs());
  system_rhs.reinit(dof_handler.n_dofs());
}

void Sapphire::Hydro::ConservationEq::assemble_system() {
  std::cout << "Assemble system" << std::endl;
}

void Sapphire::Hydro::ConservationEq::solve() {
  std::cout << "Solve" << std::endl;
}

void Sapphire::Hydro::ConservationEq::output_results() const {
  std::cout << "Output results" << std::endl;

  Vector<double> exact_solution(dof_handler.n_dofs());
  VectorTools::interpolate(dof_handler, ExactSolution(time), exact_solution);

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
  assemble_system();

  for (; time < 1; time += time_step, ++timestep_number) {
    std::cout << "  Time: " << time << std::endl;
    solve();
    output_results();
  }
}
