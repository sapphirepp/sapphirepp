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
#include <deal.II/lac/precondition_block.h>
#include <deal.II/lac/solver_richardson.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/vector.h>

#include <deal.II/meshworker/mesh_loop.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/vector_tools.h>

#include <fstream>
#include <iostream>

double Sapphire::Hydro::minmod(const std::vector<double> &values) {
  auto [min_it, max_it] = std::minmax_element(values.begin(), values.end());
  if ((*min_it) * (*max_it) < 0.0)
    return 0.0;
  else if (std::abs(*min_it) < std::abs(*max_it))
    return *min_it;
  else
    return *max_it;
}

template <int dim>
void Sapphire::Hydro::minmod(const std::vector<Tensor<1, dim>> &values,
                             const unsigned int n,
                             Tensor<1, dim> &return_value) {
  std::vector<double> component_values(n);
  for (unsigned int d = 0; d < dim; ++d) {
    for (unsigned int i = 0; i < n; ++i)
      component_values[i] = values[i][d];
    return_value[d] = minmod(component_values);
  }
}

// explicit instantiation
template void
Sapphire::Hydro::minmod<1>(const std::vector<Tensor<1, 1>> &values,
                           const unsigned int n, Tensor<1, 1> &return_value);
template void
Sapphire::Hydro::minmod<2>(const std::vector<Tensor<1, 2>> &values,
                           const unsigned int n, Tensor<1, 2> &return_value);
template void
Sapphire::Hydro::minmod<3>(const std::vector<Tensor<1, 3>> &values,
                           const unsigned int n, Tensor<1, 3> &return_value);

template <int dim>
Sapphire::Hydro::ConservationEq<dim>::ConservationEq(
    TensorFunction<1, dim, double> *beta, Function<dim> *initial_condition,
    Function<dim> *boundary_values, Function<dim> *exact_solution)
    : beta(beta), initial_condition(initial_condition),
      boundary_values(boundary_values), exact_solution(exact_solution),
      mpi_communicator(MPI_COMM_WORLD), mapping(), fe(1),
      dof_handler(triangulation), quadrature_formula(fe.tensor_degree() + 1),
      face_quadrature_formula(fe.tensor_degree() + 1), error_with_time(),
      pcout(std::cout,
            (Utilities::MPI::this_mpi_process(mpi_communicator) == 0)),
      computing_timer(mpi_communicator, pcout, TimerOutput::never,
                      TimerOutput::wall_times) {
  pcout << "Setup conservation equation" << std::endl;

  // Point<dim> x({1, 1});
  Point<dim> x;
  Tensor<1, dim> beta_value = this->beta->value(x);
  pcout << "  beta(x=[" << x << "]) = [" << beta_value << "]" << std::endl;

  time = 0.0;
  // time_step = 0.001;
  // time_step = 0.01;
  time_step = 0.05;
  timestep_number = 0;
}

template <int dim> void Sapphire::Hydro::ConservationEq<dim>::make_grid() {
  TimerOutput::Scope t(computing_timer, "Make grid");
  pcout << "Make grid" << std::endl;

  GridGenerator::hyper_cube(triangulation, -1, 1);
  triangulation.refine_global(7);
  // triangulation.refine_global(5);
  pcout << "  Number of active cells:       " << triangulation.n_active_cells()
        << std::endl;
}

template <int dim> void Sapphire::Hydro::ConservationEq<dim>::setup_system() {
  TimerOutput::Scope t(computing_timer, "Setup system");
  pcout << "Setup system" << std::endl;

  dof_handler.distribute_dofs(fe);

  DynamicSparsityPattern dsp(dof_handler.n_dofs());
  DoFTools::make_flux_sparsity_pattern(dof_handler, dsp);
  sparcity_pattern.copy_from(dsp);

  mass_matrix.reinit(sparcity_pattern);
  dg_matrix.reinit(sparcity_pattern);
  system_matrix.reinit(sparcity_pattern);
  solution.reinit(dof_handler.n_dofs());
  old_solution.reinit(dof_handler.n_dofs());
  system_rhs.reinit(dof_handler.n_dofs());

  VectorTools::interpolate(dof_handler, *initial_condition, solution);

  constraints.clear();
  constraints.close();
}

template <int dim>
void Sapphire::Hydro::ConservationEq<dim>::assemble_system() {
  TimerOutput::Scope t(computing_timer, "Assemble system");
  pcout << "Assemble system" << std::endl;

  mass_matrix = 0;
  dg_matrix = 0;
  system_matrix = 0;
  system_rhs = 0;

  beta->set_time(time);
  boundary_values->set_time(time);

  using Iterator = typename DoFHandler<dim>::active_cell_iterator;

  const auto cell_worker = [&](const Iterator &cell,
                               ScratchData<dim> &scratch_data,
                               CopyData &copy_data) {
    FEValues<dim> &fe_values = scratch_data.fe_values;

    fe_values.reinit(cell);
    const unsigned int n_dofs = fe_values.get_fe().n_dofs_per_cell();
    Tensor<1, dim> beta_value;

    copy_data.reinit(cell, n_dofs);

    for (const unsigned int q_index : fe_values.quadrature_point_indices()) {
      beta_value = beta->value(fe_values.quadrature_point(q_index));
      for (const unsigned int i : fe_values.dof_indices()) {
        for (const unsigned int j : fe_values.dof_indices()) {
          copy_data.cell_mass_matrix(i, j) +=
              (fe_values.shape_value(i, q_index) *
               fe_values.shape_value(j, q_index) * fe_values.JxW(q_index));
          copy_data.cell_dg_matrix(i, j) -=
              beta_value *
              (fe_values.shape_grad(i, q_index) *
               fe_values.shape_value(j, q_index) * fe_values.JxW(q_index));
        }
      }
    }
  };

  const auto boundary_worker = [&](const Iterator &cell,
                                   const unsigned int &face_no,
                                   ScratchData<dim> &scratch_data,
                                   CopyData &copy_data) {
    scratch_data.fe_face_values.reinit(cell, face_no);
    const FEFaceValuesBase<dim> &fe_face_values = scratch_data.fe_face_values;

    // const unsigned int n_dofs = fe_face_values.get_fe().n_dofs_per_cell();

    // BoundaryValues<dim> boundary_values(beta, time);
    // boundary_values->set_time(time);
    // std::vector<double> boundary_values_vector(
    //     fe_face_values.get_quadrature_points().size());
    // boundary_values.value_list(fe_face_values.get_quadrature_points(),
    //                            boundary_values_vector);
    Vector<double> boundary_value(1);

    for (const unsigned int q_index :
         fe_face_values.quadrature_point_indices()) {
      const double v_dot_n =
          beta->value(fe_face_values.quadrature_point(q_index)) *
          fe_face_values.normal_vector(q_index);

      if (v_dot_n > 0.0) { // outflow boundary
        for (const unsigned int i : fe_face_values.dof_indices()) {
          for (const unsigned int j : fe_face_values.dof_indices()) {

            copy_data.cell_dg_matrix(i, j) +=
                v_dot_n * fe_face_values.shape_value(i, q_index) *
                fe_face_values.shape_value(j, q_index) *
                fe_face_values.JxW(q_index);
          }
        }
      } else { // inflow boundary
        for (const unsigned int i : fe_face_values.dof_indices()) {
          boundary_values->vector_value(
              fe_face_values.quadrature_point(q_index), boundary_value);

          copy_data.cell_rhs(i) += -v_dot_n * boundary_value[0] *
                                   fe_face_values.shape_value(i, q_index) *
                                   fe_face_values.JxW(q_index);
        }
      }
    }
  };

  const auto face_worker =
      [&](const Iterator &cell, const unsigned int &face_no,
          const unsigned int &subface_no, const Iterator &neighbor_cell,
          const unsigned int &neighbor_face_no,
          const unsigned int &neighbor_subface_no,
          ScratchData<dim> &scratch_data, CopyData &copy_data) {
        // supress unused variable warning
        (void)subface_no;
        (void)neighbor_subface_no;

        FEFaceValues<dim> &fe_face_values = scratch_data.fe_face_values;
        fe_face_values.reinit(cell, face_no);

        FEFaceValues<dim> &fe_face_values_neighbor =
            scratch_data.fe_face_values_neighbor;
        fe_face_values_neighbor.reinit(neighbor_cell, neighbor_face_no);

        copy_data.face_data.emplace_back();
        CopyDataFace &copy_data_face = copy_data.face_data.back();

        const unsigned int n_dofs = fe_face_values.get_fe().n_dofs_per_cell();
        copy_data_face.reinit(cell, neighbor_cell, n_dofs);

        for (const unsigned int q_index :
             fe_face_values.quadrature_point_indices()) {
          const double v_dot_n =
              beta->value(fe_face_values.quadrature_point(q_index)) *
              fe_face_values.normal_vector(q_index);

          for (const unsigned int i : fe_face_values.dof_indices()) {
            for (const unsigned int j : fe_face_values.dof_indices()) {

              copy_data_face.cell_dg_matrix_11(i, j) +=
                  0.5 * v_dot_n *
                  (fe_face_values.shape_value(i, q_index) *
                   fe_face_values.shape_value(j, q_index) *
                   fe_face_values.JxW(q_index));

              copy_data_face.cell_dg_matrix_21(i, j) -=
                  0.5 * v_dot_n *
                  (fe_face_values_neighbor.shape_value(i, q_index) *
                   fe_face_values.shape_value(j, q_index) *
                   fe_face_values.JxW(q_index));

              copy_data_face.cell_dg_matrix_12(i, j) +=
                  0.5 * v_dot_n *
                  (fe_face_values.shape_value(i, q_index) *
                   fe_face_values_neighbor.shape_value(j, q_index) *
                   fe_face_values.JxW(q_index));

              copy_data_face.cell_dg_matrix_22(i, j) -=
                  0.5 * v_dot_n *
                  (fe_face_values_neighbor.shape_value(i, q_index) *
                   fe_face_values_neighbor.shape_value(j, q_index) *
                   fe_face_values.JxW(q_index));

              // upwind flux
              const double eta = 1.0;
              // const double eta = 0.0; // central flux

              copy_data_face.cell_dg_matrix_11(i, j) +=
                  0.5 * eta * std::abs(v_dot_n) *
                  (fe_face_values.shape_value(i, q_index) *
                   fe_face_values.shape_value(j, q_index) *
                   fe_face_values.JxW(q_index));

              copy_data_face.cell_dg_matrix_21(i, j) -=
                  0.5 * eta * std::abs(v_dot_n) *
                  (fe_face_values_neighbor.shape_value(i, q_index) *
                   fe_face_values.shape_value(j, q_index) *
                   fe_face_values.JxW(q_index));

              copy_data_face.cell_dg_matrix_12(i, j) -=
                  0.5 * eta * std::abs(v_dot_n) *
                  (fe_face_values.shape_value(i, q_index) *
                   fe_face_values_neighbor.shape_value(j, q_index) *
                   fe_face_values.JxW(q_index));

              copy_data_face.cell_dg_matrix_22(i, j) +=
                  0.5 * eta * std::abs(v_dot_n) *
                  (fe_face_values_neighbor.shape_value(i, q_index) *
                   fe_face_values_neighbor.shape_value(j, q_index) *
                   fe_face_values.JxW(q_index));
            }
          }
        }
      };

  const auto copier = [&](const CopyData &c) {
    constraints.distribute_local_to_global(c.cell_mass_matrix,
                                           c.local_dof_indices, mass_matrix);
    constraints.distribute_local_to_global(c.cell_dg_matrix,
                                           c.local_dof_indices, dg_matrix);
    constraints.distribute_local_to_global(c.cell_rhs, c.local_dof_indices,
                                           system_rhs);

    for (auto &cdf : c.face_data) {
      constraints.distribute_local_to_global(cdf.cell_dg_matrix_11,
                                             cdf.local_dof_indices,
                                             cdf.local_dof_indices, dg_matrix);
      constraints.distribute_local_to_global(cdf.cell_dg_matrix_21,
                                             cdf.local_dof_indices_neighbor,
                                             cdf.local_dof_indices, dg_matrix);
      constraints.distribute_local_to_global(
          cdf.cell_dg_matrix_12, cdf.local_dof_indices,
          cdf.local_dof_indices_neighbor, dg_matrix);
      constraints.distribute_local_to_global(
          cdf.cell_dg_matrix_22, cdf.local_dof_indices_neighbor,
          cdf.local_dof_indices_neighbor, dg_matrix);

      // for (unsigned int i = 0; i < cdf.local_dof_indices.size(); ++i) {
      //   for (unsigned int j = 0; j < cdf.local_dof_indices.size(); ++j) {
      //     dg_matrix.add(cdf.local_dof_indices[i], cdf.local_dof_indices[j],
      //                   cdf.cell_dg_matrix_11(i, j));
      //     dg_matrix.add(cdf.local_dof_indices_neighbor[i],
      //                   cdf.local_dof_indices[j], cdf.cell_dg_matrix_12(i,
      //                   j));
      //     dg_matrix.add(cdf.local_dof_indices[i],
      //                   cdf.local_dof_indices_neighbor[j],
      //                   cdf.cell_dg_matrix_21(i, j));
      //     dg_matrix.add(cdf.local_dof_indices_neighbor[i],
      //                   cdf.local_dof_indices_neighbor[j],
      //                   cdf.cell_dg_matrix_22(i, j));
      //   }
      // }
    }
  };

  ScratchData<dim> scratch_data(mapping, fe, quadrature_formula,
                                face_quadrature_formula);
  CopyData copy_data;
  pcout << "Assembling System Matrices" << std::endl;

  MeshWorker::mesh_loop(dof_handler.begin_active(), dof_handler.end(),
                        cell_worker, copier, scratch_data, copy_data,
                        MeshWorker::assemble_own_cells |
                            MeshWorker::assemble_boundary_faces |
                            MeshWorker::assemble_own_interior_faces_once,
                        boundary_worker, face_worker);
  pcout << "   Done" << std::endl;
}

template <int dim>
void Sapphire::Hydro::ConservationEq<dim>::assemble_system_old() {
  TimerOutput::Scope t(computing_timer, "Assemble system");
  pcout << "Assemble system" << std::endl;

  FEValues<dim> fe_values(mapping, fe, quadrature_formula,
                          update_values | update_gradients |
                              update_quadrature_points | update_JxW_values);

  FullMatrix<double> cell_mass_matrix(fe.dofs_per_cell, fe.dofs_per_cell);
  FullMatrix<double> cell_dg_matrix(fe.dofs_per_cell, fe.dofs_per_cell);
  Vector<double> cell_rhs(fe.dofs_per_cell);

  std::vector<types::global_dof_index> local_dof_indices(fe.dofs_per_cell);
  std::vector<types::global_dof_index> local_dof_indices_neighbor(
      fe.dofs_per_cell);

  FEFaceValues<dim> fe_face_values(fe_values.get_fe(), face_quadrature_formula,
                                   update_values | update_quadrature_points |
                                       update_normal_vectors |
                                       update_JxW_values);
  FEFaceValues<dim> fe_face_values_neighbor(
      fe_values.get_fe(), face_quadrature_formula,
      update_values | update_quadrature_points | update_normal_vectors |
          update_JxW_values);
  FullMatrix<double> cell_dg_matrix_11(fe.dofs_per_cell, fe.dofs_per_cell);
  FullMatrix<double> cell_dg_matrix_12(fe.dofs_per_cell, fe.dofs_per_cell);
  FullMatrix<double> cell_dg_matrix_21(fe.dofs_per_cell, fe.dofs_per_cell);
  FullMatrix<double> cell_dg_matrix_22(fe.dofs_per_cell, fe.dofs_per_cell);

  Tensor<1, dim> beta_value;
  beta->set_time(time);
  boundary_values->set_time(time);

  mass_matrix = 0;
  dg_matrix = 0;
  system_matrix = 0;
  system_rhs = 0;

  for (const auto &cell : dof_handler.active_cell_iterators()) {
    cell_mass_matrix = 0;
    cell_dg_matrix = 0;
    cell_rhs = 0;

    fe_values.reinit(cell);

    for (const unsigned int q_index : fe_values.quadrature_point_indices()) {
      beta_value = beta->value(fe_values.quadrature_point(q_index));
      for (const unsigned int i : fe_values.dof_indices()) {
        for (const unsigned int j : fe_values.dof_indices()) {
          cell_mass_matrix(i, j) +=
              (fe_values.shape_value(i, q_index) *
               fe_values.shape_value(j, q_index) * fe_values.JxW(q_index));
          cell_dg_matrix(i, j) -=
              beta_value *
              (fe_values.shape_grad(i, q_index) *
               fe_values.shape_value(j, q_index) * fe_values.JxW(q_index));
        }
      }
    }

    cell->get_dof_indices(local_dof_indices);
    constraints.distribute_local_to_global(cell_mass_matrix, local_dof_indices,
                                           mass_matrix);
    constraints.distribute_local_to_global(cell_dg_matrix, local_dof_indices,
                                           dg_matrix);
    constraints.distribute_local_to_global(cell_rhs, local_dof_indices,
                                           system_rhs);

    // Face integration
    for (const auto &face_no : cell->face_indices()) {
      cell_dg_matrix_11 = 0;
      cell_dg_matrix_12 = 0;
      cell_dg_matrix_21 = 0;
      cell_dg_matrix_22 = 0;

      if (cell->at_boundary(face_no)) {

        fe_face_values.reinit(cell, face_no);
        for (const unsigned int q_index :
             fe_face_values.quadrature_point_indices()) {
          beta_value = beta->value(fe_face_values.quadrature_point(q_index));

          for (const unsigned int i : fe_face_values.dof_indices()) {
            for (const unsigned int j : fe_face_values.dof_indices()) {

              cell_dg_matrix_11(i, j) +=
                  0.5 *
                  (std::abs(beta_value *
                            fe_face_values.normal_vector(q_index)) +
                   beta_value * fe_face_values.normal_vector(q_index)) *
                  (fe_face_values.shape_value(i, q_index) *
                   fe_face_values.shape_value(j, q_index) *
                   fe_face_values.JxW(q_index));
            }
          }
        }

        cell->get_dof_indices(local_dof_indices);
        constraints.distribute_local_to_global(
            cell_dg_matrix_11, local_dof_indices, local_dof_indices, dg_matrix);
      } else {

        fe_face_values.reinit(cell, face_no);
        auto neighbor = cell->neighbor(face_no);
        auto neigbor_face_no = cell->neighbor_of_neighbor(face_no);
        fe_face_values_neighbor.reinit(neighbor, neigbor_face_no);
        if (neighbor < cell) {
          continue;
        }

        for (const unsigned int q_index :
             fe_face_values.quadrature_point_indices()) {
          beta_value = beta->value(fe_face_values.quadrature_point(q_index));

          for (const unsigned int i : fe_face_values.dof_indices()) {
            for (const unsigned int j : fe_face_values.dof_indices()) {

              cell_dg_matrix_11(i, j) +=
                  0.5 * beta_value * fe_face_values.normal_vector(q_index) *
                  (fe_face_values.shape_value(i, q_index) *
                   fe_face_values.shape_value(j, q_index) *
                   fe_face_values.JxW(q_index));

              cell_dg_matrix_21(i, j) -=
                  0.5 * beta_value * fe_face_values.normal_vector(q_index) *
                  (fe_face_values_neighbor.shape_value(i, q_index) *
                   fe_face_values.shape_value(j, q_index) *
                   fe_face_values.JxW(q_index));

              cell_dg_matrix_12(i, j) +=
                  0.5 * beta_value * fe_face_values.normal_vector(q_index) *
                  (fe_face_values.shape_value(i, q_index) *
                   fe_face_values_neighbor.shape_value(j, q_index) *
                   fe_face_values.JxW(q_index));

              cell_dg_matrix_22(i, j) -=
                  0.5 * beta_value * fe_face_values.normal_vector(q_index) *
                  (fe_face_values_neighbor.shape_value(i, q_index) *
                   fe_face_values_neighbor.shape_value(j, q_index) *
                   fe_face_values.JxW(q_index));
            }
          }
        }

        cell->get_dof_indices(local_dof_indices);
        neighbor->get_dof_indices(local_dof_indices_neighbor);
        constraints.distribute_local_to_global(
            cell_dg_matrix_11, local_dof_indices, local_dof_indices, dg_matrix);
        constraints.distribute_local_to_global(cell_dg_matrix_21,
                                               local_dof_indices_neighbor,
                                               local_dof_indices, dg_matrix);
        constraints.distribute_local_to_global(
            cell_dg_matrix_12, local_dof_indices, local_dof_indices_neighbor,
            dg_matrix);
        constraints.distribute_local_to_global(
            cell_dg_matrix_22, local_dof_indices_neighbor,
            local_dof_indices_neighbor, dg_matrix);
      }
    }
  }
}

template <int dim>
void Sapphire::Hydro::ConservationEq<dim>::assemble_time_step() {
  TimerOutput::Scope t(computing_timer, "Time step");
  pcout << "Time step" << std::endl;

  Vector<double> tmp(dof_handler.n_dofs());
  // Euler step

  // tmp = 0;
  // mass_matrix.vmult(tmp, old_solution);
  // system_rhs.add(1.0, tmp);

  // tmp = 0;
  // dg_matrix.vmult(tmp, old_solution);
  // system_rhs.add(-time_step, tmp);

  // system_matrix.copy_from(mass_matrix);

  // Theta method
  // const double theta = 0.0; // Forward Euler
  const double theta = 0.5; // Cranc-Nicholson
  // const double theta = 1.0; // Backward Euler

  system_rhs *= time_step;

  mass_matrix.vmult(tmp, old_solution);
  system_rhs.add(1.0, tmp);

  dg_matrix.vmult(tmp, old_solution);
  system_rhs.add(-(1. - theta) * time_step, tmp);

  system_matrix.copy_from(mass_matrix);
  system_matrix.add(theta * time_step, dg_matrix);

  /** constant solution */
  // mass_matrix.vmult(tmp, old_solution);
  // system_rhs = tmp;
}

template <int dim> void Sapphire::Hydro::ConservationEq<dim>::solve() {
  TimerOutput::Scope t(computing_timer, "Solve");
  pcout << "Solve" << std::endl;

  SolverControl solver_control(1000, 1e-12);
  // SolverCG<Vector<double>> solver(solver_control);
  SolverRichardson<Vector<double>> solver(solver_control);

  // PreconditionBlockJacobi preconditioner;
  // preconditioner.initialize(system_matrix);

  PreconditionSSOR<SparseMatrix<double>> preconditioner;
  preconditioner.initialize(system_matrix);

  solver.solve(system_matrix, solution, system_rhs, preconditioner);

  // constraints.distribute(solution);

  pcout << "   Solver converged in " << solver_control.last_step()
        << " iterations." << std::endl;
}

template <int dim>
void Sapphire::Hydro::ConservationEq<dim>::output_results() const {
  pcout << "Output results" << std::endl;

  Vector<double> exact_solution_values(dof_handler.n_dofs());
  exact_solution->set_time(time);
  VectorTools::interpolate(dof_handler, *exact_solution, exact_solution_values);
  // VectorTools::interpolate(dof_handler, ExactSolution(a, time),
  // exact_solution);

  DataOut<dim> data_out;

  data_out.attach_dof_handler(dof_handler);
  data_out.add_data_vector(solution, "Solution");
  data_out.add_data_vector(exact_solution_values, "ExactSolution");

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

template <int dim>
void Sapphire::Hydro::ConservationEq<dim>::process_results() {
  TimerOutput::Scope t(computing_timer, "Process results");
  pcout << "Process results" << std::endl;

  Vector<float> difference_per_cell(triangulation.n_active_cells());
  // ExactSolution<dim> exact_solution_function(beta, time);
  exact_solution->set_time(time);

  // Use different quadrature for error computation
  const QTrapezoid<1> q_trapez;
  const QIterated<dim> q_iterated(q_trapez, fe.degree * 2 + 1);

  VectorTools::integrate_difference(mapping, dof_handler, solution,
                                    *exact_solution, difference_per_cell,
                                    q_iterated, VectorTools::L2_norm);
  float L2_error = VectorTools::compute_global_error(
      triangulation, difference_per_cell, VectorTools::L2_norm);
  pcout << "   L2 error:\t\t" << L2_error << std::endl;

  VectorTools::integrate_difference(mapping, dof_handler, solution,
                                    *exact_solution, difference_per_cell,
                                    q_iterated, VectorTools::Linfty_norm);
  float Linf_error = VectorTools::compute_global_error(
      triangulation, difference_per_cell, VectorTools::Linfty_norm);
  pcout << "   L-infinity error:\t" << Linf_error << std::endl;

  // pcout << "   OLD QUADATURE" << std::endl;

  // VectorTools::integrate_difference(
  //     mapping, dof_handler, solution, exact_solution,
  //     difference_per_cell, quadrature_formula, VectorTools::L2_norm);
  // float L2_error_old = VectorTools::compute_global_error(
  //     triangulation, difference_per_cell, VectorTools::L2_norm);
  // pcout << "   L2 error:\t\t" << L2_error_old << std::endl;

  // VectorTools::integrate_difference(
  //     mapping, dof_handler, solution, exact_solution,
  //     difference_per_cell, quadrature_formula, VectorTools::Linfty_norm);
  // float Linf_error_old = VectorTools::compute_global_error(
  //     triangulation, difference_per_cell, VectorTools::Linfty_norm);
  // pcout << "   L-infinity error:\t" << Linf_error_old << std::endl;

  error_with_time.grow_or_shrink(error_with_time.size() + 1);
  error_with_time[error_with_time.size() - 1] = L2_error;
}

template <int dim> void Sapphire::Hydro::ConservationEq<dim>::run() {
  pcout << "Run conservation equation" << std::endl;
  make_grid();
  setup_system();
  {
    TimerOutput::Scope t(computing_timer, "Output results");
    output_results();
  }

  for (; time < 1;) {
    time += time_step, ++timestep_number;
    pcout << "  Time: " << time << std::endl;
    old_solution = solution;
    assemble_system();
    // assemble_system_old();
    assemble_time_step();
    solve();
    {
      TimerOutput::Scope t(computing_timer, "Output results");
      output_results();
    }
    process_results();
  }

  computing_timer.print_summary();
  computing_timer.reset();

  if (pcout.is_active()) {
    pcout << "   L2 error with time:" << std::endl;
    error_with_time.print(pcout.get_stream());
  }
}

// explicit instantiation
template class Sapphire::Hydro::ConservationEq<1>;
template class Sapphire::Hydro::ConservationEq<2>;
template class Sapphire::Hydro::ConservationEq<3>;

template <int dim>
Sapphire::Hydro::BurgersEq<dim>::BurgersEq(Function<dim> *initial_condition,
                                           Function<dim> *boundary_values,
                                           Function<dim> *exact_solution)
    : initial_condition(initial_condition), boundary_values(boundary_values),
      exact_solution(exact_solution), mpi_communicator(MPI_COMM_WORLD),
      mapping(), fe(1), dof_handler(triangulation),
      quadrature_formula(fe.tensor_degree() + 1),
      face_quadrature_formula(fe.tensor_degree() + 1), error_with_time(),
      time(0.0),
      //
      // time_step(0.0005),
      time_step(0.001),
      // time_step(0.002),
      // time_step(0.1),
      timestep_number(0),
      pcout(std::cout,
            (Utilities::MPI::this_mpi_process(mpi_communicator) == 0)),
      computing_timer(mpi_communicator, pcout, TimerOutput::never,
                      TimerOutput::wall_times) {
  AssertDimension(dim, 1);
  pcout << "Setup conservation equation" << std::endl;
}

template <int dim> void Sapphire::Hydro::BurgersEq<dim>::make_grid() {
  TimerOutput::Scope t(computing_timer, "Make grid");
  pcout << "Make grid" << std::endl;

  GridGenerator::hyper_cube(triangulation, -1, 1);
  // triangulation.refine_global(5);
  triangulation.refine_global(7);
  // triangulation.refine_global(9);
  pcout << "  Number of active cells:       " << triangulation.n_active_cells()
        << std::endl;
}

template <int dim> void Sapphire::Hydro::BurgersEq<dim>::setup_system() {
  TimerOutput::Scope t(computing_timer, "Setup system");
  pcout << "Setup system" << std::endl;

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

  VectorTools::interpolate(dof_handler, *initial_condition, solution);

  constraints.clear();
  constraints.close();

  mark_for_limiter.reinit(triangulation.n_active_cells());
}

template <int dim>
void Sapphire::Hydro::BurgersEq<dim>::assemble_mass_matrix() {
  TimerOutput::Scope t(computing_timer, "Assemble mass matrix");
  pcout << "Assemble mass matrix" << std::endl;

  mass_matrix = 0;

  using Iterator = typename DoFHandler<dim>::active_cell_iterator;

  const auto cell_worker = [&](const Iterator &cell,
                               ScratchData<dim> &scratch_data,
                               CopyData &copy_data) {
    FEValues<dim> &fe_values = scratch_data.fe_values;

    fe_values.reinit(cell);
    const unsigned int n_dofs = fe_values.get_fe().n_dofs_per_cell();
    copy_data.reinit(cell, n_dofs);

    for (const unsigned int q_index : fe_values.quadrature_point_indices()) {
      for (const unsigned int i : fe_values.dof_indices()) {
        for (const unsigned int j : fe_values.dof_indices()) {
          copy_data.cell_mass_matrix(i, j) +=
              (fe_values.shape_value(i, q_index) *
               fe_values.shape_value(j, q_index) * fe_values.JxW(q_index));
        }
      }
    }
  };

  const auto copier = [&](const CopyData &c) {
    constraints.distribute_local_to_global(c.cell_mass_matrix,
                                           c.local_dof_indices, mass_matrix);
  };

  ScratchData<dim> scratch_data(mapping, fe, quadrature_formula,
                                face_quadrature_formula);
  CopyData copy_data;

  MeshWorker::mesh_loop(dof_handler.begin_active(), dof_handler.end(),
                        cell_worker, copier, scratch_data, copy_data,
                        MeshWorker::assemble_own_cells);
}

template <int dim>
double Sapphire::Hydro::BurgersEq<dim>::compute_numerical_flux(
    const Tensor<1, dim> &flux_1, const Tensor<1, dim> &flux_2,
    const Tensor<1, dim> &n, const double &value_1,
    const double &value_2) const {
  double numerical_flux = 0;
  switch (flux_type) {
  case FluxType::Central: {
    numerical_flux += 0.5 * (flux_1 + flux_2) * n;
    break;
  }

  case FluxType::Upwind: {
    const double eta = 1.0;
    numerical_flux += 0.5 * (flux_1 + flux_2) * n;
    numerical_flux += 0.5 * eta * (std::abs(flux_1 * n) - std::abs(flux_2 * n));
    break;
  }

  case FluxType::LaxFriedrich: {
    const double C = 3; // TODO_BE: Calculate C
    numerical_flux += 0.5 * (flux_1 + flux_2) * n;
    numerical_flux += 0.5 * C * (value_1 - value_2);
    break;
  }

  default:
    // TODO_BE: Assert(false) or throw?
    Assert(false, ExcNotImplemented());
    break;
  }

  return numerical_flux;
}

template <int dim> void Sapphire::Hydro::BurgersEq<dim>::assemble_dg_vector() {
  TimerOutput::Scope t(computing_timer, "Assemble DG vector");
  pcout << "    Assemble DG vector" << std::endl;

  dg_vector = 0;
  boundary_values->set_time(current_time);

  using Iterator = typename DoFHandler<dim>::active_cell_iterator;

  const auto cell_worker = [&](const Iterator &cell,
                               ScratchData<dim> &scratch_data,
                               CopyData &copy_data) {
    FEValues<dim> &fe_values = scratch_data.fe_values;

    fe_values.reinit(cell);
    const unsigned int n_dofs = fe_values.get_fe().n_dofs_per_cell();

    Tensor<1, dim> flux_value;
    std::vector<double> current_solution_values(n_dofs);
    fe_values.get_function_values(current_solution, current_solution_values);

    copy_data.reinit(cell, n_dofs);

    for (const unsigned int q_index : fe_values.quadrature_point_indices()) {
      flux_value[0] = beta * current_solution_values[q_index] *
                      current_solution_values[q_index];

      for (const unsigned int i : fe_values.dof_indices()) {
        copy_data.cell_vector(i) -= flux_value *
                                    fe_values.shape_grad(i, q_index) *
                                    fe_values.JxW(q_index);
      }
    }
  };

  const auto boundary_worker = [&](const Iterator &cell,
                                   const unsigned int &face_no,
                                   ScratchData<dim> &scratch_data,
                                   CopyData &copy_data) {
    scratch_data.fe_face_values.reinit(cell, face_no);
    const FEFaceValuesBase<dim> &fe_face_values = scratch_data.fe_face_values;

    const unsigned int n_dofs = fe_face_values.get_fe().n_dofs_per_cell();

    double boundary_value;
    std::vector<double> current_solution_values(n_dofs);
    fe_face_values.get_function_values(current_solution,
                                       current_solution_values);

    for (const unsigned int q_index :
         fe_face_values.quadrature_point_indices()) {
      const double f_dot_n = beta * current_solution_values[q_index] *
                             current_solution_values[q_index] *
                             fe_face_values.normal_vector(q_index)[0];

      if (f_dot_n > 0.0) { // outflow boundary
        for (const unsigned int i : fe_face_values.dof_indices()) {

          copy_data.cell_vector(i) += f_dot_n *
                                      fe_face_values.shape_value(i, q_index) *
                                      fe_face_values.JxW(q_index);
        }
      } else { // inflow boundary
        for (const unsigned int i : fe_face_values.dof_indices()) {
          boundary_value =
              boundary_values->value(fe_face_values.quadrature_point(q_index));
          const double boundary_f_dot_n =
              beta * boundary_value * boundary_value *
              fe_face_values.normal_vector(q_index)[0];

          copy_data.cell_vector(i) += boundary_f_dot_n *
                                      fe_face_values.shape_value(i, q_index) *
                                      fe_face_values.JxW(q_index);
        }
      }
    }
  };

  const auto face_worker =
      [&](const Iterator &cell, const unsigned int &face_no,
          const unsigned int &subface_no, const Iterator &neighbor_cell,
          const unsigned int &neighbor_face_no,
          const unsigned int &neighbor_subface_no,
          ScratchData<dim> &scratch_data, CopyData &copy_data) {
        // supress unused variable warning
        (void)subface_no;
        (void)neighbor_subface_no;

        FEFaceValues<dim> &fe_face_values = scratch_data.fe_face_values;
        fe_face_values.reinit(cell, face_no);

        FEFaceValues<dim> &fe_face_values_neighbor =
            scratch_data.fe_face_values_neighbor;
        fe_face_values_neighbor.reinit(neighbor_cell, neighbor_face_no);

        copy_data.face_data.emplace_back();
        CopyDataFace &copy_data_face = copy_data.face_data.back();

        const unsigned int n_dofs = fe_face_values.get_fe().n_dofs_per_cell();
        copy_data_face.reinit(cell, neighbor_cell, n_dofs);

        std::vector<double> current_solution_values_1(n_dofs);
        fe_face_values.get_function_values(current_solution,
                                           current_solution_values_1);
        std::vector<double> current_solution_values_2(n_dofs);
        fe_face_values_neighbor.get_function_values(current_solution,
                                                    current_solution_values_2);

        Tensor<1, dim> flux_1, flux_2;
        double flux_dot_n;

        for (const unsigned int q_index :
             fe_face_values.quadrature_point_indices()) {

          flux_1[0] = beta * current_solution_values_1[q_index] *
                      current_solution_values_1[q_index];
          flux_2[0] = beta * current_solution_values_2[q_index] *
                      current_solution_values_2[q_index];

          flux_dot_n = compute_numerical_flux(
              flux_1, flux_2, fe_face_values.normal_vector(q_index),
              current_solution_values_1[q_index],
              current_solution_values_2[q_index]);

          for (const unsigned int i : fe_face_values.dof_indices()) {
            copy_data_face.cell_vector_1(i) +=
                flux_dot_n * fe_face_values.shape_value(i, q_index) *
                fe_face_values.JxW(q_index);

            copy_data_face.cell_vector_2(i) -=
                flux_dot_n * fe_face_values_neighbor.shape_value(i, q_index) *
                fe_face_values.JxW(q_index);
          }
        }
      };

  const auto copier = [&](const CopyData &c) {
    constraints.distribute_local_to_global(c.cell_vector, c.local_dof_indices,
                                           dg_vector);

    for (auto &cdf : c.face_data) {
      constraints.distribute_local_to_global(cdf.cell_vector_1,
                                             cdf.local_dof_indices, dg_vector);
      constraints.distribute_local_to_global(
          cdf.cell_vector_2, cdf.local_dof_indices_neighbor, dg_vector);
    }
  };

  ScratchData<dim> scratch_data(mapping, fe, quadrature_formula,
                                face_quadrature_formula);
  CopyData copy_data;

  MeshWorker::mesh_loop(dof_handler.begin_active(), dof_handler.end(),
                        cell_worker, copier, scratch_data, copy_data,
                        MeshWorker::assemble_own_cells |
                            MeshWorker::assemble_boundary_faces |
                            MeshWorker::assemble_own_interior_faces_once,
                        boundary_worker, face_worker);
}

template <int dim> void Sapphire::Hydro::BurgersEq<dim>::assemble_system() {
  TimerOutput::Scope t(computing_timer, "Assemble system");
  pcout << "    Assemble system" << std::endl;

  system_matrix = 0;
  system_rhs = 0;

  /** Nothing to do here, RHS of equation is zero */
}

template <int dim>
void Sapphire::Hydro::BurgersEq<dim>::compute_limited_slope(
    const double &cell_average, const Tensor<1, dim> cell_average_grad,
    const std::vector<double> &neighbor_cell_averages,
    const std::vector<Tensor<1, dim>> &neighbor_distance,
    const unsigned int n_neighbors, Tensor<1, dim> &limited_slope) const {
  switch (limiter) {
  case SlopeLimiter::NoLimiter: {
    Assert(false, ExcMessage("Slope limiter is set to NoLimiter, so this "
                             "function should not be called"));
    limited_slope = cell_average_grad;
    break;
  }

  case SlopeLimiter::CellAverage: {
    // TODO_BE: remove test case
    limited_slope = 0;
    break;
  }

  case SlopeLimiter::LinearReconstruction: {
    // TODO_BE: remove test case
    limited_slope = cell_average_grad;
    break;
  }

  case SlopeLimiter::MinMod: {
    std::vector<Tensor<1, dim>> slopes(n_neighbors + 1);
    slopes[0] = cell_average_grad;
    for (unsigned int i = 0; i < n_neighbors; ++i) {
      slopes[i + 1] = (neighbor_cell_averages[i] - cell_average) /
                      (neighbor_distance[i].norm_square() / 2.) *
                      neighbor_distance[i];
    }
    minmod(slopes, n_neighbors + 1, limited_slope);
    break;
  }

  case SlopeLimiter::GerneralizedSlopeLimiter:
    // limited_slope[0] = 100;
    break;
    // Use MUSCL limiter
    // TODO_BE: Make generalized limiting independent of limiter

  case SlopeLimiter::MUSCL: {
    std::vector<Tensor<1, dim>> slopes(n_neighbors + 1);
    slopes[0] = cell_average_grad;
    for (unsigned int i = 0; i < n_neighbors; ++i) {
      slopes[i + 1] = (neighbor_cell_averages[i] - cell_average) /
                      neighbor_distance[i].norm_square() * neighbor_distance[i];
    }
    minmod(slopes, n_neighbors + 1, limited_slope);
    break;
  }

  default:
    Assert(false, ExcNotImplemented());
    break;
  }
}

template <int dim> void Sapphire::Hydro::BurgersEq<dim>::slope_limiter() {
  if (limiter == SlopeLimiter::NoLimiter)
    return;
  TimerOutput::Scope t(computing_timer, "Slope limiter");
  pcout << "    Slope limiter" << std::endl;

  Vector<double> limited_solution(dof_handler.n_dofs());

  // boundary_values->set_time(current_time);

  using Iterator = typename DoFHandler<dim>::active_cell_iterator;

  const auto cell_worker = [&](const Iterator &cell,
                               ScratchDataSlopeLimiter<dim> &scratch_data,
                               CopyDataSlopeLimiter &copy_data) {
    FEValues<dim> &fe_values = scratch_data.fe_values;
    FEValues<dim> &fe_values_neighbor = scratch_data.fe_values_neighbor;
    FEValues<dim> &fe_values_interpolate = scratch_data.fe_values_interpolate;

    fe_values.reinit(cell);
    fe_values_interpolate.reinit(cell);

    const unsigned int n_dofs = fe_values.get_fe().n_dofs_per_cell();
    copy_data.reinit(cell, n_dofs);

    std::vector<types::global_dof_index> &local_dof_indices =
        copy_data.local_dof_indices;
    std::vector<types::global_dof_index> &local_dof_indices_neighbor =
        copy_data.local_dof_indices_neighbor;
    const Point<dim> cell_center = cell->center();
    bool limit_cell = true;

    // Calculate the average of current cell
    double cell_average = 0;
    Tensor<1, dim> cell_average_grad;
    for (const unsigned int q_index : fe_values.quadrature_point_indices()) {
      for (const unsigned int i : fe_values.dof_indices()) {
        cell_average += solution[local_dof_indices[i]] *
                        fe_values.shape_value(i, q_index) *
                        fe_values.JxW(q_index);
        cell_average_grad += solution[local_dof_indices[i]] *
                             fe_values.shape_grad(i, q_index) *
                             fe_values.JxW(q_index);
      }
    }
    cell_average /= cell->measure();
    cell_average_grad /= cell->measure();

    // Calculate the averages of the neighbour cells
    std::vector<double> neighbor_cell_averages(cell->n_faces());
    std::vector<Tensor<1, dim>> neighbor_distance(cell->n_faces());
    std::vector<Point<dim>> face_centers(cell->n_faces());
    unsigned int n_neighbors = 0;
    for (const auto face_no : cell->face_indices()) {
      if (!cell->at_boundary(face_no) && cell->neighbor(face_no)->is_active()) {
        auto neighbor = cell->neighbor(face_no);
        // local_dof_indices_neighbor.resize(dofs_per_cell); //Not needed?
        neighbor->get_dof_indices(local_dof_indices_neighbor);
        fe_values_neighbor.reinit(neighbor);

        face_centers[n_neighbors] = cell->face(face_no)->center();
        neighbor_distance[n_neighbors] = neighbor->center() - cell->center();

        neighbor_cell_averages[n_neighbors] = 0;
        for (const unsigned int q_index :
             fe_values_neighbor.quadrature_point_indices()) {
          for (const unsigned int i : fe_values_neighbor.dof_indices()) {
            neighbor_cell_averages[n_neighbors] +=
                solution[local_dof_indices_neighbor[i]] *
                fe_values_neighbor.shape_value(i, q_index) *
                fe_values_neighbor.JxW(q_index);
          }
        }
        neighbor_cell_averages[n_neighbors] /= neighbor->measure();
        n_neighbors++;
      }
    }
    // neighbor_cell_averages.resize(n_neighbors);
    // neighbor_distance.resize(n_neighbors);
    // face_centers.resize(n_neighbors);

    // Check if the cell is meets the criteria for limiting
    if (limiter == SlopeLimiter::GerneralizedSlopeLimiter) {
      limit_cell = false;

      // TODO_BE: Check implementation - something goes wrong here :(
      double face_flux = 0;
      double face_value = 0;
      std::vector<double> tmp_minmod(n_neighbors + 1);

      // TODO_BE: Check if it is possible to create a quadrature for this
      //          independent of the cell
      Quadrature<dim> face_center_quadrature(face_centers);
      FEValues fe_vales_face_center(fe_values.get_fe(), face_center_quadrature,
                                    update_values);
      fe_vales_face_center.reinit(cell);
      std::vector<double> face_center_values(n_neighbors);
      fe_vales_face_center.get_function_values(solution, face_center_values);

      for (unsigned int i_neighbor = 0; i_neighbor < n_neighbors;
           i_neighbor++) {
        // Tensor<1, dim> direction =
        //     (face_centers[i_neighbor] - cell_center) /
        //     (face_centers[i_neighbor] - cell_center).norm();
        Tensor<1, dim> direction = neighbor_distance[i_neighbor] /
                                   neighbor_distance[i_neighbor].norm();

        face_value = face_center_values[i_neighbor];

        tmp_minmod[0] = (face_value - cell_average);
        for (unsigned int i2 = 0; i2 < n_neighbors; i2++) {
          tmp_minmod[i2 + 1] = (neighbor_cell_averages[i2] - cell_average) *
                               (neighbor_distance[i2] * direction) /
                               neighbor_distance[i2].norm();
        }
        face_flux = cell_average + minmod(tmp_minmod);

        if (std::abs(face_value - face_flux) > 1e-10) {
          limit_cell = true;
          break;
        }
      }
    }

    mark_for_limiter[cell->active_cell_index()] = limit_cell ? 1.0 : 0.0;

    if (limit_cell) {
      // Calculate the limited slope
      Tensor<1, dim> limited_slope;
      compute_limited_slope(cell_average, cell_average_grad,
                            neighbor_cell_averages, neighbor_distance,
                            n_neighbors, limited_slope);

      // To calculate the updated dof-values, we use a similar functionality
      // as VectroTools::interpolate

      // Calculate the support point values
      const std::vector<Point<dim>> &support_points =
          fe_values_interpolate.get_present_fe_values().get_quadrature_points();
      AssertDimension(support_points.size(), n_dofs);
      std::vector<Vector<double>> support_point_values(n_dofs,
                                                       Vector<double>(1));

      for (unsigned int i = 0; i < n_dofs; ++i) {
        support_point_values[i] =
            cell_average + (support_points[i] - cell_center) * limited_slope;
      }

      // Transformation in these cases needed, not implemented!
      Assert(fe.conforming_space != FiniteElementData<dim>::Hcurl &&
                 fe.conforming_space != FiniteElementData<dim>::Hdiv &&
                 fe.conforming_space != FiniteElementData<dim>::H1,
             ExcNotImplemented());

      fe.convert_generalized_support_point_values_to_dof_values(
          support_point_values, copy_data.cell_vector);
    } else {
      // No limiting, just copy the values
      for (const unsigned int i : fe_values.dof_indices()) {
        copy_data.cell_vector[i] = solution[local_dof_indices[i]];
      }
    }
  };

  const auto copier = [&](const CopyDataSlopeLimiter &c) {
    constraints.distribute_local_to_global(c.cell_vector, c.local_dof_indices,
                                           limited_solution);
  };

  Quadrature<dim> support_quadrature(fe.get_generalized_support_points());

  ScratchDataSlopeLimiter<dim> scratch_data(mapping, fe, quadrature_formula,
                                            support_quadrature);
  CopyDataSlopeLimiter copy_data;

  // TODO_BE: Till end() or end_active()?
  MeshWorker::mesh_loop(dof_handler.begin_active(), dof_handler.end(),
                        cell_worker, copier, scratch_data, copy_data,
                        MeshWorker::assemble_own_cells);

  solution = limited_solution;
}

template <int dim> void Sapphire::Hydro::BurgersEq<dim>::perform_time_step() {
  TimerOutput::Scope t(computing_timer, "Time step");
  pcout << "  Time step" << std::endl;

  old_solution = solution;
  current_solution = old_solution;
  Vector<double> tmp(dof_handler.n_dofs());

  switch (scheme) {
  case TimeSteppingScheme::ForwardEuler: {
    assemble_system();
    assemble_dg_vector();

    system_rhs.add(-1.0, dg_vector);
    system_rhs *= time_step;

    mass_matrix.vmult(tmp, current_solution);
    system_rhs.add(1.0, tmp);

    system_matrix.copy_from(mass_matrix);

    solve_linear_system();

    slope_limiter();

    break;
  }

  case TimeSteppingScheme::ExplicitRK: {
    // TODO_BE: Implement slope limiter for RK
    //  Butcher's array
    Vector<double> a({0.5, 0.5, 1.});
    Vector<double> b({1. / 6, 1. / 3, 1. / 3, 1. / 6});
    Vector<double> c({0., 0.5, 0.5, 1.});
    int i = 0;

    Vector<double> k1(dof_handler.n_dofs());
    i = 0;
    tmp = 0;
    current_time = time + c[i] * time_step;
    current_solution = old_solution;
    assemble_system();
    assemble_dg_vector();
    system_matrix.copy_from(mass_matrix);
    system_rhs.add(-1.0, dg_vector);
    solve_linear_system();
    slope_limiter();
    k1 = solution;

    Vector<double> k2(dof_handler.n_dofs());
    i = 1;
    tmp = 0;
    current_time = time + c[i] * time_step;
    current_solution = old_solution;
    current_solution.add(a[i - 1] * time_step, k1);
    assemble_system();
    assemble_dg_vector();
    system_matrix.copy_from(mass_matrix);
    system_rhs.add(-1.0, dg_vector);
    solve_linear_system();
    slope_limiter();
    k2 = solution;

    Vector<double> k3(dof_handler.n_dofs());
    i = 2;
    tmp = 0;
    current_time = time + c[i] * time_step;
    current_solution = old_solution;
    current_solution.add(a[i - 1] * time_step, k2);
    assemble_system();
    assemble_dg_vector();
    system_matrix.copy_from(mass_matrix);
    system_rhs.add(-1.0, dg_vector);
    solve_linear_system();
    slope_limiter();
    k3 = solution;

    Vector<double> k4(dof_handler.n_dofs());
    i = 3;
    tmp = 0;
    current_time = time + c[i] * time_step;
    current_solution = old_solution;
    current_solution.add(a[i - 1] * time_step, k3);
    assemble_system();
    assemble_dg_vector();
    system_matrix.copy_from(mass_matrix);
    system_rhs.add(-1.0, dg_vector);
    solve_linear_system();
    slope_limiter();
    k4 = solution;

    solution = old_solution;
    solution.add(b[0] * time_step, k1);
    solution.add(b[1] * time_step, k2);
    solution.add(b[2] * time_step, k3);
    solution.add(b[3] * time_step, k4);
    slope_limiter();
    break;
  }

  default:
    Assert(false, ExcNotImplemented());
    break;
  }

  time += time_step;
  timestep_number++;
}

template <int dim> void Sapphire::Hydro::BurgersEq<dim>::solve_linear_system() {
  TimerOutput::Scope t(computing_timer, "Solve linear system");
  pcout << "    Solve linear system" << std::endl;

  SolverControl solver_control(1000, 1e-12);
  // SolverCG<Vector<double>> solver(solver_control);
  SolverRichardson<Vector<double>> solver(solver_control);

  // PreconditionBlockJacobi preconditioner;
  // preconditioner.initialize(system_matrix);

  PreconditionSSOR<SparseMatrix<double>> preconditioner;
  preconditioner.initialize(system_matrix);

  solver.solve(system_matrix, solution, system_rhs, preconditioner);

  // constraints.distribute(solution);

  pcout << "      Solver converged in " << solver_control.last_step()
        << " iterations." << std::endl;
}

template <int dim>
void Sapphire::Hydro::BurgersEq<dim>::output_results() const {
  pcout << "  Output results" << std::endl;

  Vector<double> exact_solution_values(dof_handler.n_dofs());
  exact_solution->set_time(time);
  // TODO_BE: Also speciy mapping because why not...
  VectorTools::interpolate(dof_handler, *exact_solution, exact_solution_values);
  // VectorTools::interpolate(dof_handler, ExactSolution(a, time),
  // exact_solution);

  DataOut<dim> data_out;

  data_out.attach_dof_handler(dof_handler);
  data_out.add_data_vector(solution, "Solution");
  data_out.add_data_vector(exact_solution_values, "ExactSolution");
  data_out.add_data_vector(mark_for_limiter, "LimiterMark");

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

template <int dim> void Sapphire::Hydro::BurgersEq<dim>::process_results() {
  TimerOutput::Scope t(computing_timer, "Process results");
  pcout << "  Process results" << std::endl;

  Vector<float> difference_per_cell(triangulation.n_active_cells());
  exact_solution->set_time(time);

  // Use different quadrature for error computation
  const QTrapezoid<1> q_trapez;
  const QIterated<dim> q_iterated(q_trapez, fe.degree * 2 + 1);

  VectorTools::integrate_difference(mapping, dof_handler, solution,
                                    *exact_solution, difference_per_cell,
                                    q_iterated, VectorTools::L2_norm);
  float L2_error = VectorTools::compute_global_error(
      triangulation, difference_per_cell, VectorTools::L2_norm);
  pcout << "    L2 error:\t\t" << L2_error << std::endl;

  VectorTools::integrate_difference(mapping, dof_handler, solution,
                                    *exact_solution, difference_per_cell,
                                    q_iterated, VectorTools::Linfty_norm);
  float Linf_error = VectorTools::compute_global_error(
      triangulation, difference_per_cell, VectorTools::Linfty_norm);
  pcout << "    L-infinity error:\t" << Linf_error << std::endl;

  error_with_time.grow_or_shrink(error_with_time.size() + 1);
  error_with_time[error_with_time.size() - 1] = L2_error;
}

template <int dim> void Sapphire::Hydro::BurgersEq<dim>::run() {
  pcout << "Run conservation equation" << std::endl;
  make_grid();
  setup_system();
  assemble_mass_matrix();
  // assemble_system();
  {
    TimerOutput::Scope t(computing_timer, "Output results");
    output_results();
  }

  // const double time_end = 1.0;
  const double time_end = 0.4;
  const unsigned int n_steps = int(time_end / time_step);

  for (unsigned int i = 0; i < n_steps; ++i) {
    pcout << "Step " << timestep_number << "/" << n_steps << " (time = " << time
          << ")" << std::endl;
    perform_time_step();
    {
      TimerOutput::Scope t(computing_timer, "Output results");
      output_results();
    }
    process_results();
  }

  computing_timer.print_summary();
  computing_timer.reset();

  if (pcout.is_active()) {
    pcout << "   L2 error with time:" << std::endl;
    error_with_time.print(pcout.get_stream());
  }
}

// explicit instantiation
template class Sapphire::Hydro::BurgersEq<1>;
