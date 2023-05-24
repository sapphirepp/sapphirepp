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

#include <deal.II/meshworker/mesh_loop.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/vector_tools.h>

#include <fstream>
#include <iostream>

template <int dim>
void Sapphire::Hydro::ExactSolution<dim>::vector_value(
    const Point<dim> &p, Vector<double> &values) const {
  AssertThrow(dim == 1, ExcNotImplemented());
  AssertDimension(values.size(), dim);
  // values(0) = std::sin(numbers::PI * (p(0) - a * this->get_time()));
  // values(0) = 1.0;
  const double sigma = 0.1;
  values(0) = std::exp(-(p[0] - a * this->get_time()) *
                       (p[0] - a * this->get_time()) / (2.0 * sigma * sigma));
}

// explicit template instantiation
template class Sapphire::Hydro::ExactSolution<1>;

template <int dim>
void Sapphire::Hydro::InitialCondition<dim>::vector_value(
    const Point<dim> &p, Vector<double> &values) const {
  ExactSolution<dim>(a, 0.0).vector_value(p, values);
}

// explicit template instantiation
template class Sapphire::Hydro::InitialCondition<1>;

template <int dim>
Sapphire::Hydro::ConservationEq<dim>::ConservationEq()
    : mpi_communicator(MPI_COMM_WORLD), mapping(), fe(1),
      dof_handler(triangulation), quadrature_formula(fe.tensor_degree() + 1),
      face_quadrature_formula(fe.tensor_degree()),
      pcout(std::cout,
            (Utilities::MPI::this_mpi_process(mpi_communicator) == 0)),
      computing_timer(mpi_communicator, pcout, TimerOutput::never,
                      TimerOutput::wall_times) {
  pcout << "Setup conservation equation" << std::endl;
  AssertThrow(dim == 1, ExcNotImplemented());
  time = 0.0;
  // time_step = 0.001;
  time_step = 0.01;
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

  VectorTools::interpolate(dof_handler, InitialCondition<dim>(a), solution);

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
          copy_data.cell_dg_matrix(i, j) -=
              a * (fe_values.shape_grad(i, q_index)[0] *
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

    for (const unsigned int q_index :
         fe_face_values.quadrature_point_indices()) {

      for (const unsigned int i : fe_face_values.dof_indices()) {
        for (const unsigned int j : fe_face_values.dof_indices()) {

          copy_data.cell_dg_matrix(i, j) +=
              0.5 * a *
              (std::abs(fe_face_values.normal_vector(q_index)[0]) +
               fe_face_values.normal_vector(q_index)[0]) *
              (fe_face_values.shape_value(i, q_index) *
               fe_face_values.shape_value(j, q_index) *
               fe_face_values.JxW(q_index));
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

          for (const unsigned int i : fe_face_values.dof_indices()) {
            for (const unsigned int j : fe_face_values.dof_indices()) {

              copy_data_face.cell_dg_matrix_11(i, j) +=
                  0.5 * a * fe_face_values.normal_vector(q_index)[0] *
                  (fe_face_values.shape_value(i, q_index) *
                   fe_face_values.shape_value(j, q_index) *
                   fe_face_values.JxW(q_index));

              copy_data_face.cell_dg_matrix_21(i, j) -=
                  0.5 * a * fe_face_values.normal_vector(q_index)[0] *
                  (fe_face_values_neighbor.shape_value(i, q_index) *
                   fe_face_values.shape_value(j, q_index) *
                   fe_face_values.JxW(q_index));

              copy_data_face.cell_dg_matrix_12(i, j) +=
                  0.5 * a * fe_face_values.normal_vector(q_index)[0] *
                  (fe_face_values.shape_value(i, q_index) *
                   fe_face_values_neighbor.shape_value(j, q_index) *
                   fe_face_values.JxW(q_index));

              copy_data_face.cell_dg_matrix_22(i, j) -=
                  0.5 * a * fe_face_values.normal_vector(q_index)[0] *
                  (fe_face_values_neighbor.shape_value(i, q_index) *
                   fe_face_values_neighbor.shape_value(j, q_index) *
                   fe_face_values.JxW(q_index));

              // upwind flux
              const double eta = 1.0;

              copy_data_face.cell_dg_matrix_11(i, j) +=
                  0.5 * eta * a * fe_face_values.normal_vector(q_index)[0] *
                  (fe_face_values.shape_value(i, q_index) *
                   fe_face_values.shape_value(j, q_index) *
                   fe_face_values.JxW(q_index));

              copy_data_face.cell_dg_matrix_21(i, j) -=
                  0.5 * eta * a * fe_face_values.normal_vector(q_index)[0] *
                  (fe_face_values_neighbor.shape_value(i, q_index) *
                   fe_face_values.shape_value(j, q_index) *
                   fe_face_values.JxW(q_index));

              copy_data_face.cell_dg_matrix_12(i, j) -=
                  0.5 * eta * a * fe_face_values.normal_vector(q_index)[0] *
                  (fe_face_values.shape_value(i, q_index) *
                   fe_face_values_neighbor.shape_value(j, q_index) *
                   fe_face_values.JxW(q_index));

              copy_data_face.cell_dg_matrix_22(i, j) +=
                  0.5 * eta * a * fe_face_values.normal_vector(q_index)[0] *
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
      for (const unsigned int i : fe_values.dof_indices()) {
        for (const unsigned int j : fe_values.dof_indices()) {
          cell_mass_matrix(i, j) +=
              (fe_values.shape_value(i, q_index) *
               fe_values.shape_value(j, q_index) * fe_values.JxW(q_index));
          cell_dg_matrix(i, j) -=
              a * (fe_values.shape_grad(i, q_index)[0] *
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

          for (const unsigned int i : fe_face_values.dof_indices()) {
            for (const unsigned int j : fe_face_values.dof_indices()) {

              cell_dg_matrix_11(i, j) +=
                  0.5 * a *
                  (std::abs(fe_face_values.normal_vector(q_index)[0]) +
                   fe_face_values.normal_vector(q_index)[0]) *
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

          for (const unsigned int i : fe_face_values.dof_indices()) {
            for (const unsigned int j : fe_face_values.dof_indices()) {

              cell_dg_matrix_11(i, j) +=
                  0.5 * a * fe_face_values.normal_vector(q_index)[0] *
                  (fe_face_values.shape_value(i, q_index) *
                   fe_face_values.shape_value(j, q_index) *
                   fe_face_values.JxW(q_index));

              cell_dg_matrix_21(i, j) -=
                  0.5 * a * fe_face_values.normal_vector(q_index)[0] *
                  (fe_face_values_neighbor.shape_value(i, q_index) *
                   fe_face_values.shape_value(j, q_index) *
                   fe_face_values.JxW(q_index));

              cell_dg_matrix_12(i, j) +=
                  0.5 * a * fe_face_values.normal_vector(q_index)[0] *
                  (fe_face_values.shape_value(i, q_index) *
                   fe_face_values_neighbor.shape_value(j, q_index) *
                   fe_face_values.JxW(q_index));

              cell_dg_matrix_22(i, j) -=
                  0.5 * a * fe_face_values.normal_vector(q_index)[0] *
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

  tmp = 0;
  mass_matrix.vmult(tmp, old_solution);
  system_rhs.add(1.0, tmp);

  tmp = 0;
  dg_matrix.vmult(tmp, old_solution);
  system_rhs.add(-time_step, tmp);

  system_matrix.copy_from(mass_matrix);

  // Theta method
  // const double theta = 0.5;

  // mass_matrix.vmult(tmp, old_solution);
  // system_rhs.add(1.0, tmp);

  // dg_matrix.vmult(tmp, old_solution);
  // system_rhs.add(-(1. - theta) * time_step, tmp);

  // system_matrix.copy_from(mass_matrix);
  // system_matrix.add(theta * time_step, dg_matrix);

  /** constant solution */
  // mass_matrix.vmult(tmp, old_solution);
  // system_rhs = tmp;
}

template <int dim> void Sapphire::Hydro::ConservationEq<dim>::solve() {
  TimerOutput::Scope t(computing_timer, "Solve");
  pcout << "Solve" << std::endl;

  SolverControl solver_control(1000, 1e-12);
  SolverCG<Vector<double>> solver(solver_control);

  // PreconditionSSOR<SparseMatrix<double>> preconditioner;
  // preconditioner.initialize(system_matrix, 1.2);

  PreconditionIdentity preconditioner;

  solver.solve(system_matrix, solution, system_rhs, preconditioner);

  // constraints.distribute(solution);

  pcout << "   " << solver_control.last_step() << " CG iterations."
        << std::endl;
}

template <int dim>
void Sapphire::Hydro::ConservationEq<dim>::output_results() const {
  pcout << "Output results" << std::endl;

  Vector<double> exact_solution(dof_handler.n_dofs());
  VectorTools::interpolate(dof_handler, ExactSolution<dim>(a, time + time_step),
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

template <int dim> void Sapphire::Hydro::ConservationEq<dim>::run() {
  pcout << "Run conservation equation" << std::endl;
  make_grid();
  setup_system();
  {
    TimerOutput::Scope t(computing_timer, "Output results");
    output_results();
  }
  timestep_number++;

  for (; time < 1; time += time_step, ++timestep_number) {
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
  }

  computing_timer.print_summary();
  computing_timer.reset();
}

// explicit instantiation
template class Sapphire::Hydro::ConservationEq<1>;
