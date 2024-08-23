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
 * @file mhd-solver.cpp
 * @author Florian Schulze (florian.schulze@mpi-hd.mpg.de)
 * @brief Implement @ref sapphirepp::MHD::MHDSolver
 */

#include "mhd-solver.h"

#include <deal.II/base/data_out_base.h>
#include <deal.II/base/discrete_time.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/base/utilities.h>
#include <deal.II/base/vectorization.h>

#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_update_flags.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/grid/filtered_iterator.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools.h>

#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/petsc_precondition.h>
#include <deal.II/lac/petsc_solver.h>
#include <deal.II/lac/sparsity_tools.h>
#include <deal.II/lac/vector_operation.h>

#include <deal.II/meshworker/assemble_flags.h>
#include <deal.II/meshworker/mesh_loop.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/vector_tools_integrate_difference.h>
#include <deal.II/numerics/vector_tools_project.h>

#include <array>
#include <cmath>
#include <fstream>
#include <string>
#include <vector>

#include "sapphirepp-logstream.h"


namespace sapphirepp
{
  namespace internal
  {
    namespace MHDSolver
    {
      using namespace dealii;

      // The mesh_loop function requires helper data types
      template <unsigned int dim, unsigned int spacedim>
      class ScratchData
      {
      public:
        // Constructor
        ScratchData(
          const Mapping<dim, spacedim>       &mapping,
          const FiniteElement<dim, spacedim> &fe,
          const Quadrature<dim>              &quadrature,
          const Quadrature<dim - 1>          &quadrature_face,
          const UpdateFlags update_flags = update_values | update_gradients |
                                           update_quadrature_points |
                                           update_JxW_values,
          const UpdateFlags face_update_flags = update_values |
                                                update_quadrature_points |
                                                update_JxW_values |
                                                update_normal_vectors,
          const UpdateFlags neighbor_face_update_flags = update_values)
          : fe_values(mapping, fe, quadrature, update_flags)
          , fe_values_face(mapping, fe, quadrature_face, face_update_flags)
          , fe_values_face_neighbor(mapping,
                                    fe,
                                    quadrature_face,
                                    neighbor_face_update_flags)
          , fe_values_subface_neighbor(mapping,
                                       fe,
                                       quadrature_face,
                                       neighbor_face_update_flags)
        {}

        // Copy Constructor
        ScratchData(const ScratchData<dim, spacedim> &scratch_data)
          : fe_values(scratch_data.fe_values.get_mapping(),
                      scratch_data.fe_values.get_fe(),
                      scratch_data.fe_values.get_quadrature(),
                      scratch_data.fe_values.get_update_flags())
          , fe_values_face(scratch_data.fe_values_face.get_mapping(),
                           scratch_data.fe_values_face.get_fe(),
                           scratch_data.fe_values_face.get_quadrature(),
                           scratch_data.fe_values_face.get_update_flags())
          , fe_values_face_neighbor(
              scratch_data.fe_values_face_neighbor.get_mapping(),
              scratch_data.fe_values_face_neighbor.get_fe(),
              scratch_data.fe_values_face_neighbor.get_quadrature(),
              scratch_data.fe_values_face_neighbor.get_update_flags())
          , fe_values_subface_neighbor(
              scratch_data.fe_values_subface_neighbor.get_mapping(),
              scratch_data.fe_values_subface_neighbor.get_fe(),
              scratch_data.fe_values_subface_neighbor.get_quadrature(),
              scratch_data.fe_values_subface_neighbor.get_update_flags())
        {}

        FEValues<dim, spacedim>        fe_values;
        FEFaceValues<dim, spacedim>    fe_values_face;
        FEFaceValues<dim, spacedim>    fe_values_face_neighbor;
        FESubfaceValues<dim, spacedim> fe_values_subface_neighbor;
      };



      struct CopyDataFace
      {
        Vector<double> cell_dg_rhs_1;
        Vector<double> cell_dg_rhs_2;

        std::vector<types::global_dof_index> local_dof_indices;
        std::vector<types::global_dof_index> local_dof_indices_neighbor;

        template <typename Iterator>
        void
        reinit(const Iterator &cell,
               const Iterator &neighbor_cell,
               unsigned int    dofs_per_cell)
        {
          cell_dg_rhs_1.reinit(dofs_per_cell);
          cell_dg_rhs_2.reinit(dofs_per_cell);

          local_dof_indices.resize(dofs_per_cell);
          cell->get_dof_indices(local_dof_indices);

          local_dof_indices_neighbor.resize(dofs_per_cell);
          neighbor_cell->get_dof_indices(local_dof_indices_neighbor);
        }
      };



      struct CopyData
      {
        Vector<double>                       cell_dg_rhs;
        std::vector<types::global_dof_index> local_dof_indices;
        // std::vector<types::global_dof_index> local_dof_indices_neighbor;
        std::vector<CopyDataFace> face_data;

        template <typename Iterator>
        void
        reinit(const Iterator &cell, unsigned int dofs_per_cell)
        {
          cell_dg_rhs.reinit(dofs_per_cell);

          local_dof_indices.resize(dofs_per_cell);
          cell->get_dof_indices(local_dof_indices);
        }
      };
    } // namespace MHDSolver
  }   // namespace internal
} // namespace sapphirepp


template <unsigned int dim>
sapphirepp::MHD::MHDSolver<dim>::MHDSolver(
  const MHDParameters<dim>      &mhd_parameters,
  const PhysicalParameters      &physical_parameters,
  const Utils::OutputParameters &output_parameters)
  : mhd_parameters{mhd_parameters}
  , physical_parameters{physical_parameters}
  , output_parameters{output_parameters}
  , mhd_equations(mhd_parameters.adiabatic_index)
  , numerical_flux(MHDEquations<spacedim>(
      mhd_parameters
        .adiabatic_index)) // TODO: Change back after correcting templates
  , mpi_communicator{MPI_COMM_WORLD}
  , triangulation(mpi_communicator)
  , dof_handler(triangulation)
  , mapping()
  , fe(FE_DGQ<dim, spacedim>(mhd_parameters.polynomial_degree),
       MHDEquations<dim>::n_components)
  , quadrature(fe.tensor_degree() + 1)
  , quadrature_face(fe.tensor_degree() + 1)
  , pcout(std::cout,
          ((Utilities::MPI::this_mpi_process(mpi_communicator) == 0) &&
           (saplog.get_verbosity() >= 3)))
  , timer(mpi_communicator, pcout, TimerOutput::never, TimerOutput::wall_times)
{
  LogStream::Prefix p("MHD", saplog);
  LogStream::Prefix p2("Constructor", saplog);
  saplog << mhd_flags << std::endl;
  saplog << "dim_mhd=" << dim << std::endl;
}



template <unsigned int dim>
void
sapphirepp::MHD::MHDSolver<dim>::setup()
{
  LogStream::Prefix p("MHD", saplog);
  saplog << "Setup MHD equation solver. \t[" << Utilities::System::get_time()
         << "]" << std::endl;
  LogStream::Prefix p2("Setup", saplog);
  timer.reset();
  make_grid();
  setup_system();

  {
    TimerOutput::Scope timer_section(timer, "Assemble mass matrix - MHD");
    saplog << "Assemble mass matrix." << std::endl;
    MatrixCreator::create_mass_matrix(mapping,
                                      dof_handler,
                                      quadrature,
                                      mass_matrix);
  }

  {
    TimerOutput::Scope timer_section(timer, "Project initial condition - MHD");
    InitialConditionMHD<spacedim> initial_condition_function(
      physical_parameters);
    PETScWrappers::MPI::Vector initial_condition(locally_owned_dofs,
                                                 mpi_communicator);
    project(initial_condition_function, initial_condition);
    // Here a non ghosted vector, is copied into a ghosted vector. I think
    // that is the moment where the ghost cells are filled.
    locally_relevant_current_solution = initial_condition;
  }
}



template <unsigned int dim>
void
sapphirepp::MHD::MHDSolver<dim>::run()
{
  setup();
  LogStream::Prefix p("MHD", saplog);

  DiscreteTime discrete_time(0,
                             mhd_parameters.final_time,
                             mhd_parameters.time_step);
  for (; discrete_time.is_at_end() == false; discrete_time.advance_time())
    {
      saplog << "Time step " << std::setw(6) << std::right
             << discrete_time.get_step_number()
             << " at t = " << discrete_time.get_current_time() << " \t["
             << Utilities::System::get_time() << "]" << std::endl;

      if ((discrete_time.get_step_number() %
           output_parameters.output_frequency) == 0)
        output_results(discrete_time.get_step_number(),
                       discrete_time.get_current_time());

      switch (mhd_parameters.time_stepping_method)
        {
          case TimeSteppingMethodMHD::forward_euler:
            forward_euler_method(discrete_time.get_current_time(),
                                 discrete_time.get_next_step_size());
            break;
          default:
            AssertThrow(false, ExcNotImplemented());
        }
    }

  // Output at the final result
  output_results(discrete_time.get_step_number(),
                 discrete_time.get_current_time());

  saplog << "Simulation ended at t = " << discrete_time.get_current_time()
         << " \t[" << Utilities::System::get_time() << "]" << std::endl;

  timer.print_wall_time_statistics(mpi_communicator);
}



template <unsigned int dim>
void
sapphirepp::MHD::MHDSolver<dim>::make_grid()
{
  TimerOutput::Scope timer_section(timer, "Grid setup - MHD");
  saplog << "Create the grid" << std::endl;

  switch (mhd_parameters.grid_type)
    {
      case GridTypeMHD::hypercube:
        {
          saplog << "Create the grid from hyper rectangle" << std::endl;
          GridGenerator::subdivided_hyper_rectangle(triangulation,
                                                    mhd_parameters.n_cells,
                                                    mhd_parameters.p1,
                                                    mhd_parameters.p2,
                                                    true);
          break;
        }
      case GridTypeMHD::file:
        {
          saplog << "Read grid from file \"" << mhd_parameters.grid_file << "\""
                 << std::endl;
          GridIn<dim, spacedim> grid_in(triangulation);
          std::ifstream         input(mhd_parameters.grid_file);
          grid_in.read_ucd(input);
          Assert(triangulation.all_reference_cells_are_hyper_cube(),
                 ExcNotImplemented("The grid must consist of hypercubes."));
          Assert(triangulation.has_hanging_nodes() == false,
                 ExcNotImplemented("The grid must not have hanging nodes."));
          break;
        }
      default:
        AssertThrow(false, ExcNotImplemented());
    }

  // GridGenerator::hyper_cube(triangulation, -5., 5., colorize);
  // triangulation.refine_global(6);
  saplog << "The grid was created: "
         << "	#cells=" << triangulation.n_cells()
         << ",	#active cells=" << triangulation.n_global_active_cells()
         << std::endl;

  // Periodic boundary conditions with MeshWorker. Mailing list
  // https://groups.google.com/g/dealii/c/WlOiww5UVxc/m/mtQJDUwiBQAJ
  //
  // "If you call add_periodicity() on a Triangulation object, the
  // periodic faces are treated as internal faces in MeshWorker. This
  // means that you will not access them in a "integrate_boundary_term"
  // function but in a "integrate_face_term" function. "
  std::vector<
    GridTools::PeriodicFacePair<typename Triangulation::cell_iterator>>
    matched_pairs;
  for (unsigned int i = 0; i < dim; ++i)
    {
      if (mhd_parameters.boundary_conditions[2 * i] ==
            BoundaryConditionsMHD::periodic or
          mhd_parameters.boundary_conditions[2 * i + 1] ==
            BoundaryConditionsMHD::periodic)
        {
          AssertThrow(mhd_parameters.boundary_conditions[2 * i] ==
                          BoundaryConditionsMHD::periodic and
                        mhd_parameters.boundary_conditions[2 * i + 1] ==
                          BoundaryConditionsMHD::periodic,
                      ExcMessage(
                        "Periodic boundary conditions did not match."));
          GridTools::collect_periodic_faces(
            triangulation, 2 * i, 2 * i + 1, i, matched_pairs);
        }
    }
  triangulation.add_periodicity(matched_pairs);
}



template <unsigned int dim>
void
sapphirepp::MHD::MHDSolver<dim>::setup_system()
{
  TimerOutput::Scope timer_section(timer, "FE system - MHD");
  saplog << "Setup the finite element system" << std::endl;

  dof_handler.clear();
  dof_handler.distribute_dofs(fe);

  const unsigned int n_dofs = dof_handler.n_dofs();
  locally_owned_dofs        = dof_handler.locally_owned_dofs();
  locally_relevant_dofs = DoFTools::extract_locally_relevant_dofs(dof_handler);

  saplog << "The degrees of freedom were distributed:"
         << "	n_dofs=" << n_dofs
         << ", locally_owned_dofs=" << locally_owned_dofs.n_elements()
         << ", locally_relevant_dofs=" << locally_relevant_dofs.n_elements()
         << std::endl;

  // constraints.clear();
  // DoFTools::make_periodicity_constraints(matched_pairs, constraints);

  // Let me see if I have to initialise the constraint object
  // constraints.clear();
  // constraints.reinit(locally_relevant_dofs);
  // This is an rather obscure line. I do not know why I need it. (cf.
  // example 23) constraints.close();

  // Vectors
  locally_owned_previous_solution.reinit(locally_owned_dofs, mpi_communicator);
  locally_relevant_current_solution.reinit(locally_owned_dofs,
                                           locally_relevant_dofs,
                                           mpi_communicator);
  dg_rhs.reinit(locally_owned_dofs, mpi_communicator);
  system_rhs.reinit(locally_owned_dofs, mpi_communicator);

  DynamicSparsityPattern dsp(locally_relevant_dofs);
  // NON-PERIODIC
  DoFTools::make_flux_sparsity_pattern(dof_handler, dsp);
  SparsityTools::distribute_sparsity_pattern(dsp,
                                             locally_owned_dofs,
                                             mpi_communicator,
                                             locally_relevant_dofs);
  saplog << "Number of local nonzero matrix elements: "
         << dsp.n_nonzero_elements() << std::endl;
  // PERIODIC
  // DoFTools::make_flux_sparsity_pattern(dof_handler, dsp, constraints,
  // false);
  mass_matrix.reinit(locally_owned_dofs,
                     locally_owned_dofs,
                     dsp,
                     mpi_communicator);
  system_matrix.reinit(locally_owned_dofs,
                       locally_owned_dofs,
                       dsp,
                       mpi_communicator);
}



template <unsigned int dim>
void
sapphirepp::MHD::MHDSolver<dim>::project(
  const Function<spacedim>   &f,
  PETScWrappers::MPI::Vector &projected_function) const
{
  saplog << "Project a function onto the finite element space" << std::endl;
  LogStream::Prefix p("project", saplog);

  // Create right hand side
  PETScWrappers::MPI::Vector rhs(locally_owned_dofs, mpi_communicator);
  VectorTools::create_right_hand_side(mapping, dof_handler, quadrature, f, rhs);

  // Solve the system
  PETScWrappers::PreconditionNone preconditioner;
  preconditioner.initialize(mass_matrix);

  SolverControl           solver_control(1000, 1e-12);
  PETScWrappers::SolverCG cg(solver_control, mpi_communicator);
  cg.solve(mass_matrix, projected_function, rhs, preconditioner);
  saplog << "Solved in " << solver_control.last_step() << " iterations."
         << std::endl;
  // At the moment I am assuming, that I do not have constraints. Hence, I
  // do not need the following line.
  // constraints.distribute(projected_function);
}



template <unsigned int dim>
void
sapphirepp::MHD::MHDSolver<dim>::assemble_dg_rhs(const double time)
{
  TimerOutput::Scope timer_section(timer, "DG rhs - MHD");

  using Iterator = typename DoFHandler<dim, spacedim>::active_cell_iterator;
  using namespace sapphirepp::internal::MHDSolver;

  static_cast<void>(time); // suppress compiler warning
  dg_rhs = 0;

  // assemble cell terms
  const auto cell_worker = [&](const Iterator             &cell,
                               ScratchData<dim, spacedim> &scratch_data,
                               CopyData                   &copy_data) {
    FEValues<dim, spacedim> &fe_v = scratch_data.fe_values;
    // reinit cell
    fe_v.reinit(cell);

    const unsigned int n_dofs = fe_v.get_fe().n_dofs_per_cell();
    // reinit the cell_matrix
    copy_data.reinit(cell, n_dofs);

    const std::vector<Point<spacedim>> &q_points = fe_v.get_quadrature_points();
    const std::vector<double>          &JxW      = fe_v.get_JxW_values();

    std::vector<Vector<double>> states(
      q_points.size(), Vector<double>(MHDEquations<dim>::n_components));
    typename MHDEquations<dim>::flux_type flux_matrix;

    fe_v.get_function_values(locally_relevant_current_solution, states);

    for (const unsigned int q_index : fe_v.quadrature_point_indices())
      {
        mhd_equations.compute_flux_matrix(states[q_index], flux_matrix);

        for (unsigned int i : fe_v.dof_indices())
          {
            // const unsigned int component_i =
            //   fe_v.get_fe().system_to_component_index(i).first;

            for (unsigned int c = 0; c < MHDEquations<dim>::n_components; ++c)
              {
                // F[c] * \grad \phi[c]_i
                copy_data.cell_dg_rhs(i) +=
                  flux_matrix[c] * fe_v.shape_grad_component(i, q_index, c) *
                  JxW[q_index];

                /** @todo Add source term */
              }
          }
      }
  };


  // assemble boundary face terms
  const auto boundary_worker = [&](const Iterator             &cell,
                                   const unsigned int         &face_no,
                                   ScratchData<dim, spacedim> &scratch_data,
                                   CopyData                   &copy_data) {
    FEFaceValues<dim, spacedim> &fe_v_face = scratch_data.fe_values_face;
    fe_v_face.reinit(cell, face_no);

    // NOTE: copy_data is not reinitialized, the cell_workers contribution
    // to the cell_dg_matrix should not be deleted

    const unsigned int boundary_id = cell->face(face_no)->boundary_id();
    const BoundaryConditionsMHD boundary_condition =
      mhd_parameters.boundary_conditions[boundary_id];

    const std::vector<Point<spacedim>> &q_points =
      fe_v_face.get_quadrature_points();
    const std::vector<double>              &JxW = fe_v_face.get_JxW_values();
    const std::vector<Tensor<1, spacedim>> &normals =
      fe_v_face.get_normal_vectors();

    // Initialise state and flux vectors
    std::vector<Vector<double>> states(
      q_points.size(), Vector<double>(MHDEquations<dim>::n_components));
    typename MHDEquations<dim>::state_type virtual_neighbor_state(
      MHDEquations<dim>::n_components);
    typename MHDEquations<dim>::state_type numerical_normal_flux(
      MHDEquations<dim>::n_components);

    // Compute states
    fe_v_face.get_function_values(locally_relevant_current_solution, states);

    for (unsigned int q_index : fe_v_face.quadrature_point_indices())
      {
        // Calculate virtual neighbor state according to boundary condition
        switch (boundary_condition)
          {
            case BoundaryConditionsMHD::zero_inflow:
              {
                /** @todo Check if this BC is equivalent to zero inflow */
                virtual_neighbor_state = states[q_index];
                break;
              }
            case BoundaryConditionsMHD::periodic:
            default:
              Assert(false, ExcNotImplemented());
          }

        numerical_flux.compute_numerical_normal_flux(normals[q_index],
                                                     states[q_index],
                                                     virtual_neighbor_state,
                                                     numerical_normal_flux);


        for (unsigned int i : fe_v_face.dof_indices())
          {
            // const unsigned int component_i =
            //   fe_v_face.get_fe().system_to_component_index(i).first;

            for (unsigned int c = 0; c < MHDEquations<dim>::n_components; ++c)
              {
                // - n * F[c] * \phi[c]_{i,-}
                copy_data.cell_dg_rhs(i) +=
                  -numerical_normal_flux[c] *
                  fe_v_face.shape_value_component(i, q_index, c) * JxW[q_index];
              }
          }
      }
  };

  // assemble interior face terms
  // NOTE: The face worker assumes a grid consisting of rectangular cells
  // with the following pattern of interior face normals
  //    ^
  // *--|--*
  // |     |
  // -     -->   ^ y
  // |     |     |
  // *-----*     *--> x
  // Hence, n_x = 1. and n_y = 1.
  const auto face_worker = [&](const Iterator             &cell,
                               const unsigned int         &face_no,
                               const unsigned int         &subface_no,
                               const Iterator             &neighbor_cell,
                               const unsigned int         &neighbor_face_no,
                               const unsigned int         &neighbor_subface_no,
                               ScratchData<dim, spacedim> &scratch_data,
                               CopyData                   &copy_data) {
    // We assumes that faces are only  assembled once. And hence subface_no, can
    // be ignored.
    static_cast<void>(subface_no);
    // We are not using mesh refinement yet. Hence neighbor_subface_no can be
    // ignored.
    static_cast<void>(neighbor_subface_no);

    FEFaceValues<dim, spacedim> &fe_v_face = scratch_data.fe_values_face;
    fe_v_face.reinit(cell, face_no);

    FEFaceValues<dim, spacedim> &fe_v_face_neighbor =
      scratch_data.fe_values_face_neighbor;
    fe_v_face_neighbor.reinit(neighbor_cell, neighbor_face_no);

    // Create an element at the end of the vector containing the face data
    copy_data.face_data.emplace_back();
    CopyDataFace      &copy_data_face = copy_data.face_data.back();
    const unsigned int n_dofs         = fe_v_face.get_fe().n_dofs_per_cell();
    copy_data_face.reinit(cell, neighbor_cell, n_dofs);

    const std::vector<Point<spacedim>> &q_points =
      fe_v_face.get_quadrature_points();
    const std::vector<double>              &JxW = fe_v_face.get_JxW_values();
    const std::vector<Tensor<1, spacedim>> &normals =
      fe_v_face.get_normal_vectors();

    // Initialise state and flux vectors
    std::vector<Vector<double>> states(
      q_points.size(), Vector<double>(MHDEquations<dim>::n_components));
    std::vector<Vector<double>> states_neighbor(
      q_points.size(), Vector<double>(MHDEquations<dim>::n_components));
    typename MHDEquations<dim>::state_type numerical_normal_flux(
      MHDEquations<dim>::n_components);

    // Compute states
    fe_v_face.get_function_values(locally_relevant_current_solution, states);
    fe_v_face_neighbor.get_function_values(locally_relevant_current_solution,
                                           states_neighbor);

    for (unsigned int q_index : fe_v_face.quadrature_point_indices())
      {
        numerical_flux.compute_numerical_normal_flux(normals[q_index],
                                                     states[q_index],
                                                     states_neighbor[q_index],
                                                     numerical_normal_flux);

        // cell_dg_rhs_1
        for (unsigned int i : fe_v_face.dof_indices())
          {
            // const unsigned int component_i =
            //   fe_v_face.get_fe().system_to_component_index(i).first;

            for (unsigned int c = 0; c < MHDEquations<dim>::n_components; ++c)
              {
                // - n * F[c] * \phi[c]_{i,-}
                copy_data_face.cell_dg_rhs_1(i) +=
                  -numerical_normal_flux[c] *
                  fe_v_face.shape_value_component(i, q_index, c) * JxW[q_index];
              }
          }

        // cell_dg_rhs_2
        for (unsigned int i : fe_v_face_neighbor.dof_indices())
          {
            // const unsigned int component_i =
            //   fe_v_face_neighbor.get_fe().system_to_component_index(i).first;

            for (unsigned int c = 0; c < MHDEquations<dim>::n_components; ++c)
              {
                // + n * F[c] * \phi[c]_{i,+}
                copy_data_face.cell_dg_rhs_2(i) +=
                  numerical_normal_flux[c] *
                  fe_v_face_neighbor.shape_value_component(i, q_index, c) *
                  JxW[q_index];
              }
          }
      }
  };

  // copier for the mesh_loop function
  const auto copier = [&](const CopyData &c) {
    constraints.distribute_local_to_global(c.cell_dg_rhs,
                                           c.local_dof_indices,
                                           dg_rhs);
    for (auto &cdf : c.face_data)
      {
        constraints.distribute_local_to_global(cdf.cell_dg_rhs_1,
                                               cdf.local_dof_indices,
                                               dg_rhs);
        constraints.distribute_local_to_global(cdf.cell_dg_rhs_2,
                                               cdf.local_dof_indices_neighbor,
                                               dg_rhs);
      }
  };

  ScratchData<dim, spacedim> scratch_data(mapping,
                                          fe,
                                          quadrature,
                                          quadrature_face);
  CopyData                   copy_data;
  saplog << "Begin the assembly of the DG rhs." << std::endl;
  const auto filtered_iterator_range =
    dof_handler.active_cell_iterators() | IteratorFilters::LocallyOwnedCell();
  MeshWorker::mesh_loop(filtered_iterator_range,
                        cell_worker,
                        copier,
                        scratch_data,
                        copy_data,
                        MeshWorker::assemble_own_cells |
                          MeshWorker::assemble_boundary_faces |
                          MeshWorker::assemble_ghost_faces_once |
                          MeshWorker::assemble_own_interior_faces_once,
                        boundary_worker,
                        face_worker);
  dg_rhs.compress(VectorOperation::add);
  saplog << "The DG rhs was assembled." << std::endl;
}



template <unsigned int dim>
void
sapphirepp::MHD::MHDSolver<dim>::forward_euler_method(const double time,
                                                      const double time_step)
{
  TimerOutput::Scope timer_section(timer, "FE method - MHD");
  LogStream::Prefix  p("FE_method", saplog);
  // Equation: mass_matrix  f(time + time_step) = mass_matrix f(time) +
  // time_step * rhs(time)

  // Assemble the right hand side
  assemble_dg_rhs(time);

  locally_owned_previous_solution = locally_relevant_current_solution;

  mass_matrix.vmult(system_rhs, locally_owned_previous_solution);
  system_rhs.add(time_step, dg_rhs);

  SolverControl              solver_control(1000, 1e-6 * system_rhs.l2_norm());
  PETScWrappers::SolverGMRES solver(solver_control, mpi_communicator);

  PETScWrappers::PreconditionBlockJacobi preconditioner;
  preconditioner.initialize(mass_matrix);

  // The x vector has to be a vector of locally_owned_dofs, so we make use of
  // the locally_owned_previous_solution
  solver.solve(mass_matrix,
               locally_owned_previous_solution,
               system_rhs,
               preconditioner);

  // Update the solution
  locally_relevant_current_solution = locally_owned_previous_solution;

  saplog << "Solver converged in " << solver_control.last_step()
         << " iterations." << std::endl;
}



template <unsigned int dim>
void
sapphirepp::MHD::MHDSolver<dim>::output_results(
  const unsigned int time_step_number,
  const double       cur_time)
{
  TimerOutput::Scope     timer_section(timer, "Output - MHD");
  DataOut<dim, spacedim> data_out;
  data_out.attach_dof_handler(dof_handler);
  /**
   * @todo Use
   * @dealref{DataComponentInterpretation,namespaceDataComponentInterpretation}
   * to output u and B as a vector field. Problem: At the moment deal.II only
   * allows to output a `dim` dimensional vector field, with `dim` being the
   * dimension of the grid.
   */
  data_out.add_data_vector(locally_relevant_current_solution,
                           MHDEquations<dim>::create_component_name_list());

  // Output the partition of the mesh
  Vector<float> subdomain(triangulation.n_active_cells());
  for (unsigned int i = 0; i < subdomain.size(); ++i)
    subdomain(i) = static_cast<float>(triangulation.locally_owned_subdomain());
  data_out.add_data_vector(subdomain, "subdomain");

  // Adapt the output to the polynomial degree of the shape functions
  data_out.build_patches(mhd_parameters.polynomial_degree);
  output_parameters.write_results<dim, spacedim>(data_out,
                                                 time_step_number,
                                                 cur_time);
}



template <unsigned int dim>
double
sapphirepp::MHD::MHDSolver<dim>::compute_global_error(
  const Function<spacedim>         &exact_solution,
  const VectorTools::NormType      &cell_norm,
  const VectorTools::NormType      &global_norm,
  const Function<spacedim, double> *weight) const
{
  LogStream::Prefix p("MHD", saplog);
  saplog << "Compute the global error" << std::endl;
  LogStream::Prefix p2("Error", saplog);

  Vector<float> cell_errors(triangulation.n_locally_owned_active_cells());

  const QTrapezoid<1>  q_trapezoid;
  const QIterated<dim> q_iterated(q_trapezoid,
                                  mhd_parameters.polynomial_degree * 2 + 1);

  VectorTools::integrate_difference(mapping,
                                    dof_handler,
                                    locally_relevant_current_solution,
                                    exact_solution,
                                    cell_errors,
                                    q_iterated,
                                    cell_norm,
                                    weight);

  double global_error =
    VectorTools::compute_global_error(triangulation, cell_errors, global_norm);

  saplog << "Global error: " << global_error << std::endl;

  return global_error;
}



template <unsigned int dim>
double
sapphirepp::MHD::MHDSolver<dim>::compute_weighted_norm(
  const VectorTools::NormType      &cell_norm,
  const VectorTools::NormType      &global_norm,
  const Function<spacedim, double> *weight) const
{
  LogStream::Prefix p("MHD", saplog);
  saplog << "Compute the weighted norm" << std::endl;
  LogStream::Prefix p2("Norm", saplog);

  Vector<float> cell_norms(triangulation.n_locally_owned_active_cells());

  const QTrapezoid<1>  q_trapezoid;
  const QIterated<dim> q_iterated(q_trapezoid,
                                  mhd_parameters.polynomial_degree * 2 + 1);

  const Functions::ZeroFunction<spacedim> zero_function(
    mhd_equations.n_components);

  VectorTools::integrate_difference(mapping,
                                    dof_handler,
                                    locally_relevant_current_solution,
                                    zero_function,
                                    cell_norms,
                                    q_iterated,
                                    cell_norm,
                                    weight);

  double global_weighted_norm =
    VectorTools::compute_global_error(triangulation, cell_norms, global_norm);

  saplog << "Global weighted norm: " << global_weighted_norm << std::endl;

  return global_weighted_norm;
}



template <unsigned int dim>
const sapphirepp::MHD::MHDEquations<dim> &
sapphirepp::MHD::MHDSolver<dim>::get_mhd_equations() const
{
  return mhd_equations;
}



template <unsigned int dim>
const typename sapphirepp::MHD::MHDSolver<dim>::Triangulation &
sapphirepp::MHD::MHDSolver<dim>::get_triangulation() const
{
  return triangulation;
}



template <unsigned int dim>
const dealii::DoFHandler<dim, sapphirepp::MHD::MHDSolver<dim>::spacedim> &
sapphirepp::MHD::MHDSolver<dim>::get_dof_handler() const
{
  return dof_handler;
}



template <unsigned int dim>
const dealii::PETScWrappers::MPI::Vector &
sapphirepp::MHD::MHDSolver<dim>::get_current_solution() const
{
  return locally_relevant_current_solution;
}



template <unsigned int dim>
const dealii::TimerOutput &
sapphirepp::MHD::MHDSolver<dim>::get_timer() const
{
  return timer;
}



// Explicit instantiations
template class sapphirepp::MHD::MHDSolver<1>;
template class sapphirepp::MHD::MHDSolver<2>;
template class sapphirepp::MHD::MHDSolver<3>;
