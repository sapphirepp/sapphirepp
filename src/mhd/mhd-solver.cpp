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
#include <deal.II/base/table.h>
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
#include <deal.II/lac/full_matrix.h>
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
#include <functional>
#include <limits>
#include <string>
#include <vector>

#include "sapphirepp-logstream.h"


namespace sapphirepp
{

  /** @cond sapinternal */
  namespace sapinternal
  {
    namespace MHDSolver
    {
      using namespace dealii;

      // The mesh_loop function requires helper data types
      template <unsigned int dim>
      class ScratchData
      {
      public:
        // Constructor
        ScratchData(
          const Mapping<dim>        &mapping,
          const FiniteElement<dim>  &fe,
          const Quadrature<dim>     &quadrature,
          const Quadrature<dim - 1> &quadrature_face,
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
        ScratchData(const ScratchData<dim> &scratch_data)
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

        FEValues<dim>        fe_values;
        FEFaceValues<dim>    fe_values_face;
        FEFaceValues<dim>    fe_values_face_neighbor;
        FESubfaceValues<dim> fe_values_subface_neighbor;
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
        unsigned int                         cell_index;
        double                               min_dt;
        Vector<double>                       cell_dg_rhs;
        std::vector<types::global_dof_index> local_dof_indices;
        // std::vector<types::global_dof_index> local_dof_indices_neighbor;
        std::vector<CopyDataFace> face_data;

        template <typename Iterator>
        void
        reinit(const Iterator &cell, unsigned int dofs_per_cell)
        {
          cell_index = cell->active_cell_index();
          min_dt     = 0;

          cell_dg_rhs.reinit(dofs_per_cell);

          local_dof_indices.resize(dofs_per_cell);
          cell->get_dof_indices(local_dof_indices);
        }
      };



      // The helper data types for slope limiter mesh_loop
      template <unsigned int dim>
      class ScratchDataSlopeLimiter
      {
      public:
        // Constructor
        ScratchDataSlopeLimiter(const Mapping<dim>       &mapping,
                                const FiniteElement<dim> &fe,
                                const Quadrature<dim>    &quadrature)
          : fe_values_gradient(mapping,
                               fe,
                               quadrature,
                               update_gradients | update_JxW_values)
          , fe_values_support(mapping,
                              fe,
                              Quadrature<dim>(
                                fe.get_generalized_support_points()),
                              update_quadrature_points | update_values)
        {}

        // Copy Constructor
        ScratchDataSlopeLimiter(
          const ScratchDataSlopeLimiter<dim> &scratch_data)
          : fe_values_gradient(
              scratch_data.fe_values_gradient.get_mapping(),
              scratch_data.fe_values_gradient.get_fe(),
              scratch_data.fe_values_gradient.get_quadrature(),
              scratch_data.fe_values_gradient.get_update_flags())
          , fe_values_support(scratch_data.fe_values_support.get_mapping(),
                              scratch_data.fe_values_support.get_fe(),
                              scratch_data.fe_values_support.get_quadrature(),
                              scratch_data.fe_values_support.get_update_flags())
        {}

        FEValues<dim> fe_values_gradient;
        FEValues<dim> fe_values_support;
      };



      struct CopyDataSlopeLimiter
      {
        /**
         * @todo Use a `dealii::Vector<double>` for `cell_dof_values`. At the
         *       moment using
         *       `convert_generalized_support_point_values_to_dof_values` forces
         *       us to use `std::vector`.
         */
        std::vector<types::global_dof_index>     local_dof_indices;
        std::vector<double>                      cell_dof_values;
        ArrayView<const types::global_dof_index> cell_indices; // unused

        template <typename Iterator>
        void
        reinit(const Iterator                   &cell,
               unsigned int                      dofs_per_cell,
               const PETScWrappers::MPI::Vector &solution)
        {
          local_dof_indices.resize(dofs_per_cell);
          cell->get_dof_indices(local_dof_indices);

          cell_dof_values.resize(dofs_per_cell);
          cell->get_dof_values(solution,
                               cell_dof_values.begin(),
                               cell_dof_values.end());

          // std::vector<types::global_dof_index> indices_vector(dofs_per_cell);
          // std::iota(indices_vector.begin(), indices_vector.end(), 0);
          // cell_indices =
          //   ArrayView<types::global_dof_index>(indices_vector.data(),
          //                                      indices_vector.size());
        }
      };
    } // namespace MHDSolver
  } // namespace sapinternal
} // namespace sapphirepp
/** @endcond */


template <unsigned int dim>
sapphirepp::MHD::MHDSolver<dim>::MHDSolver(
  const MHDParameters<dim>      &mhd_parameters,
  const PhysicalParameters      &physical_parameters,
  const Utils::OutputParameters &output_parameters)
  : mhd_parameters{mhd_parameters}
  , physical_parameters{physical_parameters}
  , output_parameters{output_parameters}
  , mhd_equations(mhd_parameters.adiabatic_index)
  , numerical_flux(mhd_equations)
  , slope_limiter()
  , mpi_communicator{MPI_COMM_WORLD}
  , triangulation(mpi_communicator)
  , dof_handler(triangulation)
  , mapping()
  , fe(FE_DGQ<dim>(mhd_parameters.polynomial_degree), n_components)
  , quadrature(fe.tensor_degree() + 1)
  , quadrature_face(fe.tensor_degree() + 1)
  , pcout(saplog.to_condition_ostream(3))
  , timer(mpi_communicator, pcout, TimerOutput::never, TimerOutput::wall_times)
{
  LogStream::Prefix p("MHD", saplog);
  LogStream::Prefix p2("Constructor", saplog);
  saplog << mhd_flags << std::endl;
  saplog << "dim_mhd=" << dim << std::endl;

  AssertThrow(
    (1 <= dim) && (dim <= 3),
    ExcMessage("The dimension must be greater than or equal to one and smaller "
               "or equal to three."));

  if ((mhd_flags & MHDFlags::conserved_limiting) != MHDFlags::none)
    {
      AssertThrow((mhd_flags & MHDFlags::no_limiting) != MHDFlags::none,
                  ExcMessage("Limiting must be activated to use limiting "
                             "on conserved variables."));
    }
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
    InitialConditionMHD<dim> initial_condition_function(
      physical_parameters, mhd_equations.adiabatic_index);
    project(initial_condition_function, locally_owned_solution);
    // Here a non ghosted vector, is copied into a ghosted vector. I think
    // that is the moment where the ghost cells are filled.
    locally_relevant_current_solution = locally_owned_solution;
  }

  apply_limiter();
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
          case TimeSteppingMethodMHD::erk2:
          case TimeSteppingMethodMHD::erk4:
            explicit_runge_kutta(discrete_time.get_current_time(),
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
          GridIn<dim>   grid_in(triangulation);
          std::ifstream input(mhd_parameters.grid_file);
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
  locally_owned_solution.reinit(locally_owned_dofs, mpi_communicator);
  locally_relevant_current_solution.reinit(locally_owned_dofs,
                                           locally_relevant_dofs,
                                           mpi_communicator);
  dg_rhs.reinit(locally_owned_dofs, mpi_communicator);
  system_rhs.reinit(locally_owned_dofs, mpi_communicator);

  cell_average.resize(triangulation.n_active_cells(),
                      Vector<double>(n_components));
  shock_indicator.reinit(triangulation.n_active_cells());
  positivity_limiter_indicator.reinit(triangulation.n_active_cells());
  magnetic_divergence.reinit(triangulation.n_active_cells());
  cell_dt.reinit(triangulation.n_active_cells());

  DynamicSparsityPattern       dsp(locally_relevant_dofs);
  Table<2, DoFTools::Coupling> coupling(n_components, n_components);
  coupling.fill(DoFTools::Coupling::none);
  for (unsigned int i = 0; i < n_components; i++)
    coupling(i, i) = DoFTools::Coupling::always;

  // NON-PERIODIC
  DoFTools::make_sparsity_pattern(dof_handler, coupling, dsp);
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
}



template <unsigned int dim>
void
sapphirepp::MHD::MHDSolver<dim>::project(
  const Function<dim>        &f,
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
sapphirepp::MHD::MHDSolver<dim>::compute_cell_average()
{
  TimerOutput::Scope t(timer, "Cell average - MHD");
  saplog << "Compute cell average" << std::endl;
  AssertDimension(cell_average.size(), triangulation.n_active_cells());

  FEValues<dim> fe_v(mapping,
                     fe,
                     quadrature,
                     update_values | update_JxW_values);

  const unsigned int          n_q_points = quadrature.size();
  std::vector<Vector<double>> solution_values(n_q_points,
                                              Vector<double>(n_components));

  // Compute cell_average for locally_active and gost cells
  for (const auto &cell : dof_handler.active_cell_iterators())
    if (!cell->is_artificial())
      {
        fe_v.reinit(cell);
        const unsigned int cell_index = cell->active_cell_index();

        fe_v.get_function_values(locally_relevant_current_solution,
                                 solution_values);
        const std::vector<double> &JxW = fe_v.get_JxW_values();

        cell_average[cell_index] = 0.;

        for (const unsigned int q_index : fe_v.quadrature_point_indices())
          for (unsigned int c = 0; c < n_components; ++c)
            cell_average[cell_index][c] +=
              solution_values[q_index][c] * JxW[q_index];

        cell_average[cell_index] /= cell->measure();
      }
}



template <unsigned int dim>
void
sapphirepp::MHD::MHDSolver<dim>::compute_shock_indicator()
{
  if constexpr ((mhd_flags & MHDFlags::no_shock_indicator) != MHDFlags::none)
    {
      shock_indicator = 2.;
      return;
    }
  TimerOutput::Scope timer_section(timer, "Shock indicator - MHD");
  saplog << "Compute shock indicator" << std::endl;
  LogStream::Prefix p("ShockIndicator", saplog);

  using Iterator = typename DoFHandler<dim>::active_cell_iterator;
  using namespace sapphirepp::sapinternal::MHDSolver;

  shock_indicator = 0;

  // empty cell and boundary worker
  std::function<void(const Iterator &, ScratchData<dim> &, CopyData &)>
    empty_cell_worker;
  std::function<void(
    const Iterator &, const unsigned int &, ScratchData<dim> &, CopyData &)>
    empty_boundary_worker;

  // assemble interior face terms
  const auto face_worker = [&](const Iterator     &cell,
                               const unsigned int &face_no,
                               const unsigned int &subface_no,
                               const Iterator     &neighbor_cell,
                               const unsigned int &neighbor_face_no,
                               const unsigned int &neighbor_subface_no,
                               ScratchData<dim>   &scratch_data,
                               CopyData           &copy_data) {
    static_cast<void>(copy_data);
    Assert((cell->level() == neighbor_cell->level()) &&
             (neighbor_cell->has_children() == false),
           ExcNotImplemented("Mesh refinement not yet implemented."));
    static_cast<void>(subface_no);
    static_cast<void>(neighbor_subface_no);

    FEFaceValues<dim> &fe_v_face = scratch_data.fe_values_face;
    fe_v_face.reinit(cell, face_no);
    const unsigned int cell_index           = cell->active_cell_index();
    double             cell_shock_indicator = 0.;

    // saplog << "Indicator in cell " << cell_index << " on face " << face_no
    //        << std::endl;


    // Only integrate over inflow cell boundaries.
    // We use face by face discrimination using the cell_average flow and an
    // abitany face normal
    /**
     * @todo Should we use the MHD flux instead of momentum to discriminate
     *       outflow faces? For the energy component, there is an additional
     *       contribution pointing in direction of the B-field.
     */
    Tensor<1, dim> momentum;
    for (unsigned int d = 0; d < dim; ++d)
      momentum[d] = cell_average[cell_index][first_momentum_component + d];
    // Return in case of outflow face
    if (momentum * fe_v_face.normal_vector(0) >= 0)
      {
        // saplog << "Outflow face" << std::endl;
        return;
      }
    // Continue on inflow faces
    // saplog << "Inflow face" << std::endl;


    FEFaceValues<dim> &fe_v_face_neighbor =
      scratch_data.fe_values_face_neighbor;
    fe_v_face_neighbor.reinit(neighbor_cell, neighbor_face_no);

    // Use Extractor for indicator variable y, e.g.energy or density
    const FEValuesExtractors::Scalar variable(energy_component);
    const unsigned int               n_q_points = fe_v_face.n_quadrature_points;
    const std::vector<double>       &JxW        = fe_v_face.get_JxW_values();
    double                           face_norm  = 0.;
    std::vector<double>              face_values(n_q_points);
    std::vector<double>              face_values_neighbor(n_q_points);


    fe_v_face[variable].get_function_values(locally_relevant_current_solution,
                                            face_values);
    fe_v_face_neighbor[variable].get_function_values(
      locally_relevant_current_solution, face_values_neighbor);

    // saplog << "Values face:";
    // for (const auto &tmp : face_values)
    //   saplog << " " << tmp;
    // saplog << std::endl;
    // saplog << "Values neighbor:";
    // for (const auto &tmp : face_values_neighbor)
    //   saplog << " " << tmp;
    // saplog << std::endl;


    for (unsigned int q_index : fe_v_face.quadrature_point_indices())
      {
        // (y_j - y_nb)
        cell_shock_indicator +=
          (face_values[q_index] - face_values_neighbor[q_index]) * JxW[q_index];
        face_norm += JxW[q_index];
      }
    // saplog << "indicator value: " << cell_shock_indicator << std::endl;


    // Normalize the indicator variable
    const double dx     = cell->minimum_vertex_distance();
    const double degree = fe_v_face.get_fe().tensor_degree();
    // saplog << "degree=" << degree << ", dx=" << dx << std::endl;
    const double cell_norm =
      std::fabs(cell_average[cell_index][energy_component]);
    const double normalization =
      std::pow(dx, 0.5 * (degree + 1)) * face_norm * cell_norm;

    cell_shock_indicator = std::fabs(cell_shock_indicator) / normalization;
    // saplog << "indicator value: " << cell_shock_indicator << std::endl;
    shock_indicator[cell_index] += cell_shock_indicator;
    // saplog << "shock indicator: " << shock_indicator[cell_index] <<
    // std::endl;
  };

  /** @todo Use valid or empty copier? */
  // copier for the mesh_loop function
  const auto copier = [&](const CopyData &c) {
    static_cast<void>(c);
    // Empty
  };
  static_cast<void>(copier);
  // empty copier
  std::function<void(const CopyData &)> empty_copier;

  ScratchData<dim> scratch_data(mapping, fe, quadrature, quadrature_face);
  CopyData         copy_data;
  saplog << "Begin the assembly of the shock indicator." << std::endl;
  const auto filtered_iterator_range =
    dof_handler.active_cell_iterators() | IteratorFilters::LocallyOwnedCell();
  MeshWorker::mesh_loop(filtered_iterator_range,
                        empty_cell_worker,
                        empty_copier,
                        scratch_data,
                        copy_data,
                        MeshWorker::assemble_ghost_faces_both |
                          MeshWorker::assemble_own_interior_faces_both,
                        empty_boundary_worker,
                        face_worker);
  saplog << "The shock indicator was assembled." << std::endl;
}



template <unsigned int dim>
bool
sapphirepp::MHD::MHDSolver<dim>::indicate_positivity_limiting(
  const std::vector<state_type> &states) const
{
  if constexpr ((mhd_flags & MHDFlags::no_positivity_limiting) !=
                MHDFlags::none)
    return false;

  const double eps = 1e-10;

  for (const state_type &state : states)
    {
      const double pressure = mhd_equations.compute_pressure_unsafe(state);

      if ((state[density_component] <= eps) ||
          (state[energy_component] <= eps) || (pressure <= eps))
        return true;
    }

  return false;
}



template <unsigned int dim>
void
sapphirepp::MHD::MHDSolver<dim>::apply_limiter()
{
  if constexpr ((mhd_flags & MHDFlags::no_limiting) != MHDFlags::none)
    return;

  TimerOutput::Scope t(timer, "Limiter - MHD");
  saplog << "Limit solution" << std::endl;
  LogStream::Prefix p("Limiter", saplog);

  using Iterator = typename DoFHandler<dim>::active_cell_iterator;
  using namespace sapphirepp::sapinternal::MHDSolver;


  /**
   * To limit the solution we use the following procedure:
   *
   *   1. Precompute cell averages using @ref compute_cell_average().
   *   2. Compute the limited solution/limited gradient using neighbor cells
   *      in @ref cell_worker.
   *   3. Test if positivity limiting is needed.
   *      We do this by testing if the current/limited solution results in
   *      non-admissible states on any point of some quadrature
   *      (e.g. generalized support points).
   *   4. Update DoFs on each cell using generalized support points
   *      @dealref{convert_generalized_support_point_values_to_dof_values(),classFiniteElement,a8f2b9817bcf98ebd43239c6da46de809}.
   *      This only works for basis functions that have support points,
   *      i.e. for LagrangePolynomials, but not for Legendre polynomials.
   *   5. Distribute to the @ref locally_relevant_current_solution.
   *
   * There are some competing alternative ways to perform step 4:
   *
   *   - Compute the projection of the limited solution onto the DoFs:
   *     1. Compute the RHS for the projection by folding with the basis
   *        functions solution in @ref cell_worker.
   *     2. Distribute the RHS to the global @ref system_rhs in @ref copier.
   *     3. Project the limited solution on the @ref locally_owned_solution
   *        \f$ f \f$ using the @ref mass_matrix \f$ M \f$ and the
   *        @ref system_rhs \f$ b \f$, \f$ M f = b \f$.
   *   - Locally invert the mass matrix. Disadvantage is, that we need to
   * solve a (small) linear system on each cell.
   *   - If one has basis functions with a diagonal mass matrix, the inversion
   *     is trivial.
   *
   * One general advantage of doing this cell wise is, that we only have to
   * update the solution in cells that are actually limited.
   */

  /** Step 1: Precompute cell averages and shock_indicator. */
  compute_cell_average();
  if constexpr ((mhd_flags & MHDFlags::no_shock_indicator) == MHDFlags::none)
    compute_shock_indicator();

  AssertDimension(cell_average.size(), triangulation.n_active_cells());
  AssertDimension(shock_indicator.size(), triangulation.n_active_cells());
  AssertDimension(positivity_limiter_indicator.size(),
                  triangulation.n_active_cells());
  Assert(fe.has_generalized_support_points(),
         ExcMessage("The slope limiter uses generalized support points "
                    "to compute the limited DoFs."));

  locally_owned_solution       = 0;
  positivity_limiter_indicator = 0.;

  // assemble cell terms
  const auto cell_worker = [&](const Iterator               &cell,
                               ScratchDataSlopeLimiter<dim> &scratch_data,
                               CopyDataSlopeLimiter         &copy_data) {
    FEValues<dim>     &fe_v_support = scratch_data.fe_values_support;
    const unsigned int cell_index   = cell->active_cell_index();

    // Check if cell avg is valid state
    const state_type &cell_avg = get_cell_average(cell);
    Assert(cell_avg[density_component] > 0.,
           ExcNonAdmissibleState<dim>(cell_avg,
                                      "Invalid cell average. "
                                      "Can not perform limiting.",
                                      cell->center()));
    Assert(cell_avg[energy_component] > 0.,
           ExcNonAdmissibleState<dim>(cell_avg,
                                      "Invalid cell average. "
                                      "Can not perform limiting.",
                                      cell->center()));
    Assert(mhd_equations.compute_pressure_unsafe(cell_avg) > 0.,
           ExcNonAdmissibleState<dim>(cell_avg,
                                      "Invalid cell average. "
                                      "Can not perform limiting.",
                                      cell->center()));

    // reinit copy_data, copies unlimited DoFs
    copy_data.reinit(cell,
                     fe.n_dofs_per_cell(),
                     locally_relevant_current_solution);

    bool                        limit_cell = false;
    flux_type                   limited_gradient;
    std::vector<Vector<double>> support_point_values(
      fe_v_support.n_quadrature_points, Vector<double>(n_components));


    /**
     * Step 2: Compute limited gradient.
     * If limiting is needed, this results in:
     *   @ref limit_cell = true
     *   @ref limited_cell_gradient = limited gradient
     * Otherwise this block has no effect.
     */
    // Only limit shock indicated cells
    if ((shock_indicator[cell_index] > 1.) ||
        ((mhd_flags & MHDFlags::no_shock_indicator) != MHDFlags::none))
      {
        // reinit cell
        FEValues<dim> &fe_v_grad = scratch_data.fe_values_gradient;
        fe_v_grad.reinit(cell);

        const std::vector<double> &JxW = fe_v_grad.get_JxW_values();

        double                 diff;
        flux_type              cell_avg_gradient;
        flux_type              tmp_gradient;
        std::vector<flux_type> neighbor_gradients;
        neighbor_gradients.reserve(cell->n_faces());

        // Compute cell average gradient
        std::vector<std::vector<Tensor<1, dim>>> solution_gradients(
          fe_v_grad.n_quadrature_points,
          std::vector<Tensor<1, dim>>(n_components));
        /** @todo Add local version of `get_function_gradients` to @dealii */
        // fe_v_grad.get_function_gradients(copy_data.cell_dof_values,
        //                                  copy_data.cell_indices,
        //                                  solution_gradients);
        fe_v_grad.get_function_gradients(locally_relevant_current_solution,
                                         solution_gradients);
        for (unsigned int c = 0; c < n_components; ++c)
          {
            for (const unsigned int q_index :
                 fe_v_grad.quadrature_point_indices())
              cell_avg_gradient[c] +=
                solution_gradients[q_index][c] * JxW[q_index];
            cell_avg_gradient[c] /= cell->measure();
          }

        // Compute gradients by neighbor cells (ignore non-periodic boundary)
        for (const auto face_no : cell->face_indices())
          {
            if (!cell->at_boundary(face_no) ||
                cell->has_periodic_neighbor(face_no))
              {
                auto neighbor = cell->neighbor_or_periodic_neighbor(face_no);
                const state_type    &neighbor_avg = get_cell_average(neighbor);
                const Tensor<1, dim> distance =
                  SlopeLimiter<dim, divergence_cleaning>::
                    compute_periodic_distance_cell_neighbor(cell, face_no);

                for (unsigned int c = 0; c < n_components; ++c)
                  for (unsigned int d = 0; d < dim; ++d)
                    tmp_gradient[c][d] =
                      (neighbor_avg[c] - cell_avg[c]) / distance[d];
                neighbor_gradients.push_back(tmp_gradient);
              }
          }

        if constexpr ((mhd_flags & MHDFlags::conserved_limiting) !=
                      MHDFlags::none)
          {
            // Primitive Limiting
            diff = SlopeLimiter<dim, divergence_cleaning>::minmod_gradients(
              cell_avg_gradient,
              neighbor_gradients,
              limited_gradient,
              cell->minimum_vertex_distance());
          }
        else
          {
            // Characteristic Limiting
            /** @todo Move these variables into scratch data */
            flux_type              char_cell_avg_gradient;
            flux_type              char_limited_gradient;
            std::vector<flux_type> char_neighbor_gradients(
              neighbor_gradients.size());

            std::array<dealii::FullMatrix<double>, dim> left_matrices;
            std::array<dealii::FullMatrix<double>, dim> right_matrices;
            for (unsigned int d = 0; d < dim; ++d)
              {
                left_matrices[d]  = FullMatrix<double>(n_components);
                right_matrices[d] = FullMatrix<double>(n_components);
              }

            mhd_equations.compute_transformation_matrices(cell_avg,
                                                          left_matrices,
                                                          right_matrices);

            // Convert to characteristic variables
            mhd_equations.convert_gradient_characteristic_to_conserved(
              cell_avg_gradient, left_matrices, char_cell_avg_gradient);
            for (unsigned int i = 0; i < neighbor_gradients.size(); ++i)
              mhd_equations.convert_gradient_characteristic_to_conserved(
                neighbor_gradients[i],
                left_matrices,
                char_neighbor_gradients[i]);

            // Computed limited gradient
            diff = SlopeLimiter<dim, divergence_cleaning>::minmod_gradients(
              char_cell_avg_gradient,
              char_neighbor_gradients,
              char_limited_gradient,
              cell->minimum_vertex_distance());

            // Convert back to conserved variables
            mhd_equations.convert_gradient_characteristic_to_conserved(
              char_limited_gradient, right_matrices, limited_gradient);
          }

        // Only limit if average gradient was changed
        limit_cell = (diff > 1e-10);

        // Enforce divergence free limited B-field
        if (limit_cell)
          {
            SlopeLimiter<dim, divergence_cleaning>::
              enforce_divergence_free_limited_gradient(limited_gradient);
          }
      }


    /**
     * Step 3: Test if positivity limiting is needed
     * If positivity limiting is needed, this results in:
     *   @ref limit_cell = true
     *   ? @ref limited_cell_gradient = 0
     *   @ref positivity_limiter_indicator[cell_index] = 1
     * Otherwise this block has always these effects:
     *   @ref fe_v_support.reinit(cell)
     *   ? @ref support_point_values = current/limited values
     */
    if constexpr ((mhd_flags & MHDFlags::no_positivity_limiting) ==
                  MHDFlags::none)
      {
        fe_v_support.reinit(cell);

        // Calculate values at quadrature points or use limited values
        if (limit_cell)
          {
            const std::vector<Point<dim>> &support_points =
              fe_v_support.get_quadrature_points();

            // Compute limited solution values at support points
            for (const unsigned int q_index :
                 fe_v_support.quadrature_point_indices())
              for (unsigned int c = 0; c < n_components; ++c)
                support_point_values[q_index][c] +=
                  cell_avg[c] + limited_gradient[c] *
                                  (support_points[q_index] - cell->center());
          }
        else
          {
            /** @todo Use local version of `get_function_values` */
            // fe_v_support.get_function_gradients(copy_data.cell_dof_values,
            //                                     copy_data.cell_indices,
            //                                     support_point_values);
            fe_v_support.get_function_values(locally_relevant_current_solution,
                                             support_point_values);
          }

        const bool limit_positivity =
          indicate_positivity_limiting(support_point_values);

        if (limit_positivity)
          {
            limit_cell                               = true;
            positivity_limiter_indicator[cell_index] = 1;

            /** @todo If possible, only do one of the following */
            for (unsigned int c = 0; c < n_components; ++c)
              limited_gradient[c] = 0.;
            // Set cell values to cell average
            for (auto &support_point_value : support_point_values)
              support_point_value = cell_avg;
          }
      }


    /**
     * Step 4: Update DoF values.
     * If limiting is needed, this results in:
     *   @ref copy_data.cell_dof_values = new (limited) DoF values
     * Otherwise this block has no effect.
     */
    if (limit_cell)
      {
        /**
         * @todo This assumes the FE has generalized support points,
         *       and that @ref support_point_value is already set.
         *       Furthermore, we don't use @ref limited_gradient.
         */
        fe.convert_generalized_support_point_values_to_dof_values(
          support_point_values, copy_data.cell_dof_values);
      }
  };


  /** Step 5: Distribute to the @ref locally_relevant_current_solution. */
  // copier for the mesh_loop function
  const auto copier = [&](const CopyDataSlopeLimiter &c) {
    constraints.distribute_local_to_global(c.cell_dof_values,
                                           c.local_dof_indices,
                                           locally_owned_solution);
  };


  ScratchDataSlopeLimiter<dim> scratch_data(mapping, fe, quadrature);
  CopyDataSlopeLimiter         copy_data;
  saplog << "Begin limiting of solution." << std::endl;
  const auto filtered_iterator_range =
    dof_handler.active_cell_iterators() | IteratorFilters::LocallyOwnedCell();
  MeshWorker::mesh_loop(filtered_iterator_range,
                        cell_worker,
                        copier,
                        scratch_data,
                        copy_data,
                        MeshWorker::assemble_own_cells);
  locally_owned_solution.compress(VectorOperation::add);
  locally_relevant_current_solution = locally_owned_solution;
  saplog << "The solution was limited." << std::endl;

  /** @todo Remove this line */
  compute_magnetic_divergence();
}



template <unsigned int dim>
void
sapphirepp::MHD::MHDSolver<dim>::compute_magnetic_divergence()
{
  TimerOutput::Scope timer_section(timer, "Magnetic divergence - MHD");
  saplog << "Compute magnetic divergence" << std::endl;
  LogStream::Prefix p("MagneticDivergence", saplog);

  using Iterator = typename DoFHandler<dim>::active_cell_iterator;
  using namespace sapphirepp::sapinternal::MHDSolver;

  magnetic_divergence = 0;

  // empty boundary worker
  std::function<void(
    const Iterator &, const unsigned int &, ScratchData<dim> &, CopyData &)>
    empty_boundary_worker;


  // assemble cell terms
  const auto cell_worker = [&](const Iterator   &cell,
                               ScratchData<dim> &scratch_data,
                               CopyData         &copy_data) {
    static_cast<void>(copy_data);

    FEValues<dim> &fe_v = scratch_data.fe_values;
    // reinit cell
    fe_v.reinit(cell);
    const unsigned int cell_index               = cell->active_cell_index();
    double             cell_magnetic_divergence = 0.;

    // const unsigned int n_dofs     = fe_v.get_fe().n_dofs_per_cell();
    const unsigned int n_q_points = fe_v.n_quadrature_points;

    const FEValuesExtractors::Vector magnetic_field(first_magnetic_component);
    // const std::vector<Point<dim>> &q_points =
    // fe_v.get_quadrature_points();
    const std::vector<double> &JxW = fe_v.get_JxW_values();
    std::vector<double>        divergence_values(n_q_points);

    fe_v[magnetic_field].get_function_divergences(
      locally_relevant_current_solution, divergence_values);

    for (const unsigned int q_index : fe_v.quadrature_point_indices())
      {
        // divB
        cell_magnetic_divergence +=
          std::fabs(divergence_values[q_index]) * JxW[q_index];
      }

    magnetic_divergence[cell_index] +=
      cell_magnetic_divergence / cell->measure();
    // saplog << "Magnetic divergence in cell " << cell_index << ": "
    //        << cell_magnetic_divergence / cell->measure() << std::endl;
  };


  // assemble interior face terms
  const auto face_worker = [&](const Iterator     &cell,
                               const unsigned int &face_no,
                               const unsigned int &subface_no,
                               const Iterator     &neighbor_cell,
                               const unsigned int &neighbor_face_no,
                               const unsigned int &neighbor_subface_no,
                               ScratchData<dim>   &scratch_data,
                               CopyData           &copy_data) {
    static_cast<void>(copy_data);

    Assert((cell->level() == neighbor_cell->level()) &&
             (neighbor_cell->has_children() == false),
           ExcNotImplemented("Mesh refinement not yet implemented."));
    static_cast<void>(subface_no);
    static_cast<void>(neighbor_subface_no);

    FEFaceValues<dim> &fe_v_face = scratch_data.fe_values_face;
    fe_v_face.reinit(cell, face_no);
    const unsigned int cell_index               = cell->active_cell_index();
    double             face_magnetic_divergence = 0.;
    double             face_norm                = 0.;


    FEFaceValues<dim> &fe_v_face_neighbor =
      scratch_data.fe_values_face_neighbor;
    fe_v_face_neighbor.reinit(neighbor_cell, neighbor_face_no);

    const FEValuesExtractors::Vector magnetic_field(first_magnetic_component);
    const unsigned int               n_q_points = fe_v_face.n_quadrature_points;
    const std::vector<double>       &JxW        = fe_v_face.get_JxW_values();
    std::vector<Tensor<1, dim>>      face_values(n_q_points);
    std::vector<Tensor<1, dim>>      face_values_neighbor(n_q_points);


    fe_v_face[magnetic_field].get_function_values(
      locally_relevant_current_solution, face_values);
    fe_v_face_neighbor[magnetic_field].get_function_values(
      locally_relevant_current_solution, face_values_neighbor);

    for (unsigned int q_index : fe_v_face.quadrature_point_indices())
      {
        // (B_nb - B_j) * n
        face_magnetic_divergence +=
          std::fabs((face_values_neighbor[q_index] - face_values[q_index]) *
                    fe_v_face.normal_vector(q_index)) *
          JxW[q_index];
        face_norm += JxW[q_index];
      }


    magnetic_divergence[cell_index] += face_magnetic_divergence / face_norm;
  };

  /** @todo Use valid or empty copier? */
  // copier for the mesh_loop function
  const auto copier = [&](const CopyData &c) {
    static_cast<void>(c);
    // Empty
  };
  static_cast<void>(copier);
  // empty copier
  std::function<void(const CopyData &)> empty_copier;

  ScratchData<dim> scratch_data(mapping,
                                fe,
                                quadrature,
                                quadrature_face,
                                update_gradients | update_JxW_values);
  CopyData         copy_data;
  saplog << "Begin the assembly of the magnetic divergence." << std::endl;
  const auto filtered_iterator_range =
    dof_handler.active_cell_iterators() | IteratorFilters::LocallyOwnedCell();
  MeshWorker::mesh_loop(filtered_iterator_range,
                        cell_worker,
                        empty_copier,
                        scratch_data,
                        copy_data,
                        MeshWorker::assemble_own_cells |
                          MeshWorker::assemble_own_interior_faces_once,
                        // MeshWorker::assemble_own_interior_faces_both,
                        // MeshWorker::assemble_ghost_faces_both |
                        empty_boundary_worker,
                        face_worker);
  saplog << "The magnetic divergence was assembled." << std::endl;

  const double global_divergence =
    Utilities::MPI::sum(magnetic_divergence.l1_norm(), mpi_communicator);
  saplog << "Global divergence of magnetic field: " << global_divergence
         << std::endl;
}



template <unsigned int dim>
void
sapphirepp::MHD::MHDSolver<dim>::assemble_dg_rhs(const double time)
{
  TimerOutput::Scope timer_section(timer, "DG rhs - MHD");

  using Iterator = typename DoFHandler<dim>::active_cell_iterator;
  using namespace sapphirepp::sapinternal::MHDSolver;

  static_cast<void>(time); // suppress compiler warning
  dg_rhs    = 0;
  cell_dt   = 0;
  global_dt = std::numeric_limits<double>::max();

  // assemble cell terms
  const auto cell_worker = [&](const Iterator   &cell,
                               ScratchData<dim> &scratch_data,
                               CopyData         &copy_data) {
    FEValues<dim> &fe_v = scratch_data.fe_values;
    // reinit cell
    fe_v.reinit(cell);

    const unsigned int n_dofs = fe_v.get_fe().n_dofs_per_cell();
    // reinit the cell_matrix
    copy_data.reinit(cell, n_dofs);

    const std::vector<Point<dim>> &q_points = fe_v.get_quadrature_points();
    const std::vector<double>     &JxW      = fe_v.get_JxW_values();

    std::vector<Vector<double>> states(q_points.size(),
                                       Vector<double>(n_components));
    flux_type                   flux_matrix;
    double                      max_eigenvalue = 0.;

    fe_v.get_function_values(locally_relevant_current_solution, states);

    for (const unsigned int q_index : fe_v.quadrature_point_indices())
      {
        mhd_equations.compute_flux_matrix(states[q_index], flux_matrix);
        max_eigenvalue =
          std::max(max_eigenvalue,
                   mhd_equations.compute_maximum_eigenvalue(states[q_index]));

        for (unsigned int i : fe_v.dof_indices())
          {
            // const unsigned int component_i =
            //   fe_v.get_fe().system_to_component_index(i).first;

            for (unsigned int c = 0; c < n_components; ++c)
              {
                // F[c] * \grad \phi[c]_i
                copy_data.cell_dg_rhs(i) +=
                  flux_matrix[c] * fe_v.shape_grad_component(i, q_index, c) *
                  JxW[q_index];

                /** @todo Add source term */
              }
          }
      }

    const double dx  = cell->minimum_vertex_distance();
    copy_data.min_dt = dx / ((2. * fe.degree + 1.) * max_eigenvalue);
  };


  // assemble boundary face terms
  const auto boundary_worker = [&](const Iterator     &cell,
                                   const unsigned int &face_no,
                                   ScratchData<dim>   &scratch_data,
                                   CopyData           &copy_data) {
    FEFaceValues<dim> &fe_v_face = scratch_data.fe_values_face;
    fe_v_face.reinit(cell, face_no);

    // NOTE: copy_data is not reinitialized, the cell_workers contribution
    // to the cell_dg_matrix should not be deleted

    const unsigned int boundary_id = cell->face(face_no)->boundary_id();
    const BoundaryConditionsMHD boundary_condition =
      mhd_parameters.boundary_conditions[boundary_id];

    const std::vector<Point<dim>> &q_points = fe_v_face.get_quadrature_points();
    const std::vector<double>     &JxW      = fe_v_face.get_JxW_values();
    const std::vector<Tensor<1, dim>> &normals = fe_v_face.get_normal_vectors();

    // Initialise state and flux vectors
    std::vector<Vector<double>> states(q_points.size(),
                                       Vector<double>(n_components));
    state_type                  virtual_neighbor_state(n_components);
    state_type                  numerical_normal_flux(n_components);

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

            for (unsigned int c = 0; c < n_components; ++c)
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
  const auto face_worker = [&](const Iterator     &cell,
                               const unsigned int &face_no,
                               const unsigned int &subface_no,
                               const Iterator     &neighbor_cell,
                               const unsigned int &neighbor_face_no,
                               const unsigned int &neighbor_subface_no,
                               ScratchData<dim>   &scratch_data,
                               CopyData           &copy_data) {
    // We assumes that faces are only  assembled once. And hence subface_no,
    // can be ignored.
    static_cast<void>(subface_no);
    // We are not using mesh refinement yet. Hence neighbor_subface_no can be
    // ignored.
    static_cast<void>(neighbor_subface_no);

    FEFaceValues<dim> &fe_v_face = scratch_data.fe_values_face;
    fe_v_face.reinit(cell, face_no);

    FEFaceValues<dim> &fe_v_face_neighbor =
      scratch_data.fe_values_face_neighbor;
    fe_v_face_neighbor.reinit(neighbor_cell, neighbor_face_no);

    // Create an element at the end of the vector containing the face data
    copy_data.face_data.emplace_back();
    CopyDataFace      &copy_data_face = copy_data.face_data.back();
    const unsigned int n_dofs         = fe_v_face.get_fe().n_dofs_per_cell();
    copy_data_face.reinit(cell, neighbor_cell, n_dofs);

    const std::vector<Point<dim>> &q_points = fe_v_face.get_quadrature_points();
    const std::vector<double>     &JxW      = fe_v_face.get_JxW_values();
    const std::vector<Tensor<1, dim>> &normals = fe_v_face.get_normal_vectors();

    // Initialise state and flux vectors
    std::vector<Vector<double>> states(q_points.size(),
                                       Vector<double>(n_components));
    std::vector<Vector<double>> states_neighbor(q_points.size(),
                                                Vector<double>(n_components));
    state_type                  numerical_normal_flux(n_components);

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

            for (unsigned int c = 0; c < n_components; ++c)
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

            for (unsigned int c = 0; c < n_components; ++c)
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
    cell_dt[c.cell_index] = c.min_dt;
    global_dt             = std::min(global_dt, c.min_dt);
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

  ScratchData<dim> scratch_data(mapping, fe, quadrature, quadrature_face);
  CopyData         copy_data;
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

  global_dt = Utilities::MPI::min(global_dt, mpi_communicator);
  saplog << "Maximum CFL time step: " << global_dt << std::endl;
  Assert(global_dt > 0, ExcMessage("The time step must be greater than 0."));
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
  AssertThrow(time_step <= global_dt,
              ExcMessage(
                "Violated CFL condition: " + std::to_string(time_step) + " > " +
                std::to_string(global_dt)));

  mass_matrix.vmult(system_rhs, locally_owned_solution);
  system_rhs.add(time_step, dg_rhs);

  SolverControl              solver_control(1000, 1e-6 * system_rhs.l2_norm());
  PETScWrappers::SolverGMRES solver(solver_control, mpi_communicator);

  PETScWrappers::PreconditionBlockJacobi preconditioner;
  preconditioner.initialize(mass_matrix);

  // The x vector has to be a vector of locally_owned_dofs, so we make use of
  // the locally_owned_solution
  solver.solve(mass_matrix, locally_owned_solution, system_rhs, preconditioner);

  // Update the solution
  locally_relevant_current_solution = locally_owned_solution;

  saplog << "Solver converged in " << solver_control.last_step()
         << " iterations." << std::endl;

  apply_limiter();
}



template <unsigned int dim>
void
sapphirepp::MHD::MHDSolver<dim>::explicit_runge_kutta(const double time,
                                                      const double time_step)
{
  TimerOutput::Scope timer_section(timer, "ERK - MHD");
  LogStream::Prefix  p("ERK", saplog);

  // Number of stages: erk2 uses 2 stages, erk 4 uses 5 stages
  const unsigned int s =
    mhd_parameters.time_stepping_method == TimeSteppingMethodMHD::erk4 ? 5 : 2;
  FullMatrix<double> alpha(s, s);
  FullMatrix<double> beta(s, s);
  Vector<double>     gamma(s);

  switch (mhd_parameters.time_stepping_method)
    {
      case TimeSteppingMethodMHD::erk2:
        {
          saplog << "Define Butcher's array for ERK2 method" << std::endl;
          alpha[0][0] = 1.;
          alpha[1][0] = 0.5;
          alpha[1][1] = 0.5;
          beta[0][0]  = 1.;
          beta[1][0]  = 0.;
          beta[1][1]  = 0.5;
          gamma[0]    = 0.;
          gamma[1]    = 1.;
          break;
        }
      case TimeSteppingMethodMHD::erk4:
        {
          saplog << "Define Butcher's array for ERK4 method" << std::endl;
          alpha[0][0] = 1.;
          alpha[1][0] = 0.44437049406734;
          alpha[1][1] = 0.55562950593266;
          alpha[2][0] = 0.62010185138540;
          alpha[2][1] = 0.;
          alpha[2][2] = 0.37989814861460;
          alpha[3][0] = 0.17807995410773;
          alpha[3][1] = 0.;
          alpha[3][2] = 0.0;
          alpha[3][3] = 0.82192004589227;
          alpha[4][0] = 0.00683325884039;
          alpha[4][1] = 0.;
          alpha[4][2] = 0.51723167208978;
          alpha[4][3] = 0.12759831133288;
          alpha[4][4] = 0.34833675773694;

          beta[0][0] = 0.39175222700392;
          beta[1][0] = 0.;
          beta[1][1] = 0.36841059262959;
          beta[2][0] = 0.;
          beta[2][1] = 0.;
          beta[2][2] = 0.25189177424738;
          beta[3][0] = 0.;
          beta[3][1] = 0.;
          beta[3][2] = 0.;
          beta[3][3] = 0.54497475021237;
          beta[4][0] = 0.;
          beta[4][1] = 0.;
          beta[4][2] = 0.0;
          beta[4][3] = 0.08460416338212;
          beta[4][4] = 0.22600748319395;

          gamma[0] = 0.;
          gamma[1] = 0.39175222700392;
          gamma[2] = 0.58607968896780;
          gamma[3] = 0.47454236302687;
          gamma[4] = 0.93501063100924;
          break;
        }
      case TimeSteppingMethodMHD::forward_euler:
      default:
        {
          AssertThrow(false, ExcNotImplemented());
        }
    }

  PETScWrappers::PreconditionBlockJacobi preconditioner;
  preconditioner.initialize(mass_matrix);
  SolverControl              solver_control(1000);
  PETScWrappers::SolverGMRES solver(solver_control, mpi_communicator);
  // PETScWrappers::SolverCG solver(solver_control, mpi_communicator);

  PETScWrappers::MPI::Vector temp(locally_owned_dofs, mpi_communicator);
  std::vector<PETScWrappers::MPI::Vector> locally_owned_staged_solution(s);
  std::vector<PETScWrappers::MPI::Vector> locally_owned_staged_dg_rhs(s);
  for (unsigned int i = 0; i < s; ++i)
    {
      locally_owned_staged_solution[i].reinit(locally_owned_dofs,
                                              mpi_communicator);
      locally_owned_staged_dg_rhs[i].reinit(locally_owned_dofs,
                                            mpi_communicator);
    }

  // Stage i=0
  locally_owned_staged_solution[0] = locally_relevant_current_solution;
  assemble_dg_rhs(time + gamma[0] * time_step);
  AssertThrow(time_step <= global_dt,
              ExcMessage(
                "Violated CFL condition: " + std::to_string(time_step) + " > " +
                std::to_string(global_dt)));
  locally_owned_staged_dg_rhs[0] = dg_rhs;

  // Stage 0<i<s
  for (unsigned int i = 1; i < s; ++i)
    {
      temp = locally_owned_staged_solution[0];
      temp *= alpha[i - 1][0];
      mass_matrix.vmult(system_rhs, temp);
      system_rhs.add(beta[i - 1][0] * time_step,
                     locally_owned_staged_dg_rhs[0]);
      for (unsigned int j = 1; j < i; ++j)
        {
          temp = locally_owned_staged_solution[j];
          temp *= alpha[i - 1][j];
          mass_matrix.vmult_add(system_rhs, temp);
          system_rhs.add(beta[i - 1][j] * time_step,
                         locally_owned_staged_dg_rhs[j]);
        }

      solver_control.set_tolerance(1e-6 * system_rhs.l2_norm());
      solver.solve(mass_matrix,
                   locally_owned_staged_solution[i],
                   system_rhs,
                   preconditioner);
      saplog << "Stage i: " << i << "	Solver converged in "
             << solver_control.last_step() << " iterations." << std::endl;

      locally_relevant_current_solution = locally_owned_staged_solution[i];
      apply_limiter();
      assemble_dg_rhs(time + gamma[i] * time_step);
      AssertThrow(time_step <= global_dt,
                  ExcMessage(
                    "Violated CFL condition: " + std::to_string(time_step) +
                    " > " + std::to_string(global_dt)));
      locally_owned_staged_dg_rhs[i] = dg_rhs;
    }

  // Stage  i=s
  temp = locally_owned_staged_solution[0];
  temp *= alpha[s - 1][0];
  mass_matrix.vmult(system_rhs, temp);
  system_rhs.add(beta[s - 1][0] * time_step, locally_owned_staged_dg_rhs[0]);
  for (unsigned int j = 1; j < s; ++j)
    {
      temp = locally_owned_staged_solution[j];
      temp *= alpha[s - 1][j];
      mass_matrix.vmult_add(system_rhs, temp);
      system_rhs.add(beta[s - 1][j] * time_step,
                     locally_owned_staged_dg_rhs[j]);
    }

  solver_control.set_tolerance(1e-6 * system_rhs.l2_norm());
  solver.solve(mass_matrix, locally_owned_solution, system_rhs, preconditioner);
  saplog << "Stage i: " << s << "	Solver converged in "
         << solver_control.last_step() << " iterations." << std::endl;

  // Update the solution
  locally_relevant_current_solution = locally_owned_solution;

  apply_limiter();
}



template <unsigned int dim>
void
sapphirepp::MHD::MHDSolver<dim>::output_results(
  const unsigned int time_step_number,
  const double       cur_time)
{
  TimerOutput::Scope timer_section(timer, "Output - MHD");
  DataOut<dim>       data_out;
  data_out.attach_dof_handler(dof_handler);
  data_out.add_data_vector(
    locally_relevant_current_solution,
    MHDEquations<dim, divergence_cleaning>::create_component_name_list(),
    DataOut<dim>::type_dof_data,
    MHDEquations<dim,
                 divergence_cleaning>::create_component_interpretation_list());

  /** @todo [Remove Debug] */
  data_out.add_data_vector(get_cell_average_component(density_component),
                           "average_roh",
                           DataOut<dim>::type_cell_data);
  data_out.add_data_vector(shock_indicator,
                           "shock_indicator",
                           DataOut<dim>::type_cell_data);
  data_out.add_data_vector(positivity_limiter_indicator,
                           "positivity_limiter",
                           DataOut<dim>::type_cell_data);
  data_out.add_data_vector(magnetic_divergence,
                           "magnetic_divergence",
                           DataOut<dim>::type_cell_data);
  data_out.add_data_vector(cell_dt, "cell_dt", DataOut<dim>::type_cell_data);
  /** @todo [Remove Debug] */

  // Output the partition of the mesh
  Vector<float> subdomain(triangulation.n_active_cells());
  for (unsigned int i = 0; i < subdomain.size(); ++i)
    subdomain(i) = static_cast<float>(triangulation.locally_owned_subdomain());
  data_out.add_data_vector(subdomain, "subdomain");

  // Adapt the output to the polynomial degree of the shape functions
  data_out.build_patches(mhd_parameters.polynomial_degree);
  output_parameters.write_results<dim>(data_out, time_step_number, cur_time);
}



template <unsigned int dim>
double
sapphirepp::MHD::MHDSolver<dim>::compute_global_error(
  const Function<dim>         &exact_solution,
  const VectorTools::NormType &cell_norm,
  const VectorTools::NormType &global_norm,
  const Function<dim, double> *weight) const
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
  const VectorTools::NormType &cell_norm,
  const VectorTools::NormType &global_norm,
  const Function<dim, double> *weight) const
{
  LogStream::Prefix p("MHD", saplog);
  saplog << "Compute the weighted norm" << std::endl;
  LogStream::Prefix p2("Norm", saplog);

  Vector<float> cell_norms(triangulation.n_locally_owned_active_cells());

  const QTrapezoid<1>  q_trapezoid;
  const QIterated<dim> q_iterated(q_trapezoid,
                                  mhd_parameters.polynomial_degree * 2 + 1);

  const Functions::ZeroFunction<dim> zero_function(mhd_equations.n_components);

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
const sapphirepp::MHD::
  MHDEquations<dim, sapphirepp::MHD::MHDSolver<dim>::divergence_cleaning> &
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
const dealii::DoFHandler<dim> &
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
const std::vector<dealii::Vector<double>> &
sapphirepp::MHD::MHDSolver<dim>::get_cell_average() const
{
  return cell_average;
}



template <unsigned int dim>
dealii::Vector<double>
sapphirepp::MHD::MHDSolver<dim>::get_cell_average_component(
  unsigned int component) const
{
  AssertIndexRange(component, n_components);

  Vector<double> cell_average_component(triangulation.n_active_cells());

  for (unsigned int i = 0; i < cell_average_component.size(); ++i)
    cell_average_component[i] = cell_average[i][component];

  return cell_average_component;
}



template <unsigned int dim>
const dealii::Vector<double> &
sapphirepp::MHD::MHDSolver<dim>::get_shock_indicator() const
{
  return shock_indicator;
}



template <unsigned int dim>
const dealii::Vector<float> &
sapphirepp::MHD::MHDSolver<dim>::get_positivity_limiter_indicator() const
{
  return positivity_limiter_indicator;
}



template <unsigned int dim>
const dealii::Vector<double> &
sapphirepp::MHD::MHDSolver<dim>::get_magnetic_divergence() const
{
  return magnetic_divergence;
}



template <unsigned int dim>
const dealii::Vector<double> &
sapphirepp::MHD::MHDSolver<dim>::get_cell_dt() const
{
  return cell_dt;
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
