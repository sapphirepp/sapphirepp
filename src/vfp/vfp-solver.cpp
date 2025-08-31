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
 * @file vfp-solver.cpp
 * @author Nils Schween (nils.schween@mpi-hd.mpg.de),
 * @author Florian Schulze (florian.schulze@mpi-hd.mpg.de)
 * @brief Implement @ref sapphirepp::VFP::VFPSolver
 */

#include "vfp-solver.h"

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

#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <numeric>
#include <string>
#include <vector>

#include "config.h"
#include "particle-functions.h"
#include "probe-location.h"
#include "sapphirepp-logstream.h"
#include "vfp-flags.h"



/** @cond sapinternal */
namespace sapinternal
{
  using namespace dealii;

  // The mesh_loop function requires helper data types
  template <unsigned int dim_ps>
  class ScratchData
  {
  public:
    // Constructor
    ScratchData(const Mapping<dim_ps>        &mapping,
                const FiniteElement<dim_ps>  &fe,
                const Quadrature<dim_ps>     &quadrature,
                const Quadrature<dim_ps - 1> &quadrature_face,
                const UpdateFlags             update_flags = update_values |
                                                 update_gradients |
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
    ScratchData(const ScratchData<dim_ps> &scratch_data)
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

    FEValues<dim_ps>        fe_values;
    FEFaceValues<dim_ps>    fe_values_face;
    FEFaceValues<dim_ps>    fe_values_face_neighbor;
    FESubfaceValues<dim_ps> fe_values_subface_neighbor;
  };



  struct CopyDataFace
  {
    FullMatrix<double> cell_dg_matrix_11;
    FullMatrix<double> cell_dg_matrix_12;
    FullMatrix<double> cell_dg_matrix_21;
    FullMatrix<double> cell_dg_matrix_22;

    std::vector<types::global_dof_index> local_dof_indices;
    std::vector<types::global_dof_index> local_dof_indices_neighbor;

    template <typename Iterator>
    void
    reinit(const Iterator &cell,
           const Iterator &neighbor_cell,
           unsigned int    dofs_per_cell)
    {
      cell_dg_matrix_11.reinit(dofs_per_cell, dofs_per_cell);
      cell_dg_matrix_12.reinit(dofs_per_cell, dofs_per_cell);
      cell_dg_matrix_21.reinit(dofs_per_cell, dofs_per_cell);
      cell_dg_matrix_22.reinit(dofs_per_cell, dofs_per_cell);

      local_dof_indices.resize(dofs_per_cell);
      cell->get_dof_indices(local_dof_indices);

      local_dof_indices_neighbor.resize(dofs_per_cell);
      neighbor_cell->get_dof_indices(local_dof_indices_neighbor);
    }
  };



  struct CopyData
  {
    FullMatrix<double>                   cell_matrix;
    Vector<double>                       cell_rhs;
    std::vector<types::global_dof_index> local_dof_indices;
    std::vector<types::global_dof_index> local_dof_indices_neighbor;
    std::vector<CopyDataFace>            face_data;

    template <typename Iterator>
    void
    reinit(const Iterator &cell, unsigned int dofs_per_cell)
    {
      cell_matrix.reinit(dofs_per_cell, dofs_per_cell);
      cell_rhs.reinit(dofs_per_cell);

      local_dof_indices.resize(dofs_per_cell);
      cell->get_dof_indices(local_dof_indices);
    }
  };
} // namespace sapinternal
/** @endcond */



template <unsigned int dim>
sapphirepp::VFP::VFPSolver<dim>::VFPSolver(
  const VFPParameters<dim_ps>   &vfp_parameters,
  const PhysicalParameters      &physical_parameters,
  const Utils::OutputParameters &output_parameters)
  : vfp_parameters{vfp_parameters}
  , physical_parameters{physical_parameters}
  , output_parameters{output_parameters}
  , pde_system(vfp_parameters.expansion_order)
  , upwind_flux(pde_system, vfp_parameters, physical_parameters)
  , mpi_communicator{MPI_COMM_WORLD}
  , triangulation(mpi_communicator)
  , dof_handler(triangulation)
  , mapping()
  , fe(FE_DGQ<dim_ps>(vfp_parameters.polynomial_degree), pde_system.system_size)
  , quadrature(fe.tensor_degree() + 1)
  , quadrature_face(fe.tensor_degree() + 1)
  , probe_location(vfp_parameters, output_parameters, pde_system.lms_indices)
  , scaling_spectral_index{3.}
  , reflective_bc_signature{{std::vector<double>(pde_system.system_size),
                             std::vector<double>(pde_system.system_size),
                             std::vector<double>(pde_system.system_size)}}
  , scattering_frequency(physical_parameters)
  , source_function(physical_parameters, pde_system.system_size)
  , magnetic_field(physical_parameters)
  , background_velocity_field(physical_parameters)
  , bc_value_function(physical_parameters, pde_system.system_size)
  , debug_input_functions(scattering_frequency,
                          source_function,
                          magnetic_field,
                          background_velocity_field)
  , debug_input_functions_dof_handler(triangulation)
  , pcout(saplog.to_condition_ostream(3))
  , timer(mpi_communicator, pcout, TimerOutput::never, TimerOutput::wall_times)
{
  LogStream::Prefix prefix_vfp("VFP", saplog);
  LogStream::Prefix prefix("Constructor", saplog);
  saplog << vfp_flags << std::endl;
  saplog << "dim_ps=" << dim_ps << ", dim_cs=" << dim_cs << std::endl;
  if constexpr ((vfp_flags & VFPFlags::scaled_distribution_function) !=
                VFPFlags::none)
    saplog << "Scaling spectral index s: " << scaling_spectral_index
           << std::endl;

  // Compute sign changes of expansion coefficients for reflecting boundaries
  std::vector<std::array<unsigned int, 3>> lms_indices =
    pde_system.create_lms_indices(pde_system.system_size);

  for (unsigned int i = 0; i < pde_system.system_size; ++i)
    {
      // lower/upper x (y-z plane)
      reflective_bc_signature[0][i] =
        ((lms_indices[i][0] - lms_indices[i][1]) % 2 == 0 ? 1 : -1);
      // lower/upper y (x-z plane)
      reflective_bc_signature[1][i] =
        ((lms_indices[i][1] + lms_indices[i][2]) % 2 == 0 ? 1 : -1);
      // lower/upper z (x-y plane)
      reflective_bc_signature[2][i] = lms_indices[i][2] % 2 == 0 ? 1 : -1;
    }

  // Consistency checks for vfp_flags:
  AssertThrow(
    1 <= dim and dim <= 3,
    ExcMessage(
      "The total dimension must be greater than or equal to one and smaller or"
      " equal to three."));
  if (((vfp_flags & VFPFlags::spatial_advection) == VFPFlags::none) &&
      (dim_cs != 0))
    {
      saplog.print_warning(
        "Spatial advection is deactivated, but dim_cs > 0. \n"
        "  If the spatial advection term is deactivated, \n"
        "  the distribution function is assumed to be homogeneous \n"
        "  i.e. the dimension of the configuration space should be zero.");
    }
  if ((vfp_flags & VFPFlags::time_independent_fields) != VFPFlags::none)
    {
      AssertThrow(((vfp_flags & VFPFlags::spatial_advection) !=
                   VFPFlags::none) ||
                    ((vfp_flags & VFPFlags::rotation) != VFPFlags::none),
                  ExcMessage("Either the spatial advection or the rotation "
                             "term must be activated if the fields are time "
                             "independent."));
    }
  if ((vfp_flags & VFPFlags::time_independent_source) != VFPFlags::none)
    {
      AssertThrow((vfp_flags & VFPFlags::source) != VFPFlags::none,
                  ExcMessage(
                    "The source term must be activated if the source is time "
                    "independent."));
    }
}



template <unsigned int dim>
void
sapphirepp::VFP::VFPSolver<dim>::setup()
{
  TimerOutput::Scope timer_section_setup(timer, "VFP - Setup");
  LogStream::Prefix  prefix_vfp("VFP", saplog);
  saplog << "Setup VFP equation solver. \t[" << Utilities::System::get_time()
         << "]" << std::endl;
  LogStream::Prefix prefix("Setup", saplog);
  make_grid();
  setup_system();

  {
    if constexpr ((VFPFlags::time_evolution & vfp_flags) != VFPFlags::none)
      {
        TimerOutput::Scope timer_section(timer, "VFP - Mass matrix");
        saplog << "Assemble mass matrix." << std::endl;
        MatrixCreator::create_mass_matrix(mapping,
                                          dof_handler,
                                          quadrature,
                                          mass_matrix);
      }
  }

  {
    if constexpr ((VFPFlags::time_evolution & vfp_flags) != VFPFlags::none)
      {
        TimerOutput::Scope           timer_section(timer, "VFP - Project IC");
        InitialValueFunction<dim_ps> initial_value_function(
          physical_parameters, pde_system.system_size);
        PETScWrappers::MPI::Vector initial_condition(locally_owned_dofs,
                                                     mpi_communicator);
        project(initial_value_function, initial_condition);
        // Here a non ghosted vector, is copied into a ghosted vector. I think
        // that is the moment where the ghost cells are filled.
        locally_relevant_current_solution = initial_condition;
      }
  }

  // Assemble the dg matrix for t = 0
  assemble_dg_matrix(0);
  // Source term at t = 0;
  if constexpr ((vfp_flags & VFPFlags::source) != VFPFlags::none)
    {
      TimerOutput::Scope timer_source(timer, "VFP - Source");
      source_function.set_time(0);
      VectorTools::create_right_hand_side(mapping,
                                          dof_handler,
                                          quadrature,
                                          source_function,
                                          locally_owned_current_source);
    }
}



template <unsigned int dim>
void
sapphirepp::VFP::VFPSolver<dim>::run()
{
  setup();
  LogStream::Prefix prefix_vfp("VFP", saplog);
  if constexpr ((vfp_flags & VFPFlags::time_evolution) == VFPFlags::none)
    {
      steady_state_solve();
      output_results(0, 0);
      saplog << "Simulation ended. " << " \t\t["
             << Utilities::System::get_time() << "]" << std::endl;
    }
  else
    {
      DiscreteTime discrete_time(0,
                                 vfp_parameters.final_time,
                                 vfp_parameters.time_step);
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

          switch (vfp_parameters.time_stepping_method)
            {
              case TimeSteppingMethod::forward_euler:
              case TimeSteppingMethod::backward_euler:
              case TimeSteppingMethod::crank_nicolson:
                theta_method(discrete_time.get_current_time(),
                             discrete_time.get_next_step_size());
                break;
              case TimeSteppingMethod::erk4:
                explicit_runge_kutta(discrete_time.get_current_time(),
                                     discrete_time.get_next_step_size());
                break;
              case TimeSteppingMethod::lserk4:
                low_storage_explicit_runge_kutta(
                  discrete_time.get_current_time(),
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
    }

  {
    LogStream::Prefix              prefix("Summary", saplog);
    Utilities::System::MemoryStats memory_stats;
    Utilities::System::get_memory_stats(memory_stats);
    saplog << "Peak (local) resident memory size (HWM):    \t" //
           << memory_stats.VmHWM << " KiB"                     //
           << " \t= " << (memory_stats.VmHWM >> 10) << " MiB " //
           << " \t= " << (memory_stats.VmHWM >> 20) << " GiB " //
           << std::endl;
    saplog << "Peak (local) virtual/available memory size: \t"  //
           << memory_stats.VmPeak << " KiB"                     //
           << " \t= " << (memory_stats.VmPeak >> 10) << " MiB " //
           << " \t= " << (memory_stats.VmPeak >> 20) << " GiB " //
           << std::endl;

    timer.print_wall_time_statistics(mpi_communicator);
  }
}



template <unsigned int dim>
void
sapphirepp::VFP::VFPSolver<dim>::make_grid()
{
  TimerOutput::Scope timer_section(timer, "VFP - Grid setup");
  saplog << "Create the grid" << std::endl;

  switch (vfp_parameters.grid_type)
    {
      case GridType::hypercube:
        {
          saplog << "Create the grid from hyper rectangle" << std::endl;
          GridGenerator::subdivided_hyper_rectangle(triangulation,
                                                    vfp_parameters.n_cells,
                                                    vfp_parameters.p1,
                                                    vfp_parameters.p2,
                                                    true);
          break;
        }
      case GridType::shock:
        {
          make_shock_grid();
          break;
        }

      case GridType::file:
        {
          saplog << "Read grid from file \"" << vfp_parameters.grid_file << "\""
                 << std::endl;
          GridIn<dim_ps> grid_in(triangulation);
          std::ifstream  input(vfp_parameters.grid_file);
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
  for (unsigned int i = 0; i < dim_ps; ++i)
    {
      if (vfp_parameters.boundary_conditions[2 * i] ==
            BoundaryConditions::periodic or
          vfp_parameters.boundary_conditions[2 * i + 1] ==
            BoundaryConditions::periodic)
        {
          AssertThrow(vfp_parameters.boundary_conditions[2 * i] ==
                          BoundaryConditions::periodic and
                        vfp_parameters.boundary_conditions[2 * i + 1] ==
                          BoundaryConditions::periodic,
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
sapphirepp::VFP::VFPSolver<dim>::make_shock_grid()
{
  saplog << "Create shock grid" << std::endl;
  LogStream::Prefix("GridGeneration");
  Assert(dim_cs >= 1,
         ExcMessage("The shock grid is only implemented in x-direction, "
                    "but your don't have a x-coordinate."));
  AssertThrow(vfp_parameters.p1[0] < 0.,
              ExcMessage("x_min must be smaller than zero for the "
                         "shock grid."));
  AssertThrow(vfp_parameters.p2[0] > 0.,
              ExcMessage("x_max must be greater than zero for the "
                         "shock grid."));

  unsigned int       n_cells_upstream   = 0; // Negative x values
  unsigned int       n_cells_downstream = 0; // Positive x values
  const unsigned int n_cells_shock      = vfp_parameters.n_cells_shock;
  const double       shock_width        = vfp_parameters.shock_width;
  const double       scaling_factor     = vfp_parameters.scaling_factor_shock;

  const double step_size_shock = shock_width / n_cells_shock;
  saplog << "The shock has " << n_cells_shock << " cells with width "
         << step_size_shock << std::endl;

  std::vector<double> step_size_upstream;
  double              length_upstream = shock_width / 2.;
  for (unsigned int i = 1; length_upstream < -vfp_parameters.p1[0];
       ++i, ++n_cells_upstream)
    {
      const double step_size = std::pow(scaling_factor, i) * step_size_shock;
      step_size_upstream.push_back(step_size);
      length_upstream += step_size;
    }
  saplog << "Upstream has " << n_cells_upstream
         << " cells reaching x_min=" << std::to_string(-length_upstream)
         << " <= " << vfp_parameters.p1[0] << std::endl;

  std::vector<double> step_size_downstream;
  double              length_downstream = shock_width / 2.;
  for (unsigned int i = 1; length_downstream < vfp_parameters.p2[0];
       ++i, ++n_cells_downstream)
    {
      const double step_size = std::pow(scaling_factor, i) * step_size_shock;
      step_size_downstream.push_back(step_size);
      length_downstream += step_size;
    }
  saplog << "Downstream has " << n_cells_downstream
         << " cells reaching x_max=" << std::to_string(length_downstream)
         << " >= " << vfp_parameters.p2[0] << std::endl;

  std::vector<std::vector<double>> step_sizes;

  std::vector<double> step_sizes_x(n_cells_downstream + n_cells_upstream +
                                   n_cells_shock);
  for (unsigned int i = 0; i < n_cells_upstream; ++i)
    step_sizes_x[n_cells_upstream - 1 - i] = step_size_upstream[i];
  for (unsigned int i = 0; i < n_cells_shock; ++i)
    step_sizes_x[n_cells_upstream + i] = step_size_shock;
  for (unsigned int i = 0; i < n_cells_downstream; ++i)
    step_sizes_x[n_cells_upstream + n_cells_shock + i] =
      step_size_downstream[i];
  step_sizes.push_back(step_sizes_x);

  for (unsigned int i = 1; i < dim; ++i)
    {
      const double step_size = (vfp_parameters.p2[i] - vfp_parameters.p1[i]) /
                               vfp_parameters.n_cells[i];
      step_sizes.push_back(
        std::vector<double>(vfp_parameters.n_cells[i], step_size));
    }

  Point<dim> p1 = vfp_parameters.p1;
  Point<dim> p2 = vfp_parameters.p2;

  p1[0] = -length_upstream;
  p2[0] = +length_downstream;

  GridGenerator::subdivided_hyper_rectangle(
    triangulation, step_sizes, p1, p2, true);

  saplog << "Created shock grid from (" << p1 << ") to (" << p2 << ") with "
         << step_sizes[0].size() << " cells in x-direction." << std::endl;

  if (step_sizes[0].size() != vfp_parameters.n_cells[0])
    saplog.print_warning("The shock grid uses " +
                         Utilities::to_string(step_sizes[0].size()) +
                         " cells in x-direction. \n"
                         "The parameter file entry 'Number of cells = " +
                         Utilities::to_string(vfp_parameters.n_cells[0]) +
                         ", ...' is ignored for x.");
}



template <unsigned int dim>
void
sapphirepp::VFP::VFPSolver<dim>::setup_system()
{
  TimerOutput::Scope timer_section(timer, "VFP - Setup FE system");
  saplog << "Setup the finite element system" << std::endl;

  dof_handler.clear();
  dof_handler.distribute_dofs(fe);

  if (output_parameters.debug_input_functions)
    {
      const FESystem<dim_ps> debug_input_functions_fe(
        FE_DGQ<dim_ps>(vfp_parameters.polynomial_degree),
        debug_input_functions.n_components);

      debug_input_functions_dof_handler.clear();
      debug_input_functions_dof_handler.distribute_dofs(
        debug_input_functions_fe);
    }

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
  system_rhs.reinit(locally_owned_dofs, mpi_communicator);

  if constexpr ((vfp_flags & VFPFlags::source) != VFPFlags::none)
    locally_owned_current_source.reinit(locally_owned_dofs, mpi_communicator);
  // TODO: add if to check if there are non-homogeneous bc
  locally_owned_current_bc.reinit(locally_owned_dofs, mpi_communicator);

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
  dg_matrix.reinit(locally_owned_dofs,
                   locally_owned_dofs,
                   dsp,
                   mpi_communicator);
  // NOTE: DealII does not allow to use different sparsity patterns for
  // matrices, which you would like to add. Even though the the mass matrix
  // differs from the dg matrix.

  // For the steady state solver it is not necessary to allocate memory for the
  // mass_matrix and the system_matrix.
  if constexpr ((vfp_flags & VFPFlags::time_evolution) != VFPFlags::none)
    {
      mass_matrix.reinit(locally_owned_dofs,
                         locally_owned_dofs,
                         dsp,
                         mpi_communicator);
      system_matrix.reinit(locally_owned_dofs,
                           locally_owned_dofs,
                           dsp,
                           mpi_communicator);
    }
  probe_location.reinit(triangulation, mapping);

  Utilities::System::MemoryStats memory_stats;
  Utilities::System::get_memory_stats(memory_stats);
  saplog << "Current (local) memory usage: \t"               //
         << memory_stats.VmRSS << " KiB"                     //
         << " \t= " << (memory_stats.VmRSS >> 10) << " MiB " //
         << " \t= " << (memory_stats.VmRSS >> 20) << " GiB " //
         << std::endl;
}



template <unsigned int dim>
void
sapphirepp::VFP::VFPSolver<dim>::project(
  const Function<dim>        &f,
  PETScWrappers::MPI::Vector &projected_function) const
{
  saplog << "Project a function onto the finite element space" << std::endl;
  LogStream::Prefix prefix("project", saplog);

  // Create right hand side
  PETScWrappers::MPI::Vector rhs(locally_owned_dofs, mpi_communicator);
  VectorTools::create_right_hand_side(mapping, dof_handler, quadrature, f, rhs);

  // Solve the system
  PETScWrappers::PreconditionNone preconditioner;
  preconditioner.initialize(mass_matrix);

  SolverControl           solver_control(vfp_parameters.solver_max_iter,
                               vfp_parameters.solver_tolerance);
  PETScWrappers::SolverCG cg(solver_control);
  cg.solve(mass_matrix, projected_function, rhs, preconditioner);
  saplog << "Solved in " << solver_control.last_step() << " iterations."
         << std::endl;
  // At the moment I am assuming, that I do not have constraints. Hence, I
  // do not need the following line.
  // constraints.distribute(projected_function);
}



template <unsigned int dim>
void
sapphirepp::VFP::VFPSolver<dim>::assemble_dg_matrix(const double time)
{
  TimerOutput::Scope timer_section(timer, "VFP - DG matrix");
  /*
    What kind of loops are there ?
    1. Loop over all cells (this happens inside the mesh_loop)
    2. Loop over the degrees of freedom on each cell
    - the method system_to_component_index() returns the index of the non-zero
    component of the vector-valued shape function which corresponds to the
    indices (l,m,s)
  */
  using Iterator = typename DoFHandler<dim_ps>::active_cell_iterator;
  using namespace sapinternal;
  upwind_flux.set_time(time);

  const std::vector<LAPACKFullMatrix<double>> &advection_matrices =
    pde_system.get_advection_matrices();
  const std::vector<LAPACKFullMatrix<double>> &generator_rotation_matrices =
    pde_system.get_generator_rotation_matrices();
  const std::vector<LAPACKFullMatrix<double>> &adv_mat_products =
    pde_system.get_adv_mat_products();
  const std::vector<LAPACKFullMatrix<double>> &adv_x_gen_matrices =
    pde_system.get_adv_cross_gen();
  const std::vector<LAPACKFullMatrix<double>> &t_matrices =
    pde_system.get_t_matrices();

  // Collision term (essentially a reaction term)
  const Vector<double> &collision_matrix = pde_system.get_collision_matrix();

  scattering_frequency.set_time(time);
  magnetic_field.set_time(time);
  background_velocity_field.set_time(time);

  bc_value_function.set_time(time);

  ParticleVelocity<dim_ps, logarithmic_p> particle_velocity(
    vfp_parameters.mass);
  ParticleGamma<dim_ps, logarithmic_p> particle_gamma(vfp_parameters.mass);

  // I do not no the meaning of the following "const" specifier
  const auto cell_worker = [&](const Iterator      &cell,
                               ScratchData<dim_ps> &scratch_data,
                               CopyData            &copy_data) {
    FEValues<dim_ps> &fe_v = scratch_data.fe_values;
    // reinit cell
    fe_v.reinit(cell);

    const unsigned int n_dofs = fe_v.get_fe().n_dofs_per_cell();
    // reinit the cell_matrix
    copy_data.reinit(cell, n_dofs);

    const std::vector<Point<dim_ps>> &q_points = fe_v.get_quadrature_points();
    const std::vector<double>        &JxW      = fe_v.get_JxW_values();
    // Background Plasma
    // NOTE: Second argument constructs empty vectors with 2 values. They
    // are copied q_points.size() times.
    std::vector<Vector<double>> velocities(q_points.size(), Vector<double>(3));
    background_velocity_field.vector_value_list(q_points, velocities);
    std::vector<double> div_velocities(q_points.size());
    background_velocity_field.divergence_list(q_points, div_velocities);
    std::vector<Vector<double>> material_derivative_vel(q_points.size(),
                                                        Vector<double>(3));
    background_velocity_field.material_derivative_list(q_points,
                                                       material_derivative_vel);

    std::vector<FullMatrix<double>> jacobians_vel(q_points.size(),
                                                  FullMatrix<double>(3));
    background_velocity_field.jacobian_list(q_points, jacobians_vel);

    // Magnetic field
    std::vector<Vector<double>> magnetic_field_values(q_points.size(),
                                                      Vector<double>(3));
    magnetic_field.vector_value_list(q_points, magnetic_field_values);

    // Scattering frequency
    std::vector<double> frequencies(q_points.size());
    scattering_frequency.value_list(q_points, frequencies);

    // Particle
    std::vector<double> particle_velocities(q_points.size());
    particle_velocity.value_list(q_points, particle_velocities);

    std::vector<double> particle_gammas(q_points.size());
    particle_gamma.value_list(q_points, particle_gammas);
    for (const unsigned int q_index : fe_v.quadrature_point_indices())
      {
        for (unsigned int i : fe_v.dof_indices())
          {
            const unsigned int component_i =
              fe_v.get_fe().system_to_component_index(i).first;
            for (unsigned int j : fe_v.dof_indices())
              {
                const unsigned int component_j =
                  fe_v.get_fe().system_to_component_index(j).first;
                if constexpr ((vfp_flags & VFPFlags::collision) !=
                              VFPFlags::none)
                  {
                    if (component_i == component_j)
                      {
                        // scattering_frequency * l(l+1) * \phi_i * \phi_j
                        copy_data.cell_matrix(i, j) +=
                          frequencies[q_index] * collision_matrix[component_i] *
                          fe_v.shape_value(i, q_index) *
                          fe_v.shape_value(j, q_index) * JxW[q_index];
                      }
                  }
                // NOTE: If spatial advection term is deactivated, then
                // dim_cs must equal zero, i.e. the distribution function is
                // assumed to be homogeneous (it does not depend on x, y
                // and z). And if dim_cs = 0, there for loop is not entered.
                for (unsigned int coordinate = 0; coordinate < dim_cs;
                     ++coordinate)
                  {
                    // -[partial_x \phi_i * (u_x \delta_ij + Ax_ij)
                    //   + partial_y \phi_i * (u_y \delta_ij + Ay_ij)] *
                    //   phi_j
                    if (component_i == component_j)
                      {
                        copy_data.cell_matrix(i, j) -=
                          fe_v.shape_grad(i, q_index)[coordinate] *
                          velocities[q_index][coordinate] *
                          fe_v.shape_value(j, q_index) * JxW[q_index];

                        // - [\partial_x(u_x\delta_ij + Ax_ij) +
                        // \partial_y(u_y\delta_ij +
                        // - Ay_ij) ] \phi_i \phi_j where \partial_x/y
                        // Ax/y_ij = 0
                        copy_data.cell_matrix(i, j) -=
                          div_velocities[q_index] *
                          fe_v.shape_value(i, q_index) *
                          fe_v.shape_value(j, q_index) * JxW[q_index];
                      }
                    else
                      {
                        // NOTE: Many zeros are added here, because the
                        // matrices Ax, Ay and A_z are sparse. TODO:
                        // Performance check. If too bad, return to the
                        // strategy, which was used in v0.6.5
                        if ((vfp_flags & VFPFlags::momentum) != VFPFlags::none)
                          {
                            copy_data.cell_matrix(i, j) -=
                              fe_v.shape_grad(i, q_index)[coordinate] *
                              particle_velocities[q_index] *
                              advection_matrices[coordinate](component_i,
                                                             component_j) *
                              fe_v.shape_value(j, q_index) * JxW[q_index];
                          }
                        else
                          {
                            // fixed energy case (i.e. transport only)
                            copy_data.cell_matrix(i, j) -=
                              fe_v.shape_grad(i, q_index)[coordinate] *
                              vfp_parameters.velocity *
                              advection_matrices[coordinate](component_i,
                                                             component_j) *
                              fe_v.shape_value(j, q_index) * JxW[q_index];
                          }
                      }
                  }
                if constexpr ((vfp_flags & VFPFlags::rotation) !=
                              VFPFlags::none)
                  {
                    // NOTE: All three components of the B-Field are
                    // included no matter, which dimension of the
                    // configuration space is considered
                    for (unsigned int coordinate = 0; coordinate < 3;
                         ++coordinate)
                      {
                        if ((vfp_flags & VFPFlags::momentum) != VFPFlags::none)
                          {
                            copy_data.cell_matrix(i, j) -=
                              fe_v.shape_value(i, q_index) *
                              vfp_parameters.charge *
                              magnetic_field_values[q_index][coordinate] /
                              (particle_gammas[q_index] * vfp_parameters.mass) *
                              generator_rotation_matrices[coordinate](
                                component_i, component_j) *
                              fe_v.shape_value(j, q_index) * JxW[q_index];
                          }
                        else
                          {
                            // fixed energy case (i.e. transport only)
                            copy_data.cell_matrix(i, j) -=
                              fe_v.shape_value(i, q_index) *
                              vfp_parameters.charge *
                              magnetic_field_values[q_index][coordinate] /
                              (vfp_parameters.gamma * vfp_parameters.mass) *
                              generator_rotation_matrices[coordinate](
                                component_i, component_j) *
                              fe_v.shape_value(j, q_index) * JxW[q_index];
                          }
                      }
                  }
                if constexpr ((vfp_flags & VFPFlags::momentum) !=
                              VFPFlags::none)
                  {
                    if constexpr (logarithmic_p)
                      {
                        // Momentum part
                        for (unsigned int coordinate = 0; coordinate < 3;
                             ++coordinate)
                          {
                            // \grad_phi * 1/v du^k/ dt * A_k \phi
                            copy_data.cell_matrix(i, j) +=
                              fe_v.shape_grad(i, q_index)[dim_ps - 1] /
                              particle_velocities[q_index] *
                              material_derivative_vel[q_index][coordinate] *
                              advection_matrices[coordinate](component_i,
                                                             component_j) *
                              fe_v.shape_value(j, q_index) * JxW[q_index];
                            // -\phi m/(gamma * p) * du^k/dt * A_k * \phi
                            copy_data.cell_matrix(i, j) -=
                              fe_v.shape_value(i, q_index) *
                              vfp_parameters.mass /
                              (particle_gammas[q_index] *
                               std::exp(q_points[q_index][dim_ps - 1])) *
                              material_derivative_vel[q_index][coordinate] *
                              advection_matrices[coordinate](component_i,
                                                             component_j) *
                              fe_v.shape_value(j, q_index) * JxW[q_index];
                            // \phi 1/v * du^k \ dt (A x \Omega)_k \phi
                            copy_data.cell_matrix(i, j) +=
                              fe_v.shape_value(i, q_index) /
                              particle_velocities[q_index] *
                              material_derivative_vel[q_index][coordinate] *
                              adv_x_gen_matrices[coordinate](component_i,
                                                             component_j) *
                              fe_v.shape_value(j, q_index) * JxW[q_index];

                            if constexpr ((vfp_flags &
                                           VFPFlags::
                                             scaled_distribution_function) !=
                                          VFPFlags::none)
                              {
                                // \phi scaling_spectral_index 1/v * du^k/dt *
                                // A_k * \phi
                                copy_data.cell_matrix(i, j) +=
                                  scaling_spectral_index *
                                  fe_v.shape_value(i, q_index) /
                                  particle_velocities[q_index] *
                                  material_derivative_vel[q_index][coordinate] *
                                  advection_matrices[coordinate](component_i,
                                                                 component_j) *
                                  fe_v.shape_value(j, q_index) * JxW[q_index];
                              }
                          }
                        for (unsigned int coordinate_1 = 0; coordinate_1 < 3;
                             ++coordinate_1)
                          {
                            for (unsigned int coordinate_2 = coordinate_1;
                                 coordinate_2 < 3;
                                 ++coordinate_2)
                              {
                                if (coordinate_1 == coordinate_2)
                                  {
                                    // \grad_phi
                                    // \jacobian[coordinate_1][coordinate_2]
                                    // Ap_coordinate_1,coordinate_2 * \phi
                                    copy_data.cell_matrix(i, j) +=
                                      fe_v.shape_grad(i, q_index)[dim_ps - 1] *
                                      jacobians_vel[q_index][coordinate_1]
                                                   [coordinate_2] *
                                      adv_mat_products
                                        [3 * coordinate_1 -
                                         coordinate_1 * (coordinate_1 + 1) / 2 +
                                         coordinate_2](component_i,
                                                       component_j) *
                                      fe_v.shape_value(j, q_index) *
                                      JxW[q_index];
                                    if constexpr (
                                      (vfp_flags &
                                       VFPFlags::
                                         scaled_distribution_function) !=
                                      VFPFlags::none)
                                      {
                                        // \phi * scaling_spectral_index *
                                        // \jacobian[coordinate_1][coordinate_2]
                                        // Ap_coordinate_1,coordinate_2 * \phi
                                        copy_data.cell_matrix(i, j) +=
                                          fe_v.shape_value(i, q_index) *
                                          scaling_spectral_index *
                                          jacobians_vel[q_index][coordinate_1]
                                                       [coordinate_2] *
                                          adv_mat_products
                                            [3 * coordinate_1 -
                                             coordinate_1 * (coordinate_1 + 1) /
                                               2 +
                                             coordinate_2](component_i,
                                                           component_j) *
                                          fe_v.shape_value(j, q_index) *
                                          JxW[q_index];
                                      }
                                  }
                                else
                                  {
                                    // symmetry
                                    // component_1, component_2
                                    // \grad_phi
                                    // \jacobian[coordinate_1][coordinate_2]
                                    // Ap_coordinate_1,coordinate_2 * \phi
                                    copy_data.cell_matrix(i, j) +=
                                      fe_v.shape_grad(i, q_index)[dim_ps - 1] *
                                      jacobians_vel[q_index][coordinate_1]
                                                   [coordinate_2] *
                                      adv_mat_products
                                        [3 * coordinate_1 -
                                         coordinate_1 * (coordinate_1 + 1) / 2 +
                                         coordinate_2](component_i,
                                                       component_j) *
                                      fe_v.shape_value(j, q_index) *
                                      JxW[q_index];

                                    // component_2, component_1
                                    // \grad_phi p
                                    // \jacobian[coordinate_1][coordinate_2]
                                    // Ap_coordinate_1,coordinate_2 * \phi
                                    copy_data.cell_matrix(i, j) +=
                                      fe_v.shape_grad(i, q_index)[dim_ps - 1] *
                                      jacobians_vel[q_index][coordinate_2]
                                                   [coordinate_1] *
                                      adv_mat_products
                                        [3 * coordinate_1 -
                                         coordinate_1 * (coordinate_1 + 1) / 2 +
                                         coordinate_2](component_i,
                                                       component_j) *
                                      fe_v.shape_value(j, q_index) *
                                      JxW[q_index];

                                    if constexpr (
                                      (vfp_flags &
                                       VFPFlags::
                                         scaled_distribution_function) !=
                                      VFPFlags::none)
                                      {
                                        // component_1, component_2
                                        // \phi * scaling_spectral_index *
                                        // \jacobian[coordinate_1][coordinate_2]
                                        // Ap_coordinate_1,coordinate_2 * \phi
                                        copy_data.cell_matrix(i, j) +=
                                          fe_v.shape_value(i, q_index) *
                                          scaling_spectral_index *
                                          jacobians_vel[q_index][coordinate_1]
                                                       [coordinate_2] *
                                          adv_mat_products
                                            [3 * coordinate_1 -
                                             coordinate_1 * (coordinate_1 + 1) /
                                               2 +
                                             coordinate_2](component_i,
                                                           component_j) *
                                          fe_v.shape_value(j, q_index) *
                                          JxW[q_index];

                                        // component_2, component_1
                                        // \phi * scaling_spectral_index
                                        // \jacobian[coordinate_1][coordinate_2]
                                        // Ap_coordinate_1,coordinate_2 * \phi
                                        copy_data.cell_matrix(i, j) +=
                                          fe_v.shape_value(i, q_index) *
                                          scaling_spectral_index *
                                          jacobians_vel[q_index][coordinate_2]
                                                       [coordinate_1] *
                                          adv_mat_products
                                            [3 * coordinate_1 -
                                             coordinate_1 * (coordinate_1 + 1) /
                                               2 +
                                             coordinate_2](component_i,
                                                           component_j) *
                                          fe_v.shape_value(j, q_index) *
                                          JxW[q_index];
                                      }
                                  }
                              }
                          }
                        for (unsigned int coordinate_1 = 0; coordinate_1 < 3;
                             ++coordinate_1)
                          {
                            for (unsigned int coordinate_2 = 0;
                                 coordinate_2 < 3;
                                 ++coordinate_2)
                              {
                                // \phi *
                                // jacobian[coordinate_1][coordinate_2]
                                // * T_coordinate_2,coordinate_1 * \phi
                                copy_data.cell_matrix(i, j) +=
                                  fe_v.shape_value(i, q_index) *
                                  jacobians_vel[q_index][coordinate_1]
                                               [coordinate_2] *
                                  t_matrices[coordinate_2 * 3 + coordinate_1](
                                    component_i, component_j) *
                                  fe_v.shape_value(j, q_index) * JxW[q_index];
                              }
                          }
                      }
                    else
                      {
                        // Momentum part
                        for (unsigned int coordinate = 0; coordinate < 3;
                             ++coordinate)
                          {
                            // \grad_phi * \gamma * mass * du^k/ dt * A_k
                            // \phi
                            copy_data.cell_matrix(i, j) +=
                              fe_v.shape_grad(i, q_index)[dim_ps - 1] *
                              particle_gammas[q_index] * vfp_parameters.mass *
                              material_derivative_vel[q_index][coordinate] *
                              advection_matrices[coordinate](component_i,
                                                             component_j) *
                              fe_v.shape_value(j, q_index) * JxW[q_index];
                            // \phi v * du^k/dt * A_k * \phi
                            copy_data.cell_matrix(i, j) +=
                              fe_v.shape_value(i, q_index) *
                              particle_velocities[q_index] *
                              material_derivative_vel[q_index][coordinate] *
                              advection_matrices[coordinate](component_i,
                                                             component_j) *
                              fe_v.shape_value(j, q_index) * JxW[q_index];
                            // \phi 1/v * du^k \ dt (A x \Omega)_k \phi
                            copy_data.cell_matrix(i, j) +=
                              fe_v.shape_value(i, q_index) /
                              particle_velocities[q_index] *
                              material_derivative_vel[q_index][coordinate] *
                              adv_x_gen_matrices[coordinate](component_i,
                                                             component_j) *
                              fe_v.shape_value(j, q_index) * JxW[q_index];

                            if constexpr ((vfp_flags &
                                           VFPFlags::
                                             scaled_distribution_function) !=
                                          VFPFlags::none)
                              {
                                // \phi scaling_spectral_index 1/v * du^k/dt *
                                // A_k * \phi
                                copy_data.cell_matrix(i, j) +=
                                  scaling_spectral_index *
                                  fe_v.shape_value(i, q_index) /
                                  particle_velocities[q_index] *
                                  material_derivative_vel[q_index][coordinate] *
                                  advection_matrices[coordinate](component_i,
                                                                 component_j) *
                                  fe_v.shape_value(j, q_index) * JxW[q_index];
                              }
                          }
                        for (unsigned int coordinate_1 = 0; coordinate_1 < 3;
                             ++coordinate_1)
                          {
                            for (unsigned int coordinate_2 = coordinate_1;
                                 coordinate_2 < 3;
                                 ++coordinate_2)
                              {
                                if (coordinate_1 == coordinate_2)
                                  {
                                    // \grad_phi p
                                    // \jacobian[coordinate_1][coordinate_2]
                                    // Ap_coordinate_1,coordinate_2 * \phi
                                    copy_data.cell_matrix(i, j) +=
                                      fe_v.shape_grad(i, q_index)[dim_ps - 1] *
                                      q_points[q_index][dim_ps - 1] *
                                      jacobians_vel[q_index][coordinate_1]
                                                   [coordinate_2] *
                                      adv_mat_products
                                        [3 * coordinate_1 -
                                         coordinate_1 * (coordinate_1 + 1) / 2 +
                                         coordinate_2](component_i,
                                                       component_j) *
                                      fe_v.shape_value(j, q_index) *
                                      JxW[q_index];
                                    // \phi *
                                    // jacobian[coordinate_1][coordinate_2]
                                    // * Ap_coordinate_1, coordinate_2 *
                                    // \phi
                                    copy_data.cell_matrix(i, j) +=
                                      fe_v.shape_value(i, q_index) *
                                      jacobians_vel[q_index][coordinate_1]
                                                   [coordinate_2] *
                                      adv_mat_products
                                        [3 * coordinate_1 -
                                         coordinate_1 * (coordinate_1 + 1) / 2 +
                                         coordinate_2](component_i,
                                                       component_j) *
                                      fe_v.shape_value(j, q_index) *
                                      JxW[q_index];

                                    if constexpr (
                                      (vfp_flags &
                                       VFPFlags::
                                         scaled_distribution_function) !=
                                      VFPFlags::none)
                                      {
                                        // \phi * scaling_spectral_index *
                                        // \jacobian[coordinate_1][coordinate_2]
                                        // Ap_coordinate_1,coordinate_2 * \phi
                                        copy_data.cell_matrix(i, j) +=
                                          fe_v.shape_value(i, q_index) *
                                          scaling_spectral_index *
                                          jacobians_vel[q_index][coordinate_1]
                                                       [coordinate_2] *
                                          adv_mat_products
                                            [3 * coordinate_1 -
                                             coordinate_1 * (coordinate_1 + 1) /
                                               2 +
                                             coordinate_2](component_i,
                                                           component_j) *
                                          fe_v.shape_value(j, q_index) *
                                          JxW[q_index];
                                      }
                                  }
                                else
                                  {
                                    // symmetry
                                    // component_1, component_2
                                    // \grad_phi p
                                    // \jacobian[coordinate_1][coordinate_2]
                                    // Ap_coordinate_1,coordinate_2 * \phi
                                    copy_data.cell_matrix(i, j) +=
                                      fe_v.shape_grad(i, q_index)[dim_ps - 1] *
                                      q_points[q_index][dim_ps - 1] *
                                      jacobians_vel[q_index][coordinate_1]
                                                   [coordinate_2] *
                                      adv_mat_products
                                        [3 * coordinate_1 -
                                         coordinate_1 * (coordinate_1 + 1) / 2 +
                                         coordinate_2](component_i,
                                                       component_j) *
                                      fe_v.shape_value(j, q_index) *
                                      JxW[q_index];

                                    // \phi *
                                    // jacobian[coordinate_1][coordinate_2]
                                    // * Ap_coordinate_1, coordinate_2 *
                                    // \phi
                                    copy_data.cell_matrix(i, j) +=
                                      fe_v.shape_value(i, q_index) *
                                      jacobians_vel[q_index][coordinate_1]
                                                   [coordinate_2] *
                                      adv_mat_products
                                        [3 * coordinate_1 -
                                         coordinate_1 * (coordinate_1 + 1) / 2 +
                                         coordinate_2](component_i,
                                                       component_j) *
                                      fe_v.shape_value(j, q_index) *
                                      JxW[q_index];

                                    // component_2, component_1
                                    // \grad_phi p
                                    // \jacobian[coordinate_1][coordinate_2]
                                    // Ap_coordinate_1,coordinate_2 * \phi
                                    copy_data.cell_matrix(i, j) +=
                                      fe_v.shape_grad(i, q_index)[dim_ps - 1] *
                                      q_points[q_index][dim_ps - 1] *
                                      jacobians_vel[q_index][coordinate_2]
                                                   [coordinate_1] *
                                      adv_mat_products
                                        [3 * coordinate_1 -
                                         coordinate_1 * (coordinate_1 + 1) / 2 +
                                         coordinate_2](component_i,
                                                       component_j) *
                                      fe_v.shape_value(j, q_index) *
                                      JxW[q_index];

                                    // \phi *
                                    // jacobian[coordinate_1][coordinate_2]
                                    // * Ap_coordinate_1, coordinate_2 *
                                    // \phi
                                    copy_data.cell_matrix(i, j) +=
                                      fe_v.shape_value(i, q_index) *
                                      jacobians_vel[q_index][coordinate_2]
                                                   [coordinate_1] *
                                      adv_mat_products
                                        [3 * coordinate_1 -
                                         coordinate_1 * (coordinate_1 + 1) / 2 +
                                         coordinate_2](component_i,
                                                       component_j) *
                                      fe_v.shape_value(j, q_index) *
                                      JxW[q_index];
                                    if constexpr (
                                      (vfp_flags &
                                       VFPFlags::
                                         scaled_distribution_function) !=
                                      VFPFlags::none)
                                      {
                                        // component_1, component_2
                                        // \phi * scaling_spectral_index *
                                        // \jacobian[coordinate_1][coordinate_2]
                                        // Ap_coordinate_1,coordinate_2 * \phi
                                        copy_data.cell_matrix(i, j) +=
                                          fe_v.shape_value(i, q_index) *
                                          scaling_spectral_index *
                                          jacobians_vel[q_index][coordinate_1]
                                                       [coordinate_2] *
                                          adv_mat_products
                                            [3 * coordinate_1 -
                                             coordinate_1 * (coordinate_1 + 1) /
                                               2 +
                                             coordinate_2](component_i,
                                                           component_j) *
                                          fe_v.shape_value(j, q_index) *
                                          JxW[q_index];

                                        // component_2, component_1
                                        // \phi * scaling_spectral_index
                                        // \jacobian[coordinate_1][coordinate_2]
                                        // Ap_coordinate_1,coordinate_2 * \phi
                                        copy_data.cell_matrix(i, j) +=
                                          fe_v.shape_value(i, q_index) *
                                          scaling_spectral_index *
                                          jacobians_vel[q_index][coordinate_2]
                                                       [coordinate_1] *
                                          adv_mat_products
                                            [3 * coordinate_1 -
                                             coordinate_1 * (coordinate_1 + 1) /
                                               2 +
                                             coordinate_2](component_i,
                                                           component_j) *
                                          fe_v.shape_value(j, q_index) *
                                          JxW[q_index];
                                      }
                                  }
                              }
                          }
                        for (unsigned int coordinate_1 = 0; coordinate_1 < 3;
                             ++coordinate_1)
                          {
                            for (unsigned int coordinate_2 = 0;
                                 coordinate_2 < 3;
                                 ++coordinate_2)
                              {
                                // \phi *
                                // jacobian[coordinate_1][coordinate_2]
                                // * T_coordinate_2,coordinate_1 * \phi
                                copy_data.cell_matrix(i, j) +=
                                  fe_v.shape_value(i, q_index) *
                                  jacobians_vel[q_index][coordinate_1]
                                               [coordinate_2] *
                                  t_matrices[coordinate_2 * 3 + coordinate_1](
                                    component_i, component_j) *
                                  fe_v.shape_value(j, q_index) * JxW[q_index];
                              }
                          }
                      }
                  }
              }
          }
      }
  };
  // assemble boundary face terms
  const auto boundary_worker = [&](const Iterator      &cell,
                                   const unsigned int  &face_no,
                                   ScratchData<dim_ps> &scratch_data,
                                   CopyData            &copy_data) {
    scratch_data.fe_values_face.reinit(cell, face_no);
    const FEFaceValuesBase<dim_ps> &fe_face_v = scratch_data.fe_values_face;
    // Every shape function on the cell could contribute to the face
    // integral, hence n_facet_dofs = n_dofs_per_cell
    const unsigned int n_facet_dofs = fe_face_v.get_fe().n_dofs_per_cell();
    // NOTE: copy_data is not reinitialized, the cell_workers contribution
    // to the cell_dg_matrix should not be deleted

    const unsigned int       boundary_id = cell->face(face_no)->boundary_id();
    const BoundaryConditions boundary_condition =
      vfp_parameters.boundary_conditions[boundary_id];
    const unsigned int boundary_coordinate = boundary_id / 2;

    const std::vector<Point<dim_ps>> &q_points =
      fe_face_v.get_quadrature_points();
    const std::vector<double>            &JxW = fe_face_v.get_JxW_values();
    const std::vector<Tensor<1, dim_ps>> &normals =
      fe_face_v.get_normal_vectors();
    // NOLINTBEGIN(google-readability-casting)
    std::vector<FullMatrix<double>> positive_flux_matrices(
      q_points.size(), FullMatrix<double>(pde_system.system_size));
    std::vector<FullMatrix<double>> negative_flux_matrices(
      q_points.size(), FullMatrix<double>(pde_system.system_size));
    // NOLINTEND(google-readability-casting)

    upwind_flux.compute_upwind_fluxes(q_points,
                                      normals,
                                      positive_flux_matrices,
                                      negative_flux_matrices);

    std::vector<Vector<double>> bc_values(
      q_points.size(), Vector<double>(pde_system.system_size));
    bc_value_function.bc_vector_value_list(q_points, boundary_id, bc_values);

    for (unsigned int q_index : fe_face_v.quadrature_point_indices())
      {
        for (unsigned int i = 0; i < n_facet_dofs; ++i)
          {
            const unsigned int component_i =
              fe_face_v.get_fe().system_to_component_index(i).first;
            for (unsigned int j = 0; j < n_facet_dofs; ++j)
              {
                const unsigned int component_j =
                  fe_face_v.get_fe().system_to_component_index(j).first;

                // TODO: Removed constexpr check - why was it here?
                // if constexpr ((vfp_flags & VFPFlags::spatial_advection) !=
                //               VFPFlags::none)

                switch (boundary_condition)
                  {
                    case BoundaryConditions::continuous:
                      {
                        copy_data.cell_matrix(i, j) +=
                          fe_face_v.shape_value(i, q_index) *
                          positive_flux_matrices[q_index](component_i,
                                                          component_j) *
                          fe_face_v.shape_value(j, q_index) * JxW[q_index];
                        copy_data.cell_matrix(i, j) +=
                          fe_face_v.shape_value(i, q_index) *
                          negative_flux_matrices[q_index](component_i,
                                                          component_j) *
                          fe_face_v.shape_value(j, q_index) * JxW[q_index];
                        break;
                      }
                    case BoundaryConditions::zero_inflow:
                      {
                        copy_data.cell_matrix(i, j) +=
                          fe_face_v.shape_value(i, q_index) *
                          positive_flux_matrices[q_index](component_i,
                                                          component_j) *
                          fe_face_v.shape_value(j, q_index) * JxW[q_index];
                        break;
                      }
                    case BoundaryConditions::reflective:
                      {
                        copy_data.cell_matrix(i, j) +=
                          fe_face_v.shape_value(i, q_index) *
                          positive_flux_matrices[q_index](component_i,
                                                          component_j) *
                          fe_face_v.shape_value(j, q_index) * JxW[q_index];
                        copy_data.cell_matrix(i, j) +=
                          fe_face_v.shape_value(i, q_index) *
                          negative_flux_matrices[q_index](component_i,
                                                          component_j) *
                          reflective_bc_signature[boundary_coordinate]
                                                 [component_j] *
                          fe_face_v.shape_value(j, q_index) * JxW[q_index];
                        break;
                      }
                    case BoundaryConditions::inflow:
                      {
                        copy_data.cell_matrix(i, j) +=
                          fe_face_v.shape_value(i, q_index) *
                          positive_flux_matrices[q_index](component_i,
                                                          component_j) *
                          fe_face_v.shape_value(j, q_index) * JxW[q_index];
                        break;
                      }
                    case BoundaryConditions::periodic:
                      break;
                    default:
                      Assert(false, ExcNotImplemented());
                  }
              }
            // Inflow boundary conditions result in a right-hand side
            if (boundary_condition == BoundaryConditions::inflow)
              {
                for (unsigned int k = 0; k < pde_system.system_size; ++k)
                  copy_data.cell_rhs(i) -=
                    fe_face_v.shape_value(i, q_index) *
                    negative_flux_matrices[q_index](component_i, k) *
                    bc_values[q_index][k] * JxW[q_index];
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
  const auto face_worker = [&](const Iterator      &cell,
                               const unsigned int  &face_no,
                               const unsigned int  &subface_no,
                               const Iterator      &neighbor_cell,
                               const unsigned int  &neighbor_face_no,
                               const unsigned int  &neighbor_subface_no,
                               ScratchData<dim_ps> &scratch_data,
                               CopyData            &copy_data) {
    // NOTE: The flag MeshWorker::assemble_own_interior_faces_both will not
    // work for this face_worker. It implicitly assumes that faces are only
    // assembled once. And hence subface_no, can be ignored
    (void)subface_no; // suppress compiler warning

    FEFaceValues<dim_ps> &fe_v_face = scratch_data.fe_values_face;
    fe_v_face.reinit(cell, face_no);

    // NOTE: I do not know how to initialise this reference whose value
    // depends on the value of neighbor_subface_no. Subfaces only exists if
    // the triangulation was refined differently in different parts of the
    // domain. Since I am currently testing, I am always refining globally,
    // so no subfaces exist and I just ignore this case. Hence
    // neighbor_subface_no can be ignored
    (void)neighbor_subface_no;
    FEFaceValues<dim_ps> &fe_v_face_neighbor =
      scratch_data.fe_values_face_neighbor;
    fe_v_face_neighbor.reinit(neighbor_cell, neighbor_face_no);
    // Create an element at the end of the vector containing the face data
    copy_data.face_data.emplace_back();
    CopyDataFace      &copy_data_face = copy_data.face_data.back();
    const unsigned int n_dofs         = fe_v_face.get_fe().n_dofs_per_cell();
    copy_data_face.reinit(cell, neighbor_cell, n_dofs);

    const std::vector<Point<dim_ps>> &q_points =
      fe_v_face.get_quadrature_points();
    const std::vector<double>            &JxW = fe_v_face.get_JxW_values();
    const std::vector<Tensor<1, dim_ps>> &normals =
      fe_v_face.get_normal_vectors();
    // For every interior face there is an in- and outflow represented by
    // the corresponding flux matrices
    // NOLINTBEGIN(google-readability-casting)
    std::vector<FullMatrix<double>> positive_flux_matrices(
      q_points.size(), FullMatrix<double>(pde_system.system_size));
    std::vector<FullMatrix<double>> negative_flux_matrices(
      q_points.size(), FullMatrix<double>(pde_system.system_size));
    // NOLINTEND(google-readability-casting)

    upwind_flux.compute_upwind_fluxes(q_points,
                                      normals,
                                      positive_flux_matrices,
                                      negative_flux_matrices);

    for (unsigned int q_index : fe_v_face.quadrature_point_indices())
      {
        // cell_dg_matrix_11
        for (unsigned int i : fe_v_face.dof_indices())
          {
            const unsigned int component_i =
              fe_v_face.get_fe().system_to_component_index(i).first;
            for (unsigned int j : fe_v_face.dof_indices())
              {
                unsigned int component_j =
                  fe_v_face.get_fe().system_to_component_index(j).first;
                copy_data_face.cell_dg_matrix_11(i, j) +=
                  fe_v_face.shape_value(i, q_index) *
                  positive_flux_matrices[q_index](component_i, component_j) *
                  fe_v_face.shape_value(j, q_index) * JxW[q_index];
              }
          }
        // cell_dg_matrix_12
        for (unsigned int i : fe_v_face.dof_indices())
          {
            const unsigned int component_i =
              fe_v_face.get_fe().system_to_component_index(i).first;
            for (unsigned int j : fe_v_face_neighbor.dof_indices())
              {
                unsigned int component_j = fe_v_face_neighbor.get_fe()
                                             .system_to_component_index(j)
                                             .first;
                copy_data_face.cell_dg_matrix_12(i, j) -=
                  fe_v_face_neighbor.shape_value(i, q_index) *
                  positive_flux_matrices[q_index](component_i, component_j) *
                  fe_v_face.shape_value(j, q_index) * JxW[q_index];
              }
          }
        // cell_dg_matrix_21
        for (unsigned int i : fe_v_face_neighbor.dof_indices())
          {
            const unsigned int component_i =
              fe_v_face_neighbor.get_fe().system_to_component_index(i).first;
            for (unsigned int j : fe_v_face.dof_indices())
              {
                unsigned int component_j =
                  fe_v_face.get_fe().system_to_component_index(j).first;
                copy_data_face.cell_dg_matrix_21(i, j) +=
                  fe_v_face.shape_value(i, q_index) *
                  negative_flux_matrices[q_index](component_i, component_j) *
                  fe_v_face_neighbor.shape_value(j, q_index) * JxW[q_index];
              }
          }
        // cell_dg_matrix_22
        for (unsigned int i : fe_v_face_neighbor.dof_indices())
          {
            const unsigned int component_i =
              fe_v_face_neighbor.get_fe().system_to_component_index(i).first;
            for (unsigned int j : fe_v_face_neighbor.dof_indices())
              {
                unsigned int component_j = fe_v_face_neighbor.get_fe()
                                             .system_to_component_index(j)
                                             .first;
                copy_data_face.cell_dg_matrix_22(i, j) -=
                  fe_v_face_neighbor.shape_value(i, q_index) *
                  negative_flux_matrices[q_index](component_i, component_j) *
                  fe_v_face_neighbor.shape_value(j, q_index) * JxW[q_index];
              }
          }
      }
  };
  // copier for the mesh_loop function
  const auto copier = [&](const CopyData &c) {
    constraints.distribute_local_to_global(c.cell_matrix,
                                           c.cell_rhs,
                                           c.local_dof_indices,
                                           dg_matrix,
                                           locally_owned_current_bc);
    for (auto &cdf : c.face_data)
      {
        for (unsigned int i = 0; i < cdf.local_dof_indices.size(); ++i)
          for (unsigned int j = 0; j < cdf.local_dof_indices.size(); ++j)
            {
              dg_matrix.add(cdf.local_dof_indices[i],
                            cdf.local_dof_indices[j],
                            cdf.cell_dg_matrix_11(i, j));
              dg_matrix.add(cdf.local_dof_indices_neighbor[i],
                            cdf.local_dof_indices[j],
                            cdf.cell_dg_matrix_12(i, j));
              dg_matrix.add(cdf.local_dof_indices[i],
                            cdf.local_dof_indices_neighbor[j],
                            cdf.cell_dg_matrix_21(i, j));
              dg_matrix.add(cdf.local_dof_indices_neighbor[i],
                            cdf.local_dof_indices_neighbor[j],
                            cdf.cell_dg_matrix_22(i, j));
            }
      }
  };

  ScratchData<dim_ps> scratch_data(mapping, fe, quadrature, quadrature_face);
  CopyData            copy_data;
  saplog << "Begin the assembly of the DG matrix." << std::endl;
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
  dg_matrix.compress(VectorOperation::add);
  locally_owned_current_bc.compress(VectorOperation::add);
  saplog << "The DG matrix was assembled." << std::endl;
}



template <unsigned int dim>
void
sapphirepp::VFP::VFPSolver<dim>::steady_state_solve()
{
  TimerOutput::Scope timer_section(timer, "VFP - Steady state");
  saplog << "Steady-state solve" << std::endl;
  LogStream::Prefix prefix("SteadyState", saplog);

  SolverControl              solver_control(vfp_parameters.solver_max_iter);
  PETScWrappers::SolverGMRES solver(solver_control);

  PETScWrappers::PreconditionBlockJacobi preconditioner;
  preconditioner.initialize(dg_matrix);
  // Right-hand side is determined by the source and the non-homogeneous
  // boundary condition
  if constexpr ((vfp_flags & VFPFlags::source) != VFPFlags::none)
    system_rhs.add(1., locally_owned_current_source);
  // TODO: add check if non-homogeneous bc is set
  system_rhs.add(1., locally_owned_current_bc);
  solver_control.set_tolerance(vfp_parameters.solver_tolerance *
                               system_rhs.l2_norm());
  saplog << "Use solver_control(" << vfp_parameters.solver_max_iter << ", "
         << vfp_parameters.solver_tolerance
         << "*rhs = " << solver_control.tolerance() << ")" << std::endl;
  saplog << "Start steady-state solver" << std::endl;
  {
    TimerOutput::Scope timer_solve(timer, "VFP - Matrix solve");
    solver.solve(dg_matrix,
                 locally_owned_previous_solution,
                 system_rhs,
                 preconditioner);
  }

  locally_relevant_current_solution = locally_owned_previous_solution;

  saplog << "Solver converged in " << solver_control.last_step()
         << " iterations." << std::endl;
}



template <unsigned int dim>
void
sapphirepp::VFP::VFPSolver<dim>::theta_method(const double time,
                                              const double time_step)
{
  TimerOutput::Scope timer_section(timer, "VFP - Theta method");
  LogStream::Prefix  prefix("ThetaMethod", saplog);
  // Equation: (mass_matrix + time_step * theta * dg_matrix(time +
  // time_step)) f(time + time_step) = (mass_matrix - time_step * (1 -
  // theta) * dg_matrix(time) ) f(time) + time_step * theta * s(time +
  // time_step) + time_step * (1 - theta) * s(time)
  const double theta = vfp_parameters.theta;

  locally_owned_previous_solution = locally_relevant_current_solution;

  mass_matrix.vmult(system_rhs, locally_owned_previous_solution);
  PETScWrappers::MPI::Vector tmp(locally_owned_dofs, mpi_communicator);
  dg_matrix.vmult(tmp, locally_owned_previous_solution);
  system_rhs.add(-time_step * (1 - theta), tmp);
  // Source term
  if constexpr ((vfp_flags & VFPFlags::source) != VFPFlags::none)
    {
      if constexpr ((vfp_flags & VFPFlags::time_independent_source) ==
                    VFPFlags::none) // time dependent source
        {
          system_rhs.add((1 - theta) * time_step, locally_owned_current_source);
          // Update the source term
          TimerOutput::Scope timer_source(timer, "VFP - Source");
          source_function.set_time(time + time_step);
          VectorTools::create_right_hand_side(mapping,
                                              dof_handler,
                                              quadrature,
                                              source_function,
                                              locally_owned_current_source);
          system_rhs.add(theta * time_step, locally_owned_current_source);
        }
      else
        {
          system_rhs.add(time_step, locally_owned_current_source);
        }
    }
  // Boundary term at current time step
  system_rhs.add((1 - theta) * time_step, locally_owned_current_bc);

  // Update DG matrix
  // Since the the dg_matrix depends on the velocity field (and/or the
  // magnetic field) and the velocity field may depend on time, it needs to
  // reassembled every time step. This is not true for the mass matrix ( but
  // it may if the grid adapts after a specified amount of time steps)
  if constexpr ((vfp_flags & VFPFlags::time_independent_fields) ==
                VFPFlags::none) // time dependent fields
    {
      // the non-homogeneous boundary term is updated and, thus, must be nulled
      locally_owned_current_bc = 0;
      dg_matrix                = 0;
      assemble_dg_matrix(time + time_step);
    }

  // Boundary term at next time step
  system_rhs.add(theta * time_step, locally_owned_current_bc);

  system_matrix.copy_from(mass_matrix);
  system_matrix.add(time_step * theta, dg_matrix);

  SolverControl              solver_control(vfp_parameters.solver_max_iter,
                               vfp_parameters.solver_tolerance *
                                 system_rhs.l2_norm());
  PETScWrappers::SolverGMRES solver(solver_control);

  // PETScWrappers::PreconditionBoomerAMG preconditioner;
  // PETScWrappers::PreconditionBoomerAMG::AdditionalData data;
  // data.symmetric_operator = true;
  // preconditioner.initialize(system_matrix, data);

  // PETScWrappers::PreconditionSOR preconditioner;
  // preconditioner.initialize(system_matrix);

  PETScWrappers::PreconditionBlockJacobi preconditioner;
  preconditioner.initialize(system_matrix);

  {
    TimerOutput::Scope timer_solve(timer, "VFP - Matrix solve");
    solver.solve(system_matrix,
                 locally_owned_previous_solution,
                 system_rhs,
                 preconditioner);
  }

  // Update the solution
  locally_relevant_current_solution = locally_owned_previous_solution;

  saplog << "Solver converged in " << solver_control.last_step()
         << " iterations." << std::endl;
}



template <unsigned int dim>
void
sapphirepp::VFP::VFPSolver<dim>::explicit_runge_kutta(const double time,
                                                      const double time_step)
{
  TimerOutput::Scope timer_section(timer, "VFP - ERK");
  LogStream::Prefix  prefix("ERK", saplog);
  // ERK 4
  // \df(t)/dt = - mass_matrix_inv * (dg_matrix(t) * f(t) - s(t))
  // Butcher's array
  Vector<double> a({0.5, 0.5, 1.});
  Vector<double> b({1. / 6, 1. / 3, 1. / 3, 1. / 6});
  // NOTE: c is only necessary if the velocity field and the magnetic field
  // are time dependent.
  Vector<double> c({0., 0.5, 0.5, 1.});

  // The mass matrix needs to be "inverted" in every stage
  PETScWrappers::PreconditionNone preconditioner;
  preconditioner.initialize(mass_matrix);

  SolverControl           solver_control(vfp_parameters.solver_max_iter);
  PETScWrappers::SolverCG cg(solver_control);

  // I need the previous solution to compute k_0, k_1, k_2, k_3
  locally_owned_previous_solution = locally_relevant_current_solution;
  // locally_owned_current_solution is the result (or return) vector of the
  // ERK4 method. I cannot use locally_owned previous solution, because it
  // is required for the computation of the ks.
  PETScWrappers::MPI::Vector locally_owned_current_solution =
    locally_owned_previous_solution;
  // a temporary vector to compute f(time) + a*time_step*k_{i-1}
  PETScWrappers::MPI::Vector temp(locally_owned_dofs, mpi_communicator);

  // k_0
  PETScWrappers::MPI::Vector k_0(locally_owned_dofs, mpi_communicator);
  // dg_matrix(time)
  dg_matrix.vmult(system_rhs, locally_owned_previous_solution);
  if constexpr ((vfp_flags & VFPFlags::source) != VFPFlags::none)
    {
      system_rhs.add(-1., locally_owned_current_source);
    }
  system_rhs.add(-1., locally_owned_current_bc);
  solver_control.set_tolerance(vfp_parameters.solver_tolerance *
                               system_rhs.l2_norm());
  {
    TimerOutput::Scope timer_solve(timer, "VFP - Matrix solve");
    cg.solve(mass_matrix, k_0, system_rhs, preconditioner);
  }
  saplog << "Stage s: " << 0 << "	Solver converged in "
         << solver_control.last_step() << " iterations." << std::endl;
  k_0 *= -1.;
  locally_owned_current_solution.add(b[0] * time_step, k_0);

  // k_1
  PETScWrappers::MPI::Vector k_1(locally_owned_dofs, mpi_communicator);
  // Compute dg_matrix(time + c[1] * time_step) if the fields are time
  // dependent
  if constexpr ((vfp_flags & VFPFlags::time_independent_fields) ==
                VFPFlags::none) // time dependent fields
    {
      dg_matrix                = 0;
      locally_owned_current_bc = 0;
      assemble_dg_matrix(time + c[1] * time_step);
    }
  temp.add(1., locally_owned_previous_solution, a[0] * time_step, k_0);
  dg_matrix.vmult(system_rhs, temp);
  if constexpr ((vfp_flags & VFPFlags::source) != VFPFlags::none)
    {
      if constexpr ((vfp_flags & VFPFlags::time_independent_source) ==
                    VFPFlags::none) // time dependent source
        {
          TimerOutput::Scope timer_source(timer, "VFP - Source");
          source_function.set_time(time + c[1] * time_step);
          VectorTools::create_right_hand_side(mapping,
                                              dof_handler,
                                              quadrature,
                                              source_function,
                                              locally_owned_current_source);
          system_rhs.add(-1., locally_owned_current_source);
        }
      else
        {
          system_rhs.add(-1., locally_owned_current_source);
        }
    }
  system_rhs.add(-1., locally_owned_current_bc);
  solver_control.set_tolerance(vfp_parameters.solver_tolerance *
                               system_rhs.l2_norm());
  {
    TimerOutput::Scope timer_solve(timer, "VFP - Matrix solve");
    cg.solve(mass_matrix, k_1, system_rhs, preconditioner);
  }
  saplog << "	Stage s: " << 1 << "	Solver converged in "
         << solver_control.last_step() << " iterations." << std::endl;
  k_1 *= -1.;
  locally_owned_current_solution.add(b[1] * time_step, k_1);
  temp = 0;

  // k_2
  PETScWrappers::MPI::Vector k_2(locally_owned_dofs, mpi_communicator);
  // NOTE: For k_2 it is not necessary to reassemble the dg_matrix,
  // since c[1] = c[2]
  temp.add(1., locally_owned_previous_solution, a[1] * time_step, k_1);
  dg_matrix.vmult(system_rhs, temp);
  if constexpr ((vfp_flags & VFPFlags::source) != VFPFlags::none)
    {
      if constexpr ((vfp_flags & VFPFlags::time_independent_source) ==
                    VFPFlags::none) // time dependent source
        {
          TimerOutput::Scope timer_source(timer, "VFP - Source");
          source_function.set_time(time + c[2] * time_step);
          VectorTools::create_right_hand_side(mapping,
                                              dof_handler,
                                              quadrature,
                                              source_function,
                                              locally_owned_current_source);
          system_rhs.add(-1., locally_owned_current_source);
        }
      else
        {
          system_rhs.add(-1., locally_owned_current_source);
        }
    }
  system_rhs.add(-1., locally_owned_current_bc);
  solver_control.set_tolerance(vfp_parameters.solver_tolerance *
                               system_rhs.l2_norm());
  {
    TimerOutput::Scope timer_solve(timer, "VFP - Matrix solve");
    cg.solve(mass_matrix, k_2, system_rhs, preconditioner);
  }
  saplog << "	Stage s: " << 2 << "	Solver converged in "
         << solver_control.last_step() << " iterations." << std::endl;
  k_2 *= -1.;
  locally_owned_current_solution.add(b[2] * time_step, k_2);
  temp = 0;

  // k_3
  PETScWrappers::MPI::Vector k_3(locally_owned_dofs, mpi_communicator);
  // Compute dg_matrix(time + c[1] * time_step) if the fields are time
  // dependent
  if constexpr ((vfp_flags & VFPFlags::time_independent_fields) ==
                VFPFlags::none) // time dependent fields
    {
      dg_matrix                = 0;
      locally_owned_current_bc = 0;
      assemble_dg_matrix(time + c[3] * time_step);
    }
  temp.add(1., locally_owned_previous_solution, a[2] * time_step, k_2);
  dg_matrix.vmult(system_rhs, temp);
  if constexpr ((vfp_flags & VFPFlags::source) != VFPFlags::none)
    {
      if constexpr ((vfp_flags & VFPFlags::time_independent_source) ==
                    VFPFlags::none) // time dependent source
        {
          TimerOutput::Scope timer_source(timer, "VFP - Source");
          source_function.set_time(time + c[3] * time_step);
          VectorTools::create_right_hand_side(mapping,
                                              dof_handler,
                                              quadrature,
                                              source_function,
                                              locally_owned_current_source);
          system_rhs.add(-1., locally_owned_current_source);
        }
      else
        {
          system_rhs.add(-1., locally_owned_current_source);
        }
    }
  system_rhs.add(-1., locally_owned_current_bc);
  solver_control.set_tolerance(vfp_parameters.solver_tolerance *
                               system_rhs.l2_norm());
  {
    TimerOutput::Scope timer_solve(timer, "VFP - Matrix solve");
    cg.solve(mass_matrix, k_3, system_rhs, preconditioner);
  }
  saplog << "	Stage s: " << 3 << "	Solver converged in "
         << solver_control.last_step() << " iterations." << std::endl;
  k_3 *= -1.;
  locally_owned_current_solution.add(b[3] * time_step, k_3);
  temp = 0;

  locally_relevant_current_solution = locally_owned_current_solution;
}



template <unsigned int dim>
void
sapphirepp::VFP::VFPSolver<dim>::low_storage_explicit_runge_kutta(
  const double time,
  const double time_step)
{
  TimerOutput::Scope timer_section(timer, "VFP - LSERK");
  LogStream::Prefix  prefix("LSERK", saplog);
  // \df(t)/dt = - mass_matrix_inv * (dg_matrix(t) * f(t) - s(t))
  // see Hesthaven p.64
  Vector<double> a({0.,
                    -567301805773. / 1357537059087,
                    -2404267990393. / 2016746695238,
                    -3550918686646. / 2091501179385,
                    -1275806237668. / 842570457699});
  Vector<double> b({1432997174477. / 9575080441755,
                    5161836677717. / 13612068292357,
                    1720146321549. / 2090206949498,
                    3134564353537. / 4481467310338,
                    2277821191437. / 14882151754819});
  // NOTE: I only need c if the velocity field, the magnetic field or the
  // non-homogeneous boundary are time dependent
  Vector<double> c({0.,
                    1432997174477. / 9575080441755,
                    2526269341429. / 6820363962896,
                    2006345519317. / 3224310063776,
                    2802321613138. / 2924317926251});

  // The mass matrix needs to be "inverted" in every stage

  // PETScWrappers::PreconditionBoomerAMG preconditioner;
  // PETScWrappers::PreconditionBoomerAMG::AdditionalData data;
  // data.symmetric_operator = true;
  PETScWrappers::PreconditionNone preconditioner;
  preconditioner.initialize(mass_matrix);

  SolverControl           solver_control(vfp_parameters.solver_max_iter);
  PETScWrappers::SolverCG cg(solver_control);
  // NOTE: The locally_relevant_current_solution is a "ghosted" vector and
  // it cannot be written to. It is necessary to use a vector that does not
  // contain ghost cells. We extract the locally owned part with the equal
  // sign operator.
  locally_owned_previous_solution = locally_relevant_current_solution;

  PETScWrappers::MPI::Vector k(locally_owned_dofs, mpi_communicator);
  PETScWrappers::MPI::Vector temp(locally_owned_dofs, mpi_communicator);
  for (unsigned int s = 0; s < 5; ++s)
    {
      // only assemble the dg_matrix in every stage if the fields are time
      // dependent
      if constexpr ((vfp_flags & VFPFlags::time_independent_fields) ==
                    VFPFlags::none) // time dependent fields
        assemble_dg_matrix(time + c[s] * time_step);

      dg_matrix.vmult(system_rhs, locally_owned_previous_solution);
      // non-homogeneous boundary
      system_rhs.add(-1., locally_owned_current_bc);

      if constexpr ((vfp_flags & VFPFlags::time_independent_fields) ==
                    VFPFlags::none) // time dependent fields
        {
          dg_matrix                = 0;
          locally_owned_current_bc = 0;
        }

      if constexpr ((vfp_flags & VFPFlags::source) != VFPFlags::none)
        {
          if constexpr ((vfp_flags & VFPFlags::time_independent_source) ==
                        VFPFlags::none) // time dependent source
            {
              TimerOutput::Scope timer_source(timer, "VFP - Source");
              source_function.set_time(time + c[s] * time_step);
              VectorTools::create_right_hand_side(mapping,
                                                  dof_handler,
                                                  quadrature,
                                                  source_function,
                                                  locally_owned_current_source);
              system_rhs.add(-1., locally_owned_current_source);
            }
          else
            {
              system_rhs.add(-1., locally_owned_current_source);
            }
        }

      solver_control.set_tolerance(vfp_parameters.solver_tolerance *
                                   system_rhs.l2_norm());
      {
        TimerOutput::Scope timer_solve(timer, "VFP - Matrix solve");
        cg.solve(mass_matrix, temp, system_rhs, preconditioner);
      }
      saplog << "	Stage s: " << s << "	Solver converged in "
             << solver_control.last_step() << " iterations." << std::endl;

      k.sadd(a[s], -time_step, temp);

      locally_owned_previous_solution.add(b[s], k);
    }
  // Currently I assume that there are no constraints
  // constraints.distribute(locally_relevant_current_solution);
  locally_relevant_current_solution = locally_owned_previous_solution;
}



template <unsigned int dim>
void
sapphirepp::VFP::VFPSolver<dim>::output_results(
  const unsigned int time_step_number,
  const double       cur_time)
{
  TimerOutput::Scope timer_section(timer, "VFP - Output");
  LogStream::Prefix  prefix("Output", saplog);
  saplog << "Output results at t = " << cur_time << std::endl;
  DataOut<dim_ps> data_out;
  data_out.attach_dof_handler(dof_handler);
  data_out.add_data_vector(
    locally_relevant_current_solution,
    PDESystem::create_component_name_list(
      pde_system.system_size,
      ((vfp_flags & VFPFlags::scaled_distribution_function) != VFPFlags::none) ?
        "g_" :
        "f_"));

  // Output the partition of the mesh
  Vector<float> subdomain(triangulation.n_active_cells());
  for (unsigned int i = 0; i < subdomain.size(); ++i)
    subdomain(i) = static_cast<float>(triangulation.locally_owned_subdomain());
  data_out.add_data_vector(subdomain, "subdomain");

  if (output_parameters.debug_input_functions)
    {
      // Add input functions
      PETScWrappers::MPI::Vector temp(
        debug_input_functions_dof_handler.locally_owned_dofs(),
        mpi_communicator);
      VectorTools::interpolate(debug_input_functions_dof_handler,
                               debug_input_functions,
                               temp);
      data_out.add_data_vector(
        debug_input_functions_dof_handler,
        temp,
        debug_input_functions.create_component_name_list());
    }

  // Adapt the output to the polynomial degree of the shape functions
  data_out.build_patches(vfp_parameters.polynomial_degree);
  output_parameters.write_results<dim>(data_out, time_step_number, cur_time);

  probe_location.probe_all_points(
    dof_handler,
    mapping,
    locally_relevant_current_solution,
    time_step_number,
    cur_time,
    ((vfp_flags & VFPFlags::scaled_distribution_function) != VFPFlags::none) ?
      "g" :
      "f");
}



template <unsigned int dim>
double
sapphirepp::VFP::VFPSolver<dim>::compute_global_error(
  const Function<dim_ps>         &exact_solution,
  const VectorTools::NormType    &cell_norm,
  const VectorTools::NormType    &global_norm,
  const Function<dim_ps, double> *weight) const
{
  LogStream::Prefix prefix_vfp("VFP", saplog);
  saplog << "Compute the global error" << std::endl;
  LogStream::Prefix prefix("Error", saplog);

  Vector<float> cell_errors(triangulation.n_locally_owned_active_cells());

  const QTrapezoid<1>  q_trapezoid;
  const QIterated<dim> q_iterated(q_trapezoid,
                                  vfp_parameters.polynomial_degree * 2 + 1);

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
sapphirepp::VFP::VFPSolver<dim>::compute_weighted_norm(
  const VectorTools::NormType    &cell_norm,
  const VectorTools::NormType    &global_norm,
  const Function<dim_ps, double> *weight) const
{
  LogStream::Prefix prefix_vfp("VFP", saplog);
  saplog << "Compute the weighted norm" << std::endl;
  LogStream::Prefix prefix("Norm", saplog);

  Vector<float> cell_norms(triangulation.n_locally_owned_active_cells());

  const QTrapezoid<1>  q_trapezoid;
  const QIterated<dim> q_iterated(q_trapezoid,
                                  vfp_parameters.polynomial_degree * 2 + 1);

  const Functions::ZeroFunction<dim_ps> zero_function(pde_system.system_size);

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
const sapphirepp::VFP::PDESystem &
sapphirepp::VFP::VFPSolver<dim>::get_pde_system() const
{
  return pde_system;
}



template <unsigned int dim>
const typename sapphirepp::VFP::VFPSolver<dim>::Triangulation &
sapphirepp::VFP::VFPSolver<dim>::get_triangulation() const
{
  return triangulation;
}



template <unsigned int dim>
const dealii::DoFHandler<dim> &
sapphirepp::VFP::VFPSolver<dim>::get_dof_handler() const
{
  return dof_handler;
}



template <unsigned int dim>
const dealii::PETScWrappers::MPI::Vector &
sapphirepp::VFP::VFPSolver<dim>::get_current_solution() const
{
  return locally_relevant_current_solution;
}



template <unsigned int dim>
const dealii::TimerOutput &
sapphirepp::VFP::VFPSolver<dim>::get_timer() const
{
  return timer;
}



// Explicit instantiations
template class sapphirepp::VFP::VFPSolver<1>;
template class sapphirepp::VFP::VFPSolver<2>;
template class sapphirepp::VFP::VFPSolver<3>;
