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
#include <string>
#include <vector>

#include "sapphirepp-logstream.h"
#include "slope-limiter.h"


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



      // The helper data types for slope limiter mesh_loop
      template <unsigned int dim, unsigned int spacedim>
      class ScratchDataSlopeLimiter
      {
      public:
        // Constructor
        ScratchDataSlopeLimiter(const Mapping<dim, spacedim>       &mapping,
                                const FiniteElement<dim, spacedim> &fe,
                                const Quadrature<dim>              &quadrature,
                                const UpdateFlags update_flags =
                                  update_values | update_gradients |
                                  update_quadrature_points | update_JxW_values)
          : fe_values(mapping, fe, quadrature, update_flags)
        {}

        // Copy Constructor
        ScratchDataSlopeLimiter(
          const ScratchDataSlopeLimiter<dim, spacedim> &scratch_data)
          : fe_values(scratch_data.fe_values.get_mapping(),
                      scratch_data.fe_values.get_fe(),
                      scratch_data.fe_values.get_quadrature(),
                      scratch_data.fe_values.get_update_flags())
        {}

        FEValues<dim, spacedim> fe_values;
      };



      struct CopyDataSlopeLimiter
      {
        Vector<double>                       cell_system_rhs;
        std::vector<types::global_dof_index> local_dof_indices;

        template <typename Iterator>
        void
        reinit(const Iterator &cell, unsigned int dofs_per_cell)
        {
          cell_system_rhs.reinit(dofs_per_cell);

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
  , numerical_flux(mhd_equations)
  , mpi_communicator{MPI_COMM_WORLD}
  , triangulation(mpi_communicator)
  , dof_handler(triangulation)
  , mapping()
  , fe(FE_DGQLegendre<dim, spacedim>(mhd_parameters.polynomial_degree),
       MHDEquations::n_components)
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
  locally_owned_solution.reinit(locally_owned_dofs, mpi_communicator);
  locally_relevant_current_solution.reinit(locally_owned_dofs,
                                           locally_relevant_dofs,
                                           mpi_communicator);
  dg_rhs.reinit(locally_owned_dofs, mpi_communicator);
  system_rhs.reinit(locally_owned_dofs, mpi_communicator);

  cell_average.resize(triangulation.n_active_cells(),
                      Vector<double>(MHDEquations::n_components));

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
sapphirepp::MHD::MHDSolver<dim>::compute_cell_average()
{
  /** @todo Add compute_cell_average() function to @dealii library? */
  /**
   * @note This function assumes the @ref fe "FESystem" is made up of primitive
   *       @dealref{FE_DGQLegendre,classFE__DGQLegendre} polynomials. The cell
   *       average is then given by the 0th index of each component.
   */
  TimerOutput::Scope t(timer, "Cell average");
  saplog << "Compute cell average" << std::endl;
  AssertDimension(cell_average.size(), triangulation.n_active_cells());
  Assert(fe.is_primitive() && (fe.n_base_elements() == 1) &&
           Utilities::match_at_string_start(fe.base_element(0).get_name(),
                                            "FE_DGQLegendre"),
         ExcMessage("This function assumes that the FESystem is composed of "
                    "primitive FE_DGQLegendre elements."));

  Vector<double> local_dof_values(fe.n_dofs_per_cell());

  /**
   * @todo Use
   *       @dealref{mesh_loop,group__MeshWorker,ga76ec61fbd188fb320fe8ca166a79b322}
   *       for the cell loop?
   */
  // Compute cell_average for locally_active and gost cells
  for (const auto &cell : dof_handler.active_cell_iterators())
    if (!cell->is_artificial())
      {
        const unsigned int cell_index = cell->active_cell_index();

        // Use 0th Legendre polynomial
        const auto cell_dofs = cell->as_dof_handler_iterator(dof_handler);
        cell_dofs->get_dof_values(locally_relevant_current_solution,
                                  local_dof_values);
        for (unsigned int c = 0; c < MHDEquations::n_components; ++c)
          {
            const unsigned int i_avg    = fe.component_to_system_index(c, 0);
            cell_average[cell_index][c] = local_dof_values[i_avg];
          }
      }
}



template <unsigned int dim>
void
sapphirepp::MHD::MHDSolver<dim>::apply_limiter()
{
  TimerOutput::Scope t(timer, "Limiter");
  saplog << "Limit solution" << std::endl;
  LogStream::Prefix p("Limiter", saplog);

  using Iterator = typename DoFHandler<dim, spacedim>::active_cell_iterator;
  using namespace sapphirepp::internal::MHDSolver;


  /**
   * To limit the solution we use the following procedure:
   *
   *   1. Precompute cell averages using @ref compute_cell_average().
   *   2. Compute the limited solution/limited gradient using neighbor cells
   *      in @ref cell_worker.
   *   3. Compute the projection of the limited solution onto the DoFs:
   *      1. Compute the RHS for the projection by folding with the basis
   *         functions solution in @ref cell_worker.
   *      2. Distribute the RHS to the global @ref system_rhs in @ref copier.
   *      3. Project the limited solution on the @ref locally_owned_solution
   *         \f$ f \f$ using the @ref mass_matrix \f$ M \f$ and the
   *         @ref system_rhs \f$ b \f$, \f$ M f = b \f$.
   *   4. Distribute to the @ref locally_relevant_current_solution.
   *
   * There are some competing alternative ways to perform step 3 cell-wise:
   *
   *   - Use support points to convert point values to DoF values using
   *     @dealref{convert_generalized_support_point_values_to_dof_values(),classFiniteElement,a8f2b9817bcf98ebd43239c6da46de809}.
   *     This only works for basis functions that have support points,
   *     i.e. for LagrangePolynomials, but not for Legendre polynomials.
   *   - Locally invert the mass matrix. Disadvantage is, that we need to solve
   *     a (small) linear system on each cell.
   *   - If one has basis functions with a diagonal mass matrix, the inversion
   *     is trivial.
   *
   * One general advantage of doing this cell wise is, that we only have to
   * update the solution in cells that are actually limited.
   */

  compute_cell_average();
  // return; // Debug: no limiter

  system_rhs = 0;

  // assemble cell terms
  const auto cell_worker = [&](const Iterator &cell,
                               ScratchDataSlopeLimiter<dim, spacedim>
                                                    &scratch_data,
                               CopyDataSlopeLimiter &copy_data) {
    FEValues<dim, spacedim> &fe_v = scratch_data.fe_values;
    // reinit cell
    fe_v.reinit(cell);

    saplog << "Limit cell " << cell->active_cell_index() << std::endl;

    const unsigned int n_dofs = fe_v.get_fe().n_dofs_per_cell();
    // reinit copy_data
    copy_data.reinit(cell, n_dofs);

    const std::vector<Point<spacedim>> &q_points = fe_v.get_quadrature_points();
    const std::vector<double>          &JxW      = fe_v.get_JxW_values();

    const MHDEquations::state_type &cell_avg = get_cell_average(cell);
    saplog << "cell average: " << cell_avg << std::endl;

    MHDEquations::flux_type              cell_avg_gradient;
    MHDEquations::flux_type              limited_gradient;
    MHDEquations::flux_type              tmp_gradient;
    std::vector<MHDEquations::flux_type> neighbor_gradients;
    neighbor_gradients.reserve(cell->n_faces());

    Vector<double> local_dof_values(fe.n_dofs_per_cell());
    const auto     cell_dofs = cell->as_dof_handler_iterator(dof_handler);
    cell_dofs->get_dof_values(locally_relevant_current_solution,
                              local_dof_values);
    saplog << "local dofs: " << local_dof_values << std::endl;

    // Compute cell average gradient
    {
      TimerOutput::Scope t2(timer, "Limiter - AvgGrad");
      std::vector<std::vector<Tensor<1, spacedim>>> solution_gradients(
        fe_v.n_quadrature_points,
        std::vector<Tensor<1, spacedim>>(MHDEquations::n_components));
      fe_v.get_function_gradients(locally_relevant_current_solution,
                                  solution_gradients);
      for (unsigned int c = 0; c < MHDEquations::n_components; ++c)
        {
          for (const unsigned int q_index : fe_v.quadrature_point_indices())
            cell_avg_gradient[c] +=
              solution_gradients[q_index][c] * JxW[q_index];
          cell_avg_gradient[c] /= cell->measure();
          saplog << "cell average gradient[" << c
                 << "]: " << cell_avg_gradient[c] << std::endl;
        }

      /**
       * The calculation of the average gradient is not as simple as shown here,
       * since the derivatives of the Legendre polynomials do *not* form a
       * orthonormal basis. Instead, we could use recurrence relations for the
       * derivatives, to sum the different DoFs and find the average gradient.
       * (This is probably very hard for d>1, so it might be easiest to stay
       * with just integrating the average.)
       *
       * To set the limited solution, we could however use the reverse logic to
       * whats implemented below. In that case, we explicitly want to set all
       * higher orders to 0. Therefore, the gradient is set by the polynomial 1.
       */
      for (unsigned int c = 0; c < MHDEquations::n_components; ++c)
        {
          cell_avg_gradient[c] = 0;
          for (unsigned int d = 0; d < dim; ++d)
            {
              const unsigned int k            = fe.tensor_degree();
              const unsigned int tensor_index = std::pow(k + 1, d);
              const unsigned int index_grad =
                fe.component_to_system_index(c, tensor_index);
              // TODO: check Normalization: 4^2*sqrt(3)
              cell_avg_gradient[c][d] =
                local_dof_values[index_grad] * (16. * std::sqrt(3.));
            }
          saplog << "cell average gradient new[" << c
                 << "]: " << cell_avg_gradient[c] << std::endl;
        }
    }

    // Compute gradients by neighbor cells (ignore non-periodic boundary)
    {
      TimerOutput::Scope t2(timer, "Limiter - Neighbors");
      for (const auto face_no : cell->face_indices())
        {
          if (!cell->at_boundary(face_no) ||
              cell->has_periodic_neighbor(face_no))
            {
              auto neighbor = cell->neighbor_or_periodic_neighbor(face_no);
              const MHDEquations::state_type &neighbor_avg =
                get_cell_average(neighbor);
              const Tensor<1, spacedim> distance =
                SlopeLimiter::compute_periodic_distance_cell_neighbor(cell,
                                                                      face_no);

              saplog << "neighbor " << face_no << ": " << distance << std::endl;
              saplog << "avg: " << neighbor_avg << std::endl;

              for (unsigned int c = 0; c < MHDEquations::n_components; ++c)
                {
                  for (unsigned int d = 0; d < spacedim; ++d)
                    tmp_gradient[c][d] =
                      (neighbor_avg[c] - cell_avg[c]) / distance[d];

                  saplog << "neighbor gradient[" << c
                         << "]: " << tmp_gradient[c] << std::endl;
                }
              neighbor_gradients.push_back(tmp_gradient);
            }
        }
    }


    // Computed limited gradient
    {
      TimerOutput::Scope t2(timer, "Limiter - MinMod");
      const double diff = SlopeLimiter::minmod_gradients(cell_avg_gradient,
                                                         neighbor_gradients,
                                                         limited_gradient);
      saplog << "Difference between cell and limited gradient:" << diff
             << std::endl;
      for (unsigned int c = 0; c < MHDEquations::n_components; ++c)
        saplog << "limited_gradient[" << c << "]: " << limited_gradient[c]
               << std::endl;
    }


    {
      TimerOutput::Scope t2(timer, "Limiter - LimitSolution");
      // Compute limited solution values at q_points
      std::vector<Vector<double>> limited_solution_values(
        fe_v.n_quadrature_points, Vector<double>(MHDEquations::n_components));

      for (const unsigned int q_index : fe_v.quadrature_point_indices())
        {
          for (unsigned int c = 0; c < MHDEquations::n_components; ++c)
            {
              /** Use limited_gradient */
              limited_solution_values[q_index][c] =
                cell_avg[c] +
                limited_gradient[c] * (q_points[q_index] - cell->center());

              // /** Use cell average */
              // limited_solution_values[q_index][c] = cell_avg[c];
            }
        }
      // /** Reconstruct unlimited solution */
      // fe_v.get_function_values(locally_relevant_current_solution,
      //                          limited_solution_values);


      // Integrate with basis functions
      for (const unsigned int q_index : fe_v.quadrature_point_indices())
        for (unsigned int i : fe_v.dof_indices())
          for (unsigned int c = 0; c < MHDEquations::n_components; ++c)
            copy_data.cell_system_rhs(i) +=
              limited_solution_values[q_index][c] *
              fe_v.shape_value_component(i, q_index, c) * JxW[q_index];
    }


    saplog << "End cell" << std::endl << std::endl;
  };

  // copier for the mesh_loop function
  const auto copier = [&](const CopyDataSlopeLimiter &c) {
    constraints.distribute_local_to_global(c.cell_system_rhs,
                                           c.local_dof_indices,
                                           system_rhs);
  };

  ScratchDataSlopeLimiter<dim, spacedim> scratch_data(mapping, fe, quadrature);
  CopyDataSlopeLimiter                   copy_data;
  saplog << "Begin the assembly of the system rhs for projection." << std::endl;
  const auto filtered_iterator_range =
    dof_handler.active_cell_iterators() | IteratorFilters::LocallyOwnedCell();
  MeshWorker::mesh_loop(filtered_iterator_range,
                        cell_worker,
                        copier,
                        scratch_data,
                        copy_data,
                        MeshWorker::assemble_own_cells);
  system_rhs.compress(VectorOperation::add);
  saplog << "The system rhs was assembled." << std::endl;

  {
    TimerOutput::Scope t2(timer, "Limiter - Project");

    // Globally project limited solution onto locallY_owned_solution
    SolverControl solver_control(1000, 1e-6 * system_rhs.l2_norm());
    PETScWrappers::SolverGMRES solver(solver_control, mpi_communicator);

    PETScWrappers::PreconditionBlockJacobi preconditioner;
    preconditioner.initialize(mass_matrix);

    // The x vector has to be a vector of locally_owned_dofs, so we make use of
    // the locally_owned_solution
    solver.solve(mass_matrix,
                 locally_owned_solution,
                 system_rhs,
                 preconditioner);

    // Update the solution
    locally_relevant_current_solution = locally_owned_solution;

    saplog << "Solver converged in " << solver_control.last_step()
           << " iterations." << std::endl;
  }
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
      q_points.size(), Vector<double>(MHDEquations::n_components));
    typename MHDEquations::flux_type flux_matrix;

    fe_v.get_function_values(locally_relevant_current_solution, states);

    for (const unsigned int q_index : fe_v.quadrature_point_indices())
      {
        mhd_equations.compute_flux_matrix(states[q_index], flux_matrix);

        for (unsigned int i : fe_v.dof_indices())
          {
            // const unsigned int component_i =
            //   fe_v.get_fe().system_to_component_index(i).first;

            for (unsigned int c = 0; c < MHDEquations::n_components; ++c)
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
      q_points.size(), Vector<double>(MHDEquations::n_components));
    typename MHDEquations::state_type virtual_neighbor_state(
      MHDEquations::n_components);
    typename MHDEquations::state_type numerical_normal_flux(
      MHDEquations::n_components);

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

            for (unsigned int c = 0; c < MHDEquations::n_components; ++c)
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
      q_points.size(), Vector<double>(MHDEquations::n_components));
    std::vector<Vector<double>> states_neighbor(
      q_points.size(), Vector<double>(MHDEquations::n_components));
    typename MHDEquations::state_type numerical_normal_flux(
      MHDEquations::n_components);

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

            for (unsigned int c = 0; c < MHDEquations::n_components; ++c)
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

            for (unsigned int c = 0; c < MHDEquations::n_components; ++c)
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
      assemble_dg_rhs(time + gamma[i] * time_step);
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
  data_out.add_data_vector(
    locally_relevant_current_solution,
    MHDEquations::create_component_name_list(),
    DataOut<dim, spacedim>::type_dof_data,
    MHDEquations::create_component_interpretation_list());

  /** @todo [Remove Debug] */
  // Get one component of cell_averages
  Vector<double> cell_average_component(triangulation.n_active_cells());
  for (unsigned int i = 0; i < cell_average_component.size(); ++i)
    cell_average_component[i] = cell_average[i][0];

  data_out.add_data_vector(cell_average_component,
                           "average_roh",
                           DataOut<dim, spacedim>::type_cell_data);
  /** @todo [Remove Debug] */

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
const sapphirepp::MHD::MHDEquations &
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
