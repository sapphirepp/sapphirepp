#include "vfp-equation-solver.h"

#include <deal.II/grid/grid_in.h>

namespace Sapphire
{
  namespace VFP
  {
    using namespace dealii;

    // The mesh_loop function requires helper data types
    template <int dim_ps>
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
                  const UpdateFlags face_update_flags =
                    update_values | update_quadrature_points |
                    update_JxW_values | update_normal_vectors,
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
  } // namespace VFP
} // namespace Sapphire

Sapphire::VFP::VFPEquationSolver::VFPEquationSolver(
  const VFPSolverControl            &vfp_solver_control,
  const PhysicalProperties          &physical_properties,
  const Utils::OutputModule<dim_ps> &output_module)
  : vfp_solver_control(vfp_solver_control)
  , physical_properties(physical_properties)
  , mpi_communicator(MPI_COMM_WORLD)
  , n_mpi_procs(Utilities::MPI::n_mpi_processes(mpi_communicator))
  , rank(Utilities::MPI::this_mpi_process(mpi_communicator))
  , pcout(std::cout, (rank == 0))
  , output_module(output_module)
  , triangulation(mpi_communicator)
  , dof_handler(triangulation)
  , mapping()
  , fe(FE_DGQ<dim_ps>(vfp_solver_control.polynomial_degree),
       (vfp_solver_control.expansion_order + 1) *
         (vfp_solver_control.expansion_order + 1))
  , quadrature(fe.tensor_degree() + 1)
  , quadrature_face(fe.tensor_degree() + 1)
  , pde_system(vfp_solver_control.expansion_order)
  , upwind_flux(pde_system, vfp_solver_control, physical_properties)
  , timer(mpi_communicator, pcout, TimerOutput::never, TimerOutput::wall_times)
{
  // It should be checked if the results folder exists and if not it should
  // be tried to create it Only one processor needs to check and create the
  // folder if (rank == 0)
  //   if (!std::filesystem::exists(vfp_solver_control.results_path + "/" +
  //                                vfp_solver_control.simulation_id))
  //     std::filesystem::create_directories(vfp_solver_control.results_path
  //     +
  //                                         "/" +
  //                                         vfp_solver_control.simulation_id);
}

void
Sapphire::VFP::VFPEquationSolver::run()
{
  LogStream::Prefix p("VFP", saplog);
  make_grid();
  setup_system();
  assemble_mass_matrix();

  // Project the initial values
  InitialValueFunction<dim_ps> initial_value_function(physical_properties,
                                                      expansion_order);
  PETScWrappers::MPI::Vector   initial_condition(locally_owned_dofs,
                                               mpi_communicator);
  project(initial_value_function, initial_condition);
  // Here a non ghosted vector, is copied into a ghosted vector. I think
  // that is the moment where the ghost cells are filled.
  locally_relevant_current_solution = initial_condition;
  // Output t = 0
  std::vector<XDMFEntry> xdmf_entries;
  output_results(0);

  // Assemble the dg matrix for t = 0
  assemble_dg_matrix(0);
  // Source term at t = 0;
  if constexpr ((flags & TermFlags::source) != TermFlags::none)
    {
      Source<dim_ps> source_function(physical_properties, expansion_order);
      source_function.set_time(0);
      compute_source_term(source_function);
    }

  double time_step  = vfp_solver_control.time_step;
  double final_time = vfp_solver_control.final_time;
  saplog << "The time stepping loop is entered:" << std::endl;
  unsigned int time_step_number = 1;
  for (double time = 0.; time < final_time;
       time += time_step, ++time_step_number)
    {
      saplog << "Time step " << time_step_number << " at t = " << time
             << std::endl;
      // Time stepping method
      if (vfp_solver_control.time_stepping_method ==
            Utils::TimeSteppingMethod::forward_euler ||
          vfp_solver_control.time_stepping_method ==
            Utils::TimeSteppingMethod::backward_euler ||
          vfp_solver_control.time_stepping_method ==
            Utils::TimeSteppingMethod::crank_nicolson)
        theta_method(time, time_step);
      else if (vfp_solver_control.time_stepping_method ==
               Utils::TimeSteppingMethod::erk4)
        explicit_runge_kutta(time, time_step);
      else if (vfp_solver_control.time_stepping_method ==
               Utils::TimeSteppingMethod::lserk)
        low_storage_explicit_runge_kutta(time, time_step);

      output_results(time_step_number);
    }
  saplog << "The simulation ended." << std::endl;

  timer.print_summary();
  timer.reset();
}

void
Sapphire::VFP::VFPEquationSolver::make_grid()
{
  TimerOutput::Scope timer_section(timer, "Grid setup");
  saplog << "Create the grid" << std::endl;
  // Colorise = true means to set boundary ids (default for 1D)
  bool colorise = vfp_solver_control.periodicity[0] ||
                  vfp_solver_control.periodicity[1] ||
                  vfp_solver_control.periodicity[2];

  switch (vfp_solver_control.grid_type)
    {
      case Utils::GridType::hypercube:
        {
          saplog << "Create the grid from hyper rectangle" << std::endl;
          GridGenerator::subdivided_hyper_rectangle(triangulation,
                                                    vfp_solver_control.n_cells,
                                                    vfp_solver_control.p1,
                                                    vfp_solver_control.p2,
                                                    colorise);
          break;
        }
      case Utils::GridType::file:
        {
          saplog << "Read grid from file \"" << vfp_solver_control.grid_file
                 << "\"" << std::endl;
          GridIn<dim_ps> grid_in(triangulation);
          std::ifstream  input(vfp_solver_control.grid_file);
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

  // GridGenerator::hyper_cube(triangulation, -5., 5., colorise);
  // triangulation.refine_global(6);
  saplog << "The grid was created: "
         << "	#cells=" << triangulation.n_cells()
         << ",	#active cells=" << triangulation.n_global_active_cells()
         << std::endl;

  if (colorise)
    {
      saplog << "Set up periodic boundary conditions" << std::endl;
      // Periodic boundary conditions with MeshWorker. Mailinglist
      // https://groups.google.com/g/dealii/c/WlOiww5UVxc/m/mtQJDUwiBQAJ
      //
      // "If you call add_periodicity() on a Triangulation object, the
      // periodic faces are treated as internal faces in MeshWorker. This
      // means that you will not access them in a "integrate_boundary_term"
      // function but in a "integrate_face_term" function. "
      std::vector<GridTools::PeriodicFacePair<
        typename Triangulation<dim_ps>::cell_iterator>>
        matched_pairs;

      // Fill the matched_pairs vector manually if dim_ps = 1. At the moment
      // there is no instance of the template collect_period_faces for
      // MeshType = parallel::shared::Triangulation<1,1>.
      // https://github.com/dealii/dealii/issues/14879
      if constexpr (dim_ps == 1)
        {
          if (vfp_solver_control.periodicity[0])
            {
              matched_pairs.resize(1);
              matched_pairs[0].cell[0]     = triangulation.begin();
              matched_pairs[0].cell[1]     = triangulation.last();
              matched_pairs[0].face_idx[0] = 0;
              matched_pairs[0].face_idx[1] = 1;
              std::bitset<3> temp_bitset;
              temp_bitset[0]               = true;
              matched_pairs[0].orientation = temp_bitset;
            }
        }
      else
        {
          for (unsigned int i = 0; i < dim_ps; ++i)
            {
              if (vfp_solver_control.periodicity[i])
                GridTools::collect_periodic_faces(
                  triangulation, 2 * i, 2 * i + 1, i, matched_pairs);
            }
        }
      triangulation.add_periodicity(matched_pairs);
    }
}

void
Sapphire::VFP::VFPEquationSolver::setup_system()
{
  TimerOutput::Scope timer_section(timer, "FE system");
  saplog << "Setup the finite element system" << std::endl;

  dof_handler.clear();
  dof_handler.distribute_dofs(fe);

  const unsigned int n_dofs = dof_handler.n_dofs();
  saplog << "The degrees of freedom were distributed:"
         << "	n_dofs=" << n_dofs << std::endl;

  locally_owned_dofs    = dof_handler.locally_owned_dofs();
  locally_relevant_dofs = DoFTools::extract_locally_relevant_dofs(dof_handler);

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

  DynamicSparsityPattern dsp(locally_relevant_dofs);
  // NON-PERIODIC
  DoFTools::make_flux_sparsity_pattern(dof_handler, dsp);
  SparsityTools::distribute_sparsity_pattern(dsp,
                                             locally_owned_dofs,
                                             mpi_communicator,
                                             locally_relevant_dofs);
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
  mass_matrix.reinit(locally_owned_dofs,
                     locally_owned_dofs,
                     dsp,
                     mpi_communicator);
  system_matrix.reinit(locally_owned_dofs,
                       locally_owned_dofs,
                       dsp,
                       mpi_communicator);
}

void
Sapphire::VFP::VFPEquationSolver::assemble_mass_matrix()
{
  TimerOutput::Scope timer_section(timer, "Mass matrix");

  using Iterator = typename DoFHandler<dim_ps>::active_cell_iterator;

  const auto cell_worker = [](const Iterator      &cell,
                              ScratchData<dim_ps> &scratch_data,
                              CopyData            &copy_data) {
    FEValues<dim_ps> &fe_v = scratch_data.fe_values;
    fe_v.reinit(cell);
    const unsigned int n_dofs = fe_v.get_fe().n_dofs_per_cell();
    // reinit the cell_matrix
    copy_data.reinit(cell, n_dofs);
    const std::vector<double> &JxW = fe_v.get_JxW_values();

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
                // mass matrix
                copy_data.cell_matrix(i, j) +=
                  (component_i == component_j ? fe_v.shape_value(i, q_index) *
                                                  fe_v.shape_value(j, q_index) :
                                                0) *
                  JxW[q_index];
              }
          }
      }
  };

  const auto copier = [&](const CopyData &c) {
    constraints.distribute_local_to_global(c.cell_matrix,
                                           c.local_dof_indices,
                                           mass_matrix);
  };
  ScratchData<dim_ps> scratch_data(mapping, fe, quadrature, quadrature_face);
  CopyData            copy_data;
  // Perfom the integration loop only over the locally owned cells
  const auto filtered_iterator_range =
    dof_handler.active_cell_iterators() | IteratorFilters::LocallyOwnedCell();

  saplog << "Begin the assembly of the mass matrix." << std::endl;
  MeshWorker::mesh_loop(filtered_iterator_range,
                        cell_worker,
                        copier,
                        scratch_data,
                        copy_data,
                        MeshWorker::assemble_own_cells);
  mass_matrix.compress(VectorOperation::add);
  saplog << "The mass matrix was assembled." << std::endl;
}

void
Sapphire::VFP::VFPEquationSolver::assemble_dg_matrix(const double time)
{
  TimerOutput::Scope timer_section(timer, "DG matrix");
  /*
    What kind of loops are there ?
    1. Loop over all cells (this happens inside the mesh_loop)
    2. Loop over the degrees of freedom on each cell
    - the metod system_to_componet_index() returns the index of the non-zero
    component of the vector-valued shape function which corresponds to the
    indices (l,m,s)
  */
  using Iterator = typename DoFHandler<dim_ps>::active_cell_iterator;
  upwind_flux.set_time(time);

  const std::vector<LAPACKFullMatrix<double>> &advection_matrices =
    pde_system.get_advection_matrices();
  const std::vector<LAPACKFullMatrix<double>> &generator_rotation_matrices =
    pde_system.get_generator_matrices();
  const std::vector<LAPACKFullMatrix<double>> &adv_mat_products =
    pde_system.get_adv_mat_products();
  const std::vector<LAPACKFullMatrix<double>> &adv_x_gen_matrices =
    pde_system.get_adv_cross_gen();
  const std::vector<LAPACKFullMatrix<double>> &t_matrices =
    pde_system.get_t_matrices();

  // Collision term (essentially a reaction term)
  const Vector<double> &collision_matrix = pde_system.get_collision_matrix();

  BackgroundVelocityField<dim_ps> background_velocity_field(
    physical_properties);
  background_velocity_field.set_time(time);
  MagneticField<dim_ps> magnetic_field(physical_properties);
  magnetic_field.set_time(time);
  ScatteringFrequency<dim_ps> scattering_frequency(physical_properties);
  scattering_frequency.set_time(time);

  ParticleVelocity<dim_ps> particle_velocity;
  ParticleGamma<dim_ps>    particle_gamma;

  // For the transport only case, the energy, the Lorentz factor and the
  // velocity are defined in TransportOnly struct
  TransportOnly transport_only;

  ParticleProperties particle_properties;

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
    std::vector<std::vector<Vector<double>>> jacobians_vel(
      q_points.size(), std::vector<Vector<double>>(3, Vector<double>(3)));
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
                if constexpr ((flags & TermFlags::collision) != TermFlags::none)
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
                // dim_cs must euqal zero, i.e. the distributin function is
                // assumend to be homogeneous (it does not depend on x, y
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
                        // NOTE: Many zerso are added here, because the
                        // matrices Ax, Ay and A_z are sparse. TODO:
                        // Performance check. If too bad, return to the
                        // strategy, which was used in v0.6.5
                        if ((flags & TermFlags::momentum) != TermFlags::none)
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
                              transport_only.velocity *
                              advection_matrices[coordinate](component_i,
                                                             component_j) *
                              fe_v.shape_value(j, q_index) * JxW[q_index];
                          }
                      }
                  }
                if constexpr ((flags & TermFlags::magnetic) != TermFlags::none)
                  {
                    // NOTE: All three components of the B-Field are
                    // included no matter, which dimension of the
                    // configuration space is considered
                    for (unsigned int coordinate = 0; coordinate < 3;
                         ++coordinate)
                      {
                        if ((flags & TermFlags::momentum) != TermFlags::none)
                          {
                            copy_data.cell_matrix(i, j) -=
                              fe_v.shape_value(i, q_index) *
                              particle_properties.charge *
                              magnetic_field_values[q_index][coordinate] /
                              (particle_gammas[q_index] *
                               particle_properties.mass) *
                              generator_rotation_matrices[coordinate](
                                component_i, component_j) *
                              fe_v.shape_value(j, q_index) * JxW[q_index];
                          }
                        else
                          {
                            // fixed energy case (i.e. transport only)
                            copy_data.cell_matrix(i, j) -=
                              fe_v.shape_value(i, q_index) *
                              particle_properties.charge *
                              magnetic_field_values[q_index][coordinate] /
                              (transport_only.gamma *
                               particle_properties.mass) *
                              generator_rotation_matrices[coordinate](
                                component_i, component_j) *
                              fe_v.shape_value(j, q_index) * JxW[q_index];
                          }
                      }
                  }
                if constexpr ((flags & TermFlags::momentum) != TermFlags::none)
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
                              particle_properties.mass /
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
                              particle_gammas[q_index] *
                              particle_properties.mass *
                              material_derivative_vel[q_index][coordinate] *
                              advection_matrices[coordinate](component_i,
                                                             component_j) *
                              fe_v.shape_value(j, q_index) * JxW[q_index];
                            // \phi m * v * du^k/dt * A_k * \phi
                            copy_data.cell_matrix(i, j) +=
                              fe_v.shape_value(i, q_index) *
                              particle_properties.mass *
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
    // NOTE: copy_data is not reinitialised, the cell_workers contribution
    // to the cell_dg_matrix should not be deleted

    const std::vector<Point<dim_ps>> &q_points =
      fe_face_v.get_quadrature_points();
    const std::vector<double>            &JxW = fe_face_v.get_JxW_values();
    const std::vector<Tensor<1, dim_ps>> &normals =
      fe_face_v.get_normal_vectors();
    std::vector<FullMatrix<double>> positive_flux_matrices(
      q_points.size(), FullMatrix<double>(num_exp_coefficients));
    std::vector<FullMatrix<double>> negative_flux_matrices(
      q_points.size(), FullMatrix<double>(num_exp_coefficients));

    upwind_flux.compute_upwind_fluxes(q_points,
                                      normals,
                                      positive_flux_matrices,
                                      negative_flux_matrices);
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
                if constexpr ((flags & TermFlags::spatial_advection) !=
                              TermFlags::none)
                  {
                    // Outflow boundary: Everyhing with a positive flux
                    // along the direction of the normal leaves the boundary
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
                  }
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
    // Create an element at the end of the vector containig the face data
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
    std::vector<FullMatrix<double>> positive_flux_matrices(
      q_points.size(), FullMatrix<double>(num_exp_coefficients));
    std::vector<FullMatrix<double>> negative_flux_matrices(
      q_points.size(), FullMatrix<double>(num_exp_coefficients));

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
                                           c.local_dof_indices,
                                           dg_matrix);
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
  saplog << "The DG matrix was assembled." << std::endl;
}

template <int dim>
void
Sapphire::VFP::VFPEquationSolver::project(
  const Function<dim>        &f,
  PETScWrappers::MPI::Vector &projected_function)
{
  TimerOutput::Scope timer_section(timer, "Project f onto the FEM space");
  saplog << "Project a function onto the finite element space" << std::endl;
  LogStream::Prefix p("project", saplog);
  // Create right hand side
  FEValues<dim_ps> fe_v(mapping,
                        fe,
                        quadrature,
                        update_values | update_quadrature_points |
                          update_JxW_values);

  const unsigned int n_dofs = fe.n_dofs_per_cell();
  Vector<double>     cell_rhs(n_dofs);

  std::vector<types::global_dof_index> local_dof_indices(n_dofs);

  PETScWrappers::MPI::Vector rhs(locally_owned_dofs, mpi_communicator);
  for (const auto &cell : dof_handler.active_cell_iterators())
    {
      if (cell->is_locally_owned())
        {
          cell_rhs = 0;
          fe_v.reinit(cell);

          const std::vector<Point<dim_ps>> &q_points =
            fe_v.get_quadrature_points();
          const std::vector<double> &JxW = fe_v.get_JxW_values();

          std::vector<Vector<double>> function_values(
            q_points.size(), Vector<double>(f.n_components));
          f.vector_value_list(q_points, function_values);

          for (const unsigned int q_index : fe_v.quadrature_point_indices())
            {
              for (unsigned int i : fe_v.dof_indices())
                {
                  const unsigned int component_i =
                    fe.system_to_component_index(i).first;
                  cell_rhs(i) += fe_v.shape_value(i, q_index) *
                                 function_values[q_index][component_i] *
                                 JxW[q_index];
                }
            }
          cell->get_dof_indices(local_dof_indices);

          constraints.distribute_local_to_global(cell_rhs,
                                                 local_dof_indices,
                                                 rhs);
        }
    }
  rhs.compress(VectorOperation::add);

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

template <int dim>
void
Sapphire::VFP::VFPEquationSolver::compute_source_term(
  const Function<dim> &source_function)
{
  FEValues<dim_ps> fe_v(mapping,
                        fe,
                        quadrature,
                        update_values | update_quadrature_points |
                          update_JxW_values);

  const unsigned int n_dofs = fe.n_dofs_per_cell();
  Vector<double>     cell_rhs(n_dofs);

  std::vector<types::global_dof_index> local_dof_indices(n_dofs);
  // reinitialis the source term
  locally_owned_current_source.reinit(locally_owned_dofs, mpi_communicator);
  for (const auto &cell : dof_handler.active_cell_iterators())
    {
      if (cell->is_locally_owned())
        {
          cell_rhs = 0;
          fe_v.reinit(cell);

          const std::vector<Point<dim_ps>> &q_points =
            fe_v.get_quadrature_points();
          const std::vector<double> &JxW = fe_v.get_JxW_values();

          std::vector<Vector<double>> source_values(
            q_points.size(), Vector<double>(num_exp_coefficients));
          source_function.vector_value_list(q_points, source_values);

          for (const unsigned int q_index : fe_v.quadrature_point_indices())
            {
              for (unsigned int i : fe_v.dof_indices())
                {
                  const unsigned int component_i =
                    fe.system_to_component_index(i).first;
                  cell_rhs(i) += fe_v.shape_value(i, q_index) *
                                 source_values[q_index][component_i] *
                                 JxW[q_index];
                }
            }
          cell->get_dof_indices(local_dof_indices);

          constraints.distribute_local_to_global(cell_rhs,
                                                 local_dof_indices,
                                                 locally_owned_current_source);
        }
    }
  locally_owned_current_source.compress(VectorOperation::add);
}

void
Sapphire::VFP::VFPEquationSolver::theta_method(const double time,
                                               const double time_step)
{
  TimerOutput::Scope timer_section(timer, "Theta method");
  LogStream::Prefix  p("theta_method", saplog);
  // Equation: (mass_matrix + time_step * theta * dg_matrix(time +
  // time_step)) f(time + time_step) = (mass_matrix - time_step * (1 -
  // theta) * dg_matrix(time) ) f(time) + time_step * theta * s(time +
  // time_step) + time_step * (1 - theta) * s(time)
  const double theta = vfp_solver_control.theta;

  locally_owned_previous_solution = locally_relevant_current_solution;

  mass_matrix.vmult(system_rhs, locally_owned_previous_solution);
  PETScWrappers::MPI::Vector tmp(locally_owned_dofs, mpi_communicator);
  dg_matrix.vmult(tmp, locally_owned_previous_solution);
  system_rhs.add(-time_step * (1 - theta), tmp);
  // Source term
  if constexpr ((flags & TermFlags::source) != TermFlags::none)
    {
      if constexpr (time_dependent_source)
        {
          system_rhs.add((1 - theta) * time_step, locally_owned_current_source);
          // Update the source term
          Source<dim_ps> source_function(physical_properties, expansion_order);
          source_function.set_time(time + time_step);
          compute_source_term(source_function);
          system_rhs.add(theta * time_step, locally_owned_current_source);
        }
      else
        {
          system_rhs.add(time_step, locally_owned_current_source);
        }
    }
  // Since the the dg_matrix depends on the velocity field (and/or the
  // magnetice field) and the velocity field may depend on time, it needs to
  // reassembled every time step. This is not true for the mass matrix ( but
  // it may if the grid adapts after a specified amount of time steps)
  if constexpr (time_dependent_fields)
    {
      dg_matrix = 0;
      assemble_dg_matrix(time + time_step);
    }

  system_matrix.copy_from(mass_matrix);
  system_matrix.add(time_step * theta, dg_matrix);

  SolverControl                   solver_control(1000, 1e-12);
  PETScWrappers::SolverRichardson solver(solver_control, mpi_communicator);

  // PETScWrappers::PreconditionBoomerAMG preconditioner;
  // PETScWrappers::PreconditionBoomerAMG::AdditionalData data;
  // data.symmetric_operator = true;
  // preconditioner.initialize(system_matrix, data);

  // PETScWrappers::PreconditionSOR preconditioner;
  // preconditioner.initialize(system_matrix);

  PETScWrappers::PreconditionBlockJacobi preconditioner;
  preconditioner.initialize(system_matrix);

  solver.solve(system_matrix,
               locally_owned_previous_solution,
               system_rhs,
               preconditioner);

  // Update the solution
  locally_relevant_current_solution = locally_owned_previous_solution;

  saplog << "Solver converged in " << solver_control.last_step()
         << " iterations." << std::endl;
}

void
Sapphire::VFP::VFPEquationSolver::explicit_runge_kutta(const double time,
                                                       const double time_step)
{
  TimerOutput::Scope timer_section(timer, "ERK4");
  LogStream::Prefix  p("ERK4", saplog);
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

  SolverControl           solver_control(1000, 1e-12);
  PETScWrappers::SolverCG cg(solver_control, mpi_communicator);

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
  if constexpr ((flags & TermFlags::source) != TermFlags::none)
    {
      system_rhs.add(-1., locally_owned_current_source);
    }
  cg.solve(mass_matrix, k_0, system_rhs, preconditioner);
  saplog << "Stage s: " << 0 << "	Solver converged in "
         << solver_control.last_step() << " iterations." << std::endl;
  k_0 *= -1.;
  locally_owned_current_solution.add(b[0] * time_step, k_0);

  // k_1
  PETScWrappers::MPI::Vector k_1(locally_owned_dofs, mpi_communicator);
  // Comute dg_matrix(time + c[1] * time_step) if the fields are time
  // dependent
  if constexpr (time_dependent_fields)
    {
      dg_matrix = 0;
      assemble_dg_matrix(time + c[1] * time_step);
    }
  temp.add(1., locally_owned_previous_solution, a[0] * time_step, k_0);
  dg_matrix.vmult(system_rhs, temp);
  if constexpr ((flags & TermFlags::source) != TermFlags::none)
    {
      if constexpr (time_dependent_source)
        {
          Source<dim_ps> source_function(physical_properties, expansion_order);
          source_function.set_time(time + c[1] * time_step);
          compute_source_term(source_function);
          system_rhs.add(-1., locally_owned_current_source);
        }
      else
        {
          system_rhs.add(-1., locally_owned_current_source);
        }
    }
  cg.solve(mass_matrix, k_1, system_rhs, preconditioner);
  saplog << "	Stage s: " << 1 << "	Solver converged in "
         << solver_control.last_step() << " iterations." << std::endl;
  k_1 *= -1.;
  locally_owned_current_solution.add(b[1] * time_step, k_1);
  temp = 0;

  // k_2
  PETScWrappers::MPI::Vector k_2(locally_owned_dofs, mpi_communicator);
  // NOTE: For k_2 it is not necessary to reassamble the dg_matrix,
  // since c[1] = c[2]
  temp.add(1., locally_owned_previous_solution, a[1] * time_step, k_1);
  dg_matrix.vmult(system_rhs, temp);
  if constexpr ((flags & TermFlags::source) != TermFlags::none)
    {
      if constexpr (time_dependent_source)
        {
          Source<dim_ps> source_function(physical_properties, expansion_order);
          source_function.set_time(time + c[2] * time_step);
          compute_source_term(source_function);
          system_rhs.add(-1., locally_owned_current_source);
        }
      else
        {
          system_rhs.add(-1., locally_owned_current_source);
        }
    }
  cg.solve(mass_matrix, k_2, system_rhs, preconditioner);
  saplog << "	Stage s: " << 2 << "	Solver converged in "
         << solver_control.last_step() << " iterations." << std::endl;
  k_2 *= -1.;
  locally_owned_current_solution.add(b[2] * time_step, k_2);
  temp = 0;

  // k_3
  PETScWrappers::MPI::Vector k_3(locally_owned_dofs, mpi_communicator);
  // Comute dg_matrix(time + c[1] * time_step) if the fields are time
  // dependent
  if constexpr (time_dependent_fields)
    {
      dg_matrix = 0;
      assemble_dg_matrix(time + c[3] * time_step);
    }
  temp.add(1., locally_owned_previous_solution, a[2] * time_step, k_2);
  dg_matrix.vmult(system_rhs, temp);
  if constexpr ((flags & TermFlags::source) != TermFlags::none)
    {
      if constexpr (time_dependent_source)
        {
          Source<dim_ps> source_function(physical_properties, expansion_order);
          source_function.set_time(time + c[3] * time_step);
          compute_source_term(source_function);
          system_rhs.add(-1., locally_owned_current_source);
        }
      else
        {
          system_rhs.add(-1., locally_owned_current_source);
        }
    }
  cg.solve(mass_matrix, k_3, system_rhs, preconditioner);
  saplog << "	Stage s: " << 3 << "	Solver converged in "
         << solver_control.last_step() << " iterations." << std::endl;
  k_3 *= -1.;
  locally_owned_current_solution.add(b[3] * time_step, k_3);
  temp = 0;

  locally_relevant_current_solution = locally_owned_current_solution;
}

void
Sapphire::VFP::VFPEquationSolver::low_storage_explicit_runge_kutta(
  const double time,
  const double time_step)
{
  TimerOutput::Scope timer_section(timer, "LSERK");
  LogStream::Prefix  p("LSERK", saplog);
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
  // NOTE: I only need c if the velocity field and the magnetic field are
  // time dependent
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

  SolverControl           solver_control(1000, 1e-12);
  PETScWrappers::SolverCG cg(solver_control, mpi_communicator);
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
      if constexpr (time_dependent_fields)
        assemble_dg_matrix(time + c[s] * time_step);

      dg_matrix.vmult(system_rhs, locally_owned_previous_solution);
      if constexpr (time_dependent_fields)
        dg_matrix = 0;

      if constexpr ((flags & TermFlags::source) != TermFlags::none)
        {
          if constexpr (time_dependent_source)
            {
              Source<dim_ps> source_function(physical_properties,
                                             expansion_order);
              source_function.set_time(time + c[s] * time_step);
              compute_source_term(source_function);
              system_rhs.add(-1., locally_owned_current_source);
            }
          else
            {
              system_rhs.add(-1., locally_owned_current_source);
            }
        }

      cg.solve(mass_matrix, temp, system_rhs, preconditioner);
      saplog << "	Stage s: " << s << "	Solver converged in "
             << solver_control.last_step() << " iterations." << std::endl;

      k.sadd(a[s], -time_step, temp);

      locally_owned_previous_solution.add(b[s], k);
    }
  // Currently I assume that there are no constraints
  // constraints.distribute(locally_relevant_current_solution);
  locally_relevant_current_solution = locally_owned_previous_solution;
}

void
Sapphire::VFP::VFPEquationSolver::output_results(
  const unsigned int time_step_number)
{
  TimerOutput::Scope timer_section(timer, "Output");
  DataOut<dim_ps>    data_out;
  data_out.attach_dof_handler(dof_handler);
  // Create a vector of strings with names for the components of the
  // solution
  std::vector<std::string> component_names(num_exp_coefficients);
  const std::vector<std::array<unsigned int, 3>> &lms_indices =
    pde_system.get_lms_indices();

  for (unsigned int i = 0; i < num_exp_coefficients; ++i)
    {
      const std::array<unsigned int, 3> &lms = lms_indices[i];
      component_names[i]                     = "f_" + std::to_string(lms[0]) +
                           std::to_string(lms[1]) + std::to_string(lms[2]);
    }

  data_out.add_data_vector(locally_relevant_current_solution, component_names);

  // Output the partition of the mesh
  Vector<float> subdomain(triangulation.n_active_cells());
  for (unsigned int i = 0; i < subdomain.size(); ++i)
    subdomain(i) = triangulation.locally_owned_subdomain();
  data_out.add_data_vector(subdomain, "subdomain");

  // Adapt the output to the polynomial degree of the shape functions
  data_out.build_patches(vfp_solver_control.polynomial_degree);
  output_module.write_results(data_out, time_step_number);
}
