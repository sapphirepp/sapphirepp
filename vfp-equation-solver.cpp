#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/base/function.h>
#include <deal.II/base/index_set.h>
#include <deal.II/base/mpi.h>
#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/patterns.h>
#include <deal.II/base/quadrature.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/tensor_function.h>
#include <deal.II/base/timer.h>
#include <deal.II/base/types.h>
#include <deal.II/base/utilities.h>
#include <deal.II/base/vectorization.h>
#include <deal.II/fe/fe_update_flags.h>
#include <deal.II/grid/filtered_iterator.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools.h>  // neeed for periodic boundary conditions
                                      // (collect_periodic_faces) and
                                      // partitioning of the distributed mesh
// #include <deal.II/grid/grid_refinement.h>
// #include <deal.II/grid/tria.h>
#include <deal.II/distributed/tria.h>
// #include <deal.II/distributed/grid_refinement.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/lapack_full_matrix.h>
// PetscWrappers
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/lac/petsc_precondition.h>
#include <deal.II/lac/petsc_solver.h>
#include <deal.II/lac/petsc_sparse_matrix.h>
#include <deal.II/lac/petsc_vector.h>
#include <deal.II/lac/sparsity_tools.h>  // need for distribute_sparsity_pattern
                                         // to fit into a global numbering
#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_q1.h>
#include <deal.II/lac/vector_operation.h>
#include <deal.II/meshworker/assemble_flags.h>
#include <deal.II/meshworker/mesh_loop.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/vector_tools_project.h>
#include <mpi.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iostream>
#include <numeric>
#include <string>
#include <vector>

// own header files
#include "compile-time-flags.h"
#include "particle-functions.h"
#include "pde-system.h"
#include "physical-setup.h"
#include "upwind-flux.h"
#include "vfp-solver-control.h"

namespace VFPEquation {
using namespace dealii;

// Initial values
template <int dim_cs, bool momentum>
class InitialValueFunction : public Function<dim_cs + momentum> {
 public:
  InitialValueFunction(unsigned int exp_order) : expansion_order{exp_order} {}
  virtual void vector_value(const Point<dim_cs + momentum> &p,
                            Vector<double> &values) const override {
    Assert(dim_cs <= 3, ExcNotImplemented());
    Assert(values.size() == (expansion_order + 1) * (expansion_order + 1),
           ExcDimensionMismatch(values.size(),
                                (expansion_order + 1) * (expansion_order + 1)));
    // The zeroth component of values corresponds to f_000, the first component
    // to f_110 etc.
    if constexpr (dim_cs == 1) {
      // values[0] = 1. * std::exp(-(std::pow(p[0], 2)));
      // if (std::abs(p[0]) < 1.) values[0] = 1.;
      // values[0] = 1. + 0.2 * p[0];
      // constant
      values[0] = 1.;
    }
    // values[0] = std::sin((1. * 3.14159265359) / 2 * p[0]) + 1.;

    if constexpr (dim_cs == 2) {
      // Gaussian
      values[0] = 1. * std::exp(-((std::pow(p[0], 2) + std::pow(p[1], 2))));

      // constant
      // values[0] = 1.;

      // constant disc
      // if (p.norm() <= 1.) values[0] = 1.;

      // Rigid rotator
      // if (std::abs(p[0]) <= 3 && std::abs(p[1]) <= 0.5) values[0] = 1.;
      // if (std::abs(p[0]) <= 0.5 && std::abs(p[1]) <= 3) values[0] = 1.;
    }
    if constexpr (dim_cs == 3)
      values[0] = 1. * std::exp(-((std::pow(p[0], 2) + std::pow(p[1], 2) +
                                   std::pow(p[2], 2))));

    // Fill all components with the same values
    // std::fill(
    //     values.begin(), values.end(),
    //     1. * std::exp(-((std::pow(p[0] - 0.5, 2) + std::pow(p[1] - 0.5, 2)) /
    //                     0.01)));
    if constexpr (momentum)
      values[0] *=
          std::exp(-(std::pow(p[dim_cs + momentum - 1] - 5.5, 2) / 0.5));
  }

 private:
  const unsigned int expansion_order = 1;
};

// The mesh_loop function requires helper data types
template <int dim_ps>
class ScratchData {
 public:
  // Constructor
  ScratchData(const Mapping<dim_ps> &mapping, const FiniteElement<dim_ps> &fe,
              const Quadrature<dim_ps> &quadrature,
              const Quadrature<dim_ps - 1> &quadrature_face,
              const UpdateFlags update_flags = update_values |
                                               update_gradients |
                                               update_quadrature_points |
                                               update_JxW_values,
              const UpdateFlags face_update_flags = update_values |
                                                    update_quadrature_points |
                                                    update_JxW_values |
                                                    update_normal_vectors,
              const UpdateFlags neighbor_face_update_flags = update_values)
      : fe_values(mapping, fe, quadrature, update_flags),
        fe_values_face(mapping, fe, quadrature_face, face_update_flags),
        fe_values_face_neighbor(mapping, fe, quadrature_face,
                                neighbor_face_update_flags),
        fe_values_subface_neighbor(mapping, fe, quadrature_face,
                                   neighbor_face_update_flags) {}

  // Copy Constructor
  ScratchData(const ScratchData<dim_ps> &scratch_data)
      : fe_values(scratch_data.fe_values.get_mapping(),
                  scratch_data.fe_values.get_fe(),
                  scratch_data.fe_values.get_quadrature(),
                  scratch_data.fe_values.get_update_flags()),
        fe_values_face(scratch_data.fe_values_face.get_mapping(),
                       scratch_data.fe_values_face.get_fe(),
                       scratch_data.fe_values_face.get_quadrature(),
                       scratch_data.fe_values_face.get_update_flags()),
        fe_values_face_neighbor(
            scratch_data.fe_values_face_neighbor.get_mapping(),
            scratch_data.fe_values_face_neighbor.get_fe(),
            scratch_data.fe_values_face_neighbor.get_quadrature(),
            scratch_data.fe_values_face_neighbor.get_update_flags()),
        fe_values_subface_neighbor(
            scratch_data.fe_values_subface_neighbor.get_mapping(),
            scratch_data.fe_values_subface_neighbor.get_fe(),
            scratch_data.fe_values_subface_neighbor.get_quadrature(),
            scratch_data.fe_values_subface_neighbor.get_update_flags()) {}

  FEValues<dim_ps> fe_values;
  FEFaceValues<dim_ps> fe_values_face;
  FEFaceValues<dim_ps> fe_values_face_neighbor;
  FESubfaceValues<dim_ps> fe_values_subface_neighbor;
};

struct CopyDataFace {
  FullMatrix<double> cell_dg_matrix_11;
  FullMatrix<double> cell_dg_matrix_12;
  FullMatrix<double> cell_dg_matrix_21;
  FullMatrix<double> cell_dg_matrix_22;

  std::vector<types::global_dof_index> local_dof_indices;
  std::vector<types::global_dof_index> local_dof_indices_neighbor;

  template <typename Iterator>
  void reinit(const Iterator &cell, const Iterator &neighbor_cell,
              unsigned int dofs_per_cell) {
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

struct CopyData {
  FullMatrix<double> cell_matrix;
  Vector<double> cell_rhs;
  std::vector<types::global_dof_index> local_dof_indices;
  std::vector<types::global_dof_index> local_dof_indices_neighbor;
  std::vector<CopyDataFace> face_data;

  template <typename Iterator>
  void reinit(const Iterator &cell, unsigned int dofs_per_cell) {
    cell_matrix.reinit(dofs_per_cell, dofs_per_cell);
    cell_rhs.reinit(dofs_per_cell);

    local_dof_indices.resize(dofs_per_cell);
    cell->get_dof_indices(local_dof_indices);
  }
};

class VFPEquationSolver {
 public:
  VFPEquationSolver(const VFPSolverControl &control);
  void run();

 private:
  // VFPSolverControl
  const VFPSolverControl &vfp_solver_control;

  static constexpr int dim_ps = VFPSolverControl::dim;
  static constexpr int dim_cs = VFPSolverControl::dim_configuration_space;
  static constexpr TermFlags flags = VFPSolverControl::terms;
  static constexpr bool time_dependent_fields =
      VFPSolverControl::time_dependent_fields;
  // ((flags & TermFlags::momentum) != TermFlags::none) ? dim_cs + 1 : dim_cs;

  // Triangulation
  void make_grid();
  // Setup data structures for the linear system
  void setup_system();
  // Matrix assembly
  void assemble_mass_matrix();
  void assemble_dg_matrix(const double time);
  // Time stepping methods
  void theta_method_solve_system();
  void theta_method(double theta);
  void explicit_runge_kutta();
  void low_storage_explicit_runge_kutta(const double time,
                                        const double time_step);
  // Output
  void output_results(const unsigned int time_step_number) const;

  // auxiliary functions
  template <int dim>
  void project(const Function<dim> &f,
               PETScWrappers::MPI::Vector &projected_function);

  MPI_Comm mpi_communicator;
  const unsigned int n_mpi_procs;
  const unsigned int rank;

  ConditionalOStream pcout;

  // NOTE: Parallel distribution does not allow 1D triangulation. This excludes
  // the 1D transport (i.e. no momentum terms) only case. To reimplement this
  // case, take a look at Step-16, Step-17 how to deal with copys of
  // triangulations and dof handler on every mpi process.
  parallel::distributed::Triangulation<dim_ps> triangulation;
  DoFHandler<dim_ps> dof_handler;

  IndexSet locally_owned_dofs;
  IndexSet locally_relevant_dofs;

  // NOTE: The explicit use of a mapping is most likely related to the usage
  // of mesh_loop as well
  const MappingQ1<dim_ps> mapping;
  const FESystem<dim_ps> fe;  // TODO: const is probably wrong

  // NOTE: Quadratures are members of this class (and not e.g. part of the
  // assemble_system method), because I am using the mesh_loop function
  const QGauss<dim_ps> quadrature;
  const QGauss<dim_ps - 1> quadrature_face;

  // The constraints object is used for the distribute_local_to_global()
  // function, (maybe for the periodic boundary conditions) and, in the future,
  // for AMR
  const AffineConstraints<double> constraints;

  // PDE System
  PDESystem pde_system;
  UpwindFlux<dim_ps> upwind_flux;

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

  // Upwind flux matrices
  // Eigenvectors
  std::vector<FullMatrix<double>> eigenvectors_advection_matrices;
  std::vector<FullMatrix<double>> eigenvectors_adv_mat_prod_matrices;

  // Eigenvalues
  Vector<double> eigenvalues_adv;
  std::vector<Vector<double>> eigenvalues_adv_mat_prod;

  SparsityPattern sparsity_pattern;
  PETScWrappers::MPI::SparseMatrix mass_matrix;
  PETScWrappers::MPI::SparseMatrix dg_matrix;
  PETScWrappers::MPI::SparseMatrix system_matrix;

  PETScWrappers::MPI::Vector system_rhs;
  PETScWrappers::MPI::Vector locally_owned_previous_solution;
  PETScWrappers::MPI::Vector locally_relevant_current_solution;

  const int expansion_order = vfp_solver_control.expansion_order;
  const unsigned int num_exp_coefficients =
      static_cast<unsigned int>((expansion_order + 1) * (expansion_order + 1));

  // particle
  ParticleProperties particle_properties;

  TimerOutput timer;
};

VFPEquationSolver::VFPEquationSolver(const VFPSolverControl &control)
    : vfp_solver_control(control),
      mpi_communicator(MPI_COMM_WORLD),
      n_mpi_procs(Utilities::MPI::n_mpi_processes(mpi_communicator)),
      rank(Utilities::MPI::this_mpi_process(mpi_communicator)),
      pcout(std::cout, (rank == 0)),
      triangulation(mpi_communicator),
      dof_handler(triangulation),
      mapping(),
      fe(FE_DGQ<dim_ps>(vfp_solver_control.polynomial_degree),
         (vfp_solver_control.expansion_order + 1) *
             (vfp_solver_control.expansion_order + 1)),
      quadrature(fe.tensor_degree() + 1),
      quadrature_face(fe.tensor_degree() + 1),
      pde_system(vfp_solver_control.expansion_order),
      upwind_flux(pde_system, VFPSolverControl::momentum),
      timer(mpi_communicator, pcout, TimerOutput::never,
            TimerOutput::wall_times) {}

void VFPEquationSolver::run() {
  make_grid();
  setup_system();
  assemble_mass_matrix();

  // Project the initial values
  InitialValueFunction<dim_cs, VFPSolverControl::momentum> iv(expansion_order);
  project(iv, locally_relevant_current_solution);

  // parameters of the time stepping method
  double time_step = vfp_solver_control.time_step;
  double time = 0.;
  double final_time = vfp_solver_control.final_time;
  unsigned int time_step_number = 0;

  // Output time step zero
  output_results(time_step_number);

  // if the fields are time independent the dg matrix is not assembled inside
  // the time stepping methods. But it needs to be assembled once.
  if constexpr (!time_dependent_fields) assemble_dg_matrix(0.);

  time += time_step;
  ++time_step_number;

  pcout << "The time stepping loop is entered: \n";
  for (; time <= final_time; time += time_step, ++time_step_number) {
    pcout << "Time step " << time_step_number << " at t = " << time << "\n";
    // Time stepping method
    // theta_method(0.5);
    // explicit_runge_kutta();
    low_storage_explicit_runge_kutta(time, time_step);
    // NOTE: I cannot create TimerOutput::Scope inside output_results(),
    // because it is declared const.
    {
      TimerOutput::Scope timer_section(timer, "Output");
      output_results(time_step_number);
    }
  }
  pcout << "The simulation ended. \n";
}

void VFPEquationSolver::make_grid() {
  TimerOutput::Scope timer_section(timer, "Grid setup");
  // Number of refinements
  unsigned int num_refinements = vfp_solver_control.num_refinements;
  if constexpr ((flags & TermFlags::momentum) != TermFlags::none) {
    if constexpr (dim_cs == 1) {
      unsigned int n_cells = 1 << num_refinements;
      std::vector<unsigned int> repititions{n_cells, n_cells};
      Point<dim_ps> p1{-5., 1.};
      Point<dim_ps> p2{5., 10.};
      // // Colorize = true means to set boundary ids (default for 1D)
      GridGenerator::subdivided_hyper_rectangle(triangulation, repititions, p1,
                                                p2, true);
      // Periodic boundary conditions with MeshWorker. Mailinglist
    // https://groups.google.com/g/dealii/c/WlOiww5UVxc/m/mtQJDUwiBQAJ
    //
    // "If you call add_periodicity() on a Triangulation object, the periodic
    // faces are treated as internal faces in MeshWorker. This means that you
    // will not access them in a "integrate_boundary_term" function but in a
    // "integrate_face_term" function. "
    std::vector<GridTools::PeriodicFacePair<
        typename Triangulation<dim_ps>::cell_iterator>>
        matched_pairs;
    GridTools::collect_periodic_faces(triangulation, 0, 1, 0, matched_pairs);
    // GridTools::collect_periodic_faces(triangulation, 2, 3, 1, matched_pairs);
    triangulation.add_periodicity(matched_pairs);

      // GridGenerator::hyper_cube(triangulation, 1., 6., true);
      // triangulation.refine_global(num_refinements);
    }
    if constexpr (dim_cs == 2) {
      unsigned int n_cells = 1 << num_refinements;
      std::vector<unsigned int> repititions{n_cells, n_cells, n_cells};
      Point<dim_ps> p1{-5., -5., 1.};
      Point<dim_ps> p2{5., 5., 10.};
      GridGenerator::subdivided_hyper_rectangle(triangulation, repititions, p1,
                                                p2, true);
    }
  } else {
    GridGenerator::hyper_cube(triangulation, -5., 5., true);

    triangulation.refine_global(num_refinements);
  }

  // std::ofstream out("grid.vtk");
  // GridOut grid_out;
  // grid_out.write_vtk(triangulation, out);
  // std::cout << "	Grid written to grid.vtk"
  //           << "\n";
  pcout << "The grid was created: \n"
        << "	Number of active cells: "
        << triangulation.n_global_active_cells() << "\n";
}

void VFPEquationSolver::setup_system() {
  TimerOutput::Scope timer_section(timer, "FE system");

  dof_handler.distribute_dofs(fe);

  const unsigned int n_dofs = dof_handler.n_dofs();
  pcout << "The degrees of freedom were distributed: \n"
        << "	Number of degrees of freedom: " << n_dofs << "\n";

  locally_owned_dofs = dof_handler.locally_owned_dofs();
  locally_relevant_dofs = DoFTools::extract_locally_relevant_dofs(dof_handler);

  // constraints.clear();
  // DoFTools::make_periodicity_constraints(matched_pairs, constraints);

  // Let me see if I have to initialise the constraint object
  // constraints.clear();
  // constraints.reinit(locally_relevant_dofs);
  // This is an rather obscure line. I do not know why I need it. (cf. example
  // 23)
  // constraints.close();

  // Vectors
  locally_owned_previous_solution.reinit(locally_owned_dofs, mpi_communicator);
  locally_relevant_current_solution.reinit(
      locally_owned_dofs, locally_relevant_dofs, mpi_communicator);
  system_rhs.reinit(locally_owned_dofs, mpi_communicator);

  DynamicSparsityPattern dsp(locally_relevant_dofs);
  // NON-PERIODIC
  DoFTools::make_flux_sparsity_pattern(dof_handler, dsp);
  SparsityTools::distribute_sparsity_pattern(
      dsp, locally_owned_dofs, mpi_communicator, locally_relevant_dofs);
  // PERIODIC
  // DoFTools::make_flux_sparsity_pattern(dof_handler, dsp, constraints, false);
  dg_matrix.reinit(locally_owned_dofs, locally_owned_dofs, dsp,
                   mpi_communicator);
  // NOTE: DealII does not allow to use different sparsity patterns for
  // matrices, which you would like to add. Even though the the mass matrix
  // differs from the dg matrix.
  mass_matrix.reinit(locally_owned_dofs, locally_owned_dofs, dsp,
                     mpi_communicator);
  system_matrix.reinit(locally_owned_dofs, locally_owned_dofs, dsp,
                       mpi_communicator);
}

void VFPEquationSolver::assemble_mass_matrix() {
  TimerOutput::Scope timer_section(timer, "Mass matrix");

  using Iterator = typename DoFHandler<dim_ps>::active_cell_iterator;

  const auto cell_worker = [](const Iterator &cell,
                              ScratchData<dim_ps> &scratch_data,
                              CopyData &copy_data) {
    FEValues<dim_ps> &fe_v = scratch_data.fe_values;
    fe_v.reinit(cell);
    const unsigned int n_dofs = fe_v.get_fe().n_dofs_per_cell();
    // reinit the cell_matrix
    copy_data.reinit(cell, n_dofs);
    const std::vector<double> &JxW = fe_v.get_JxW_values();

    for (const unsigned int q_index : fe_v.quadrature_point_indices()) {
      for (unsigned int i : fe_v.dof_indices()) {
        const unsigned int component_i =
            fe_v.get_fe().system_to_component_index(i).first;
        for (unsigned int j : fe_v.dof_indices()) {
          const unsigned int component_j =
              fe_v.get_fe().system_to_component_index(j).first;
          // mass matrix
          copy_data.cell_matrix(i, j) +=
              (component_i == component_j
                   ? fe_v.shape_value(i, q_index) * fe_v.shape_value(j, q_index)
                   : 0) *
              JxW[q_index];
        }
      }
    }
  };

  const auto copier = [&](const CopyData &c) {
    constraints.distribute_local_to_global(c.cell_matrix, c.local_dof_indices,
                                           mass_matrix);
  };
  ScratchData<dim_ps> scratch_data(mapping, fe, quadrature, quadrature_face);
  CopyData copy_data;
  // Perfom the integration loop only over the locally owned cells
  const auto filtered_iterator_range =
      dof_handler.active_cell_iterators() | IteratorFilters::LocallyOwnedCell();

  pcout << "Begin the assembly of the mass matrix. \n";
  MeshWorker::mesh_loop(filtered_iterator_range, cell_worker, copier,
                        scratch_data, copy_data,
                        MeshWorker::assemble_own_cells);
  mass_matrix.compress(VectorOperation::add);
  pcout << "The mass matrix was assembled. \n";
}

void VFPEquationSolver::assemble_dg_matrix(const double time) {
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

  BackgroundVelocityField<dim_ps> background_velocity_field;
  background_velocity_field.set_time(time);
  MagneticField<dim_ps> magnetic_field;
  magnetic_field.set_time(time);

  ParticleVelocity<dim_ps> particle_velocity;
  ParticleGamma<dim_ps> particle_gamma;

  // I do not no the meaning of the following "const" specifier
  const auto cell_worker = [&](const Iterator &cell,
                               ScratchData<dim_ps> &scratch_data,
                               CopyData &copy_data) {
    FEValues<dim_ps> &fe_v = scratch_data.fe_values;
    // reinit cell
    fe_v.reinit(cell);

    const unsigned int n_dofs = fe_v.get_fe().n_dofs_per_cell();
    // reinit the cell_matrix
    copy_data.reinit(cell, n_dofs);

    const std::vector<Point<dim_ps>> &q_points = fe_v.get_quadrature_points();
    const std::vector<double> &JxW = fe_v.get_JxW_values();
    // Background Plasma
    // NOTE: Second argument constructs empty vectors with 2 values. They are
    // copied q_points.size() times.
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

    // Particle
    std::vector<double> particle_velocities(q_points.size());
    particle_velocity.value_list(q_points, particle_velocities);

    std::vector<double> particle_gammas(q_points.size());
    particle_gamma.value_list(q_points, particle_gammas);

    for (unsigned int i : fe_v.dof_indices()) {
      const unsigned int component_i =
          fe_v.get_fe().system_to_component_index(i).first;
      for (unsigned int j : fe_v.dof_indices()) {
        const unsigned int component_j =
            fe_v.get_fe().system_to_component_index(j).first;
        for (const unsigned int q_index : fe_v.quadrature_point_indices()) {
          if constexpr ((flags & TermFlags::collision) != TermFlags::none) {
            if (component_i == component_j) {
              // scattering_frequency * l(l+1) * \phi_i * \phi_j
              copy_data.cell_matrix(i, j) +=
                  scattering_frequency * collision_matrix[component_i] *
                  fe_v.shape_value(i, q_index) * fe_v.shape_value(j, q_index) *
                  JxW[q_index];
            }
          }
          // NOTE: If spatial advection term is deactivated, then dim_cs must
          // euqal zero, i.e. the distributin function is assumend to be
          // homogeneous (it does not depend on x, y and z). And if dim_cs = 0,
          // there for loop is not entered.
          for (unsigned int coordinate = 0; coordinate < dim_cs; ++coordinate) {
            // -[partial_x \phi_i * (u_x \delta_ij + Ax_ij)
            //   + partial_y \phi_i * (u_y \delta_ij + Ay_ij)] * phi_j
            if (component_i == component_j) {
              copy_data.cell_matrix(i, j) -=
                  fe_v.shape_grad(i, q_index)[coordinate] *
                  velocities[q_index][coordinate] *
                  fe_v.shape_value(j, q_index) * JxW[q_index];

              // - [\partial_x(u_x\delta_ij + Ax_ij) + \partial_y(u_y\delta_ij +
              // - Ay_ij) ] \phi_i \phi_j where \partial_x/y Ax/y_ij = 0
              copy_data.cell_matrix(i, j) -=
                  div_velocities[q_index] * fe_v.shape_value(i, q_index) *
                  fe_v.shape_value(j, q_index) * JxW[q_index];

            } else {
              // NOTE: Many zerso are added here, because the matrices Ax, Ay
              // and A_z are sparse. TODO: Performance check. If too bad, return
              // to the strategy, which was used in v0.6.5
              if ((flags & TermFlags::momentum) != TermFlags::none) {
                copy_data.cell_matrix(i, j) -=
                    fe_v.shape_grad(i, q_index)[coordinate] *
                    particle_velocities[q_index] *
                    advection_matrices[coordinate](component_i, component_j) *
                    fe_v.shape_value(j, q_index) * JxW[q_index];

              } else {
                // fixed energy case (i.e. transport only)
                copy_data.cell_matrix(i, j) -=
                    fe_v.shape_grad(i, q_index)[coordinate] *
                    particle_properties.velocity *
                    advection_matrices[coordinate](component_i, component_j) *
                    fe_v.shape_value(j, q_index) * JxW[q_index];
              }
            }
          }
          if constexpr ((flags & TermFlags::magnetic) != TermFlags::none) {
            // NOTE: All three components of the B-Field are included no
            // matter, which dimension of the configuration space is considered
            for (unsigned int coordinate = 0; coordinate < 3; ++coordinate) {
              if ((flags & TermFlags::momentum) != TermFlags::none) {
                copy_data.cell_matrix(i, j) -=
                    fe_v.shape_value(i, q_index) * particle_properties.charge *
                    magnetic_field_values[q_index][coordinate] /
                    particle_gammas[q_index] *
                    generator_rotation_matrices[coordinate](component_i,
                                                            component_j) *
                    fe_v.shape_value(j, q_index) * JxW[q_index];

              } else {
                // fixed energy case (i.e. transport only)
                copy_data.cell_matrix(i, j) -=
                    fe_v.shape_value(i, q_index) * particle_properties.charge *
                    magnetic_field_values[q_index][coordinate] /
                    particle_properties.gamma *
                    generator_rotation_matrices[coordinate](component_i,
                                                            component_j) *
                    fe_v.shape_value(j, q_index) * JxW[q_index];
              }
            }
          }
          if constexpr ((flags & TermFlags::momentum) != TermFlags::none) {
            // Momentum part
            for (unsigned int coordinate = 0; coordinate < 3; ++coordinate) {
              // \grad_phi * \gamma du^k/ dt * A_k \phi
              copy_data.cell_matrix(i, j) +=
                  fe_v.shape_grad(i, q_index)[dim_ps - 1] *
                  particle_gammas[q_index] *
                  material_derivative_vel[q_index][coordinate] *
                  advection_matrices[coordinate](component_i, component_j) *
                  fe_v.shape_value(j, q_index) * JxW[q_index];
              // \phi v * du^k/dt * A_k * \phi
              copy_data.cell_matrix(i, j) +=
                  fe_v.shape_value(i, q_index) * particle_velocities[q_index] *
                  material_derivative_vel[q_index][coordinate] *
                  advection_matrices[coordinate](component_i, component_j) *
                  fe_v.shape_value(j, q_index) * JxW[q_index];
              // \phi 1/v * du^k \ dt (A x \Omega)_k \phi
              copy_data.cell_matrix(i, j) +=
                  fe_v.shape_value(i, q_index) / particle_velocities[q_index] *
                  material_derivative_vel[q_index][coordinate] *
                  adv_x_gen_matrices[coordinate](component_i, component_j) *
                  fe_v.shape_value(j, q_index) * JxW[q_index];
            }
            for (unsigned int coordinate_1 = 0; coordinate_1 < 3;
                 ++coordinate_1) {
              for (unsigned int coordinate_2 = coordinate_1; coordinate_2 < 3;
                   ++coordinate_2) {
                if (coordinate_1 == coordinate_2) {
                  // \grad_phi p \jacobian[coordinate_1][coordinate_2]
                  // Ap_coordinate_1,coordinate_2 * \phi
                  copy_data.cell_matrix(i, j) +=
                      fe_v.shape_grad(i, q_index)[dim_ps - 1] *
                      q_points[q_index][dim_ps - 1] *
                      jacobians_vel[q_index][coordinate_1][coordinate_2] *
                      adv_mat_products[3 * coordinate_1 -
                                       coordinate_1 * (coordinate_1 + 1) / 2 +
                                       coordinate_2](component_i, component_j) *
                      fe_v.shape_value(j, q_index) * JxW[q_index];
                  // \phi * jacobian[coordinate_1][coordinate_2] *
                  // Ap_coordinate_1, coordinate_2 * \phi
                  copy_data.cell_matrix(i, j) +=
                      fe_v.shape_value(i, q_index) *
                      jacobians_vel[q_index][coordinate_1][coordinate_2] *
                      adv_mat_products[3 * coordinate_1 -
                                       coordinate_1 * (coordinate_1 + 1) / 2 +
                                       coordinate_2](component_i, component_j) *
                      fe_v.shape_value(j, q_index) * JxW[q_index];
                } else {
                  // symmetry
                  // component_1, component_2
                  // \grad_phi p \jacobian[coordinate_1][coordinate_2]
                  // Ap_coordinate_1,coordinate_2 * \phi
                  copy_data.cell_matrix(i, j) +=
                      fe_v.shape_grad(i, q_index)[dim_ps - 1] *
                      q_points[q_index][dim_ps - 1] *
                      jacobians_vel[q_index][coordinate_1][coordinate_2] *
                      adv_mat_products[3 * coordinate_1 -
                                       coordinate_1 * (coordinate_1 + 1) / 2 +
                                       coordinate_2](component_i, component_j) *
                      fe_v.shape_value(j, q_index) * JxW[q_index];

                  // \phi * jacobian[coordinate_1][coordinate_2] *
                  // Ap_coordinate_1, coordinate_2 * \phi
                  copy_data.cell_matrix(i, j) +=
                      fe_v.shape_value(i, q_index) *
                      jacobians_vel[q_index][coordinate_1][coordinate_2] *
                      adv_mat_products[3 * coordinate_1 -
                                       coordinate_1 * (coordinate_1 + 1) / 2 +
                                       coordinate_2](component_i, component_j) *
                      fe_v.shape_value(j, q_index) * JxW[q_index];

                  // component_2, component_1
                  // \grad_phi p \jacobian[coordinate_1][coordinate_2]
                  // Ap_coordinate_1,coordinate_2 * \phi
                  copy_data.cell_matrix(i, j) +=
                      fe_v.shape_grad(i, q_index)[dim_ps - 1] *
                      q_points[q_index][dim_ps - 1] *
                      jacobians_vel[q_index][coordinate_2][coordinate_1] *
                      adv_mat_products[3 * coordinate_1 -
                                       coordinate_1 * (coordinate_1 + 1) / 2 +
                                       coordinate_2](component_i, component_j) *
                      fe_v.shape_value(j, q_index) * JxW[q_index];

                  // \phi * jacobian[coordinate_1][coordinate_2] *
                  // Ap_coordinate_1, coordinate_2 * \phi
                  copy_data.cell_matrix(i, j) +=
                      fe_v.shape_value(i, q_index) *
                      jacobians_vel[q_index][coordinate_2][coordinate_1] *
                      adv_mat_products[3 * coordinate_1 -
                                       coordinate_1 * (coordinate_1 + 1) / 2 +
                                       coordinate_2](component_i, component_j) *
                      fe_v.shape_value(j, q_index) * JxW[q_index];
                }
              }
            }
            for (unsigned int coordinate_1 = 0; coordinate_1 < 3;
                 ++coordinate_1) {
              for (unsigned int coordinate_2 = 0; coordinate_2 < 3;
                   ++coordinate_2) {
                // \phi * jacobian[coordinate_1][coordinate_2] *
                // T_coordinate_2,coordinate_1 * \phi
                copy_data.cell_matrix(i, j) +=
                    fe_v.shape_value(i, q_index) *
                    jacobians_vel[q_index][coordinate_1][coordinate_2] *
                    t_matrices[coordinate_2 * 3 + coordinate_1](component_i,
                                                                component_j) *
                    fe_v.shape_value(j, q_index) * JxW[q_index];
              }
            }
          }
        }
      }
    }
  };
  // assemble boundary face terms
  const auto boundary_worker =
      [&](const Iterator &cell, const unsigned int &face_no,
          ScratchData<dim_ps> &scratch_data, CopyData &copy_data) {
        scratch_data.fe_values_face.reinit(cell, face_no);
        const FEFaceValuesBase<dim_ps> &fe_face_v = scratch_data.fe_values_face;
        // Every shape function on the cell could contribute to the face
        // integral, hence n_facet_dofs = n_dofs_per_cell
        const unsigned int n_facet_dofs = fe_face_v.get_fe().n_dofs_per_cell();
        // NOTE: copy_data is not reinitialised, the cell_workers contribution
        // to the cell_dg_matrix should not be deleted

        const std::vector<Point<dim_ps>> &q_points =
            fe_face_v.get_quadrature_points();
        const std::vector<double> &JxW = fe_face_v.get_JxW_values();
        const std::vector<Tensor<1, dim_ps>> &normals =
            fe_face_v.get_normal_vectors();
        std::vector<FullMatrix<double>> positive_flux_matrices(
            q_points.size(), FullMatrix<double>(num_exp_coefficients));
        std::vector<FullMatrix<double>> negative_flux_matrices(
            q_points.size(), FullMatrix<double>(num_exp_coefficients));

        upwind_flux.compute_upwind_fluxes(
            q_points, normals, positive_flux_matrices, negative_flux_matrices);

        for (unsigned int i = 0; i < n_facet_dofs; ++i) {
          const unsigned int component_i =
              fe_face_v.get_fe().system_to_component_index(i).first;
          for (unsigned int j = 0; j < n_facet_dofs; ++j) {
            const unsigned int component_j =
                fe_face_v.get_fe().system_to_component_index(j).first;
            for (unsigned int q_index : fe_face_v.quadrature_point_indices()) {
              if constexpr ((flags & TermFlags::spatial_advection) !=
                            TermFlags::none) {
                // Outflow boundary: Everyhing with a positive flux along the
                // direction of the normal leaves the boundary
                copy_data.cell_matrix(i, j) +=
                    fe_face_v.shape_value(i, q_index) *
                    positive_flux_matrices[q_index](component_i, component_j) *
                    fe_face_v.shape_value(j, q_index) * JxW[q_index];
                copy_data.cell_matrix(i, j) +=
                    fe_face_v.shape_value(i, q_index) *
                    negative_flux_matrices[q_index](component_i, component_j) *
                    fe_face_v.shape_value(j, q_index) * JxW[q_index];
              }
            }
          }
        }
      };

  // assemble interior face terms
  // NOTE: The face worker assumes a grid consisting of rectangular cells with
  // the following pattern of interior face normals
  //    ^
  // *--|--*
  // |     |
  // -     -->   ^ y
  // |     |     |
  // *-----*     *--> x
  // Hence, n_x = 1. and n_y = 1.
  const auto face_worker = [&](const Iterator &cell,
                               const unsigned int &face_no,
                               const unsigned int &subface_no,
                               const Iterator &neighbor_cell,
                               const unsigned int &neighbor_face_no,
                               const unsigned int &neighbor_subface_no,
                               ScratchData<dim_ps> &scratch_data,
                               CopyData &copy_data) {
    // NOTE: The flag MeshWorker::assemble_own_interior_faces_both will not
    // work for this face_worker. It implicitly assumes that faces are only
    // assembled once. And hence subface_no, can be ignored
    (void)subface_no;  // suppress compiler warning

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
    CopyDataFace &copy_data_face = copy_data.face_data.back();
    const unsigned int n_dofs = fe_v_face.get_fe().n_dofs_per_cell();
    copy_data_face.reinit(cell, neighbor_cell, n_dofs);

    const std::vector<Point<dim_ps>> &q_points =
        fe_v_face.get_quadrature_points();
    const std::vector<double> &JxW = fe_v_face.get_JxW_values();
    const std::vector<Tensor<1, dim_ps>> &normals =
        fe_v_face.get_normal_vectors();
    // For every interior face there is an in- and outflow represented by the
    // corresponding flux matrices
    std::vector<FullMatrix<double>> positive_flux_matrices(
        q_points.size(), FullMatrix<double>(num_exp_coefficients));
    std::vector<FullMatrix<double>> negative_flux_matrices(
        q_points.size(), FullMatrix<double>(num_exp_coefficients));

    upwind_flux.compute_upwind_fluxes(q_points, normals, positive_flux_matrices,
                                      negative_flux_matrices);

    for (unsigned int q_index : fe_v_face.quadrature_point_indices()) {
      // cell_dg_matrix_11
      for (unsigned int i : fe_v_face.dof_indices()) {
        const unsigned int component_i =
            fe_v_face.get_fe().system_to_component_index(i).first;
        for (unsigned int j : fe_v_face.dof_indices()) {
          unsigned int component_j =
              fe_v_face.get_fe().system_to_component_index(j).first;
          copy_data_face.cell_dg_matrix_11(i, j) +=
              fe_v_face.shape_value(i, q_index) *
              positive_flux_matrices[q_index](component_i, component_j) *
              fe_v_face.shape_value(j, q_index) * JxW[q_index];
        }
      }
      // cell_dg_matrix_12
      for (unsigned int i : fe_v_face.dof_indices()) {
        const unsigned int component_i =
            fe_v_face.get_fe().system_to_component_index(i).first;
        for (unsigned int j : fe_v_face_neighbor.dof_indices()) {
          unsigned int component_j =
              fe_v_face_neighbor.get_fe().system_to_component_index(j).first;
          copy_data_face.cell_dg_matrix_12(i, j) -=
              fe_v_face_neighbor.shape_value(i, q_index) *
              positive_flux_matrices[q_index](component_i, component_j) *
              fe_v_face.shape_value(j, q_index) * JxW[q_index];
        }
      }
      // cell_dg_matrix_21
      for (unsigned int i : fe_v_face_neighbor.dof_indices()) {
        const unsigned int component_i =
            fe_v_face_neighbor.get_fe().system_to_component_index(i).first;
        for (unsigned int j : fe_v_face.dof_indices()) {
          unsigned int component_j =
              fe_v_face.get_fe().system_to_component_index(j).first;
          copy_data_face.cell_dg_matrix_21(i, j) +=
              fe_v_face.shape_value(i, q_index) *
              negative_flux_matrices[q_index](component_i, component_j) *
              fe_v_face_neighbor.shape_value(j, q_index) * JxW[q_index];
        }
      }
      // cell_dg_matrix_22
      for (unsigned int i : fe_v_face_neighbor.dof_indices()) {
        const unsigned int component_i =
            fe_v_face_neighbor.get_fe().system_to_component_index(i).first;
        for (unsigned int j : fe_v_face_neighbor.dof_indices()) {
          unsigned int component_j =
              fe_v_face_neighbor.get_fe().system_to_component_index(j).first;
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
    constraints.distribute_local_to_global(c.cell_matrix, c.local_dof_indices,
                                           dg_matrix);
    for (auto &cdf : c.face_data) {
      for (unsigned int i = 0; i < cdf.local_dof_indices.size(); ++i)
        for (unsigned int j = 0; j < cdf.local_dof_indices.size(); ++j) {
          dg_matrix.add(cdf.local_dof_indices[i], cdf.local_dof_indices[j],
                        cdf.cell_dg_matrix_11(i, j));
          dg_matrix.add(cdf.local_dof_indices_neighbor[i],
                        cdf.local_dof_indices[j], cdf.cell_dg_matrix_12(i, j));
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
  CopyData copy_data;
  pcout << "	Begin the assembly of the DG matrix. \n";
  const auto filtered_iterator_range =
      dof_handler.active_cell_iterators() | IteratorFilters::LocallyOwnedCell();
  MeshWorker::mesh_loop(
      filtered_iterator_range, cell_worker, copier, scratch_data, copy_data,
      MeshWorker::assemble_own_cells | MeshWorker::assemble_boundary_faces |
          MeshWorker::assemble_ghost_faces_once |
          MeshWorker::assemble_own_interior_faces_once,
      boundary_worker, face_worker);
  dg_matrix.compress(VectorOperation::add);
  pcout << "	The DG matrix was assembled. \n";
}

template <int dim>
void VFPEquationSolver::project(
    const Function<dim> &f, PETScWrappers::MPI::Vector &projected_function) {
  TimerOutput::Scope timer_section(timer, "Project f onto the FEM space");
  pcout << "Project a function onto the finite element space \n";
  // Create right hand side
  FEValues<dim_ps> fe_v(
      mapping, fe, quadrature,
      update_values | update_quadrature_points | update_JxW_values);

  const unsigned int n_dofs = fe.n_dofs_per_cell();
  Vector<double> cell_rhs(n_dofs);

  std::vector<types::global_dof_index> local_dof_indices(n_dofs);

  for (const auto &cell : dof_handler.active_cell_iterators()) {
    if (cell->is_locally_owned()) {
      cell_rhs = 0;
      fe_v.reinit(cell);

      const std::vector<Point<dim_ps>> &q_points = fe_v.get_quadrature_points();
      const std::vector<double> &JxW = fe_v.get_JxW_values();

      // Initial values
      std::vector<Vector<double>> function_values(
          q_points.size(), Vector<double>(num_exp_coefficients));
      f.vector_value_list(q_points, function_values);

      for (const unsigned int q_index : fe_v.quadrature_point_indices()) {
        for (unsigned int i : fe_v.dof_indices()) {
          const unsigned int component_i =
              fe.system_to_component_index(i).first;
          cell_rhs(i) += fe_v.shape_value(i, q_index) *
                         function_values[q_index][component_i] * JxW[q_index];
        }
      }
      cell->get_dof_indices(local_dof_indices);

      constraints.distribute_local_to_global(cell_rhs, local_dof_indices,
                                             system_rhs);
    }
  }
  system_rhs.compress(VectorOperation::add);

  // Solve the system
  PETScWrappers::PreconditionNone preconditioner;
  preconditioner.initialize(mass_matrix);
  // NOTE: It is not possible to directly write into a ghosted vector. And the
  // the second argument of project is meant to be a ghosted vector. Hence, we
  // create an unghosted vector and copy it later on.
  PETScWrappers::MPI::Vector temporary_non_ghosted_vector(locally_owned_dofs,
                                                          mpi_communicator);
  SolverControl solver_control(1000, 1e-12);
  PETScWrappers::SolverCG cg(solver_control, mpi_communicator);
  cg.solve(mass_matrix, temporary_non_ghosted_vector, system_rhs,
           preconditioner);
  pcout << "	Solved in " << solver_control.last_step() << " iterations."
        << std::endl;
  // At the moment I am assuming, that I do not have constraints. Hence, I do
  // not need the following line.
  // constraints.distribute(completely_distributed_solution);
  projected_function = temporary_non_ghosted_vector;

  // Reset system RHS
  system_rhs = 0;
}

// void VFPEquationSolver::theta_method_solve_system() {
//   SolverControl solver_control(1000, 1e-12);
//   SolverRichardson<Vector<double>> solver(solver_control);

//   PreconditionBlockSSOR<SparseMatrix<double>> preconditioner;
//   preconditioner.initialize(system_matrix, fe.n_dofs_per_cell());
//   solver.solve(system_matrix, current_solution, system_rhs, preconditioner);

//   std::cout << "	Solver converged in " << solver_control.last_step()
//             << " iterations."
//             << "\n";
// }

// template <TermFlags flags, int dim_cs>
// void VFPEquationSolver<flags, dim_cs>::theta_method(double theta) {
//   TimerOutput::Scope timer_section(timer, "Theta method");

//   Vector<double> tmp(locally_relevant_current_solution.size());

//   mass_matrix.vmult(system_rhs, locally_relevant_previous_solution);
//   dg_matrix.vmult(tmp, locally_relevant_previous_solution);
//   system_rhs.add(-time_step * (1 - theta), tmp);

//   // Since the the dg_matrix depends on the velocity field and the the
//   // velocity field may depend on time, it needs to reassembled every time
//   // step. This is not true for the mass matrix ( but it may if the grid
//   // adapts after a specified amount of time steps)

//   // mass_matrix = 0; dg_matrix = 0;
//   // assemble_system(time + time_step);

//   // NOTE: Currently that is not necessary, because the velocity field is
//   // constant. And I should check if I cannot update the system matrix
//   without
//   // reassembling in each time step.

//   system_matrix.copy_from(mass_matrix);
//   system_matrix.add(time_step * theta, dg_matrix);
//   theta_method_solve_system();
// }

// template <TermFlags flags, int dim_cs>
// void VFPEquationSolver<flags, dim_cs>::explicit_runge_kutta() {
//   TimerOutput::Scope timer_section(timer, "ERK4");

//   // ERK 4
//   // Butcher's array
//   Vector<double> a({0.5, 0.5, 1.});
//   Vector<double> b({1. / 6, 1. / 3, 1. / 3, 1. / 6});
//   // NOTE: c is only necessary if the velocity field and the magnetic field
//   // are time dependent.
//   // Vector<double> c({0., 0.5, 0.5, 1.});

//   // Allocate storage for the four stages
//   std::vector<Vector<double>> k(4,
//   Vector<double>(locally_relevant_current_solution.size()));
//   // The mass matrix needs to be "inverted" in every stage
//   SolverControl solver_control(1000, 1e-12);
//   SolverCG<Vector<double>> cg(solver_control);
//   // k_0
//   dg_matrix.vmult(system_rhs, locally_relevant_previous_solution);
//   cg.solve(mass_matrix, k[0], system_rhs, PreconditionIdentity());
//   std::cout << "	Stage s: " << 0 << "	Solver converged in "
//             << solver_control.last_step() << " iterations."
//             << "\n";
//   locally_relevant_current_solution.add(b[0] * time_step, k[0]);

//   Vector<double> temp(locally_relevant_current_solution.size());
//   for (unsigned int s = 1; s < 4; ++s) {
//     // NOTE: It would be nessary to reassemble the dg matrix twice for
//     // (for half a time step and a whole time step) if the velocity field and
//     // the magnetic field were time dependent.
//     temp.add(-1., locally_relevant_previous_solution, -a[s - 1] * time_step,
//     k[s - 1]);
//     // assemble_system(time + c[s]*time_step);
//     dg_matrix.vmult(system_rhs, temp);
//     cg.solve(mass_matrix, k[s], system_rhs, PreconditionIdentity());
//     std::cout << "	Stage s: " << s << "	Solver converged in "
//               << solver_control.last_step() << " iterations."
//               << "\n";

//     locally_relevant_current_solution.add(b[s] * time_step, k[s]);

//     // empty temp vector
//     temp = 0;
//   }
//   // NOTE: It is not necessary to assemble the matrix again, because the last
//   // element of the c vector is 1. (see low_storage_erk())
// }

void VFPEquationSolver::low_storage_explicit_runge_kutta(
    const double time, const double time_step) {
  TimerOutput::Scope timer_section(timer, "LSERK");

  // see Hesthaven p.64
  Vector<double> a(
      {0., -567301805773. / 1357537059087, -2404267990393. / 2016746695238,
       -3550918686646. / 2091501179385, -1275806237668. / 842570457699});
  Vector<double> b(
      {1432997174477. / 9575080441755, 5161836677717. / 13612068292357,
       1720146321549. / 2090206949498, 3134564353537. / 4481467310338,
       2277821191437. / 14882151754819});
  // NOTE: I only need c if the velocity field and the magnetic field are time
  // dependent
  Vector<double> c(
      {0., 1432997174477. / 9575080441755, 2526269341429. / 6820363962896,
       2006345519317. / 3224310063776, 2802321613138. / 2924317926251});

  // The mass matrix needs to be "inverted" in every stage

  // PETScWrappers::PreconditionBoomerAMG preconditioner;
  // PETScWrappers::PreconditionBoomerAMG::AdditionalData data;
  // data.symmetric_operator = true;
  PETScWrappers::PreconditionNone preconditioner;
  preconditioner.initialize(mass_matrix);

  SolverControl solver_control(1000, 1e-12);
  PETScWrappers::SolverCG cg(solver_control, mpi_communicator);
  // NOTE: The locally_relevant_current_solution is a "ghosted" vector and it
  // cannot be written to. It is necessary to use a vector does not contain
  // ghost cells. We extract the locally owned part with the equal sign
  // operator.
  locally_owned_previous_solution = locally_relevant_current_solution;

  PETScWrappers::MPI::Vector k(locally_owned_dofs, mpi_communicator);
  PETScWrappers::MPI::Vector temp(locally_owned_dofs, mpi_communicator);
  for (unsigned int s = 0; s < 5; ++s) {
    // only assemble the dg_matrix in every stage if the fields are time
    // dependent
    if constexpr (time_dependent_fields)
      assemble_dg_matrix(time + c[s] * time_step);

    dg_matrix.vmult(system_rhs, locally_owned_previous_solution);

    if constexpr (time_dependent_fields) dg_matrix = 0;

    cg.solve(mass_matrix, temp, system_rhs, preconditioner);
    pcout << "	Stage s: " << s << "	Solver converged in "
          << solver_control.last_step() << " iterations."
          << "\n";

    k.sadd(a[s], -time_step, temp);
    locally_owned_previous_solution.add(b[s], k);
  }
  // Currently I assume that there are no constraints
  // constraints.distribute(locally_relevant_current_solution);
  // std::cout << "Rank: " << rank << "\n";
  // std::cout << "Locally owned current solution: \n";
  // locally_owned_current_solution.print(std::cout);

  locally_relevant_current_solution = locally_owned_previous_solution;
  // std::cout << "Locally relevant current solution: \n";
  // locally_relevant_current_solution.print(std::cout);
}

void VFPEquationSolver::output_results(
    const unsigned int time_step_number) const {
  DataOut<dim_ps> data_out;
  data_out.attach_dof_handler(dof_handler);
  // Create a vector of strings with names for the components of the solution
  std::vector<std::string> component_names(num_exp_coefficients);
  const std::vector<std::array<unsigned int, 3>> &lms_indices =
      pde_system.get_lms_indices();

  for (unsigned int i = 0; i < num_exp_coefficients; ++i) {
    const std::array<unsigned int, 3> &lms = lms_indices[i];
    component_names[i] = "f_" + std::to_string(lms[0]) +
                         std::to_string(lms[1]) + std::to_string(lms[2]);
  }

  data_out.add_data_vector(locally_relevant_current_solution, component_names);

  // Output the partition of the mesh
  Vector<float> subdomain(triangulation.n_active_cells());
  for (unsigned int i = 0; i < subdomain.size(); ++i)
    subdomain(i) = triangulation.locally_owned_subdomain();
  data_out.add_data_vector(subdomain, "subdomain");

  data_out.build_patches();
  // const std::string file_path = "results/solution" +
  //                               Utilities::int_to_string(time_step_number, 3)
  //                               +
  //                               ".vtu";

  DataOutBase::VtkFlags vtk_flags;
  vtk_flags.compression_level = DataOutBase::VtkFlags::best_speed;
  data_out.set_flags(vtk_flags);

  // std::ofstream output(file_path);
  // data_out.write_vtu(output);
  data_out.write_vtu_with_pvtu_record("./results/", "solution",
                                      time_step_number, mpi_communicator, 3, 8);
}

}  // namespace VFPEquation

int main(int argc, char *argv[]) {
  try {
    using namespace VFPEquation;
    Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

    VFPSolverControl vfp_solver_control("vfp-equation.prm");
    if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
      vfp_solver_control.print_settings(std::cout);

    VFPEquationSolver vfp_equation_solver(vfp_solver_control);
    vfp_equation_solver.run();

  } catch (std::exception &exc) {
    std::cerr << std::endl
              << std::endl
              << "----------------------------------------------------"
              << std::endl;
    std::cerr << "Exception on processing: " << std::endl
              << exc.what() << std::endl
              << "Aborting!" << std::endl
              << "----------------------------------------------------"
              << std::endl;
    return 1;
  } catch (...) {
    std::cerr << std::endl
              << std::endl
              << "----------------------------------------------------"
              << std::endl;
    std::cerr << "Unknown exception!" << std::endl
              << "Aborting!" << std::endl
              << "----------------------------------------------------"
              << std::endl;
    return 1;
  }
  return 0;
}
