#include <deal.II/base/exceptions.h>
#include <deal.II/base/function.h>
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
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/tria.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/lapack_full_matrix.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/precondition_block.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/solver_richardson.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/vector.h>
// #include <deal.II/grid/grid_refinement.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_interface_values.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_q1.h>
#include <deal.II/meshworker/assemble_flags.h>
#include <deal.II/meshworker/mesh_loop.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/vector_tools_project.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iostream>
#include <numeric>
#include <string>
#include <vector>

namespace vfp_equation_solver {
using namespace dealii;
// Input data
class ParameterReader {
 public:
  ParameterReader(ParameterHandler &prm);
  void read_parameters(const std::string &input_file);

 private:
  void declare_parameters();
  ParameterHandler &parameter_handler;
};

ParameterReader::ParameterReader(ParameterHandler &prm)
    : parameter_handler(prm) {}

void ParameterReader::declare_parameters() {
  parameter_handler.enter_subsection("Mesh");
  {  // NOTE: This is a very strange syntax
    parameter_handler.declare_entry("Number of refinements", "6",
                                    Patterns::Integer(0),
                                    "Number of global mesh refinement steps");
  }
  parameter_handler.leave_subsection();

  parameter_handler.enter_subsection("Time stepping");
  {
    parameter_handler.declare_entry("Time step size", "7.8125e-3",
                                    Patterns::Double(),
                                    "Duration of the simulation.");
    parameter_handler.declare_entry("Final time", "0.4", Patterns::Double(),
                                    "Duration of the simulation.");
  }
  parameter_handler.leave_subsection();

  parameter_handler.enter_subsection("Expansion");
  {
    parameter_handler.declare_entry(
        "Expansion order", "0", Patterns::Integer(0),
        "The order of the expansion of the particle distribution function.");
  }
  parameter_handler.leave_subsection();

  parameter_handler.enter_subsection("Finite element");
  {
    parameter_handler.declare_entry("Polynomial degree", "1",
                                    Patterns::Integer(0),
                                    "The degree of the shape functions (i.e. "
                                    "the polynomials) of the finite element.");
  }
  parameter_handler.leave_subsection();

  parameter_handler.enter_subsection("Physical parameters");
  {
    parameter_handler.declare_entry(
        "Scattering frequency", "1.", Patterns::Double(),
        "Frequency at which energetic particles are scattered.");
    parameter_handler.declare_entry(
        "Particle energy", "1.", Patterns::Double(),
        "The energy of the particle in units of 10 Gev");
    // parameter_handler.declare_entry(
    //     "Plasma velocity x-component", ".1", Patterns::Double(),
    //     "The x-component of the background plasma's velocity");
    // parameter_handler.declare_entry(
    //     "Plasma velocity y-component", ".1", Patterns::Double(),
    //     "The y-component of the background plasma's velocity");
    parameter_handler.declare_entry("Magnetic field x-component", "0.",
                                    Patterns::Double(),
                                    "The x-component of the magnetic field");
    parameter_handler.declare_entry("Magnetic field y-component", "0.",
                                    Patterns::Double(),
                                    "The y-component of the magnetic field");
    parameter_handler.declare_entry("Magnetic field z-component", "1.",
                                    Patterns::Double(),
                                    "The z-component of the magnetic field");
  }
  parameter_handler.leave_subsection();
}

void ParameterReader::read_parameters(const std::string &input_file) {
  declare_parameters();
  parameter_handler.parse_input(input_file);
}

struct ParticleProperties {
  ParticleProperties(ParameterHandler &prm) : parameter_handler(prm) {
    parameter_handler.enter_subsection("Physical parameters");
    {
      particle_energy = parameter_handler.get_double("Particle energy");
      particle_velocity = 1. / (beta_0 * gamma_0) *
                          std::sqrt(gamma_0 * gamma_0 -
                                    1. / (particle_energy * particle_energy));
    }
    parameter_handler.leave_subsection();
  }
  // dimensionless values
  double particle_energy;
  double particle_velocity;
  // reference values
  const double energy_scale = 10;        // GeV
  const double proton_mass = 0.9382721;  // GeV/c^2
  const double gamma_0 = energy_scale / proton_mass;
  const double beta_0 = std::sqrt(1 - 1. / (gamma_0 * gamma_0));

 private:
  ParameterHandler &parameter_handler;
};

std::ostream &operator<<(std::ostream &os,
                         const ParticleProperties &particle_properties) {
  os << "Particle properties: \n"
     << " Dimensionless units: \n"
     << "	Particle energy: " << particle_properties.particle_energy
     << "\n"
     << "	Particle velocity: " << particle_properties.particle_velocity
     << "\n"
     << " Reference values: \n"
     << "	Energy scale (E_0): " << particle_properties.energy_scale
     << " GeV \n"
     << "	Proton_mass (m_p): " << particle_properties.proton_mass
     << " GeV/c^2 \n"
     << "	beta_0: " << particle_properties.beta_0 << "\n"
     << "	gamma_0: " << particle_properties.gamma_0 << "\n\n";
  return os;
}

// the background velocity field
template <int dim>
class BackgroundVelocityField : public Function<dim> {
 public:
  virtual void vector_value(const Point<dim> &point,
                            Vector<double> &value) const override {
    Assert(dim <= 2, ExcNotImplemented());
    // constant velocity field
    (void)point;
    if constexpr (dim == 1) value[0] = 0.2;
    if constexpr (dim == 2) {
      // constant velocity
      value[0] = 0.1;
      value[1] = 0.1;
      // rigid rotator
      // value[0] = -point[1];
      // value[1] = point[0];
    }
    if constexpr (dim == 3) {
      // constant velocity
      value[0] = 0.;
      value[1] = 0.;
      value[2] = 0.;
      // rigid rotator
      // value[0] = -point[1];
      // value[1] = point[0];
    }

  }

  void divergence_list(const std::vector<Point<dim>> &points,
                       std::vector<double> &values) {
    Assert(dim <= 2, ExcNotImplemented());
    Assert(values.size() == points.size(),
           ExcDimensionMismatch(values.size(), points.size()));
    std::fill(values.begin(), values.end(), 0.);
    // for (unsigned int q_index = 0; q_index < points.size(); ++q_index) {
    //   values[q_index] = points[q_index][0]*0.5;
    // }
  }
  // TODO: Time derivative and total derivative
};

// the magnetic field
template <int dim>
class MagneticField : public Function<dim> {
 public:
  MagneticField(ParameterHandler &prm) : parameter_handler(prm) {
    parameter_handler.enter_subsection("Physical parameters");
    {
      B_x = parameter_handler.get_double("Magnetic field x-component");
      B_y = parameter_handler.get_double("Magnetic field y-component");
      B_z = parameter_handler.get_double("Magnetic field z-component");
    }
    parameter_handler.leave_subsection();
  }
  virtual double value(const Point<dim> &point,
                       const unsigned int component) const override {
    // Check if a component was picked
    Assert(component == 0 || component == 1 || component == 2,
           ExcNotImplemented());
    // constant magnetic field
    (void)point;  // (suppresses the compiler warning unused variable)
    // constant magnetic field
    switch (component) {
      case 0:
        return B_x;
        break;
      case 1:
        return B_y;
        break;
      case 2:
        return B_z;
        break;
        // NOTE: This is very bad style, because it does not give an error, if
        // no component was chosen
      default:
        return B_x;
    }
  }

  virtual void vector_value(const Point<dim> &point,
                            Vector<double> &value) const override {
    (void)point;
    value[0] = B_x;
    value[1] = B_y;
    value[2] = B_z;
  }

 private:
  ParameterHandler &parameter_handler;
  double B_x = 0.;
  double B_y = 0.;
  double B_z = 1.;
};

enum class TermFlags {
  none = 0,
  advection = 1 << 0,
  reaction = 1 << 1,
  magnetic = 1 << 2
};

constexpr TermFlags operator|(TermFlags f1, TermFlags f2) {
  return static_cast<TermFlags>(static_cast<int>(f1) | static_cast<int>(f2));
}

constexpr TermFlags operator&(TermFlags f1, TermFlags f2) {
  return static_cast<TermFlags>(static_cast<int>(f1) & static_cast<int>(f2));
}

template <typename StreamType>
inline StreamType &operator<<(StreamType &s, TermFlags f) {
  s << "Term flags: \n";
  if ((f & TermFlags::advection) != TermFlags::none) s << "	 - Advection\n";
  if ((f & TermFlags::reaction) != TermFlags::none) s << "	 - Reaction\n";
  if ((f & TermFlags::magnetic) != TermFlags::none) s << "	 - Magnetic\n";
  return s;
}

// Initial values
template <int dim>
class InitialValueFunction : public Function<dim> {
 public:
  InitialValueFunction(unsigned int exp_order) : expansion_order{exp_order} {}
  virtual void vector_value(const Point<dim> &p,
                            Vector<double> &values) const override {
    Assert(dim <= 2, ExcNotImplemented());
    Assert(values.size() == (expansion_order + 1) * (expansion_order + 1),
           ExcDimensionMismatch(values.size(),
                                (expansion_order + 1) * (expansion_order + 1)));
    // The zeroth component of values corresponds to f_000, the first component
    // to f_110 etc.
    if constexpr (dim == 1) values[0] = 1. * std::exp(-(std::pow(p[0], 2)));
    // values[0] = std::sin((1. * 3.14159265359) / 2 * p[0]) + 1.;

    if constexpr (dim == 2) {
      // Gaussian
      values[0] = 1. * std::exp(-((std::pow(p[0], 2) + std::pow(p[1], 2))));

      // constant disc
      // if (p.norm() <= 1.) values[0] = 1.;

      // Rigid rotator
      // if (std::abs(p[0]) <= 3 && std::abs(p[1]) <= 0.5) values[0] = 1.;
      // if (std::abs(p[0]) <= 0.5 && std::abs(p[1]) <= 3) values[0] = 1.;
    }
    if constexpr (dim == 3) 
      values[0] = 1. * std::exp(-((std::pow(p[0], 2) + std::pow(p[1], 2) + std::pow(p[0], 2))));

    // Fill all components with the same values
    // std::fill(
    //     values.begin(), values.end(),
    //     1. * std::exp(-((std::pow(p[0] - 0.5, 2) + std::pow(p[1] - 0.5, 2)) /
    //                     0.01)));
  }

 private:
  const unsigned int expansion_order = 1;
};

// The mesh_loop function requires helper data types
template <int dim>
class ScratchData {
 public:
  // Constructor
  ScratchData(const Mapping<dim> &mapping, const FiniteElement<dim> &fe,
              const Quadrature<dim> &quadrature,
              const Quadrature<dim - 1> &quadrature_face,
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
  ScratchData(const ScratchData<dim> &scratch_data)
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

  FEValues<dim> fe_values;
  FEFaceValues<dim> fe_values_face;
  FEFaceValues<dim> fe_values_face_neighbor;
  FESubfaceValues<dim> fe_values_subface_neighbor;
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
  FullMatrix<double> cell_dg_matrix;
  FullMatrix<double> cell_mass_matrix;
  Vector<double> cell_rhs;
  std::vector<types::global_dof_index> local_dof_indices;
  std::vector<types::global_dof_index> local_dof_indices_neighbor;
  std::vector<CopyDataFace> face_data;

  template <typename Iterator>
  void reinit(const Iterator &cell, unsigned int dofs_per_cell) {
    cell_dg_matrix.reinit(dofs_per_cell, dofs_per_cell);
    cell_mass_matrix.reinit(dofs_per_cell, dofs_per_cell);
    cell_rhs.reinit(dofs_per_cell);

    local_dof_indices.resize(dofs_per_cell);
    cell->get_dof_indices(local_dof_indices);
  }
};

template <TermFlags flags, int dim>
class VFPEquationSolver {
 public:
  VFPEquationSolver(ParameterHandler &prm, unsigned int polynomial_degree,
                    int order_expansion);
  void run();

 private:
  void make_grid();
  void setup_pde_system();
  void prepare_upwind_fluxes();
  // NOTE: The LAPACKFullMatrix class does not implement a move constructor,
  // this is why the output argument is an input argument
  void compute_upwind_fluxes(
      const std::vector<Point<dim>> &q_points,
      const std::vector<Tensor<1, dim>> &normals,
      std::vector<FullMatrix<double>> &positive_flux_matrices,
      std::vector<FullMatrix<double>> &negative_flux_matrices);
  void setup_system();
  void assemble_mass_matrix();
  void assemble_dg_matrix(unsigned int evaluation_time);
  void project_initial_condition();
  void theta_method_solve_system();
  void theta_method(double theta);
  void explicit_runge_kutta();
  void low_storage_explicit_runge_kutta();
  void output_results() const;
  void output_compile_time_parameters() const;
  void output_index_order() const;

  ParameterHandler &parameter_handler;
  Triangulation<dim> triangulation;
  DoFHandler<dim> dof_handler;

  // NOTE: The explicit use of a mapping is most likely related to the usage
  // of mesh_loop as well
  const MappingQ1<dim> mapping;
  const FESystem<dim> fe;  // TODO: const is probably wrong

  // NOTE: Quadratures are members of this class (and not e.g. part of the
  // assemble_system method), because I am using the mesh_loop function
  const QGauss<dim> quadrature;
  const QGauss<dim - 1> quadrature_face;

  AffineConstraints<double> constraints;  // used in the L_2 projection of the
  // intial condition onto the finite
  // element space

  // PDE System data
  // Spatial advection/(magnitude) p advection
  std::vector<LAPACKFullMatrix<double>> advection_matrices;

  // Rotataion matrices (due to the magnetic field)
  std::vector<LAPACKFullMatrix<double>> generator_rotation_matrices;
  // (magnitude) p advection
  LAPACKFullMatrix<double> Ap_xx;
  LAPACKFullMatrix<double> Ap_xy;
  LAPACKFullMatrix<double> Ap_xz;
  LAPACKFullMatrix<double> Ap_yy;
  LAPACKFullMatrix<double> Ap_yz;
  LAPACKFullMatrix<double> Ap_zz;

  // Collision term (essentially a reaction term)
  Vector<double> collision_matrix;

  // Upwind flux matrices
  // Eigenvectors
  std::vector<FullMatrix<double>> eigenvectors_advection_matrices;
  // Eigenvalues
  Vector<double> eigenvalues;

  SparsityPattern sparsity_pattern;
  SparseMatrix<double> mass_matrix;
  SparseMatrix<double> dg_matrix;
  SparseMatrix<double> system_matrix;

  Vector<double> current_solution;
  Vector<double> previous_solution;

  Vector<double> system_rhs;
  const int expansion_order = 0;
  const unsigned int num_exp_coefficients;

  // parameters of the time stepping method
  double time_step = 1. / 128;
  double time = 0.;
  double final_time = .4;
  unsigned int time_step_number = 0;

  // scattering frequency
  double scattering_frequency = 1.;
  // particle
  ParticleProperties particle_properties;

  // Number of refinements
  unsigned int num_refinements = 5;
  // For the moment, I will use a quadratic domain
  // Rectangular domain
  // Point<dim> left_bottom;
  // Point<dim> right_top;

  // Map between i and l,m,s (implemented in constructor)
  std::vector<std::array<unsigned int, 3>> lms_indices;
  TimerOutput timer;
};

template <TermFlags flags, int dim>
VFPEquationSolver<flags, dim>::VFPEquationSolver(ParameterHandler &prm,
                                                 unsigned int polynomial_degree,
                                                 int order)
    : parameter_handler(prm),
      dof_handler(triangulation),
      mapping(),
      fe(FE_DGQ<dim>(polynomial_degree), (order + 1) * (order + 1)),
      quadrature(fe.tensor_degree() + 1),
      quadrature_face(fe.tensor_degree() + 1),
      expansion_order{order},
      num_exp_coefficients{
          static_cast<unsigned int>((order + 1) * (order + 1))},
      particle_properties(prm),
      lms_indices((order + 1) * (order + 1)),
      timer(std::cout, TimerOutput::summary, TimerOutput::wall_times) {
  // Create the index order
  for (int s = 0, idx = 0; s <= 1; ++s) {
    for (int l = 0; l <= expansion_order; ++l) {
      for (int m = l; m >= s; --m) {
        idx = l * (l + 1) - (s ? -1. : 1.) * m;
        lms_indices[idx][0] = l;
        lms_indices[idx][1] = m;
        lms_indices[idx][2] = s;
      }
    }
  }
  parameter_handler.enter_subsection("Mesh");
  { num_refinements = parameter_handler.get_integer("Number of refinements"); }
  parameter_handler.leave_subsection();
  parameter_handler.enter_subsection("Time stepping");
  {
    time_step = parameter_handler.get_double("Time step size");
    final_time = parameter_handler.get_double("Final time");
  }
  parameter_handler.leave_subsection();
  parameter_handler.enter_subsection("Physical parameters");
  {
    scattering_frequency = parameter_handler.get_double("Scattering frequency");
  }
  parameter_handler.leave_subsection();
}

template <TermFlags flags, int dim>
void VFPEquationSolver<flags, dim>::run() {
  // std::cout << particle_properties;
  output_compile_time_parameters();
  output_index_order();
  make_grid();
  setup_pde_system();
  setup_system();
  prepare_upwind_fluxes();
  assemble_mass_matrix();
  assemble_dg_matrix(time);
  project_initial_condition();

  current_solution = previous_solution;
  output_results();

  time += time_step;
  ++time_step_number;

  std::cout << "The time stepping loop is entered: \n";
  for (; time <= final_time; time += time_step, ++time_step_number) {
    std::cout << "	Time step " << time_step_number << " at t = " << time
              << "\n";
    // Time stepping method
    // theta_method(0.5);
    // explicit_runge_kutta();
    low_storage_explicit_runge_kutta();
    // NOTE: I cannot create TimerOutput::Scope inside output_results(),
    // because it is declared const.
    {
      TimerOutput::Scope timer_section(timer, "Output");
      output_results();
    }
    // Update solution
    previous_solution = current_solution;
  }
  std::cout << "The simulation ended. \n";
}

template <TermFlags flags, int dim>
void VFPEquationSolver<flags, dim>::make_grid() {
  TimerOutput::Scope timer_section(timer, "Grid setup");

  GridGenerator::hyper_cube(triangulation, -5., 5.);
  triangulation.refine_global(num_refinements);

  // std::ofstream out("grid.vtk");
  // GridOut grid_out;
  // grid_out.write_vtk(triangulation, out);
  // std::cout << "	Grid written to grid.vtk"
  //           << "\n";
  std::cout << "The grid was created: \n"
            << "	Number of active cells: "
            << triangulation.n_active_cells() << "\n";
}

template <TermFlags flags, int dim>
void VFPEquationSolver<flags, dim>::setup_pde_system() {
  TimerOutput::Scope timer_section(timer, "PDE system");

  // The matrix products (e.g. A_x * A_y) only yield the correct system if we
  // compute the matrices for expansion_order + 1 and later shrink them to
  // expansion_order
  unsigned int matrix_size = (expansion_order + 2) * (expansion_order + 2);

  advection_matrices.resize(3);
  for (auto &advection_matrix : advection_matrices)
    advection_matrix.reinit(matrix_size);

  generator_rotation_matrices.resize(3);
  for (auto &generator_matrix : generator_rotation_matrices)
    generator_matrix.reinit(matrix_size);

  collision_matrix.reinit(matrix_size);
  for (int s = 0; s <= 1; ++s) {
    for (int l = 0, i = 0; l <= expansion_order + 1; ++l) {
      for (int m = l; m >= s; --m) {
        i = l * (l + 1) - (s ? -1. : 1.) * m;  // (-1)^s
        for (int s_prime = 0; s_prime <= 1; ++s_prime) {
          for (int l_prime = 0, j = 0; l_prime <= expansion_order + 1;
               ++l_prime) {
            for (int m_prime = l_prime; m_prime >= s_prime; --m_prime) {
              j = l_prime * (l_prime + 1) - (s_prime ? -1. : 1.) * m_prime;
              // Ax
              if (l + 1 == l_prime && m == m_prime && s == s_prime)
                advection_matrices[0].set(
                    i, j,
                    std::sqrt(((l - m + 1.) * (l + m + 1.)) /
                              ((2. * l + 3.) * (2 * l + 1.))));
              if (l - 1 == l_prime && m == m_prime && s == s_prime)
                advection_matrices[0].set(
                    i, j,
                    std::sqrt(((l - m) * (l + m)) /
                              ((2. * l + 1.) * (2. * l - 1.))));
              // Ay
              if ((l + 1) == l_prime && (m + 1) == m_prime && s == s_prime)
                advection_matrices[1].set(
                    i, j,
                    -0.5 * std::sqrt(((l + m + 1.) * (l + m + 2.)) /
                                     ((2. * l + 3.) * (2. * l + 1.))));
              if ((l - 1) == l_prime && m + 1 == m_prime && s == s_prime)
                advection_matrices[1].set(
                    i, j,
                    0.5 * std::sqrt(((l - m - 1.) * (l - m)) /
                                    ((2. * l + 1.) * (2. * l - 1.))));
              if ((l + 1) == l_prime && (m - 1) == m_prime && s == s_prime)
                advection_matrices[1].set(
                    i, j,
                    0.5 * std::sqrt(((l - m + 1.) * (l - m + 2.)) /
                                    ((2. * l + 3.) * (2 * l + 1.))));
              if ((l - 1) == l_prime && (m - 1) == m_prime && s == s_prime)
                advection_matrices[1].set(
                    i, j,
                    -0.5 * std::sqrt(((l + m - 1.) * (l + m)) /
                                     ((2. * l + 1.) * (2. * l - 1.))));
              // Az
              // l + 1 = l_prime and s = 1 and s_prime = 0
              if (l + 1 == l_prime && m + 1 == m_prime && s == 1 &&
                  s_prime == 0)
                advection_matrices[2].set(
                    i, j,
                    0.5 * (m_prime == 0 ? 1 / std::sqrt(2) : 1.) *
                        std::sqrt(((l + m + 1.) * (l + m + 2.)) /
                                  ((2 * l + 3.) * (2 * l + 1.))));
              if (l + 1 == l_prime && m - 1 == m_prime && s == 1 &&
                  s_prime == 0)
                advection_matrices[2].set(
                    i, j,
                    0.5 * (m_prime == 0 ? 1 / std::sqrt(2) : 1.) *
                        std::sqrt(((l - m + 1.) * (l - m + 2.)) /
                                  ((2 * l + 3.) * (2 * l + 1.))));
              // NOTE: The correction are now directly included (-> new way to
              // compute the real matrices)
              if (l + 1 == l_prime && m == 0 && m_prime == 1 && s == 1 &&
                  s_prime == 0)
                advection_matrices[2](i, j) -=
                    0.5 * std::sqrt(((l - m + 1.) * (l - m + 2.)) /
                                    ((2 * l + 3.) * (2 * l + 1.)));
              if (l + 1 == l_prime && m == 1 && m_prime == 0 && s == 1 &&
                  s_prime == 0)
                advection_matrices[2](i, j) +=
                    0.5 * 1. / std::sqrt(2) *
                    std::sqrt(((l - m + 1.) * (l - m + 2.)) /
                              ((2 * l + 3.) * (2 * l + 1.)));

              // l+ 1 = l_prime and s = 0 and s_prime = 1
              if (l + 1 == l_prime && m + 1 == m_prime && s == 0 &&
                  s_prime == 1)
                advection_matrices[2].set(
                    i, j,
                    -0.5 * (m == 0 ? 1 / std::sqrt(2) : 1.) *
                        std::sqrt(((l + m + 1.) * (l + m + 2.)) /
                                  ((2 * l + 3.) * (2 * l + 1.))));
              if (l + 1 == l_prime && m - 1 == m_prime && s == 0 &&
                  s_prime == 1)
                advection_matrices[2].set(
                    i, j,
                    -0.5 * (m == 0 ? 1 / std::sqrt(2) : 1.) *
                        std::sqrt(((l - m + 1.) * (l - m + 2.)) /
                                  ((2 * l + 3.) * (2 * l + 1.))));
              // Corrections
              if (l + 1 == l_prime && m == 0 && m_prime == 1 && s == 0 &&
                  s_prime == 1)
                advection_matrices[2](i, j) -=
                    0.5 * 1. / std::sqrt(2) *
                    std::sqrt(((l - m + 1.) * (l - m + 2.)) /
                              ((2 * l + 3.) * (2 * l + 1.)));
              if (l + 1 == l_prime && m == 1 && m_prime == 0 && s == 0 &&
                  s_prime == 1)
                advection_matrices[2](i, j) +=
                    0.5 * std::sqrt(((l - m + 1.) * (l - m + 2.)) /
                                    ((2 * l + 3.) * (2 * l + 1.)));

              // l - 1 = l_prime and s = 1 and s_prime = 0
              if (l - 1 == l_prime && m + 1 == m_prime && s == 1 &&
                  s_prime == 0)
                advection_matrices[2].set(
                    i, j,
                    -0.5 * (m_prime == 0 ? 1 / std::sqrt(2) : 1.) *
                        std::sqrt(((l - m - 1.) * (l - m)) /
                                  ((2 * l + 1.) * (2 * l - 1.))));
              if (l - 1 == l_prime && m - 1 == m_prime && s == 1 &&
                  s_prime == 0)
                advection_matrices[2].set(
                    i, j,
                    -0.5 * (m_prime == 0 ? 1 / std::sqrt(2) : 1.) *
                        std::sqrt(((l + m - 1.) * (l + m)) /
                                  ((2 * l + 1.) * (2 * l - 1.))));
              // Corrections
              if (l - 1 == l_prime && m == 0 && m_prime == 1 && s == 1 &&
                  s_prime == 0)
                advection_matrices[2](i, j) +=
                    0.5 * std::sqrt(((l + m - 1.) * (l + m)) /
                                    ((2 * l + 1.) * (2 * l - 1.)));
              if (l - 1 == l_prime && m == 1 && m_prime == 0 && s == 1 &&
                  s_prime == 0)
                advection_matrices[2](i, j) -=
                    0.5 * 1. / std::sqrt(2) *
                    std::sqrt(((l + m - 1.) * (l + m)) /
                              ((2 * l + 1.) * (2 * l - 1.)));

              // l - 1 = l_prime and s = 0 and s_prime = 1
              if (l - 1 == l_prime && m + 1 == m_prime && s == 0 &&
                  s_prime == 1)
                advection_matrices[2].set(
                    i, j,
                    0.5 * (m == 0 ? 1 / std::sqrt(2) : 1.) *
                        std::sqrt(((l - m - 1.) * (l - m)) /
                                  ((2 * l + 1.) * (2 * l - 1.))));
              if (l - 1 == l_prime && m - 1 == m_prime && s == 0 &&
                  s_prime == 1)
                advection_matrices[2].set(
                    i, j,
                    0.5 * (m == 0 ? 1 / std::sqrt(2) : 1.) *
                        std::sqrt(((l + m - 1.) * (l + m)) /
                                  ((2 * l + 1.) * (2 * l - 1.))));
              // Corrections
              if (l - 1 == l_prime && m == 0 && m_prime == 1 && s == 0 &&
                  s_prime == 1)
                advection_matrices[2](i, j) +=
                    0.5 * 1. / std::sqrt(2) *
                    std::sqrt(((l + m - 1.) * (l + m)) /
                              ((2 * l + 1.) * (2 * l - 1.)));
              if (l - 1 == l_prime && m == 1 && m_prime == 0 && s == 0 &&
                  s_prime == 1)
                advection_matrices[2](i, j) -=
                    0.5 * std::sqrt(((l + m - 1.) * (l + m)) /
                                    ((2 * l + 1.) * (2 * l - 1.)));

              // Omega_x
              if (l == l_prime && m == m_prime && s == 0 && s_prime == 1) {
                generator_rotation_matrices[0].set(i, j, 1. * m);
                generator_rotation_matrices[0].set(
                    j, i,
                    -1. * m);  // Omega matrices are anti-symmetric
              }
              // Omega_y
              if (l == l_prime && (m + 1) == m_prime && s == 0 &&
                  s_prime == 1) {
                generator_rotation_matrices[1].set(
                    i, j, 0.5 * std::sqrt((l + m + 1.) * (l - m)));
                generator_rotation_matrices[1].set(
                    j, i, -0.5 * std::sqrt((l + m + 1.) * (l - m)));
              }
              if (l == l_prime && (m - 1) == m_prime && s == 0 &&
                  s_prime == 1) {
                generator_rotation_matrices[1].set(
                    i, j, 0.5 * std::sqrt((l - m + 1.) * (l + m)));
                generator_rotation_matrices[1].set(
                    j, i, -0.5 * std::sqrt((l - m + 1.) * (l + m)));
              }
              // Omega_z
              if (l == l_prime && (m + 1) == m_prime && s == s_prime) {
                generator_rotation_matrices[2].set(
                    i, j, -0.5 * std::sqrt((l + m + 1.) * (l - m)));
              }
              if (l == l_prime && (m - 1) == m_prime && s == s_prime) {
                generator_rotation_matrices[2].set(
                    i, j, 0.5 * std::sqrt((l - m + 1.) * (l + m)));
              }
              // C
              if (l == l_prime && m == m_prime && s == s_prime) {
                collision_matrix[i] = 0.5 * scattering_frequency * l * (l + 1.);
              }
            }
          }
        }
      }
    }
  }
  // TODO: Remove the correction part and derive formulas as you did for Az

  // Edit entries of the matrices around m=0 and s=0. After the transformation
  // they fall out of the pattern of the matrix elements, namely they differ by
  // a factor of 2^(1/2).
  //
  // NOTE: These corrections were not included in the above loops, because some
  // of them were overwritten. This is an effect of the s loops being outside
  // the l and m loops.
  if (expansion_order > 0) {
    for (int l = 0; l <= expansion_order + 1; ++l) {
      // Special cases for Omega
      if (l > 0) {
        // l == l_prime, m = 0, s = 0 and m_prime = 1 and s_prime = 1
        generator_rotation_matrices[1](l * (l + 1), l * (l + 1) + 1) =
            std::sqrt(2) *
            generator_rotation_matrices[1](l * (l + 1), l * (l + 1) + 1);
        // l == l_prime, m = 1, s = 1 and m_prime = 0 and s_prime = 0
        generator_rotation_matrices[1](l * (l + 1) + 1, l * (l + 1)) =
            std::sqrt(2) *
            generator_rotation_matrices[1](l * (l + 1) + 1, l * (l + 1));

        // l == l_prime, m = 0, s = 0 and m_prime = 1 and s_prime = 0
        generator_rotation_matrices[2](l * (l + 1), l * (l + 1) - 1) =
            std::sqrt(2) *
            generator_rotation_matrices[2](l * (l + 1), l * (l + 1) - 1);
        // // l == l_prime, m = 1, s = 0 and m_prime = 0 and s_prime = 0
        generator_rotation_matrices[2](l * (l + 1) - 1, l * (l + 1)) =
            std::sqrt(2) *
            generator_rotation_matrices[2](l * (l + 1) - 1, l * (l + 1));
      }
      // Special cases for A_y (necessary for every value of l)
      // Above the diagonal
      // l + 1 = l_prime, m = 0, s = 0, and m_prime = 1, s_prime = 0
      if (l != expansion_order + 1)
        advection_matrices[1](l * (l + 1), (l + 1) * (l + 2) - 1) =
            std::sqrt(2) *
            advection_matrices[1](l * (l + 1), (l + 2) * (l + 1) - 1);
      // l + 1 = l_prime, m = 1, s = 0, and m_prime = 0, s_prime = 0
      if (l != 0 && l != expansion_order + 1)
        advection_matrices[1](l * (l + 1) - 1, (l + 1) * (l + 2)) =
            std::sqrt(2) *
            advection_matrices[1](l * (l + 1) - 1, (l + 2) * (l + 1));
      // Below the diagonal
      // l - 1 = l_prime, m = 0, s = 0, and m_prime = 1, s_prime = 0
      if (l > 1)
        advection_matrices[1](l * (l + 1), l * (l - 1) - 1) =
            std::sqrt(2) * advection_matrices[1](l * (l + 1), l * (l - 1) - 1);
      // l - 1 = l_prime , m = 1, s = 0 and m_prime = 0, s_prime = 0
      if (l != 0)
        advection_matrices[1](l * (l + 1) - 1, l * (l - 1)) =
            std::sqrt(2) * advection_matrices[1](l * (l + 1) - 1, l * (l - 1));
    }
  }

  // Compute the matrix products
  Ap_xx.reinit(matrix_size);
  advection_matrices[0].mmult(Ap_xx, advection_matrices[0]);
  Ap_xy.reinit(matrix_size);
  advection_matrices[0].mmult(Ap_xy, advection_matrices[1]);
  Ap_xz.reinit(matrix_size);
  advection_matrices[0].mmult(Ap_xz, advection_matrices[2]);

  Ap_yy.reinit(matrix_size);
  advection_matrices[1].mmult(Ap_yy, advection_matrices[1]);
  Ap_yz.reinit(matrix_size);
  advection_matrices[1].mmult(Ap_yz, advection_matrices[2]);

  Ap_zz.reinit(matrix_size);
  advection_matrices[2].mmult(Ap_zz, advection_matrices[2]);

  // Shrink the matrices such that they agree with order of the expansion
  for (auto &advection_matrix : advection_matrices)
    advection_matrix.grow_or_shrink(num_exp_coefficients);

  for (auto &generator_matrix : generator_rotation_matrices)
    generator_matrix.grow_or_shrink(num_exp_coefficients);

  collision_matrix.grow_or_shrink(num_exp_coefficients);

  Ap_xx.grow_or_shrink(num_exp_coefficients);
  Ap_xy.grow_or_shrink(num_exp_coefficients);
  Ap_xz.grow_or_shrink(num_exp_coefficients);

  Ap_yy.grow_or_shrink(num_exp_coefficients);
  Ap_yz.grow_or_shrink(num_exp_coefficients);

  Ap_zz.grow_or_shrink(num_exp_coefficients);

  // std::cout << "A_x: \n";
  // advection_matrices[0].print_formatted(std::cout);

  // std::cout << "A_y: \n";
  // advection_matrices[1].print_formatted(std::cout);

  // std::cout << "A_z: \n";
  // advection_matrices[2].print_formatted(std::cout);

  // std::cout << "Omega_x: \n";
  // omega_matrices[0].print_formatted(std::cout);

  // std::cout << "Omega_y: \n";
  // omega_matrices[1].print_formatted(std::cout);

  // std::cout << "Omega_z: \n";
  // omega_matrices[2].print_formatted(std::cout);

  // std::cout << "Ap_xx: \n";
  // Ap_xx.print_formatted(std::cout);

  // std::cout << "Ap_xy: \n";
  // Ap_xy.print_formatted(std::cout);

  // std::cout << "Ap_xz: \n";
  // Ap_xz.print_formatted(std::cout);

  // std::cout << "Ap_yy: \n";
  // Ap_yy.print_formatted(std::cout);

  // std::cout << "Ap_yz: \n";
  // Ap_yz.print_formatted(std::cout);

  // std::cout << "Ap_zz: \n";
  // Ap_zz.print_formatted(std::cout);

  // std::cout << "R: \n";
  // R.print(std::cout);
}

template <TermFlags flags, int dim>
void VFPEquationSolver<flags, dim>::prepare_upwind_fluxes() {
  TimerOutput::Scope timer_section(timer, "Prepare Upwind flux");
  // NOTE: The matrix gets destroyed, when the eigenvalues are computed. This
  // requires to copy it
  LAPACKFullMatrix<double> matrix_copy;
  double tolerance = 1.e-15;
  eigenvectors_advection_matrices.resize(3);

  // The eigenvalues of A_x are also the eigenvalues of A_y and A_z. The
  // eigenvalues of A_x are the roots of the associated legendre polynomials.
  // Since the associated Legendre Polynomials are orthogonal polynomials, it
  // can be shown that all roots are contained in in the intervall [-1., 1.].
  // Physically this makes sense, because c = 1 and the eigenvalues encode the
  // speed of 'information' transport.
  matrix_copy = advection_matrices[0];  // A_x
  matrix_copy.compute_eigenvalues_symmetric(-1.2, 1.2, tolerance, eigenvalues,
                                            eigenvectors_advection_matrices[0]);
  // NOTE: The eigenvalues of A_x can computed very efficiently because A_x is
  // (when ordered correctly) a tridiagonal symmetric matrix, i.e. it is in
  // Hessenberg-Form. It seems that Lapack has function for this implemented but
  // there is no dealii interface to it.

  matrix_copy = advection_matrices[1];
  matrix_copy.compute_eigenvalues_symmetric(-1.2, 1.2, tolerance, eigenvalues,
                                            eigenvectors_advection_matrices[1]);
  matrix_copy = advection_matrices[2];
  matrix_copy.compute_eigenvalues_symmetric(-1.2, 1.2, tolerance, eigenvalues,
                                            eigenvectors_advection_matrices[2]);
  // TODO: Rotate the eigenvectors of A_x to get the eigenvectors of A_y and
  // A_z, i.e. implement the rotation matrices e^{-i\Omega_x pi/2} and e^{-i
  // \Omega_z pi/2} For now we will three times compute the same eigenvalues to
  // get the eigenvectors of A_y and A_z

  // NOTE: I actually thinkt that the next clean up step is not necessary
  // Remove eigenvalues, which are smaller than tolerance
  auto smaller_than_tolerance = [&tolerance](double value) {
    return (std::abs(value) < tolerance ? true : false);
  };
  std::replace_if(eigenvalues.begin(), eigenvalues.end(),
                  smaller_than_tolerance, 0.);
}

template <TermFlags flags, int dim>
void VFPEquationSolver<flags, dim>::compute_upwind_fluxes(
    const std::vector<Point<dim>> &q_points,
    const std::vector<Tensor<1, dim>> &normals,
    std::vector<FullMatrix<double>> &positive_flux_matrices,
    std::vector<FullMatrix<double>> &negative_flux_matrices) {
  // Determine if we are computing the flux of an interior face whose normal
  // points into the x,y,z or p direction. Since we are using a rectangular
  // grid, the normal is the same for all quadrature points and it points either
  // in the x, y, z or p direction. Hence it can determined outside the loop
  // over the quadrature points
  enum Coordinate { x, y, z, p, none = 1000 };
  Coordinate component = none;  // Produces an error if not overwritten
  if (std::abs(normals[0][x]) == 1.)
    component = x;  // std::abs is nessecary
                    // because of the boundary
                    // faces, i.e. <-| n_k = -1
                    // for some faces.
  if constexpr (dim >= 2) {
    if (std::abs(normals[0][y]) == 1.) component = y;
  }
  if constexpr (dim == 3) {
    if (std::abs(normals[0][z]) == 1.) component = z;
  }
  // TODO: Insert a if constexpr for the p case

  // Get the value of the velocity field at every quadrature point
  BackgroundVelocityField<dim> velocity_field;
  std::vector<Vector<double>> velocities(q_points.size(), Vector<double>(dim));
  velocity_field.vector_value_list(q_points, velocities);

  for (unsigned int q_index = 0; q_index < q_points.size(); ++q_index) {
    // Create a copy of the eigenvalues to compute the positive and negative
    // fluxes at interface point q at time t
    Vector<double> lambda(eigenvalues);
    // Multiply the eigenvalues of A_i with the particle velocity
    lambda *= particle_properties.particle_velocity;
    // TODO: Write an if constexpr for the transport only case (fixed velocity),
    // else compute the particle velocity from the Point p in phase space

    // Add the velocities to the eigenvalues
    lambda.add(velocities[q_index][component]);
    // std::cout << "eigenvalues plus velocity: "
    //           << "\n";
    // lambda.print(std::cout);
    //
    lambda *= normals[q_index][component];
    Vector<double> positive_lambda(eigenvalues.size());
    Vector<double> negative_lambda(eigenvalues.size());

    std::replace_copy_if(
        lambda.begin(), lambda.end(), positive_lambda.begin(),
        std::bind(std::less<double>(), std::placeholders::_1, 0.), 0.);
    std::replace_copy_if(
        lambda.begin(), lambda.end(), negative_lambda.begin(),
        std::bind(std::greater<double>(), std::placeholders::_1, 0.), 0.);

    FullMatrix<double> lambda_plus_matrix(num_exp_coefficients);
    FullMatrix<double> lambda_minus_matrix(num_exp_coefficients);

    for (unsigned int i = 0; i < num_exp_coefficients; ++i) {
      lambda_plus_matrix(i, i) = positive_lambda[i];
      lambda_minus_matrix(i, i) = negative_lambda[i];
    }
    // Return
    positive_flux_matrices[q_index].triple_product(
        lambda_plus_matrix, eigenvectors_advection_matrices[component],
        eigenvectors_advection_matrices[component], false, true);

    negative_flux_matrices[q_index].triple_product(
        lambda_minus_matrix, eigenvectors_advection_matrices[component],
        eigenvectors_advection_matrices[component], false, true);
  }
}

template <TermFlags flags, int dim>
void VFPEquationSolver<flags, dim>::setup_system() {
  TimerOutput::Scope timer_section(timer, "FE system");

  dof_handler.distribute_dofs(fe);
  const unsigned int n_dofs = dof_handler.n_dofs();
  std::cout << "The degrees of freedom were distributed: \n"
            << "	Number of degrees of freedom: " << n_dofs << "\n";

  DynamicSparsityPattern dsp(n_dofs);
  DoFTools::make_flux_sparsity_pattern(dof_handler, dsp);
  sparsity_pattern.copy_from(dsp);

  dg_matrix.reinit(sparsity_pattern);
  // NOTE: DealII does not allow to use different sparsity patterns for
  // matrices, which you would like to add. Even though the the mass matrix
  // differs from the dg matrix.
  mass_matrix.reinit(sparsity_pattern);
  system_matrix.reinit(sparsity_pattern);

  // Vectors
  current_solution.reinit(n_dofs);
  previous_solution.reinit(n_dofs);
  system_rhs.reinit(n_dofs);

  // This is an rather obscure line. I do not know why I need it. (cf. example
  // 23)
  constraints.close();
}

template <TermFlags flags, int dim>
void VFPEquationSolver<flags, dim>::assemble_mass_matrix() {
  TimerOutput::Scope timer_section(timer, "Mass matrix");

  std::cout << "The mass matrix is assembled. \n";
  FEValues<dim> fe_v(
      mapping, fe, quadrature,
      update_values | update_quadrature_points | update_JxW_values);

  const unsigned int n_dofs = fe.n_dofs_per_cell();

  FullMatrix<double> cell_mass_matrix(n_dofs, n_dofs);

  std::vector<types::global_dof_index> local_dof_indices(n_dofs);

  for (const auto &cell : dof_handler.active_cell_iterators()) {
    cell_mass_matrix = 0;
    fe_v.reinit(cell);

    const std::vector<double> &JxW = fe_v.get_JxW_values();

    for (const unsigned int q_index : fe_v.quadrature_point_indices()) {
      for (unsigned int i : fe_v.dof_indices()) {
        const unsigned int component_i = fe.system_to_component_index(i).first;
        for (unsigned int j : fe_v.dof_indices()) {
          const unsigned int component_j =
              fe.system_to_component_index(j).first;
          // mass matrix
          cell_mass_matrix(i, j) +=
              (component_i == component_j
                   ? fe_v.shape_value(i, q_index) * fe_v.shape_value(j, q_index)
                   : 0) *
              JxW[q_index];
        }
      }
    }
    cell->get_dof_indices(local_dof_indices);

    constraints.distribute_local_to_global(cell_mass_matrix, local_dof_indices,
                                           mass_matrix);
  }
}
template <TermFlags flags, int dim>
void VFPEquationSolver<flags, dim>::assemble_dg_matrix(
    unsigned int evaluation_time) {
  TimerOutput::Scope timer_section(timer, "DG matrix");

  std::cout << "The DG matrix is assembled. \n";
  /*
    What kind of loops are there ?
    1. Loop over all cells (this happens inside the mesh_loop)
    2. Loop over the degrees of freedom on each cell
    - the metod system_to_componet_index() returns the index of the non-zero
    component of the vector-valued shape function which corresponds to the
    indices (l,m,s)
  */
  using Iterator = typename DoFHandler<dim>::active_cell_iterator;
  BackgroundVelocityField<dim> background_velocity_field;
  background_velocity_field.set_time(evaluation_time);
  MagneticField<dim> magnetic_field(parameter_handler);
  magnetic_field.set_time(evaluation_time);
  // I do not no the meaning of the following "const" specifier
  const auto cell_worker = [&](const Iterator &cell,
                               ScratchData<dim> &scratch_data,
                               CopyData &copy_data) {
    FEValues<dim> &fe_v = scratch_data.fe_values;
    // reinit cell
    fe_v.reinit(cell);

    const unsigned int n_dofs = fe_v.get_fe().n_dofs_per_cell();
    // reinit the matrices cell_matrix_dg, cell_mass_matrix
    copy_data.reinit(cell, n_dofs);

    const std::vector<Point<dim>> &q_points = fe_v.get_quadrature_points();
    const std::vector<double> &JxW = fe_v.get_JxW_values();

    // NOTE: Second argument constructs empty vectors with 2 values. They are
    // copied q_points.size() times.
    std::vector<Vector<double>> velocities(q_points.size(), Vector<double>(dim));
    background_velocity_field.vector_value_list(q_points, velocities);
    std::vector<double> div_velocities(q_points.size());
    background_velocity_field.divergence_list(q_points, div_velocities);

    std::vector<Vector<double>> magnetic_field_values(q_points.size(),
                                                      Vector<double>(3));
    magnetic_field.vector_value_list(q_points, magnetic_field_values);

    for (unsigned int i : fe_v.dof_indices()) {
      const unsigned int component_i =
          fe_v.get_fe().system_to_component_index(i).first;
      for (unsigned int j : fe_v.dof_indices()) {
        const unsigned int component_j =
            fe_v.get_fe().system_to_component_index(j).first;
        for (const unsigned int q_index : fe_v.quadrature_point_indices()) {
          if (component_i == component_j) {
            if constexpr ((flags & TermFlags::reaction) != TermFlags::none) {
              // 0.5 * scattering_frequency * l(l+1) * \phi_i * \phi_j
              copy_data.cell_dg_matrix(i, j) +=
                  collision_matrix[component_i] * fe_v.shape_value(i, q_index) *
                  fe_v.shape_value(j, q_index) * JxW[q_index];
            }
            if constexpr ((flags & TermFlags::advection) != TermFlags::none) {
              // - [\partial_x(u_x\delta_ij + Ax_ij) + \partial_y(u_y\delta_ij +
              // - Ay_ij) ] \phi_i \phi_j where \partial_x/y Ax/y_ij = 0
              copy_data.cell_dg_matrix(i, j) -=
                  div_velocities[q_index] * fe_v.shape_value(i, q_index) *
                  fe_v.shape_value(j, q_index) * JxW[q_index];
            }
          }
          if constexpr ((flags & TermFlags::advection) != TermFlags::none) {
            // -[partial_x \phi_i * (u_x \delta_ij + Ax_ij)
            //   + partial_y \phi_i * (u_y \delta_ij + Ay_ij)] * phi_j
            if (component_i == component_j) {
              for (unsigned int coordinate = 0; coordinate < dim;
                   ++coordinate) {
                copy_data.cell_dg_matrix(i, j) -=
                    fe_v.shape_grad(i, q_index)[coordinate] *
                    velocities[q_index][coordinate] *
                    fe_v.shape_value(j, q_index) * JxW[q_index];
              }

            } else {
              // NOTE: Many zerso are added here, because the matrices Ax, Ay
              // and A_z are sparse. TODO: Performance check. If too bad, return
              // to the strategy, which was used in v0.6.5
              for (unsigned int coordinate = 0; coordinate < dim; ++coordinate)
                copy_data.cell_dg_matrix(i, j) -=
                    fe_v.shape_grad(i, q_index)[coordinate] *
                    particle_properties.particle_velocity *
                    advection_matrices[coordinate](component_i, component_j) *
                    fe_v.shape_value(j, q_index) * JxW[q_index];
            }
          }
          if constexpr ((flags & TermFlags::magnetic) != TermFlags::none) {
            // NOTE: All three components of the B-Field are included not
            // matter, which dimension of the configuration space is considered
            for (unsigned int coordinate = 0; coordinate < 3; ++coordinate)
              copy_data.cell_dg_matrix(i, j) -=
                  fe_v.shape_value(i, q_index) *
                  magnetic_field_values[q_index][coordinate] *
                  generator_rotation_matrices[coordinate](component_i,
                                                          component_j) *
                  fe_v.shape_value(j, q_index) * JxW[q_index];
          }
        }
      }
    }
  };
  // assemble boundary face terms
  const auto boundary_worker = [&](const Iterator &cell,
                                   const unsigned int &face_no,
                                   ScratchData<dim> &scratch_data,
                                   CopyData &copy_data) {
    scratch_data.fe_values_face.reinit(cell, face_no);
    const FEFaceValuesBase<dim> &fe_face_v = scratch_data.fe_values_face;
    // Every shape function on the cell could contribute to the face integral,
    // hence n_facet_dofs = n_dofs_per_cell
    const unsigned int n_facet_dofs = fe_face_v.get_fe().n_dofs_per_cell();
    // NOTE: copy_data is not reinitialised, the cell_workers contribution to
    // the cell_dg_matrix should not be deleted

    // NOTE: Currently I do not the velocity field here. The upwind flux may
    // depend on the quadrature point in the future.
    // const std::vector<Point<dim>> &q_points =
    // fe_face_v.get_quadrature_points(); std::vector<Vector<double>>
    // velocities(q_points.size(), Vector<double>(2));
    // background_velocity_field.vector_value_list(q_points, velocities);

    const std::vector<Point<dim>> &q_points = fe_face_v.get_quadrature_points();
    const std::vector<double> &JxW = fe_face_v.get_JxW_values();
    const std::vector<Tensor<1, dim>> &normals = fe_face_v.get_normal_vectors();
    std::vector<FullMatrix<double>> positive_flux_matrices(
        q_points.size(), FullMatrix<double>(num_exp_coefficients));
    std::vector<FullMatrix<double>> negative_flux_matrices(
        q_points.size(), FullMatrix<double>(num_exp_coefficients));

    compute_upwind_fluxes(q_points, normals, positive_flux_matrices,
                          negative_flux_matrices);
    for (unsigned int i = 0; i < n_facet_dofs; ++i) {
      const unsigned int component_i =
          fe_face_v.get_fe().system_to_component_index(i).first;
      for (unsigned int j = 0; j < n_facet_dofs; ++j) {
        const unsigned int component_j =
            fe_face_v.get_fe().system_to_component_index(j).first;
        for (unsigned int q_index : fe_face_v.quadrature_point_indices()) {
          if constexpr ((flags & TermFlags::advection) != TermFlags::none) {
            // Outflow boundary: Everyhing with a positive flux along the
            // direction of the normal leaves the boundary
            copy_data.cell_dg_matrix(i, j) +=
                fe_face_v.shape_value(i, q_index) *
                positive_flux_matrices[q_index](component_i, component_j) *
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
                               ScratchData<dim> &scratch_data,
                               CopyData &copy_data) {
    // NOTE: The flag MeshWorker::assemble_own_interior_faces_both will not
    // work for this face_worker. It implicitly assumes that faces are only
    // assembled once. And hence subface_no, can be ignored
    (void)subface_no;  // suppress compiler warning

    FEFaceValues<dim> &fe_v_face = scratch_data.fe_values_face;
    fe_v_face.reinit(cell, face_no);

    // NOTE: I do not know how to initialise this reference whose value
    // depends on the value of neighbor_subface_no. Subfaces only exists if
    // the triangulation was refined differently in different parts of the
    // domain. Since I am currently testing, I am always refining globally,
    // so no subfaces exist and I just ignore this case. Hence
    // neighbor_subface_no can be ignored
    (void)neighbor_subface_no;
    FEFaceValues<dim> &fe_v_face_neighbor =
        scratch_data.fe_values_face_neighbor;
    fe_v_face_neighbor.reinit(neighbor_cell, neighbor_face_no);
    // Create an element at the end of the vector containig the face data
    copy_data.face_data.emplace_back();
    CopyDataFace &copy_data_face = copy_data.face_data.back();
    const unsigned int n_dofs = fe_v_face.get_fe().n_dofs_per_cell();
    copy_data_face.reinit(cell, neighbor_cell, n_dofs);

    const std::vector<Point<dim>> &q_points = fe_v_face.get_quadrature_points();
    const std::vector<double> &JxW = fe_v_face.get_JxW_values();
    const std::vector<Tensor<1, dim>> &normals = fe_v_face.get_normal_vectors();
    // For every interior face there is an in- and outflow represented by the
    // corresponding flux matrices
    std::vector<FullMatrix<double>> positive_flux_matrices(
        q_points.size(), FullMatrix<double>(num_exp_coefficients));
    std::vector<FullMatrix<double>> negative_flux_matrices(
        q_points.size(), FullMatrix<double>(num_exp_coefficients));

    compute_upwind_fluxes(q_points, normals, positive_flux_matrices,
                          negative_flux_matrices);

    for (unsigned int q_index : fe_v_face.quadrature_point_indices()) {
      // cell_dg_matrix_11
      for (unsigned int i : fe_v_face.dof_indices()) {
        const unsigned int component_i =
            fe_v_face.get_fe().system_to_component_index(i).first;
        for (unsigned int j : fe_v_face.dof_indices()) {
          unsigned int component_j =
              fe_v_face.get_fe().system_to_component_index(j).first;
          if constexpr ((flags & TermFlags::advection) != TermFlags::none) {
            copy_data_face.cell_dg_matrix_11(i, j) +=
                fe_v_face.shape_value(i, q_index) *
                positive_flux_matrices[q_index](component_i, component_j) *
                fe_v_face.shape_value(j, q_index) * JxW[q_index];
          }
        }
      }
      // cell_dg_matrix_12
      for (unsigned int i : fe_v_face.dof_indices()) {
        const unsigned int component_i =
            fe_v_face.get_fe().system_to_component_index(i).first;
        for (unsigned int j : fe_v_face_neighbor.dof_indices()) {
          unsigned int component_j =
              fe_v_face_neighbor.get_fe().system_to_component_index(j).first;
          if constexpr ((flags & TermFlags::advection) != TermFlags::none) {
            copy_data_face.cell_dg_matrix_12(i, j) -=
                fe_v_face_neighbor.shape_value(i, q_index) *
                positive_flux_matrices[q_index](component_i, component_j) *
                fe_v_face.shape_value(j, q_index) * JxW[q_index];
          }
        }
      }
      // cell_dg_matrix_21
      for (unsigned int i : fe_v_face_neighbor.dof_indices()) {
        const unsigned int component_i =
            fe_v_face_neighbor.get_fe().system_to_component_index(i).first;
        for (unsigned int j : fe_v_face.dof_indices()) {
          unsigned int component_j =
              fe_v_face.get_fe().system_to_component_index(j).first;
          if constexpr ((flags & TermFlags::advection) != TermFlags::none) {
            copy_data_face.cell_dg_matrix_21(i, j) +=
                fe_v_face.shape_value(i, q_index) *
                negative_flux_matrices[q_index](component_i, component_j) *
                fe_v_face_neighbor.shape_value(j, q_index) * JxW[q_index];
          }
        }
      }
      // cell_dg_matrix_22
      for (unsigned int i : fe_v_face_neighbor.dof_indices()) {
        const unsigned int component_i =
            fe_v_face_neighbor.get_fe().system_to_component_index(i).first;
        for (unsigned int j : fe_v_face_neighbor.dof_indices()) {
          unsigned int component_j =
              fe_v_face_neighbor.get_fe().system_to_component_index(j).first;
          if constexpr ((flags & TermFlags::advection) != TermFlags::none) {
            copy_data_face.cell_dg_matrix_22(i, j) -=
                fe_v_face_neighbor.shape_value(i, q_index) *
                negative_flux_matrices[q_index](component_i, component_j) *
                fe_v_face_neighbor.shape_value(j, q_index) * JxW[q_index];
          }
        }
      }
    }
  };
  // copier for the mesh_loop function
  const auto copier = [&](const CopyData &c) {
    constraints.distribute_local_to_global(c.cell_dg_matrix,
                                           c.local_dof_indices, dg_matrix);
    constraints.distribute_local_to_global(c.cell_mass_matrix,
                                           c.local_dof_indices, mass_matrix);
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

  ScratchData<dim> scratch_data(mapping, fe, quadrature, quadrature_face);
  CopyData copy_data;

  MeshWorker::mesh_loop(dof_handler.begin_active(), dof_handler.end(),
                        cell_worker, copier, scratch_data, copy_data,
                        MeshWorker::assemble_own_cells |
                            MeshWorker::assemble_boundary_faces |
                            MeshWorker::assemble_own_interior_faces_once,
                        boundary_worker, face_worker);
}

template <TermFlags flags, int dim>
void VFPEquationSolver<flags, dim>::project_initial_condition() {
  TimerOutput::Scope timer_section(timer, "Intial condition");

  // Create right hand side
  FEValues<dim> fe_v(
      mapping, fe, quadrature,
      update_values | update_quadrature_points | update_JxW_values);

  const unsigned int n_dofs = fe.n_dofs_per_cell();
  Vector<double> cell_rhs(n_dofs);

  std::vector<types::global_dof_index> local_dof_indices(n_dofs);

  for (const auto &cell : dof_handler.active_cell_iterators()) {
    cell_rhs = 0;
    fe_v.reinit(cell);

    const std::vector<Point<dim>> &q_points = fe_v.get_quadrature_points();
    const std::vector<double> &JxW = fe_v.get_JxW_values();

    // Initial values
    InitialValueFunction<dim> initial_value_function(expansion_order);
    std::vector<Vector<double>> initial_values(
        q_points.size(), Vector<double>(num_exp_coefficients));
    initial_value_function.vector_value_list(q_points, initial_values);

    for (const unsigned int q_index : fe_v.quadrature_point_indices()) {
      for (unsigned int i : fe_v.dof_indices()) {
        const unsigned int component_i = fe.system_to_component_index(i).first;
        cell_rhs(i) += fe_v.shape_value(i, q_index) *
                       initial_values[q_index][component_i] * JxW[q_index];
      }
    }
    cell->get_dof_indices(local_dof_indices);

    constraints.distribute_local_to_global(cell_rhs, local_dof_indices,
                                           system_rhs);
  }

  // Solve the system
  SolverControl solver_control(1000, 1e-12);
  SolverCG<Vector<double>> cg(solver_control);
  cg.solve(mass_matrix, previous_solution, system_rhs, PreconditionIdentity());

  // Reset system RHS
  system_rhs = 0;
}

template <TermFlags flags, int dim>
void VFPEquationSolver<flags, dim>::theta_method_solve_system() {
  SolverControl solver_control(1000, 1e-12);
  SolverRichardson<Vector<double>> solver(solver_control);

  PreconditionBlockSSOR<SparseMatrix<double>> preconditioner;
  preconditioner.initialize(system_matrix, fe.n_dofs_per_cell());
  solver.solve(system_matrix, current_solution, system_rhs, preconditioner);

  std::cout << "	Solver converged in " << solver_control.last_step()
            << " iterations."
            << "\n";
}

template <TermFlags flags, int dim>
void VFPEquationSolver<flags, dim>::theta_method(double theta) {
  TimerOutput::Scope timer_section(timer, "Theta method");

  Vector<double> tmp(current_solution.size());

  mass_matrix.vmult(system_rhs, previous_solution);
  dg_matrix.vmult(tmp, previous_solution);
  system_rhs.add(-time_step * (1 - theta), tmp);

  // Since the the dg_matrix depends on the velocity field and the the
  // velocity field may depend on time, it needs to reassembled every time
  // step. This is not true for the mass matrix ( but it may if the grid
  // adapts after a specified amount of time steps)

  // mass_matrix = 0; dg_matrix = 0;
  // assemble_system(time + time_step);

  // NOTE: Currently that is not necessary, because the velocity field is
  // constant. And I should check if I cannot update the system matrix without
  // reassembling in each time step.

  system_matrix.copy_from(mass_matrix);
  system_matrix.add(time_step * theta, dg_matrix);
  theta_method_solve_system();
}

template <TermFlags flags, int dim>
void VFPEquationSolver<flags, dim>::explicit_runge_kutta() {
  TimerOutput::Scope timer_section(timer, "ERK4");

  // ERK 4
  // Butcher's array
  Vector<double> a({0.5, 0.5, 1.});
  Vector<double> b({1. / 6, 1. / 3, 1. / 3, 1. / 6});
  // NOTE: c is only necessary if the velocity field and the magnetic field
  // are time dependent.
  // Vector<double> c({0., 0.5, 0.5, 1.});

  // Allocate storage for the four stages
  std::vector<Vector<double>> k(4, Vector<double>(current_solution.size()));
  // The mass matrix needs to be "inverted" in every stage
  SolverControl solver_control(1000, 1e-12);
  SolverCG<Vector<double>> cg(solver_control);
  // k_0
  dg_matrix.vmult(system_rhs, previous_solution);
  cg.solve(mass_matrix, k[0], system_rhs, PreconditionIdentity());
  std::cout << "	Stage s: " << 0 << "	Solver converged in "
            << solver_control.last_step() << " iterations."
            << "\n";
  current_solution.add(b[0] * time_step, k[0]);

  Vector<double> temp(current_solution.size());
  for (unsigned int s = 1; s < 4; ++s) {
    // NOTE: It would be nessary to reassemble the dg matrix twice for
    // (for half a time step and a whole time step) if the velocity field and
    // the magnetic field were time dependent.
    temp.add(-1., previous_solution, -a[s - 1] * time_step, k[s - 1]);
    // assemble_system(time + c[s]*time_step);
    dg_matrix.vmult(system_rhs, temp);
    cg.solve(mass_matrix, k[s], system_rhs, PreconditionIdentity());
    std::cout << "	Stage s: " << s << "	Solver converged in "
              << solver_control.last_step() << " iterations."
              << "\n";

    current_solution.add(b[s] * time_step, k[s]);

    // empty temp vector
    temp = 0;
  }
  // NOTE: It is not necessary to assemble the matrix again, because the last
  // element of the c vector is 1. (see low_storage_erk())
}

template <TermFlags flags, int dim>
void VFPEquationSolver<flags, dim>::low_storage_explicit_runge_kutta() {
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
  // Vector<double> c(
  //     {0., 1432997174477. / 9575080441755, 2526269341429. / 6820363962896,
  //      2006345519317. / 3224310063776, 2802321613138. / 2924317926251});

  // The mass matrix needs to be "inverted" in every stage
  SolverControl solver_control(1000, 1e-12);
  SolverCG<Vector<double>> cg(solver_control);

  Vector<double> k(current_solution.size());
  Vector<double> temp(current_solution.size());
  for (unsigned int s = 0; s < 5; ++s) {
    // assemble_system(time + c[s]*time_step);
    dg_matrix.vmult(system_rhs, current_solution);
    cg.solve(mass_matrix, temp, system_rhs, PreconditionIdentity());
    std::cout << "	Stage s: " << s << "	Solver converged in "
              << solver_control.last_step() << " iterations."
              << "\n";

    k.sadd(a[s], -time_step, temp);
    current_solution.add(b[s], k);
  }
  // assemble_system(time + time_step)
}

template <TermFlags flags, int dim>
void VFPEquationSolver<flags, dim>::output_results() const {
  DataOut<dim> data_out;
  data_out.attach_dof_handler(dof_handler);
  // Create a vector of strings with names for the components of the solution
  std::vector<std::string> component_names(num_exp_coefficients);

  for (unsigned int i = 0; i < num_exp_coefficients; ++i) {
    const std::array<unsigned int, 3> lms = lms_indices[i];
    component_names[i] = "f_" + std::to_string(lms[0]) +
                         std::to_string(lms[1]) + std::to_string(lms[2]);
  }

  data_out.add_data_vector(current_solution, component_names);

  data_out.build_patches();
  const std::string file_path = "results/solution" +
                                Utilities::int_to_string(time_step_number, 3) +
                                ".vtu";

  DataOutBase::VtkFlags vtk_flags;
  vtk_flags.compression_level = DataOutBase::VtkFlags::best_speed;
  data_out.set_flags(vtk_flags);

  std::ofstream output(file_path);
  data_out.write_vtu(output);
}

template <TermFlags flags, int dim>
void VFPEquationSolver<flags, dim>::output_compile_time_parameters() const {
  std::cout << "Compile-time parameters: "
            << "\n";
  std::cout << "	" << flags;
  std::cout << "	Dimension: " << dim << "\n\n";
}

template <TermFlags flags, int dim>
void VFPEquationSolver<flags, dim>::output_index_order() const {
  std::cout << "Ordering of the lms indices: "
            << "\n";

  for (std::array<unsigned int, 3> lms : lms_indices)
    std::cout << "	" << lms[0] << lms[1] << lms[2] << "\n";
  std::cout << "\n";
}
}  // namespace vfp_equation_solver

int main() {
  try {
    using namespace vfp_equation_solver;

    constexpr TermFlags flags =
        TermFlags::advection | TermFlags::magnetic | TermFlags::reaction;

    ParameterHandler parameter_handler;
    ParameterReader parameter_reader(parameter_handler);
    parameter_reader.read_parameters("vfp-equation.prm");
    std::cout << "Run-time parameters: "
              << "\n";
    parameter_handler.print_parameters(std::cout, ParameterHandler::ShortPRM);
    int expansion_order;
    parameter_handler.enter_subsection("Expansion");
    { expansion_order = parameter_handler.get_integer("Expansion order"); }
    parameter_handler.leave_subsection();
    unsigned int polynomial_degree;
    parameter_handler.enter_subsection("Finite element");
    { polynomial_degree = parameter_handler.get_integer("Polynomial degree"); }
    parameter_handler.leave_subsection();

    VFPEquationSolver<flags, 3> vfp_equation_solver(
        parameter_handler, polynomial_degree, expansion_order);
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
