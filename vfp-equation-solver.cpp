#include <deal.II/base/exceptions.h>
#include <deal.II/base/function.h>
#include <deal.II/base/quadrature.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/tensor_function.h>
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
#include <string>
#include <vector>

namespace vfp_equation_solver {
using namespace dealii;
// Input data
// the velocity field field
template <int dim>
class VelocityField : public TensorFunction<1, dim> {
 public:
  virtual Tensor<1, dim> value(const Point<dim> &point) const override {
    Assert(dim == 2, ExcNotImplemented());
    (void)point;  // constant velocity field (suppresses the
                  // compiler warning unused variable)
    // constant velocity
    Tensor<1, dim> velocity({u_x, u_y});
    return velocity;
  }

  void divergence_list(const std::vector<Point<dim>> &points,
                       std::vector<double> &values) {
    Assert(dim == 2, ExcNotImplemented());
    Assert(values.size() == points.size(),
           ExcDimensionMismatch(values.size(), points.size()));
    std::fill(values.begin(), values.end(), 0.);
  }

 private:
  double u_x = 0.;
  double u_y = 0.;
};

// the magnetic field
template <int dim>
class MagneticField : public Function<dim> {
 public:
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
  double B_x = 0.;
  double B_y = 0.;
  double B_z = 1.;
};

enum TermFlags { advection = 1 << 0, reaction = 1 << 1, magnetic = 1 << 2 };

constexpr TermFlags operator|(TermFlags f1, TermFlags f2) {
  return static_cast<TermFlags>(static_cast<int>(f1) | static_cast<int>(f2));
}

constexpr TermFlags operator&(TermFlags f1, TermFlags f2) {
  return static_cast<TermFlags>(static_cast<int>(f1) & static_cast<int>(f2));
}

template <typename StreamType>
inline StreamType &operator<<(StreamType &s, TermFlags f) {
  s << "Term flags: \n";
  if (f & advection) s << "	 - Advection\n";
  if (f & reaction) s << "	 - Reaction\n";
  if (f & magnetic) s << "	 - Magnetic\n";
  return s;
}

// Initial values
template <int max_degree, int dim>
class InitialValueFunction : public Function<dim> {
 public:
  virtual void vector_value(const Point<dim> &p,
                            Vector<double> &values) const override {
    Assert(dim == 2, ExcNotImplemented());
    Assert(values.size() == (max_degree + 1) * (max_degree + 1),
           ExcDimensionMismatch(values.size(),
                                (max_degree + 1) * (max_degree + 1)));
    // The zeroth component of values corresponds to f_000, the first component
    // to f_110 etc.
    values[0] =
        1. *
        std::exp(-((std::pow(p[0] - 0.5, 2) + std::pow(p[1] - 0.5, 2)) / 0.01));

    // Fill all components with the same values
    // std::fill(
    //     values.begin(), values.end(),
    //     1. * std::exp(-((std::pow(p[0] - 0.5, 2) + std::pow(p[1] - 0.5, 2)) /
    //                     0.01)));
  }
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

template <TermFlags flags, int max_degree, int dim>
class VFPEquationSolver {
 public:
  VFPEquationSolver();
  void run();

 private:
  void make_grid();
  void setup_pde_system();
  // Enum to determine if positive or negative upwind flux should be computed
  enum class FluxDirection { positive, negative };
  enum class Coordinate { x, y };
  void prepare_upwind_fluxes(const Point<dim> &p, const Coordinate coordinate,
                             const FluxDirection flux_direction);
  void setup_system();
  void project_initial_condition();
  void assemble_system();
  void solve_system();
  void output_results() const;
  void output_parameters() const;
  void output_index_order() const;

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
  // Spatial advection
  LAPACKFullMatrix<double> Ax;
  LAPACKFullMatrix<double> Ay;
  // Magnetic field
  LAPACKFullMatrix<double> Omega_x;
  LAPACKFullMatrix<double> Omega_y;
  LAPACKFullMatrix<double> Omega_z;
  // Reaction term
  Vector<double> R;

  // Upwind flux matrices
  FullMatrix<double> pi_x_positive;
  FullMatrix<double> pi_x_negative;
  FullMatrix<double> pi_y_positive;
  FullMatrix<double> pi_y_negative;

  SparsityPattern sparsity_pattern;
  SparseMatrix<double> mass_matrix;
  SparseMatrix<double> dg_matrix;
  SparseMatrix<double> system_matrix;

  Vector<double> current_solution;
  Vector<double> previous_solution;

  Vector<double> system_rhs;

  // number of expansion coefficients
  const unsigned int num_modes = (max_degree + 1) * (max_degree + 1);

  const double time_step;
  double time;
  unsigned int time_step_number;
  const double theta;
  // penalty parameter ( for a scalar equation eta = 1 -> upwinding)
  /* const double eta; */

  // scattering frequency
  const double scattering_frequency;
  // Number of refinements
  const unsigned int num_refinements;
  // Fort the moment, I will use a quadratic domain
  // Rectangular domain
  // Point<dim> left_bottom;
  // Point<dim> right_top;

  // Map between i and l,m,s (implemented in constructor)
  std::array<std::array<unsigned int, 3>, (max_degree + 1) * (max_degree + 1)>
      lms_indices;
};

template <TermFlags flags, int max_degree, int dim>
VFPEquationSolver<flags, max_degree, dim>::VFPEquationSolver()
    : dof_handler(triangulation),
      mapping(),
      fe(FE_DGQ<dim>(1), (max_degree + 1) * (max_degree + 1)),
      quadrature(fe.tensor_degree() + 1),
      quadrature_face(fe.tensor_degree() + 1),
      time_step(1. / 128),
      time(0.),
      time_step_number(0),
      theta(.5),
      /* eta(1.), */
      scattering_frequency(1.),
      num_refinements(5) {
  for (int s = 0, idx = 0; s <= 1; ++s) {
    for (int l = 0; l <= max_degree; ++l) {
      for (int m = l; m >= s; --m) {
        idx = l * (l + 1) - (s ? -1. : 1.) * m;
        lms_indices[idx][0] = l;
        lms_indices[idx][1] = m;
        lms_indices[idx][2] = s;
      }
    }
  }
}

template <TermFlags flags, int max_degree, int dim>
void VFPEquationSolver<flags, max_degree, dim>::run() {
  std::cout << "Start the simulation: "
            << "\n\n";
  output_parameters();
  output_index_order();
  make_grid();
  setup_pde_system();
  // Compute upwind fluxes
  prepare_upwind_fluxes({0.1, 0.1}, Coordinate::x, FluxDirection::positive);
  prepare_upwind_fluxes({0.1, 0.1}, Coordinate::x, FluxDirection::negative);
  prepare_upwind_fluxes({0.1, 0.1}, Coordinate::y, FluxDirection::positive);
  prepare_upwind_fluxes({0.1, 0.1}, Coordinate::y, FluxDirection::negative);
  setup_system();
  project_initial_condition();
  current_solution = previous_solution;
  output_results();

  time += time_step;
  ++time_step_number;

  assemble_system();

  Vector<double> tmp(current_solution.size());
  for (; time <= .8; time += time_step, ++time_step_number) {
    std::cout << "	Time step " << time_step_number << " at t = " << time
              << "\n";
    mass_matrix.vmult(system_rhs, previous_solution);
    dg_matrix.vmult(tmp, previous_solution);
    system_rhs.add(-time_step * (1 - theta), tmp);

    // Since the the dg_matrix depends on the velocity field and the the
    // velocity field may depend on time, it needs to reassembled every time
    // step. This is not true for the mass matrix ( but it may if the grid
    // adapts after a specified amount of time steps)

    // mass_matrix = 0; dg_matrix = 0;
    // assemble_system();

    // NOTE: Currently that is not necessary, because the velocity field is
    // constant. And I should check if I cannot update the system matrix without
    // reassembling in each time step.
    
    system_matrix.copy_from(mass_matrix);
    system_matrix.add(time_step * theta, dg_matrix);
    solve_system();

    output_results();
    previous_solution = current_solution;
  }
}

template <TermFlags flags, int max_degree, int dim>
void VFPEquationSolver<flags, max_degree, dim>::make_grid() {
  GridGenerator::hyper_cube(triangulation);
  triangulation.refine_global(num_refinements);

  // std::ofstream out("grid.vtk");
  // GridOut grid_out;
  // grid_out.write_vtk(triangulation, out);
  // std::cout << "	Grid written to grid.vtk"
  //           << "\n";
  std::cout << "	Number of active cells: "
            << triangulation.n_active_cells() << "\n";
}

template <TermFlags flags, int max_degree, int dim>
void VFPEquationSolver<flags, max_degree, dim>::setup_pde_system() {
  Ax.reinit(num_modes);
  Ay.reinit(num_modes);
  Omega_x.reinit(num_modes);
  Omega_y.reinit(num_modes);
  Omega_z.reinit(num_modes);
  R.reinit(num_modes);
  for (int s = 0; s <= 1; ++s) {
    for (int l = 0, i = 0; l <= max_degree; ++l) {
      for (int m = l; m >= s; --m) {
        i = l * (l + 1) - (s ? -1. : 1.) * m;  // (-1)^s
        for (int s_prime = 0; s_prime <= 1; ++s_prime) {
          for (int l_prime = 0, j = 0; l_prime <= max_degree; ++l_prime) {
            for (int m_prime = l_prime; m_prime >= s_prime; --m_prime) {
              j = l_prime * (l_prime + 1) - (s_prime ? -1. : 1.) * m_prime;
              // Ax
              if (l + 1 == l_prime && m == m_prime && s == s_prime)
                Ax.set(i, j,
                       std::sqrt(((l - m + 1.) * (l + m + 1.)) /
                                 ((2. * l + 3.) * (2 * l + 1.))));
              if (l - 1 == l_prime && m == m_prime && s == s_prime)
                Ax.set(i, j,
                       std::sqrt(((l - m) * (l + m)) /
                                 ((2. * l + 1.) * (2. * l - 1.))));
              // Ay
              if ((l + 1) == l_prime && (m + 1) == m_prime && s == s_prime)
                Ay.set(i, j,
                       -0.5 * std::sqrt(((l + m + 1.) * (l + m + 2.)) /
                                        ((2. * l + 3.) * (2. * l + 1.))));
              if ((l - 1) == l_prime && m + 1 == m_prime && s == s_prime)
                Ay.set(i, j,
                       0.5 * std::sqrt(((l - m - 1.) * (l - m)) /
                                       ((2. * l + 1.) * (2. * l - 1.))));
              if ((l + 1) == l_prime && (m - 1) == m_prime && s == s_prime)
                Ay.set(i, j,
                       0.5 * std::sqrt(((l - m + 1.) * (l - m + 2.)) /
                                       ((2. * l + 3.) * (2 * l + 1.))));
              if ((l - 1) == l_prime && (m - 1) == m_prime && s == s_prime)
                Ay.set(i, j,
                       -0.5 * std::sqrt(((l + m - 1.) * (l + m)) /
                                        ((2. * l + 1.) * (2. * l - 1.))));
              // Omega_x
              if (l == l_prime && m == m_prime && s == 0 && s_prime == 1) {
                Omega_x.set(i, j, -1. * m);
                Omega_x.set(j, i, 1. * m);  // Omega matrices are anti-symmetric
              }
              // Omega_y
              if (l == l_prime && (m + 1) == m_prime && s == 0 &&
                  s_prime == 1) {
                Omega_y.set(i, j, -0.5 * std::sqrt((l + m + 1.) * (l - m)));
                Omega_y.set(j, i, 0.5 * std::sqrt((l + m + 1.) * (l - m)));
              }
              if (l == l_prime && (m - 1) == m_prime && s == 0 &&
                  s_prime == 1) {
                Omega_y.set(i, j, -0.5 * std::sqrt((l - m + 1.) * (l + m)));
                Omega_y.set(j, i, 0.5 * std::sqrt((l - m + 1.) * (l + m)));
              }
              // Omega_z
              if (l == l_prime && (m + 1) == m_prime && s == s_prime) {
                Omega_z.set(i, j, 0.5 * std::sqrt((l + m + 1.) * (l - m)));
              }
              if (l == l_prime && (m - 1) == m_prime && s == s_prime) {
                Omega_z.set(i, j, -0.5 * std::sqrt((l - m + 1.) * (l + m)));
              }
              // R
              if (l == l_prime && m == m_prime && s == s_prime) {
                R[i] = 0.5 * scattering_frequency * l * (l + 1.);
              }
            }
          }
        }
      }
    }
  }
  // Edit entries of the matrices around m=0 and s=0. After the transformation
  // they fall out of the pattern of the matrix elements, namely they differ by
  // a factor of 2^(1/2).
  //
  // NOTE: These corrections were not included in the above loops, because some
  // of them were overwritten. This is an effect of the s loops being outside
  // the l and m loops.
  if (max_degree > 0) {
    for (unsigned int l = 0; l <= max_degree; ++l) {
      // Special cases for Omega
      if (l > 0) {
        // l == l_prime, m = 0, s = 0 and m_prime = 1 and s_prime = 1
        Omega_y(l * (l + 1), l * (l + 1) + 1) =
            std::sqrt(2) * Omega_y(l * (l + 1), l * (l + 1) + 1);
        // l == l_prime, m = 1, s = 1 and m_prime = 0 and s_prime = 0
        Omega_y(l * (l + 1) + 1, l * (l + 1)) =
            std::sqrt(2) * Omega_y(l * (l + 1) + 1, l * (l + 1));

        // l == l_prime, m = 0, s = 0 and m_prime = 1 and s_prime = 0
        Omega_z(l * (l + 1), l * (l + 1) - 1) =
            std::sqrt(2) * Omega_z(l * (l + 1), l * (l + 1) - 1);
        // // l == l_prime, m = 1, s = 0 and m_prime = 0 and s_prime = 0
        Omega_z(l * (l + 1) - 1, l * (l + 1)) =
            std::sqrt(2) * Omega_z(l * (l + 1) - 1, l * (l + 1));
      }
      // Special cases for A_y (necessary for every value of l)
      // Above the diagonal
      // l + 1 = l_prime, m = 0, s = 0, and m_prime = 1, s_prime = 0
      if (l != max_degree)
        Ay(l * (l + 1), (l + 1) * (l + 2) - 1) =
            std::sqrt(2) * Ay(l * (l + 1), (l + 2) * (l + 1) - 1);
      // l + 1 = l_prime, m = 1, s = 0, and m_prime = 0, s_prime = 0
      if (l != 0 && l != max_degree)
        Ay(l * (l + 1) - 1, (l + 1) * (l + 2)) =
            std::sqrt(2) * Ay(l * (l + 1) - 1, (l + 2) * (l + 1));
      // Below the diagonal
      // l - 1 = l_prime, m = 0, s = 0, and m_prime = 1, s_prime = 0
      if (l > 1)
        Ay(l * (l + 1), l * (l - 1) - 1) =
            std::sqrt(2) * Ay(l * (l + 1), l * (l - 1) - 1);
      // l - 1 = l_prime , m = 1, s = 0 and m_prime = 0, s_prime = 0
      if (l != 0)
        Ay(l * (l + 1) - 1, l * (l - 1)) =
            std::sqrt(2) * Ay(l * (l + 1) - 1, l * (l - 1));
    }
  }
  // std::cout << "Ax: " << "\n";
  // Ax.print_formatted(std::cout);
  // std::cout << "Ay: "
  //           << "\n";
  // Ay.print_formatted(std::cout);
  std::cout << "Omega_x: "
            << "\n";
  Omega_x.print_formatted(std::cout);
  std::cout << "Omega_y: "
            << "\n";
  Omega_y.print_formatted(std::cout);
  std::cout << "Omega_z: "
            << "\n";
  Omega_z.print_formatted(std::cout);

  // std::cout << "R: "
  //           << "\n";
  // R.print(std::cout);
}

template <TermFlags flags, int max_degree, int dim>
void VFPEquationSolver<flags, max_degree, dim>::prepare_upwind_fluxes(
    const Point<dim> &p, const Coordinate coordinate,
    const FluxDirection flux_direction) {
  Vector<double> eigenvalues(num_modes);
  FullMatrix<double> eigenvectors(num_modes);
  // The matrix gets destroyed, when the eigenvalues are computed. This requires
  // to copy it
  LAPACKFullMatrix<double> CopyA;
  switch (coordinate) {
    case Coordinate::x:
      CopyA = Ax;
      break;
    case Coordinate::y:
      CopyA = Ay;
      break;
  }

  double tolerance = 1.e-8;
  CopyA.compute_eigenvalues_symmetric(-2., 2., tolerance, eigenvalues,
                                      eigenvectors);
  // std::cout << "Eigenvectors: " << "\n"
  // eigenvectors.print_formatted(std::cout);
  // Remove eigenvalues, which are smaller than tolerance
  auto smaller_than_tolerance = [&tolerance](double value) {
    return (std::abs(value) < tolerance ? true : false);
  };
  std::replace_if(eigenvalues.begin(), eigenvalues.end(),
                  smaller_than_tolerance, 0.);
  // Add the velocities to the eigenvalues (u_k \delta_ij + a_k,ij)
  VelocityField<dim> velocity_field;
  double u = velocity_field.value(p)[static_cast<unsigned int>(coordinate)];
  eigenvalues.add(u);
  // std::cout << "eigenvalues plus velocity: "
  //           << "\n";
  // eigenvalues.print(std::cout);

  // Depending on the flux direction extract either the positive or negative
  // eigenvalues
  switch (flux_direction) {
      // Replace all values smaller than zero with 0
    case FluxDirection::positive:
      std::replace_if(eigenvalues.begin(), eigenvalues.end(),
                      std::bind(std::less<double>(), std::placeholders::_1, 0.),
                      0.);
      break;
      // Replace all values greater than zero with 0
    case FluxDirection::negative:
      std::replace_if(
          eigenvalues.begin(), eigenvalues.end(),
          std::bind(std::greater<double>(), std::placeholders::_1, 0.), 0.);
      break;
  }
  // std::cout << "positive (or negative) eigenvalues of the flux matrix pi_x/y:
  //           " << "\n";
  // eigenvalues.print(std::cout);

  // Create Lambda_+/-
  FullMatrix<double> lambda(num_modes, num_modes);
  // Fill diagonal
  for (unsigned int i = 0; i < num_modes; ++i) {
    lambda(i, i) = eigenvalues[i];
  }
  // std::cout << "positive (or negative) lambda matrix: " << "\n";
  // lambda.print_formatted(std::cout);

  // Compute the flux matrices
  switch (coordinate) {
    case Coordinate::x:
      switch (flux_direction) {
        case FluxDirection::positive:
          pi_x_positive.reinit(num_modes, num_modes);
          pi_x_positive.triple_product(lambda, eigenvectors, eigenvectors,
                                       false, true);
          // std::cout << "pi_x_positive: "
          //           << "\n";
          // pi_x_positive.print_formatted(std::cout);
          break;
        case FluxDirection::negative:
          pi_x_negative.reinit(num_modes, num_modes);
          pi_x_negative.triple_product(lambda, eigenvectors, eigenvectors,
                                       false, true);
          // std::cout << "pi_x_negative: "
          //           << "\n";
          // pi_x_negative.print_formatted(std::cout);
          break;
      }
      break;
    case Coordinate::y:
      switch (flux_direction) {
        case FluxDirection::positive:
          pi_y_positive.reinit(num_modes, num_modes);
          pi_y_positive.triple_product(lambda, eigenvectors, eigenvectors,
                                       false, true);
          std::cout << "pi_y_positive: "
                    << "\n";
          pi_y_positive.print_formatted(std::cout);
          break;
        case FluxDirection::negative:
          pi_y_negative.reinit(num_modes, num_modes);
          pi_y_negative.triple_product(lambda, eigenvectors, eigenvectors,
                                       false, true);
          std::cout << "pi_y_negative: "
                    << "\n";
          pi_y_negative.print_formatted(std::cout);
          break;
      }
      break;
  }
}

template <TermFlags flags, int max_degree, int dim>
void VFPEquationSolver<flags, max_degree, dim>::setup_system() {
  dof_handler.distribute_dofs(fe);
  const unsigned int n_dofs = dof_handler.n_dofs();
  std::cout << "	Number of degrees of freedom: " << n_dofs << "\n";

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

template <TermFlags flags, int max_degree, int dim>
void VFPEquationSolver<flags, max_degree, dim>::project_initial_condition() {
  FEValues<dim> fe_v(
      mapping, fe, quadrature,
      update_values | update_quadrature_points | update_JxW_values);

  const unsigned int n_dofs = fe.n_dofs_per_cell();

  FullMatrix<double> cell_mass_matrix(n_dofs, n_dofs);
  Vector<double> cell_rhs(n_dofs);

  std::vector<types::global_dof_index> local_dof_indices(n_dofs);

  for (const auto &cell : dof_handler.active_cell_iterators()) {
    cell_mass_matrix = 0;
    cell_rhs = 0;
    fe_v.reinit(cell);

    const std::vector<Point<dim>> &q_points = fe_v.get_quadrature_points();
    const std::vector<double> &JxW = fe_v.get_JxW_values();

    // Initial values
    InitialValueFunction<max_degree, dim> initial_value_function;
    // The Vector<double> inside the standard vector will not have the
    // correct size when the function vector_value is called by the function
    // vector_value_list(). Thus I have to create an empty vector of the
    // correct size and copy it into the standard vector.
    Vector<double> empty_vector_correct_size(num_modes);
    std::vector<Vector<double>> initial_values(q_points.size(),
                                               empty_vector_correct_size);
    initial_value_function.vector_value_list(q_points, initial_values);

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
        cell_rhs(i) += fe_v.shape_value(i, q_index) *
                       initial_values[q_index][component_i] * JxW[q_index];
      }
    }
    cell->get_dof_indices(local_dof_indices);

    constraints.distribute_local_to_global(
        cell_mass_matrix, cell_rhs, local_dof_indices, mass_matrix, system_rhs);
  }

  // Solve the system
  SolverControl solver_control(1000, 1e-12);
  SolverCG<Vector<double>> cg(solver_control);
  // PreconditionSSOR<SparseMatrix<double>> preconditioner;
  // preconditioner.initialize(mass_matrix, 1.2);
  cg.solve(mass_matrix, previous_solution, system_rhs, PreconditionIdentity());

  // DataOut<dim> data_out;
  // data_out.attach_dof_handler(dof_handler);
  // data_out.add_data_vector(previous_solution, "projection");
  // data_out.build_patches();
  // std::ofstream output("projection.vtu");
  // data_out.write_vtu(output);

  // Reset mass matrix and system RHS
  mass_matrix = 0;
  system_rhs = 0;
}

template <TermFlags flags, int max_degree, int dim>
void VFPEquationSolver<flags, max_degree, dim>::assemble_system() {
  /*
    What kind of loops are there ?
    1. Loop over all cells (this happens inside the mesh_loop)
    2. Loop over the degrees of freedom on each cell
    - the metod system_to_componet_index() returns the index of the non-zero
    component of the vector-valued shape function which corresponds to the
    indices (l,m,s)
  */
  using Iterator = typename DoFHandler<dim>::active_cell_iterator;
  VelocityField<dim> beta;
  beta.set_time(time);
  MagneticField<dim> magnetic_field;
  magnetic_field.set_time(time);

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

    std::vector<Tensor<1, dim>> velocities(q_points.size());
    beta.value_list(q_points, velocities);
    std::vector<double> div_velocities(q_points.size());
    beta.divergence_list(q_points, div_velocities);
    // Second argument constructs empty vectors with 3 values. They are copied
    // q_points.size() times
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
          // mass matrix
          copy_data.cell_mass_matrix(i, j) +=
              (component_i == component_j
                   ? fe_v.shape_value(i, q_index) * fe_v.shape_value(j, q_index)
                   : 0) *
              JxW[q_index];
          // dg matrix
          if (flags & TermFlags::reaction) {
            if (component_i == component_j) {
              // 0.5 * scattering_frequency * l(l+1) * \phi_i * \phi_j
              copy_data.cell_dg_matrix(i, j) +=
                  R[component_i] * fe_v.shape_value(i, q_index) *
                  fe_v.shape_value(j, q_index) * JxW[q_index];
              // - [\partial_x(u_x\delta_ij + Ax_ij) + \partial_y(u_y\delta_ij +
              // - Ay_ij) ] \phi_i \phi_j where \partial_x/y Ax/y_ij = 0
              copy_data.cell_dg_matrix(i, j) -=
                  div_velocities[q_index] * fe_v.shape_value(i, q_index) *
                  fe_v.shape_value(j, q_index) * JxW[q_index];
            }
          }
          if (flags & TermFlags::advection) {
            // -[partial_x \phi_i * (u_x \delta_ij + Ax_ij)
            //   + partial_y \phi_i * (u_y \delta_ij + Ay_ij)] * phi_j
            if (component_i == component_j) {
              copy_data.cell_dg_matrix(i, j) -=
                  fe_v.shape_grad(i, q_index)[0] * velocities[q_index][0] *
                  fe_v.shape_value(j, q_index) * JxW[q_index];
              copy_data.cell_dg_matrix(i, j) -=
                  fe_v.shape_grad(i, q_index)[1] * velocities[q_index][1] *
                  fe_v.shape_value(j, q_index) * JxW[q_index];

            } else {
              // NOTE: Many zerso are added here, because the matrices Ax, Ay
              // are sparse. TODO: Performance check. If too bad, return to the
              // strategy, which was used in v0.6.5
              // Ax
              copy_data.cell_dg_matrix(i, j) -= fe_v.shape_grad(i, q_index)[0] *
                                                Ax(component_i, component_j) *
                                                fe_v.shape_value(j, q_index) *
                                                JxW[q_index];
              // Ay
              copy_data.cell_dg_matrix(i, j) -= fe_v.shape_grad(i, q_index)[1] *
                                                Ay(component_i, component_j) *
                                                fe_v.shape_value(j, q_index) *
                                                JxW[q_index];
            }
          }
          if (flags & TermFlags::magnetic) {
            // Omega_x
            // std::cout << "B_x: " << magnetic_field_values[q_index][0] <<
            // "\n";
            copy_data.cell_dg_matrix(i, j) +=
                fe_v.shape_value(i, q_index) *
                magnetic_field_values[q_index][0] *
                Omega_x(component_i, component_j) *
                fe_v.shape_value(j, q_index) * JxW[q_index];
            // Omega_y
            // std::cout << "B_y: " << magnetic_field_values[q_index][1] <<
            // "\n";
            copy_data.cell_dg_matrix(i, j) +=
                fe_v.shape_value(i, q_index) *
                magnetic_field_values[q_index][1] *
                Omega_y(component_i, component_j) *
                fe_v.shape_value(j, q_index) * JxW[q_index];
            // Omega_z
            // std::cout << "B_z: " << magnetic_field_values[q_index][2] <<
            // "\n";
            copy_data.cell_dg_matrix(i, j) +=
                fe_v.shape_value(i, q_index) *
                magnetic_field_values[q_index][2] *
                Omega_z(component_i, component_j) *
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
    const std::vector<Point<dim>> &q_points = fe_face_v.get_quadrature_points();
    const std::vector<double> &JxW = fe_face_v.get_JxW_values();
    const std::vector<Tensor<1, dim>> &normals = fe_face_v.get_normal_vectors();

    std::vector<Tensor<1, dim>> velocities(q_points.size());
    beta.value_list(q_points, velocities);

    for (unsigned int i = 0; i < n_facet_dofs; ++i) {
      const unsigned int component_i =
          fe_face_v.get_fe().system_to_component_index(i).first;
      for (unsigned int j = 0; j < n_facet_dofs; ++j) {
        const unsigned int component_j =
            fe_face_v.get_fe().system_to_component_index(j).first;
        for (unsigned int q_index : fe_face_v.quadrature_point_indices()) {
          // for outflow boundary: if n_k > 0 , then n_k * \phi_i *
          // pi_k_positive_ij \phi_j, else n_k * \phi_i * pi_k_negative_ij
          // \phi_j
          // Ax
          if (normals[q_index][0] == 1.) {
            copy_data.cell_dg_matrix(i, j) +=
                normals[q_index][0] * fe_face_v.shape_value(i, q_index) *
                pi_x_positive(component_i, component_j) *
                fe_face_v.shape_value(j, q_index) * JxW[q_index];
          }
          if (normals[q_index][0] == -1.) {
            copy_data.cell_dg_matrix(i, j) +=
                normals[q_index][0] * fe_face_v.shape_value(i, q_index) *
                pi_x_negative(component_i, component_j) *
                fe_face_v.shape_value(j, q_index) * JxW[q_index];
          }
          // Ay
          if (normals[q_index][1] == 1.) {
            copy_data.cell_dg_matrix(i, j) +=
                normals[q_index][1] * fe_face_v.shape_value(i, q_index) *
                pi_y_positive(component_i, component_j) *
                fe_face_v.shape_value(j, q_index) * JxW[q_index];
          }
          if (normals[q_index][1] == -1.) {
            copy_data.cell_dg_matrix(i, j) +=
                normals[q_index][1] * fe_face_v.shape_value(i, q_index) *
                pi_y_negative(component_i, component_j) *
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
  // -     -->
  // |     |
  // *-----*
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

    const std::vector<double> &JxW = fe_v_face.get_JxW_values();
    const std::vector<Tensor<1, dim>> &normals = fe_v_face.get_normal_vectors();
    const std::vector<Point<dim>> &q_points = fe_v_face.get_quadrature_points();
    std::vector<Tensor<1, dim>> velocities(q_points.size());
    beta.value_list(q_points, velocities);

    for (unsigned int q_index : fe_v_face.quadrature_point_indices()) {
      // cell_dg_matrix_11
      for (unsigned int i : fe_v_face.dof_indices()) {
        const unsigned int component_i =
            fe_v_face.get_fe().system_to_component_index(i).first;
        for (unsigned int j : fe_v_face.dof_indices()) {
          unsigned int component_j =
              fe_v_face.get_fe().system_to_component_index(j).first;
          /* std::cout << normals[q_index] << std::endl; */
          // Ax
          if (normals[q_index][0] == 1.) {
            copy_data_face.cell_dg_matrix_11(i, j) +=
                normals[q_index][0] * fe_v_face.shape_value(i, q_index) *
                pi_x_positive(component_i, component_j) *
                fe_v_face.shape_value(j, q_index) * JxW[q_index];
          }
          /* if (normals[q_index][0] == -1.) { */
          /*   /\* std::cout << "test" *\/ */
          /*   /\*           << "\n"; *\/ */
          /*   copy_data_face.cell_dg_matrix_11(i, j) += */
          /*       normals[q_index][0] * fe_v_face.shape_value(i, q_index) * */
          /*       pi_x_negative(component_i, component_j) * */
          /*       fe_v_face.shape_value(j, q_index) * JxW[q_index]; */
          /* } */

          // Ay
          if (normals[q_index][1] == 1.) {
            copy_data_face.cell_dg_matrix_11(i, j) +=
                normals[q_index][1] * fe_v_face.shape_value(i, q_index) *
                pi_y_positive(component_i, component_j) *
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
          // Ax
          if (normals[q_index][0] == 1.) {
            copy_data_face.cell_dg_matrix_12(i, j) -=
                normals[q_index][0] *
                fe_v_face_neighbor.shape_value(i, q_index) *
                pi_x_positive(component_i, component_j) *
                fe_v_face.shape_value(j, q_index) * JxW[q_index];
          }
          // if (normals[q_index][0] == -1.) {
          //   copy_data_face.cell_dg_matrix_12(i, j) -=
          //       normals[q_index][0] *
          //       fe_v_face_neighbor.shape_value(i, q_index) *
          //       pi_x_negative(component_i, component_j) *
          //       fe_v_face.shape_value(j, q_index) * JxW[q_index];
          // }
          // Ay
          if (normals[q_index][1] == 1.) {
            copy_data_face.cell_dg_matrix_12(i, j) -=
                normals[q_index][1] *
                fe_v_face_neighbor.shape_value(i, q_index) *
                pi_y_positive(component_i, component_j) *
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
          if (normals[q_index][0] == 1.) {
            copy_data_face.cell_dg_matrix_21(i, j) +=
                normals[q_index][0] * fe_v_face.shape_value(i, q_index) *
                pi_x_negative(component_i, component_j) *
                fe_v_face_neighbor.shape_value(j, q_index) * JxW[q_index];
          }
          // if (normals[q_index][0] == -1.) {
          //   copy_data_face.cell_dg_matrix_21(i, j) +=
          //       normals[q_index][0] * fe_v_face.shape_value(i, q_index) *
          //       pi_x_positive(component_i, component_j) *
          //       fe_v_face_neighbor.shape_value(j, q_index) * JxW[q_index];
          // }
          if (normals[q_index][1] == 1.) {
            copy_data_face.cell_dg_matrix_21(i, j) +=
                normals[q_index][1] * fe_v_face.shape_value(i, q_index) *
                pi_y_negative(component_i, component_j) *
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
          if (normals[q_index][0] == 1.) {
            copy_data_face.cell_dg_matrix_22(i, j) -=
                normals[q_index][0] *
                fe_v_face_neighbor.shape_value(i, q_index) *
                pi_x_negative(component_i, component_j) *
                fe_v_face_neighbor.shape_value(j, q_index) * JxW[q_index];
          }
          // if (normals[q_index][0] == -1.) {
          //   copy_data_face.cell_dg_matrix_22(i, j) -=
          //       normals[q_index][0] *
          //       fe_v_face_neighbor.shape_value(i, q_index) *
          //       pi_x_positive(component_i, component_j) *
          //       fe_v_face_neighbor.shape_value(j, q_index) * JxW[q_index];
          // }
          if (normals[q_index][1] == 1.) {
            copy_data_face.cell_dg_matrix_22(i, j) -=
                normals[q_index][1] *
                fe_v_face_neighbor.shape_value(i, q_index) *
                pi_y_negative(component_i, component_j) *
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

template <TermFlags flags, int max_degree, int dim>
void VFPEquationSolver<flags, max_degree, dim>::solve_system() {
  SolverControl solver_control(1000, 1e-12);
  SolverRichardson<Vector<double>> solver(solver_control);

  PreconditionBlockSSOR<SparseMatrix<double>> preconditioner;
  preconditioner.initialize(system_matrix, fe.n_dofs_per_cell());
  solver.solve(system_matrix, current_solution, system_rhs, preconditioner);

  std::cout << "	Solver converged in " << solver_control.last_step()
            << " iterations."
            << "\n";
}

template <TermFlags flags, int max_degree, int dim>
void VFPEquationSolver<flags, max_degree, dim>::output_results() const {
  DataOut<dim> data_out;
  data_out.attach_dof_handler(dof_handler);
  // Create a vector of strings with names for the components of the solution
  std::vector<std::string> component_names((max_degree + 1) * (max_degree + 1));

  for (unsigned int i = 0; i < (max_degree + 1) * (max_degree + 1); ++i) {
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

template <TermFlags flags, int max_degree, int dim>
void VFPEquationSolver<flags, max_degree, dim>::output_parameters() const {
  std::cout << "	Time step: " << time_step << "\n";
  std::cout << "	Theta: " << theta << "\n";
  /* std::cout << "	Eta: " << eta << "\n"; */
  std::cout << "	" << flags;
  std::cout << "	Dimension: " << dim << "\n";
  std::cout << "	Scattering frequency: " << scattering_frequency << "\n";
  std::cout << "	Number of modes: " << num_modes << "\n";
  std::cout << "	Number of global refinements: " << num_refinements
            << "\n";
}

template <TermFlags flags, int max_degree, int dim>
void VFPEquationSolver<flags, max_degree, dim>::output_index_order() const {
  std::cout << "Ordering of the lms indices: "
            << "\n";

  for (std::array<unsigned int, 3> lms : lms_indices)
    std::cout << "	" << lms[0] << lms[1] << lms[2] << "\n";
}
}  // namespace vfp_equation_solver

int main() {
  try {
    using namespace vfp_equation_solver;

    constexpr TermFlags flags = TermFlags::advection | TermFlags::reaction;
    VFPEquationSolver<flags, 1, 2> pure_advection;
    pure_advection.run();

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
