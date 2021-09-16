#include <deal.II/base/function.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/tensor_function.h>
#include <deal.II/base/types.h>
#include <deal.II/base/vectorization.h>
#include <deal.II/fe/fe_update_flags.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/tria.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/precondition_block.h>
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
#include <deal.II/meshworker/mesh_loop.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/vector_tools.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iostream>

namespace pure_advection_system {
using namespace dealii;
// Input data
// the velocity field field
template <int dim>
class VelocityField : public TensorFunction<1, dim> {
 public:
  virtual void value_list(const std::vector<Point<dim>> &points,
                          std::vector<Tensor<1, dim>> &values) const override {
    Assert(dim == 1, ExcNotImplemented());
    Assert(values.size() == points.size(),
           ExcDimensionMismatch(values.size(), points.size()));

    Tensor<1, dim> value({.5});
    // constant velocity field
    for (unsigned int i = 0; i < points.size(); ++i) {
      values[i] = value;
    }
  }

  void divergence_list(const std::vector<Point<dim>> &points,
                       std::vector<double> &values) {
    Assert(dim == 1, ExcNotImplemented());
    Assert(values.size() == points.size(),
           ExcDimensionMismatch(values.size(), points.size()));
    std::fill(values.begin(), values.end(), 0.);
  }
};

enum TermFlags { advection = 1 << 0, reaction = 1 << 1 };

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
  return s;
}

// the mesh_loop function requires helper data types
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
              const UpdateFlags interface_update_flags =
                  update_values | update_quadrature_points | update_JxW_values |
                  update_normal_vectors)
      : fe_values(mapping, fe, quadrature, update_flags),
        fe_interface_values(mapping, fe, quadrature_face,
                            interface_update_flags) {}

  // Copy Constructor
  ScratchData(const ScratchData<dim> &scratch_data)
      : fe_values(scratch_data.fe_values.get_mapping(),
                  scratch_data.fe_values.get_fe(),
                  scratch_data.fe_values.get_quadrature(),
                  scratch_data.fe_values.get_update_flags()),
        fe_interface_values(
            scratch_data.fe_interface_values.get_mapping(),
            scratch_data.fe_interface_values.get_fe(),
            scratch_data.fe_interface_values.get_quadratutre(),
            scratch_data.fe_interface_values.get_update_flags()) {}

  FEValues<dim> fe_values;
  FEInterfaceValues<dim> fe_interface_values;
};

struct CopyDataFace {
  FullMatrix<double> cell_matrix;
  std::vector<types::global_dof_index> joint_dof_indices;
};

struct CopyData {
  FullMatrix<double> cell_dg_matrix;
  FullMatrix<double> cell_mass_matrix;
  Vector<double> cell_rhs;
  std::vector<types::global_dof_index> local_dof_indices;
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
class PureAdvection {
 public:
  PureAdvection();
  void run();

 private:
  void make_grid();
  void setup_system();
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
  const double eta;

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
PureAdvection<flags, max_degree, dim>::PureAdvection()
    : dof_handler(triangulation),
      mapping(),
      fe(FE_DGQ<dim>(1), (max_degree + 1) * (max_degree + 1)),
      quadrature(fe.tensor_degree() + 1),
      quadrature_face(fe.tensor_degree() + 1),
      time_step(3. / 64),
      time(0.),
      time_step_number(0),
      theta(0.5),
      eta(1.),
      num_refinements(4) {
  for (unsigned int j = 0.5 * (max_degree + 1) * (max_degree + 2), i = 0, l = 0;
       l <= max_degree; ++l) {
    // a_lm indices
    for (unsigned int m = 0; m <= l; ++m, ++i) {
      lms_indices[i][0] = l;
      lms_indices[i][1] = m;
      lms_indices[i][2] = 0;
      // b_lm indices
      if (m != 0) {
        lms_indices[j][0] = l;
        lms_indices[j][1] = m;
        lms_indices[j][2] = 1;
        ++j;
      }
    }
  }
}

template <TermFlags flags, int max_degree, int dim>
void PureAdvection<flags, max_degree, dim>::run() {
  std::cout << "Start the simulation: "
            << "\n\n";
  output_parameters();
  make_grid();
  setup_system();
  assemble_system();
  output_index_order();
}

template <TermFlags flags, int max_degree, int dim>
void PureAdvection<flags, max_degree, dim>::make_grid() {
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
void PureAdvection<flags, max_degree, dim>::setup_system() {
  dof_handler.distribute_dofs(fe);
  const unsigned int n_dofs = dof_handler.n_dofs();
  std::cout << "	Number of degrees of freedom: " << n_dofs << "\n";

  DynamicSparsityPattern dsp(n_dofs);
  DoFTools::make_flux_sparsity_pattern(dof_handler, dsp);
  sparsity_pattern.copy_from(dsp);

  dg_matrix.reinit(sparsity_pattern);
  // TODO: Check if you can work with two different sparsity patterns.
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
void PureAdvection<flags, max_degree, dim>::assemble_system() {
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

    for (unsigned int i : fe_v.dof_indices()) {
      const unsigned int component_i =
          fe_v.get_fe().system_to_component_index(i).first;
      const std::array<unsigned int, 3> i_lms = lms_indices[component_i];

      for (unsigned int j : fe_v.dof_indices()) {
        const unsigned int component_j =
            fe_v.get_fe().system_to_component_index(j).first;
        const std::array<unsigned int, 3> j_lms = lms_indices[component_j];

        for (const unsigned int q_index : fe_v.quadrature_point_indices()) {
          // mass matrix
          copy_data.cell_mass_matrix(i, j) +=
              (component_i == component_j
                   ? fe_v.shape_value(i, q_index) * fe_v.shape_value(j, q_index)
                   : 0) *
              JxW[q_index];
          // dg matrix
          if (flags & TermFlags::reaction) {
            std::cout << "reaction"
                      << "\n";
            // l(l+1) * \phi_i * \phi_j
            if (component_i == component_j) {
              copy_data.cell_dg_matrix(i, j) +=
                  (i_lms[0] * (i_lms[0] + 1) * fe_v.shape_value(i, q_index) *
                   fe_v.shape_value(j, q_index) * JxW[q_index]);
            }
          }
          if (flags & TermFlags::advection) {
            std::cout << "advection"
                      << "\n";
            // - \div \vec{a}_ij \phi_i \phi_j where \div \vec{a}_ij =
            // 0 for i
            // != j and \div \vec{a}_ii = \div \vec{u}
            if (component_i == component_j) {
              copy_data.cell_dg_matrix(i, j) +=
                  -div_velocities[q_index] * fe_v.shape_value(i, q_index) *
                  fe_v.shape_value(j, q_index) * JxW[q_index];
            }
            // -\grad \phi_i * \vec{a}_i_j * phi_j
            if (component_i == component_j) {
              std::cout << "component_i: " << component_i << "\n";
              std::cout << "component_j: " << component_j << "\n";
              std::cout << "l, m, s: " << i_lms[0] << i_lms[1] << i_lms[2]
                        << "\n";
              std::cout << "l_prime, m_prime, s_prime: " << j_lms[0] << j_lms[1]
                        << j_lms[2] << "\n";
              std::cout << "a: " << velocities[q_index] << "\n";

              copy_data.cell_dg_matrix(i, j) +=
                  -fe_v.shape_grad(i, q_index) * velocities[q_index] *
                  fe_v.shape_value(j, q_index) * JxW[q_index];
            }
            if (i_lms[0] - 1 == j_lms[0] && i_lms[1] == j_lms[1] &&
                i_lms[2] == j_lms[2]) {
              std::cout << "component_i: " << component_i << "\n";
              std::cout << "component_j: " << component_j << "\n";
              std::cout << "l, m, s: " << i_lms[0] << i_lms[1] << i_lms[2]
                        << "\n";
              std::cout << "l_prime, m_prime, s_prime: " << j_lms[0] << j_lms[1]
                        << j_lms[2] << "\n";

              Tensor<1, dim> a;
              a[0] = (i_lms[0] - i_lms[1]) / (2. * i_lms[0] - 1.);
              std::cout << "a: " << a << "\n";

              copy_data.cell_dg_matrix(i, j) +=
                  -fe_v.shape_grad(i, q_index) * a *
                  fe_v.shape_value(j, q_index) * JxW[q_index];
            }
            if (i_lms[0] + 1 == j_lms[0] && i_lms[1] == j_lms[1] &&
                i_lms[2] == j_lms[2]) {
              std::cout << "component_i: " << component_i << "\n";
              std::cout << "component_j: " << component_j << "\n";
              std::cout << "l, m, s: " << i_lms[0] << i_lms[1] << i_lms[2]
                        << "\n";
              std::cout << "l_prime, m_prime, s_prime: " << j_lms[0] << j_lms[1]
                        << j_lms[2] << "\n";

              Tensor<1, dim> a;
              a[0] = (i_lms[0] + i_lms[1] + 1.) / (2. * i_lms[0] + 3.);
              std::cout << "a: " << a << "\n";
              copy_data.cell_dg_matrix(i, j) +=
                  -fe_v.shape_grad(i, q_index) * a *
                  fe_v.shape_value(j, q_index) * JxW[q_index];
            }
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
    scratch_data.fe_interface_values.reinit(cell, face_no);
    // Return a reference to the FEFaceValues or FESubfaceValues object
    // of the specified cell of the interface
    const FEFaceValuesBase<dim> &fe_face_v =
        scratch_data.fe_interface_values.get_fe_face_values(0);
    // The next line is unclear to me; Are the number of DOFs on an
    // interface the same as the number of DOFs on a cell when using DG?
    const unsigned int n_facet_dofs = fe_face_v.get_fe().n_dofs_per_cell();

    const std::vector<Point<dim>> &q_points = fe_face_v.get_quadrature_points();
    const std::vector<double> &JxW = fe_face_v.get_JxW_values();
    const std::vector<Tensor<1, dim>> &normals = fe_face_v.get_normal_vectors();

    std::vector<Tensor<1, dim>> velocities(q_points.size());
    beta.value_list(q_points, velocities);

    for (unsigned int i = 0; i < n_facet_dofs; ++i) {
      const unsigned int component_i =
          fe_face_v.get_fe().system_to_component_index(i).first;
      const std::array<unsigned int, 3> i_lms = lms_indices[i];

      for (unsigned int j = 0; j < n_facet_dofs; ++j) {
        const unsigned int component_j =
            fe_face_v.get_fe().system_to_component_index(j).first;
        const std::array<unsigned int, 3> j_lms = lms_indices[j];

        for (unsigned int q_index : fe_face_v.quadrature_point_indices()) {
          // for outflow boundary \vec{a}_ij * n \phi_i \phi_j
          if (component_i == component_j) {
            if (velocities[q_index] * normals[q_index] > 0) {
              copy_data.cell_dg_matrix(i, j) +=
                  fe_face_v.shape_value(i, q_index) * velocities[q_index] *
                  normals[q_index] * fe_face_v.shape_value(j, q_index) *
                  JxW[q_index];
            }
          }
          if (i_lms[0] - 1 == j_lms[0] && i_lms[1] == j_lms[1] &&
              i_lms[2] == j_lms[2]) {
            Tensor<1, dim> a;
            a[0] = (i_lms[0] - i_lms[1]) / (2. * i_lms[0] - 1.);

            if (a * normals[q_index] > 0) {
              copy_data.cell_dg_matrix(i, j) +=
                  fe_face_v.shape_value(i, q_index) * a * normals[q_index] *
                  fe_face_v.shape_value(j, q_index) * JxW[q_index];
            }
          }
          if (i_lms[0] + 1 == j_lms[0] && i_lms[1] == j_lms[1] &&
              i_lms[2] == j_lms[2]) {
            Tensor<1, dim> a;
            a[0] = (i_lms[0] + i_lms[1] + 1.) / (2. * i_lms[0] + 3.);
            if (a * normals[q_index] > 0) {
              copy_data.cell_dg_matrix(i, j) +=
                  fe_face_v.shape_value(i, q_index) * a * normals[q_index] *
                  fe_face_v.shape_value(j, q_index) * JxW[q_index];
            }
          }
        }
      }
    }
  };

  // assemble interior face terms
  const auto face_worker = [&](const Iterator &cell,
                               const unsigned int &face_no,
                               const unsigned int &subface_no,
                               const Iterator &neighbor_cell,
                               const unsigned int &neighbor_face_no,
                               const unsigned int &neighbor_subface_no,
                               ScratchData<dim> &scratch_data,
                               CopyData &copy_data) {
    FEInterfaceValues<dim> &fe_interface_v = scratch_data.fe_interface_values;
    fe_interface_v.reinit(cell, face_no, subface_no, neighbor_cell,
                          neighbor_face_no, neighbor_subface_no);

    const std::vector<Point<dim>> &q_points =
        fe_interface_v.get_quadrature_points();
    // Create an element at the end of the vector containig the face data
    copy_data.face_data.emplace_back();
    CopyDataFace &copy_data_face = copy_data.face_data.back();
    const unsigned int n_interface_dofs =
        fe_interface_v.n_current_interface_dofs();
    copy_data_face.cell_matrix.reinit(n_interface_dofs, n_interface_dofs);
    copy_data_face.joint_dof_indices =
        fe_interface_v.get_interface_dof_indices();

    const std::vector<double> &JxW = fe_interface_v.get_JxW_values();
    const std::vector<Tensor<1, dim>> &normals =
        fe_interface_v.get_normal_vectors();

    std::vector<Tensor<1, dim>> velocities(q_points.size());
    beta.value_list(q_points, velocities);

    for (unsigned int i = 0; i < n_interface_dofs; ++i) {
      unsigned int component_i =
          fe_interface_v.get_fe().system_to_component_index(i).first;
      std::array<unsigned int, 3> i_lms = lms_indices[component_i];
      for (unsigned int j = 0; j < n_interface_dofs; ++j) {
        unsigned int component_j =
            fe_interface_v.get_fe().system_to_component_index(j).first;
        std::array<unsigned int, 3> j_lms = lms_indices[component_j];
        for (unsigned int q_index = 0; q_index < q_points.size(); ++q_index) {
          if (component_i == component_j) {
            // centered flux \vec{a}_ij * n_F *
            // average_of_shape_values(\phi) * jump_in_shape_values(\phi_j)
            copy_data.cell_dg_matrix(i, j) +=
                velocities[q_index] * normals[q_index] *
                fe_interface_v.average(i, q_index) *
                fe_interface_v.jump(j, q_index) * JxW[q_index];
            // updwinding eta/2 * abs(\vec{a}_ij * n_F) *
            // jump_in_shape_values(\phi_i) * jump_in_shape_values(phi_j)
            copy_data.cell_dg_matrix(i, j) +=
                eta / 2 * std::abs(velocities[q_index] * normals[q_index]) *
                fe_interface_v.jump(i, q_index) *
                fe_interface_v.jump(j, q_index) * JxW[q_index];
          }
          if (i_lms[0] - 1 == j_lms[0] && i_lms[1] == j_lms[1] &&
              i_lms[2] == j_lms[2]) {
            Tensor<1, dim> a;
            a[0] = (i_lms[0] - i_lms[1]) / (2. * i_lms[0] - 1.);
            // centered fluxes ( no updwinding in off diagonal elements,
            // since it should increase coercivity of the bilinear form)
            copy_data.cell_dg_matrix(i, j) +=
                a * normals[q_index] * fe_interface_v.average(i, q_index) *
                fe_interface_v.jump(j, q_index) * JxW[q_index];
          }
          if (i_lms[0] + 1 == j_lms[0] && i_lms[1] == j_lms[1] &&
              i_lms[2] == j_lms[2]) {
            Tensor<1, dim> a;
            a[0] = (i_lms[0] + i_lms[1] + 1.) / (2. * i_lms[0] + 3.);
            copy_data.cell_dg_matrix(i, j) +=
                a * normals[q_index] * fe_interface_v.average(i, q_index) *
                fe_interface_v.jump(j, q_index) * JxW[q_index];
          }
        }
      }
    }
  };
  
  for (const auto &cell : dof_handler.active_cell_iterators()) {
    CopyData copy_data;
    ScratchData<dim> scratch_data{mapping, fe, quadrature, quadrature_face};
    cell_worker(cell, scratch_data, copy_data);
  }
}

template <TermFlags flags, int max_degree, int dim>
void PureAdvection<flags, max_degree, dim>::output_parameters() const {
  std::cout << "	Time step: " << time_step << "\n";
  std::cout << "	Theta: " << theta << "\n";
  std::cout << "	Eta: " << eta << "\n";
  std::cout << "	" << flags;
  std::cout << "	Number of modes: " << num_modes << "\n";
  std::cout << "	Number of global refinements: " << num_refinements
            << "\n";
}

template <TermFlags flags, int max_degree, int dim>
void PureAdvection<flags, max_degree, dim>::output_index_order() const {
  std::cout << "Ordering of the lms indices: "
            << "\n";

  for (std::array<unsigned int, 3> lms : lms_indices)
    std::cout << "	" << lms[0] << lms[1] << lms[2] << "\n";
}
}  // namespace pure_advection_system

int main() {
  using namespace pure_advection_system;
  // Test 1D velocity field
  VelocityField<1> beta;
  Point<1> point_01{0.3};
  Point<1> point_02{0.1};
  std::vector<Point<1>> points{point_01, point_02};
  // std::vector<Tensor<1, 1>> values(2);
  std::vector<double> values(2);
  beta.divergence_list(points, values);
  for (auto value : values) {
    std::cout << value << "\n";
  }

  constexpr TermFlags flags = TermFlags::advection;
  PureAdvection<flags, 2, 1> pure_advection;
  pure_advection.run();

  return 0;
}
