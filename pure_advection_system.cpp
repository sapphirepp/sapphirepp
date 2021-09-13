#include <deal.II/base/function.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/tensor_function.h>
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

#include <array>
#include <fstream>
#include <iostream>

namespace pure_advection_system {
using namespace dealii;
// Input data
// the velocity field field
template <int dim>
class VelocityField : public TensorFunction<1, dim> {
 public:
  virtual void value_list(const std::vector<Point<dim>>& points,
                          std::vector<Tensor<1, dim>>& values) const override {
    Assert(dim == 1, ExcNotImplemented());
    Assert(values.size() == points.size(),
           ExcDimensionMismatch(values.size(), points.size()));

    Tensor<1, dim> value({.5});
    // constant velocity field
    for (unsigned int i = 0; i < points.size(); ++i) {
      values[i] = value;
    }
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
inline StreamType& operator<<(StreamType& s, TermFlags f) {
  s << "Term flags: \n";
  if (f & advection) s << "	 - Advection\n";
  if (f & reaction) s << "	 - Reaction\n";
  return s;
}

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

  // NOTE: The explicit use of a mapping is most likely related to the usage of
  // mesh_loop as well
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
  std::vector<Tensor<1, 1>> values(2);
  beta.value_list(points, values);
  for (auto value : values) {
    std::cout << value << "\n";
  }

  constexpr TermFlags flags = advection;
  PureAdvection<flags, 2, 1> pure_advection;
  pure_advection.run();

  return 0;
}
