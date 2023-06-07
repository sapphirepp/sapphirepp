/**
 * @file conservation-eq.h
 * @author Florian Schulze (florian.schulze@mpi-hd.mpg.de)
 * @brief Sovle the conservation equation.
 * @version 0.1
 * @date 2023-05-17
 *
 * @copyright Copyright (c) 2023
 *
 * We consider the conservation equation
 * \f$ \frac{\partial u}{\partial t} + \nabla \cdot \mathbf{f}(u) = 0 \f$
 * where \f$ u \f$ is the solution and \f$ \mathbf{f}(u) \f$ is the flux
 * function.
 *
 */

#ifndef HYDROSOLVER_CONSERVATIONEQ_H
#define HYDROSOLVER_CONSERVATIONEQ_H

#include <deal.II/base/function.h>
#include <deal.II/base/mpi.h>
#include <deal.II/base/smartpointer.h>
#include <deal.II/base/tensor_function.h>
#include <deal.II/base/timer.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_q1.h>
#include <deal.II/grid/tria.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/vector.h>

#include <mpi.h>

namespace Sapphire {
namespace Hydro {
using namespace dealii;

template <int dim> class ScratchData {
public:
  // Constructor
  ScratchData(const Mapping<dim> &mapping, const FiniteElement<dim> &fe,
              const Quadrature<dim> &quadrature,
              const Quadrature<dim - 1> &face_quadrature,
              const UpdateFlags update_flags = update_values |
                                               update_gradients |
                                               update_quadrature_points |
                                               update_JxW_values,
              const UpdateFlags face_update_flags = update_values |
                                                    update_quadrature_points |
                                                    update_JxW_values |
                                                    update_normal_vectors |
                                                    update_JxW_values,
              const UpdateFlags neighbor_face_update_flags = update_values)
      : fe_values(mapping, fe, quadrature, update_flags),
        fe_face_values(mapping, fe, face_quadrature, face_update_flags),
        fe_face_values_neighbor(mapping, fe, face_quadrature,
                                neighbor_face_update_flags) {}

  // Copy constructor
  ScratchData(const ScratchData &scratch_data)
      : fe_values(scratch_data.fe_values.get_mapping(),
                  scratch_data.fe_values.get_fe(),
                  scratch_data.fe_values.get_quadrature(),
                  scratch_data.fe_values.get_update_flags()),
        fe_face_values(scratch_data.fe_face_values.get_mapping(),
                       scratch_data.fe_face_values.get_fe(),
                       scratch_data.fe_face_values.get_quadrature(),
                       scratch_data.fe_face_values.get_update_flags()),
        fe_face_values_neighbor(
            scratch_data.fe_face_values_neighbor.get_mapping(),
            scratch_data.fe_face_values_neighbor.get_fe(),
            scratch_data.fe_face_values_neighbor.get_quadrature(),
            scratch_data.fe_face_values_neighbor.get_update_flags()) {}

  FEValues<dim> fe_values;
  FEFaceValues<dim> fe_face_values;
  FEFaceValues<dim> fe_face_values_neighbor;
};

struct CopyDataFace {
  FullMatrix<double> cell_dg_matrix_11;
  FullMatrix<double> cell_dg_matrix_12;
  FullMatrix<double> cell_dg_matrix_21;
  FullMatrix<double> cell_dg_matrix_22;

  Vector<double> cell_vector_1;
  Vector<double> cell_vector_2;

  std::vector<types::global_dof_index> local_dof_indices;
  std::vector<types::global_dof_index> local_dof_indices_neighbor;

  template <typename Iterator>
  void reinit(const Iterator &cell, const Iterator &neighbor_cell,
              unsigned int dofs_per_cell) {
    cell_dg_matrix_11.reinit(dofs_per_cell, dofs_per_cell);
    cell_dg_matrix_12.reinit(dofs_per_cell, dofs_per_cell);
    cell_dg_matrix_21.reinit(dofs_per_cell, dofs_per_cell);
    cell_dg_matrix_22.reinit(dofs_per_cell, dofs_per_cell);

    cell_vector_1.reinit(dofs_per_cell);
    cell_vector_2.reinit(dofs_per_cell);

    local_dof_indices.resize(dofs_per_cell);
    cell->get_dof_indices(local_dof_indices);

    local_dof_indices_neighbor.resize(dofs_per_cell);
    neighbor_cell->get_dof_indices(local_dof_indices_neighbor);
  }
};

struct CopyData {
  // TODO_BE: optimize memory layout
  FullMatrix<double> cell_mass_matrix;
  FullMatrix<double> cell_dg_matrix;
  Vector<double> cell_rhs;
  Vector<double> cell_vector;
  std::vector<types::global_dof_index> local_dof_indices;
  std::vector<types::global_dof_index> local_dof_indices_neighbor;
  std::vector<CopyDataFace> face_data;

  template <typename Iterator>
  void reinit(const Iterator &cell, unsigned int dofs_per_cell) {
    cell_mass_matrix.reinit(dofs_per_cell, dofs_per_cell);
    cell_dg_matrix.reinit(dofs_per_cell, dofs_per_cell);
    cell_rhs.reinit(dofs_per_cell);
    cell_vector.reinit(dofs_per_cell);

    local_dof_indices.resize(dofs_per_cell);
    cell->get_dof_indices(local_dof_indices);
  }
};

/**
 * @brief Solve the linear advection equation.
 *
 * This class solves the linear advection equation
 * \f$ \frac{\partial u}{\partial t} + \nabla \cdot (\mathbf{\beta}(\mathbf{x})
 * u) = 0 \f$
 * where \f$ u(\mathbf{x}, t) \f$ is the solution and the flux is \f$
 * \mathbf{f}(u) = \mathbf{\beta}(\mathbf{x}) u \f$. The initial condition is
 * given by \f$ u(\mathbf{x}, 0) = u_0(\mathbf{x}) \f$ and the boundary
 * condition is given by \f$ u(\mathbf{x}, t) = u_b(\mathbf{x}, t) \f$.
 *
 * \tparam dim dimension of the problem
 */
template <int dim> class ConservationEq {
public:
  /**
   * @brief Construct a new Conservation Eq object
   *
   *
   * @param beta wind vector field \f$ \mathbf{\beta}(\mathbf{x}) \f$
   * @param initial_condition initial condition \f$ u_0(\mathbf{x}) \f$
   * @param boundary_values boundary values \f$ u_b(\mathbf{x}, t) \f$
   * @param exact_solution exact solution for comparison \f$ u(\mathbf{x}, t)
   * \f$
   */
  ConservationEq(TensorFunction<1, dim, double> *beta,
                 Function<dim> *initial_condition,
                 Function<dim> *boundary_values, Function<dim> *exact_solution);
  void run();

private:
  void make_grid();
  void setup_system();
  void assemble_system();
  void assemble_system_old();
  void assemble_time_step();
  void solve();
  void output_results() const;
  void process_results();

  const SmartPointer<TensorFunction<1, dim, double>> beta;
  const SmartPointer<Function<dim>> initial_condition;
  const SmartPointer<Function<dim>> boundary_values;
  const SmartPointer<Function<dim>> exact_solution;

  MPI_Comm mpi_communicator;

  Triangulation<dim> triangulation;
  const MappingQ1<dim> mapping;

  FE_DGQ<dim> fe;
  DoFHandler<dim> dof_handler;
  const QGauss<dim> quadrature_formula;
  const QGauss<dim - 1> face_quadrature_formula;

  AffineConstraints<double> constraints;

  SparsityPattern sparcity_pattern;
  SparseMatrix<double> mass_matrix;
  SparseMatrix<double> dg_matrix;
  SparseMatrix<double> system_matrix;

  Vector<double> solution;
  Vector<double> old_solution;
  Vector<double> system_rhs;

  Vector<float> error_with_time;

  double time;
  double time_step;
  unsigned int timestep_number;

  ConditionalOStream pcout;
  TimerOutput computing_timer;
};

/**
 * @brief Solve Burgers' equation.
 *
 * This class solves the (1d) Burgers' equation
 * \f$ \frac{\partial u}{\partial t} + \frac{1}{2} \, \frac{\partial
 * u^2}{\partial x} = 0 \f$
 * where \f$ u(\mathbf{x}, t) \f$ is the solution and
 * the flux is \f$ \mathbf{f}(u) = \frac{1}{2} u^2 \f$. The initial
 * condition is given by \f$ u(\mathbf{x}, 0) = u_0(\mathbf{x}) \f$ and the
 * boundary condition is given by \f$ u(\mathbf{x}, t) = u_b(\mathbf{x}, t) \f$.
 *
 * \tparam dim dimension of the problem (must be 1)
 */
template <int dim> class BurgersEq {
public:
  /**
   * @brief Construct a new Burgers Eq object
   *
   *
   * @param initial_condition initial condition \f$ u_0(\mathbf{x}) \f$
   * @param boundary_values boundary values \f$ u_b(\mathbf{x}, t) \f$
   * @param exact_solution exact solution for comparison \f$ u(\mathbf{x}, t)
   * \f$
   */
  BurgersEq(Function<dim> *initial_condition, Function<dim> *boundary_values,
            Function<dim> *exact_solution);
  void run();

private:
  void make_grid();
  void setup_system();
  void assemble_mass_matrix();
  void assemble_dg_vector();
  void assemble_system();
  /**
   * @brief Solve the linear system
   *
   * This function solves the linear system \f$ M \mathbf{u} = \mathbf{b} \f$
   * where \f$ M \f$ is the mass matrix, \f$ \mathbf{u} \f$ is the solution
   * vector and \f$ \mathbf{b} \f$ is the right hand side vector.
   */
  void solve_linear_system();
  void perform_time_step();
  void output_results() const;
  void process_results();

  const SmartPointer<Function<dim>> initial_condition;
  const SmartPointer<Function<dim>> boundary_values;
  const SmartPointer<Function<dim>> exact_solution;

  MPI_Comm mpi_communicator;

  Triangulation<dim> triangulation;
  const MappingQ1<dim> mapping;

  FE_DGQ<dim> fe;
  DoFHandler<dim> dof_handler;
  const QGauss<dim> quadrature_formula;
  const QGauss<dim - 1> face_quadrature_formula;

  AffineConstraints<double> constraints;

  SparsityPattern sparcity_pattern;
  SparseMatrix<double> mass_matrix;
  SparseMatrix<double> dg_matrix;
  SparseMatrix<double> system_matrix;

  Vector<double> solution;

  // solution at the current intermediate time step, used to calculate the flux
  Vector<double> current_solution;

  // solution of the last time step
  Vector<double> old_solution;

  // vector of the fluxes at the current intermediate time step
  Vector<double> dg_vector;

  // right hand side of the linear system \f$ \mathbf{b} \f$
  Vector<double> system_rhs;

  Vector<float> error_with_time;

  double time;
  double time_step;
  unsigned int timestep_number;

  ConditionalOStream pcout;
  TimerOutput computing_timer;
};

} // namespace Hydro
} // namespace Sapphire
#endif
