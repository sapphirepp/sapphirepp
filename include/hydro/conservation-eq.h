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
 * \( \frac{\partial u}{\partial t} + \nabla \cdot \mathbf{f}(u) = 0 \)
 * where \( u \) is the solution and \( \mathbf{f}(u) \) is the flux function.
 * Here the flux is given by \( \mathbf{f}(u) = a u \) with \( a \) a constant.
 *
 */

#ifndef HYDROSOLVER_CONSERVATIONEQ_H
#define HYDROSOLVER_CONSERVATIONEQ_H

#include <deal.II/base/function.h>
#include <deal.II/base/mpi.h>
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

/**
 * @brief Exact analytical solution of the conservation equation.
 *
 *  \( u(x, t) = u_0(x - a \cdot t) \)
 */
template <int dim> class ExactSolution : public Function<dim> {
public:
  ExactSolution(const Tensor<1, dim> &beta, const double time = 0.0)
      : Function<dim>(1, time), beta(beta) {}
  void vector_value(const Point<dim> &p, Vector<double> &values) const override;

private:
  const Tensor<1, dim> beta;
};

/**
 * @brief Iniitial condition for the conservation equation.
 *
 * \( u_0(x) = sin(x) \)
 */
template <int dim> class InitialCondition : public Function<dim> {
public:
  InitialCondition(const Tensor<1, dim> &beta) : Function<dim>(1), beta(beta) {}
  void vector_value(const Point<dim> &p, Vector<double> &values) const override;

private:
  const Tensor<1, dim> beta;
};

template <int dim> class BoundaryValues : public Function<dim> {
public:
  BoundaryValues(const Tensor<1, dim> &beta, const double time = 0.0)
      : Function<dim>(1, time), beta(beta) {}
  void vector_value(const Point<dim> &p, Vector<double> &values) const override;

private:
  const Tensor<1, dim> beta;
};

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
  FullMatrix<double> cell_mass_matrix;
  FullMatrix<double> cell_dg_matrix;
  Vector<double> cell_rhs;
  std::vector<types::global_dof_index> local_dof_indices;
  std::vector<types::global_dof_index> local_dof_indices_neighbor;
  std::vector<CopyDataFace> face_data;

  template <typename Iterator>
  void reinit(const Iterator &cell, unsigned int dofs_per_cell) {
    cell_mass_matrix.reinit(dofs_per_cell, dofs_per_cell);
    cell_dg_matrix.reinit(dofs_per_cell, dofs_per_cell);
    cell_rhs.reinit(dofs_per_cell);

    local_dof_indices.resize(dofs_per_cell);
    cell->get_dof_indices(local_dof_indices);
  }
};

/**
 * @brief Solve the simple conservation equation.
 *
 * This class solves the conservation equation
 * \( \frac{\partial u}{\partial t} + \nabla \cdot \mathbf{f}(u) = 0 \)
 * where \( u \) is the solution and \( \mathbf{f}(u) \) is the flux function.
 * Here the flux is given by \( \mathbf{f}(u) = a u \) with \( a \) a constant.
 */
template <int dim> class ConservationEq {
public:
  ConservationEq(const Tensor<1, dim> &beta);
  void run();

private:
  void make_grid();
  void setup_system();
  void assemble_system();
  void assemble_system_old();
  void assemble_time_step();
  void solve();
  void output_results() const;

  const Tensor<1, dim> beta;

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

  double time;
  double time_step;
  unsigned int timestep_number;

  ConditionalOStream pcout;
  TimerOutput computing_timer;
};

} // namespace Hydro
} // namespace Sapphire
#endif
