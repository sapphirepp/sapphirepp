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
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/mapping_q1.h>
#include <deal.II/grid/tria.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/vector.h>

namespace Sapphire {
namespace Hydro {
using namespace dealii;

/**
 * @brief Exact analytical solution of the conservation equation.
 *
 *  \( u(x, t) = u_0(x - a \cdot t) \)
 */
class ExactSolution : public Function<1> {
public:
  ExactSolution(double time = 0.0) : Function<1>(1, time) {}
  void vector_value(const Point<1> &p, Vector<double> &values) const override;

private:
  const double a = 1.0;
};

/**
 * @brief Iniitial condition for the conservation equation.
 *
 * \( u_0(x) = sin(x) \)
 */
class InitialCondition : public Function<1> {
public:
  InitialCondition() = default;
  void vector_value(const Point<1> &p, Vector<double> &values) const override;
};

/**
 * @brief Solve the simple conservation equation.
 *
 * This class solves the conservation equation
 * \( \frac{\partial u}{\partial t} + \nabla \cdot \mathbf{f}(u) = 0 \)
 * where \( u \) is the solution and \( \mathbf{f}(u) \) is the flux function.
 * Here the flux is given by \( \mathbf{f}(u) = a u \) with \( a \) a constant.
 */
class ConservationEq {
public:
  ConservationEq();
  void run();

private:
  void make_grid();
  void setup_system();
  void assemble_system();
  void solve();
  void output_results() const;

  const static unsigned int dim = 1; // TODO: make this a template parameter

  Triangulation<dim> triangulation;
  const MappingQ1<dim> mapping;

  FE_DGQ<dim> fe;
  DoFHandler<dim> dof_handler;
  const QGauss<dim> quadrature_formula;
  const QGauss<dim - 1> face_quadrature_formula;

  SparsityPattern sparcity_pattern;
  SparseMatrix<double> system_matrix;

  Vector<double> solution;
  Vector<double> old_solution;
  Vector<double> system_rhs;

  double time;
  double time_step;
  unsigned int timestep_number;
};

} // namespace Hydro
} // namespace Sapphire
#endif
