/**
 * @file conservation-eq.h
 * @author Florian Schulze (florian.schulze@mpi-hd.mpg.de)
 * @brief Solve the conservation equation
 * @version 0.1
 * @date 2023-05-17
 *
 *
 * Implement the class ConservationEq that solves the conservation equation
 * \f$ \frac{\partial u}{\partial t} + \nabla \cdot (\mathbf{\beta}(\mathbf{x})
 * u) = 0 \f.
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

namespace Sapphire
{
  namespace Hydro
  {
    using namespace dealii;

    /**
     * @brief Solve the linear advection equation.
     *
     * This class solves the linear advection equation
     * \f$ \frac{\partial u}{\partial t} + \nabla \cdot
     * (\mathbf{\beta}(\mathbf{x}) u) = 0 \f$ where \f$ u(\mathbf{x}, t) \f$ is
     * the solution and the flux is \f$ \mathbf{f}(u) =
     * \mathbf{\beta}(\mathbf{x}) u \f$. The initial condition is given by \f$
     * u(\mathbf{x}, 0) = u_0(\mathbf{x}) \f$ and the boundary condition is
     * given by \f$ u(\mathbf{x}, t) = u_b(\mathbf{x}, t) \f$.
     *
     * \tparam dim dimension of the problem
     */
    template <int dim>
    class ConservationEq
    {
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
                     Function<dim>                  *initial_condition,
                     Function<dim>                  *boundary_values,
                     Function<dim>                  *exact_solution);
      void
      run();

    private:
      void
      make_grid();
      void
      setup_system();
      void
      assemble_system();
      void
      assemble_system_old();
      void
      assemble_time_step();
      void
      solve();
      void
      output_results() const;
      void
      process_results();

      const SmartPointer<TensorFunction<1, dim, double>> beta;
      const SmartPointer<Function<dim>>                  initial_condition;
      const SmartPointer<Function<dim>>                  boundary_values;
      const SmartPointer<Function<dim>>                  exact_solution;

      MPI_Comm mpi_communicator;

      Triangulation<dim>   triangulation;
      const MappingQ1<dim> mapping;

      FE_DGQ<dim>           fe;
      DoFHandler<dim>       dof_handler;
      const QGauss<dim>     quadrature_formula;
      const QGauss<dim - 1> face_quadrature_formula;

      AffineConstraints<double> constraints;

      SparsityPattern      sparcity_pattern;
      SparseMatrix<double> mass_matrix;
      SparseMatrix<double> dg_matrix;
      SparseMatrix<double> system_matrix;

      Vector<double> solution;
      Vector<double> old_solution;
      Vector<double> system_rhs;

      Vector<float> error_with_time;

      double       time;
      double       time_step;
      unsigned int timestep_number;

      ConditionalOStream pcout;
      TimerOutput        computing_timer;
    };

  } // namespace Hydro
} // namespace Sapphire
#endif
