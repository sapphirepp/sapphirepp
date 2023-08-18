/**
 * @file euler-operator.h
 * @author Florian Schulze (florian.schulze@mpi-hd.mpg.de)
 * @brief Implement EulerOperator class
 * @version 0.1
 * @date 2023-08-18
 *
 * @copyright Copyright (c) 2023
 *
 */

#ifndef HYDROSOLVER_EULEROPERATOR_H
#define HYDROSOLVER_EULEROPERATOR_H

#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/function.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/time_stepping.h>
#include <deal.II/base/timer.h>
#include <deal.II/base/utilities.h>
#include <deal.II/base/vectorization.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_system.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/la_parallel_vector.h>

#include <deal.II/matrix_free/fe_evaluation.h>
#include <deal.II/matrix_free/matrix_free.h>
#include <deal.II/matrix_free/operators.h>

#include <deal.II/numerics/data_out.h>

#include <fstream>
#include <iomanip>
#include <iostream>

#include "config.h"
#include "flux.h"
#include "output-module.h"
#include "parameter-parser.h"
#include "sapphire-logstream.h"
#include "time-stepping.h"


namespace Sapphire
{
  namespace Hydro
  {
    using namespace dealii;

    template <int dim, int degree, int n_points_1d>
    class EulerOperator
    {
    public:
      static constexpr unsigned int n_quadrature_points_1d = n_points_1d;

      EulerOperator(TimerOutput &timer_output);

      void
      reinit(const Mapping<dim> &mapping, const DoFHandler<dim> &dof_handler);

      void
      set_inflow_boundary(const types::boundary_id       boundary_id,
                          std::unique_ptr<Function<dim>> inflow_function);

      void
      set_subsonic_outflow_boundary(
        const types::boundary_id       boundary_id,
        std::unique_ptr<Function<dim>> outflow_energy);

      void
      set_wall_boundary(const types::boundary_id boundary_id);

      void
      set_body_force(std::unique_ptr<Function<dim>> body_force);

      void
      apply(const double                                      current_time,
            const LinearAlgebra::distributed::Vector<Number> &src,
            LinearAlgebra::distributed::Vector<Number>       &dst) const;

      void
      perform_stage(
        const Number                                      cur_time,
        const Number                                      factor_solution,
        const Number                                      factor_ai,
        const LinearAlgebra::distributed::Vector<Number> &current_ri,
        LinearAlgebra::distributed::Vector<Number>       &vec_ki,
        LinearAlgebra::distributed::Vector<Number>       &solution,
        LinearAlgebra::distributed::Vector<Number>       &next_ri) const;

      void
      project(const Function<dim>                        &function,
              LinearAlgebra::distributed::Vector<Number> &solution) const;

      std::array<double, 3>
      compute_errors(
        const Function<dim>                              &function,
        const LinearAlgebra::distributed::Vector<Number> &solution) const;

      double
      compute_cell_transport_speed(
        const LinearAlgebra::distributed::Vector<Number> &solution) const;

      void
      initialize_vector(
        LinearAlgebra::distributed::Vector<Number> &vector) const;

    private:
      MatrixFree<dim, Number> data;

      TimerOutput &timer;

      std::map<types::boundary_id, std::unique_ptr<Function<dim>>>
        inflow_boundaries;
      std::map<types::boundary_id, std::unique_ptr<Function<dim>>>
                                     subsonic_outflow_boundaries;
      std::set<types::boundary_id>   wall_boundaries;
      std::unique_ptr<Function<dim>> body_force;

      void
      local_apply_inverse_mass_matrix(
        const MatrixFree<dim, Number>                    &data,
        LinearAlgebra::distributed::Vector<Number>       &dst,
        const LinearAlgebra::distributed::Vector<Number> &src,
        const std::pair<unsigned int, unsigned int>      &cell_range) const;

      void
      local_apply_cell(
        const MatrixFree<dim, Number>                    &data,
        LinearAlgebra::distributed::Vector<Number>       &dst,
        const LinearAlgebra::distributed::Vector<Number> &src,
        const std::pair<unsigned int, unsigned int>      &cell_range) const;

      void
      local_apply_face(
        const MatrixFree<dim, Number>                    &data,
        LinearAlgebra::distributed::Vector<Number>       &dst,
        const LinearAlgebra::distributed::Vector<Number> &src,
        const std::pair<unsigned int, unsigned int>      &face_range) const;

      void
      local_apply_boundary_face(
        const MatrixFree<dim, Number>                    &data,
        LinearAlgebra::distributed::Vector<Number>       &dst,
        const LinearAlgebra::distributed::Vector<Number> &src,
        const std::pair<unsigned int, unsigned int>      &face_range) const;
    };
  } // namespace Hydro
} // namespace Sapphire

#endif