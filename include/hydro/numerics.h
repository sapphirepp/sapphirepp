/**
 * @file numerics.h
 * @author Florian Schulze (florian.schulze@mpi-hd.mpg.de)
 * @brief Implement numerical functions to solve the hydrodynamics equations
 * @version 0.1
 * @date 2023-07-12
 */

#ifndef HYDROSOLVER_NUMERICS_H
#define HYDROSOLVER_NUMERICS_H

#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_control.h>
#include <deal.II/lac/vector.h>

#include "parameter-parser.h"

namespace Sapphire
{
  namespace Hydro
  {
    using namespace dealii;

    class HDSolverControl
    {
    public:
      HDSolverControl(const Sapphire::Utils::ParameterParser &prm);

      const TimeSteppingScheme    scheme;
      const SlopeLimiter          limiter;
      const SlopeLimiterCriterion limiter_criterion;

      const unsigned int fe_degree;
      const double       time_step;
      const double       end_time;
      const unsigned int refinement_level;

      const unsigned int max_iterations;
      const double       tolerance;
    };

    double
    minmod(const std::vector<double> &values);

    template <int dim>
    void
    minmod(const std::vector<Tensor<1, dim>> &values,
           const unsigned int                 n,
           Tensor<1, dim>                    &return_value);

    template <int dim>
    void
    compute_limited_slope(const double              &cell_average,
                          const Tensor<1, dim>       cell_average_grad,
                          const std::vector<double> &neighbor_cell_averages,
                          const std::vector<Tensor<1, dim>> &neighbor_distance,
                          const unsigned int                 n_neighbors,
                          Tensor<1, dim>                    &limited_slope,
                          const HDSolverControl             &hd_solver_control);

  } // namespace Hydro
} // namespace Sapphire
#endif