/**
 * @file slope-limiter.h
 * @author Florian Schulze (florian.schulze@mpi-hd.mpg.de)
 * @brief
 * @version 0.1
 * @date 2023-07-28
 */


#ifndef HYDROSOLVER_SLOPELIMITER_H
#define HYDROSOLVER_SLOPELIMITER_H

#include <deal.II/lac/vector.h>

#include "parameter-parser.h"

namespace Sapphire
{
  namespace Hydro
  {
    using namespace dealii;

    template <int dim>
    class SlopeLimiter
    {
    public:
      SlopeLimiter(const Sapphire::Utils::ParameterParser &prm);

      double
      minmod(const std::vector<double> &values);

      void
      minmod(const std::vector<Tensor<1, dim>> &values,
             const unsigned int                 n,
             Tensor<1, dim>                    &return_value);

      void
      compute_limited_slope(
        const double                      &cell_average,
        const Tensor<1, dim>              &cell_grad,
        const std::vector<double>         &neighbor_cell_averages,
        const std::vector<Tensor<1, dim>> &neighbor_distance,
        const unsigned int                 n_neighbors,
        Tensor<1, dim>                    &limited_slope);

      const SlopeLimiterType limiter_type;
      const SlopeLimiterCriterion limiter_criterion;
    };

  } // namespace Hydro
} // namespace Sapphire

#endif