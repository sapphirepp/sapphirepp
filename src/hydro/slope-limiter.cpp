#include "slope-limiter.h"

#include <iostream>

template <int dim>
Sapphire::Hydro::SlopeLimiter<dim>::SlopeLimiter(
  const Sapphire::Utils::ParameterParser &prm)
  : limiter_type(prm.hdsolver_limiter)
  , limiter_criterion(prm.hdsolver_limiter_criterion)
{}

template <int dim>
double
Sapphire::Hydro::SlopeLimiter<dim>::minmod(const std::vector<double> &values)
{
  auto [min_it, max_it] = std::minmax_element(values.begin(), values.end());
  if ((*min_it) * (*max_it) < 0.0)
    return 0.0;
  else if (std::abs(*min_it) < std::abs(*max_it))
    return *min_it;
  else
    return *max_it;
}

template <int dim>
void
Sapphire::Hydro::SlopeLimiter<dim>::minmod(
  const std::vector<Tensor<1, dim>> &values,
  const unsigned int                 n,
  Tensor<1, dim>                    &return_value)
{
  std::vector<double> component_values(n);
  for (unsigned int d = 0; d < dim; ++d)
    {
      for (unsigned int i = 0; i < n; ++i)
        component_values[i] = values[i][d];
      return_value[d] = minmod(component_values);
    }
}

template <int dim>
void
Sapphire::Hydro::SlopeLimiter<dim>::compute_limited_slope(
  const double                      &cell_average,
  const Tensor<1, dim>              &cell_average_grad,
  const std::vector<double>         &neighbor_cell_averages,
  const std::vector<Tensor<1, dim>> &neighbor_distance,
  const unsigned int                 n_neighbors,
  Tensor<1, dim>                    &limited_slope)
{
  switch (limiter_type)
    {
      case SlopeLimiterType::NoLimiter:
        {
          Assert(false,
                 ExcMessage("Slope limiter is set to NoLimiter, so this "
                            "function should not be called"));
        }

      case SlopeLimiterType::LinearReconstruction:
        {
          // This is a test case, and is not limiting the slope
          limited_slope = cell_average_grad;
          break;
        }

      case SlopeLimiterType::MinMod:
        {
          std::vector<Tensor<1, dim>> slopes(n_neighbors + 1);
          slopes[0] = cell_average_grad;
          for (unsigned int i = 0; i < n_neighbors; ++i)
            {
              slopes[i + 1] = (neighbor_cell_averages[i] - cell_average) /
                              (neighbor_distance[i].norm_square() / 2.) *
                              neighbor_distance[i];
            }
          minmod(slopes, n_neighbors + 1, limited_slope);
          break;
        }

      case SlopeLimiterType::MUSCL:
        {
          std::vector<Tensor<1, dim>> slopes(n_neighbors + 1);
          slopes[0] = cell_average_grad;
          for (unsigned int i = 0; i < n_neighbors; ++i)
            {
              slopes[i + 1] = (neighbor_cell_averages[i] - cell_average) /
                              neighbor_distance[i].norm_square() *
                              neighbor_distance[i];
            }
          minmod(slopes, n_neighbors + 1, limited_slope);
          break;
        }

      default:
        Assert(false, ExcNotImplemented());
        break;
    }
}

// explicit instantiation
template class Sapphire::Hydro::SlopeLimiter<1>;
template class Sapphire::Hydro::SlopeLimiter<2>;
template class Sapphire::Hydro::SlopeLimiter<3>;
