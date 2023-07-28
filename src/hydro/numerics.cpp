#include "numerics.h"

#include <iostream>

Sapphire::Hydro::HDSolverControl::HDSolverControl(
  const Sapphire::Utils::ParameterParser &prm)
  : scheme(prm.hdsolver_scheme)
  , limiter(prm.hdsolver_limiter)
  , limiter_criterion(prm.hdsolver_limiter_criterion)
  , fe_degree(prm.hdsolver_fe_degree)
  , time_step(prm.hdsolver_time_step)
  , end_time(prm.hdsolver_end_time)
  , refinement_level(prm.hdsolver_refinement_level)
  , max_iterations(prm.hdsolver_max_iterations)
  , tolerance(prm.hdsolver_tolerance)
{}

double
Sapphire::Hydro::minmod(const std::vector<double> &values)
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
Sapphire::Hydro::minmod(const std::vector<Tensor<1, dim>> &values,
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

// explicit instantiation
template void
Sapphire::Hydro::minmod<1>(const std::vector<Tensor<1, 1>> &values,
                           const unsigned int               n,
                           Tensor<1, 1>                    &return_value);
template void
Sapphire::Hydro::minmod<2>(const std::vector<Tensor<1, 2>> &values,
                           const unsigned int               n,
                           Tensor<1, 2>                    &return_value);
template void
Sapphire::Hydro::minmod<3>(const std::vector<Tensor<1, 3>> &values,
                           const unsigned int               n,
                           Tensor<1, 3>                    &return_value);

template <int dim>
void
Sapphire::Hydro::compute_limited_slope(
  const double                      &cell_average,
  const Tensor<1, dim>               cell_average_grad,
  const std::vector<double>         &neighbor_cell_averages,
  const std::vector<Tensor<1, dim>> &neighbor_distance,
  const unsigned int                 n_neighbors,
  Tensor<1, dim>                    &limited_slope,
  const HDSolverControl             &hd_solver_control)
{
  switch (hd_solver_control.limiter)
    {
      case SlopeLimiter::NoLimiter:
        {
          Assert(false,
                 ExcMessage("Slope limiter is set to NoLimiter, so this "
                            "function should not be called"));
        }

      case SlopeLimiter::LinearReconstruction:
        {
          // This is a test case, and is not limiting the slope
          limited_slope = cell_average_grad;
          break;
        }

      case SlopeLimiter::MinMod:
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

      case SlopeLimiter::MUSCL:
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
template void
Sapphire::Hydro::compute_limited_slope<1>(
  const double                    &cell_average,
  const Tensor<1, 1>               cell_average_grad,
  const std::vector<double>       &neighbor_cell_averages,
  const std::vector<Tensor<1, 1>> &neighbor_distance,
  const unsigned int               n_neighbors,
  Tensor<1, 1>                    &limited_slope,
  const HDSolverControl           &hd_solver_control);
template void
Sapphire::Hydro::compute_limited_slope<2>(
  const double                    &cell_average,
  const Tensor<1, 2>               cell_average_grad,
  const std::vector<double>       &neighbor_cell_averages,
  const std::vector<Tensor<1, 2>> &neighbor_distance,
  const unsigned int               n_neighbors,
  Tensor<1, 2>                    &limited_slope,
  const HDSolverControl           &hd_solver_control);
template void
Sapphire::Hydro::compute_limited_slope<3>(
  const double                    &cell_average,
  const Tensor<1, 3>               cell_average_grad,
  const std::vector<double>       &neighbor_cell_averages,
  const std::vector<Tensor<1, 3>> &neighbor_distance,
  const unsigned int               n_neighbors,
  Tensor<1, 3>                    &limited_slope,
  const HDSolverControl           &hd_solver_control);
