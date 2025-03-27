// -----------------------------------------------------------------------------
//
// Copyright (C) 2023 by the Sapphire++ authors
//
// This file is part of Sapphire++.
//
// Sapphire++ is free software: you can redistribute it and/or modify it under
// the terms of the GNU Lesser General Public License as published by the Free
// Software Foundation, either version 3 of the License, or (at your option) any
// later version.
//
// Sapphire++ is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
// A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more
// details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with Sapphire++. If not, see <https://www.gnu.org/licenses/>.
//
// -----------------------------------------------------------------------------

/**
 * @file slope-limiter.cpp
 * @author Florian Schulze (florian.schulze@mpi-hd.mpg.de)
 * @brief Implement utility functions related to slope limiting
 */

#include "slope-limiter.h"

#include <deal.II/base/exceptions.h>

#include <cmath>



template <unsigned int dim, bool divergence_cleaning>
double
sapphirepp::MHD::SlopeLimiter<dim, divergence_cleaning>::minmod(
  const std::vector<double> &values)
{
  auto [min_it, max_it] = std::minmax_element(values.begin(), values.end());
  if ((*min_it) * (*max_it) < 0.0)
    return 0.0;
  else if (std::abs(*min_it) < std::abs(*max_it))
    return *min_it;
  else
    return *max_it;
}



template <unsigned int dim, bool divergence_cleaning>
double
sapphirepp::MHD::SlopeLimiter<dim, divergence_cleaning>::minmod_gradients(
  const typename MHDEquations<dim, divergence_cleaning>::flux_type
    &cell_gradient,
  const std::vector<typename MHDEquations<dim, divergence_cleaning>::flux_type>
    &neighbor_gradients,
  typename MHDEquations<dim, divergence_cleaning>::flux_type &limited_gradient,
  const double                                                dx)
{
  constexpr unsigned int n_components =
    MHDEquations<dim, divergence_cleaning>::n_components;
  const double beta = 2.;
  const double M    = 0.;

  double              difference = 0.;
  std::vector<double> values;
  values.reserve(neighbor_gradients.size() + 1);
  for (unsigned int c = 0; c < n_components; ++c)
    {
      for (unsigned int d = 0; d < dim; ++d)
        {
          AssertIsFinite(cell_gradient[c][d]);
          if (std::abs(cell_gradient[c][d]) < M * dx * dx)
            {
              limited_gradient[c][d] = cell_gradient[c][d];
              continue;
            }

          values.push_back(cell_gradient[c][d]);

          for (const auto &tmp : neighbor_gradients)
            if (std::isfinite(tmp[c][d]))
              values.push_back(beta * tmp[c][d]);

          limited_gradient[c][d] = minmod(values);
          AssertIsFinite(limited_gradient[c][d]);

          difference += std::fabs(cell_gradient[c][d] - limited_gradient[c][d]);

          values.clear();
        }
    }

  difference /= static_cast<double>(n_components);
  return difference;
}



template <unsigned int dim, bool divergence_cleaning>
void
sapphirepp::MHD::SlopeLimiter<dim, divergence_cleaning>::
  enforce_divergence_free_limited_gradient(
    typename MHDEquations<dim, divergence_cleaning>::flux_type
      &limited_gradient)
{
  constexpr unsigned int first_magnetic_component =
    MHDEquations<dim, divergence_cleaning>::first_magnetic_component;
  double delta_p = 0;
  double delta_m = 0;

  for (unsigned int d = 0; d < dim; ++d)
    {
      delta_p +=
        std::max(limited_gradient[first_magnetic_component + d][d], 0.);
      delta_m +=
        std::max(-limited_gradient[first_magnetic_component + d][d], 0.);
    }

  const double delta = delta_p - delta_m;

  if (delta > 0)
    {
      for (unsigned int d = 0; d < dim; ++d)
        if (limited_gradient[first_magnetic_component + d][d] > 0)
          limited_gradient[first_magnetic_component + d][d] *=
            delta_m / delta_p;
    }
  else
    {
      for (unsigned int d = 0; d < dim; ++d)
        if (limited_gradient[first_magnetic_component + d][d] < 0)
          limited_gradient[first_magnetic_component + d][d] *=
            delta_p / delta_m;
    }
}



// Explicit instantiations
template class sapphirepp::MHD::SlopeLimiter<1, false>;
template class sapphirepp::MHD::SlopeLimiter<1, true>;
template class sapphirepp::MHD::SlopeLimiter<2, false>;
template class sapphirepp::MHD::SlopeLimiter<2, true>;
template class sapphirepp::MHD::SlopeLimiter<3, false>;
template class sapphirepp::MHD::SlopeLimiter<3, true>;
