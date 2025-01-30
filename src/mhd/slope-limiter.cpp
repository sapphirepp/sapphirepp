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



double
sapphirepp::MHD::SlopeLimiter::minmod(const std::vector<double> &values)
{
  auto [min_it, max_it] = std::minmax_element(values.begin(), values.end());
  if ((*min_it) * (*max_it) < 0.0)
    return 0.0;
  else if (std::abs(*min_it) < std::abs(*max_it))
    return *min_it;
  else
    return *max_it;
}



template <unsigned int dim>
double
sapphirepp::MHD::SlopeLimiter::minmod_gradients(
  const typename MHDEquations<dim>::flux_type              &cell_gradient,
  const std::vector<typename MHDEquations<dim>::flux_type> &neighbor_gradients,
  typename MHDEquations<dim>::flux_type                    &limited_gradient,
  const double                                              dx)
{
  const double beta = 2.;
  const double M    = 0.;

  double              difference = 0.;
  std::vector<double> values;
  values.reserve(neighbor_gradients.size() + 1);
  for (unsigned int c = 0; c < MHDEquations<dim>::n_components; ++c)
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

  difference /= static_cast<double>(MHDEquations<dim>::n_components);
  return difference;
}



/** @cond sapinternal */
// Explicit instantiations
template double
sapphirepp::MHD::SlopeLimiter::minmod_gradients<1>(
  const typename MHDEquations<1>::flux_type &,
  const std::vector<typename MHDEquations<1>::flux_type> &,
  typename MHDEquations<1>::flux_type &,
  const double);
template double
sapphirepp::MHD::SlopeLimiter::minmod_gradients<2>(
  const typename MHDEquations<2>::flux_type &,
  const std::vector<typename MHDEquations<2>::flux_type> &,
  typename MHDEquations<2>::flux_type &,
  const double);
template double
sapphirepp::MHD::SlopeLimiter::minmod_gradients<3>(
  const typename MHDEquations<3>::flux_type &,
  const std::vector<typename MHDEquations<3>::flux_type> &,
  typename MHDEquations<3>::flux_type &,
  const double);
/** @endcond */
