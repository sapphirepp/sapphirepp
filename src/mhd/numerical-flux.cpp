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
 * @file numerical-flux.cpp
 * @author Florian Schulze (florian.schulze@mpi-hd.mpg.de)
 * @brief Implement @ref sapphirepp::MHD::NumericalFlux
 */

#include "numerical-flux.h"

#include <deal.II/base/exceptions.h>

#include <cmath>



template <unsigned int spacedim>
sapphirepp::MHD::NumericalFlux<spacedim>::NumericalFlux(
  const MHDEquations &mhd_equations)
  : mhd_equations{mhd_equations}
{}



template <unsigned int spacedim>
void
sapphirepp::MHD::NumericalFlux<spacedim>::compute_numerical_normal_flux(
  const dealii::Tensor<1, spacedim>       &normal,
  const typename MHDEquations::state_type &state_1,
  const typename MHDEquations::state_type &state_2,
  typename MHDEquations::state_type       &numerical_normal_flux) const
{
  /** @todo Make calculation for vector of points */

  AssertDimension(state_1.size(), MHDEquations::n_components);
  AssertDimension(state_2.size(), MHDEquations::n_components);
  AssertDimension(numerical_normal_flux.size(), MHDEquations::n_components);

  typename MHDEquations::flux_type flux_matrix_1, flux_matrix_2;

  mhd_equations.compute_flux_matrix(state_1, flux_matrix_1);
  mhd_equations.compute_flux_matrix(state_2, flux_matrix_2);

  const double max_eigenvalue_1 =
    mhd_equations.compute_maximum_normal_eigenvalue(state_1, normal);
  const double max_eigenvalue_2 =
    mhd_equations.compute_maximum_normal_eigenvalue(state_2, normal);
  const double max_eigenvalue = std::fmax(max_eigenvalue_1, max_eigenvalue_2);

  for (unsigned int c = 0; c < MHDEquations::n_components; ++c)
    {
      /** @todo Implement other numerical fluxes */

      // /** [Central Flux] */
      // numerical_normal_flux[c] =
      //   0.5 * (normal * flux_matrix_1[c] + normal * flux_matrix_2[c]);
      // /** [Central Flux] */

      /** [local Lax-Friedrichs Flux] */
      numerical_normal_flux[c] =
        0.5 * (normal * flux_matrix_1[c] + normal * flux_matrix_2[c]) -
        0.5 * max_eigenvalue * (state_2[c] - state_1[c]);
      /** [local Lax-Friedrichs Flux] */
    }
}



// explicit instantiation
template class sapphirepp::MHD::NumericalFlux<
  sapphirepp::MHD::MHDEquations::spacedim>;
