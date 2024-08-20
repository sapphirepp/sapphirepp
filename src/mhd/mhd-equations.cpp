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
 * @file mhd-equations.cpp
 * @author Florian Schulze (florian.schulze@mpi-hd.mpg.de)
 * @brief Implement @ref sapphirepp::MHD::MHDEquations
 */

#include "mhd-equations.h"

#include <deal.II/base/exceptions.h>


template <unsigned int dim>
sapphirepp::MHD::MHDEquations<dim>::MHDEquations() = default;



template <unsigned int dim>
std::vector<std::string>
sapphirepp::MHD::MHDEquations<dim>::create_component_name_list(
  const std::string &prefix)
{
  std::vector<std::string> component_names(n_components);

  component_names[density_component] = prefix + "rho";
  component_names[energy_component]  = prefix + "E";
  for (unsigned int d = 0; d < dim_uB; ++d)
    {
      component_names[first_momentum_component + d] =
        prefix + "p_" + std::to_string(d + 1);
      component_names[first_magnetic_component + d] =
        prefix + "B_" + std::to_string(d + 1);
    }

  return component_names;
}



template <unsigned int dim>
void
sapphirepp::MHD::MHDEquations<dim>::compute_flux_matrix(
  const state_type &state,
  flux_type        &flux_matrix) const
{
  /** @todo Implement this function. So far this is only an empty implementation
   * to test the code */

  AssertDimension(state.size(), n_components);

  for (unsigned int c = 0; c < n_components; ++c)
    {
      for (unsigned int d = 0; d < dim; ++d)
        {
          flux_matrix[c][d] = 0.;
        }
      // Test case: advect each component in x direction
      flux_matrix[c][0] = 1. * state[c];
    }
}



template <unsigned int dim>
double
sapphirepp::MHD::MHDEquations<dim>::compute_maximal_eigenvalue_normal(
  const state_type             &state,
  const dealii::Tensor<1, dim> &normal) const
{
  /** @todo Implement this eigenvalue calculation */
  static_cast<void>(state);
  static_cast<void>(normal);

  return 1.0;
}



// Explicit instantiations
template class sapphirepp::MHD::MHDEquations<1>;
template class sapphirepp::MHD::MHDEquations<2>;
template class sapphirepp::MHD::MHDEquations<3>;