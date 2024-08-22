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
sapphirepp::MHD::MHDEquations<dim>::MHDEquations(const double adiabatic_index)
  : adiabatic_index{adiabatic_index}
{
  AssertThrow(adiabatic_index > 1.0,
              dealii::ExcMessage("Adiabatic index must be larger than 1.0."));
}



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
  AssertDimension(state.size(), n_components);

  const double pressure = compute_pressure(state);
  double       B2       = 0.;
  double       pB       = 0.;
  for (unsigned int d = 0; d < dim_uB; ++d)
    {
      B2 += state[first_magnetic_component + d] *
            state[first_magnetic_component + d];
      pB += state[first_momentum_component + d] *
            state[first_magnetic_component + d];
    }

  for (unsigned int j = 0; j < dim; ++j)
    {
      flux_matrix[density_component][j] = state[first_momentum_component + j];

      flux_matrix[energy_component][j] =
        state[first_momentum_component + j] / state[density_component] *
          (state[energy_component] + pressure + 1. / (8. * M_PI) * B2) -
        1. / (4. * M_PI * state[density_component]) * pB *
          state[first_magnetic_component + j];

      for (unsigned int i = 0; i < dim_uB; ++i)
        {
          flux_matrix[first_momentum_component + i][j] =
            state[first_momentum_component + j] *
              state[first_momentum_component + i] / state[density_component] -
            1. / (4. * M_PI) * state[first_magnetic_component + j] *
              state[first_magnetic_component + i];

          flux_matrix[first_magnetic_component + i][j] =
            (state[first_momentum_component + j] *
               state[first_magnetic_component + i] -
             state[first_momentum_component + i] *
               state[first_magnetic_component + j]) /
            state[density_component];
        }

      flux_matrix[first_momentum_component + j][j] +=
        pressure + 1. / (8. * M_PI) * B2;
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



template <unsigned int dim>
double
sapphirepp::MHD::MHDEquations<dim>::compute_pressure(
  const state_type &state) const
{
  AssertDimension(state.size(), n_components);

  double p2 = 0.;
  double B2 = 0.;
  for (unsigned int d = 0; d < dim_uB; ++d)
    {
      p2 += state[first_momentum_component + d] *
            state[first_momentum_component + d];
      B2 += state[first_magnetic_component + d] *
            state[first_magnetic_component + d];
    }


  return (adiabatic_index - 1.) *
         (state[energy_component] - 1. / (2. * state[density_component]) * p2 -
          1 / (8. * M_PI) * B2);
}



// Explicit instantiations
template class sapphirepp::MHD::MHDEquations<1>;
template class sapphirepp::MHD::MHDEquations<2>;
template class sapphirepp::MHD::MHDEquations<3>;