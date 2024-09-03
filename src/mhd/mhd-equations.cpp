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

#include <algorithm>
#include <cmath>



sapphirepp::MHD::MHDEquations::MHDEquations(const double adiabatic_index)
  : adiabatic_index{adiabatic_index}
{
  AssertThrow(adiabatic_index > 1.0,
              dealii::ExcMessage("Adiabatic index must be larger than 1.0."));
}



std::vector<std::string>
sapphirepp::MHD::MHDEquations::create_component_name_list(
  const std::string &prefix)
{
  std::vector<std::string> component_names(n_components);

  component_names[density_component] = prefix + "rho";
  component_names[energy_component]  = prefix + "E";
  for (unsigned int d = 0; d < spacedim; ++d)
    {
      component_names[first_momentum_component + d] = prefix + "p";
      component_names[first_magnetic_component + d] = prefix + "b";
    }

  return component_names;
}



std::vector<dealii::DataComponentInterpretation::DataComponentInterpretation>
sapphirepp::MHD::MHDEquations::create_component_interpretation_list()
{
  std::vector<dealii::DataComponentInterpretation::DataComponentInterpretation>
    data_component_interpretation(n_components);

  data_component_interpretation[density_component] =
    dealii::DataComponentInterpretation::component_is_scalar;
  data_component_interpretation[energy_component] =
    dealii::DataComponentInterpretation::component_is_scalar;
  for (unsigned int d = 0; d < spacedim; ++d)
    {
      data_component_interpretation[first_momentum_component + d] =
        dealii::DataComponentInterpretation::component_is_part_of_vector;
      data_component_interpretation[first_magnetic_component + d] =
        dealii::DataComponentInterpretation::component_is_part_of_vector;
    }

  return data_component_interpretation;
}



void
sapphirepp::MHD::MHDEquations::compute_flux_matrix(const state_type &state,
                                                   flux_type &flux_matrix) const
{
  AssertDimension(state.size(), n_components);

  const double pressure = compute_pressure(state);
  double       b2       = 0.;
  double       pb       = 0.;
  for (unsigned int d = 0; d < spacedim; ++d)
    {
      b2 += state[first_magnetic_component + d] *
            state[first_magnetic_component + d];
      pb += state[first_momentum_component + d] *
            state[first_magnetic_component + d];
    }

  for (unsigned int j = 0; j < spacedim; ++j)
    {
      flux_matrix[density_component][j] = state[first_momentum_component + j];

      flux_matrix[energy_component][j] =
        state[first_momentum_component + j] / state[density_component] *
          (state[energy_component] + pressure + 0.5 * b2) -
        1. / (state[density_component]) * pb *
          state[first_magnetic_component + j];

      for (unsigned int i = 0; i < spacedim; ++i)
        {
          flux_matrix[first_momentum_component + i][j] =
            state[first_momentum_component + j] *
              state[first_momentum_component + i] / state[density_component] -
            state[first_magnetic_component + j] *
              state[first_magnetic_component + i];

          flux_matrix[first_magnetic_component + i][j] =
            (state[first_momentum_component + j] *
               state[first_magnetic_component + i] -
             state[first_momentum_component + i] *
               state[first_magnetic_component + j]) /
            state[density_component];
        }

      flux_matrix[first_momentum_component + j][j] += pressure + 0.5 * b2;
    }
}



double
sapphirepp::MHD::MHDEquations::compute_maximum_normal_eigenvalue(
  const state_type                  &state,
  const dealii::Tensor<1, spacedim> &normal) const
{
  AssertDimension(state.size(), n_components);

  /** @todo Optimize this function */

  dealii::Vector<double> eigenvalues(8);
  compute_normale_eigenvalues(state, normal, eigenvalues);

  double max_eigenvalue =
    *std::max_element(eigenvalues.begin(), eigenvalues.end());
  double min_eigenvalue =
    *std::min_element(eigenvalues.begin(), eigenvalues.end());

  return std::max(std::fabs(max_eigenvalue), std::fabs(min_eigenvalue));
}



double
sapphirepp::MHD::MHDEquations::compute_pressure(const state_type &state) const
{
  AssertDimension(state.size(), n_components);

  double p2 = 0.;
  double b2 = 0.;
  for (unsigned int d = 0; d < spacedim; ++d)
    {
      p2 += state[first_momentum_component + d] *
            state[first_momentum_component + d];
      b2 += state[first_magnetic_component + d] *
            state[first_magnetic_component + d];
    }


  return (adiabatic_index - 1.) *
         (state[energy_component] - 1. / (2. * state[density_component]) * p2 -
          0.5 * b2);
}



void
sapphirepp::MHD::MHDEquations::compute_normale_eigenvalues(
  const state_type                  &state,
  const dealii::Tensor<1, spacedim> &normal,
  dealii::Vector<double>            &eigenvalues) const
{
  AssertDimension(state.size(), n_components);
  AssertDimension(eigenvalues.size(), 8);

  double b2 = 0.;
  for (unsigned int d = 0; d < spacedim; ++d)
    {
      b2 += state[first_magnetic_component + d] *
            state[first_magnetic_component + d];
    }
  double nu = 0.;
  double nb = 0.;
  for (unsigned int d = 0; d < spacedim; ++d)
    {
      nu += normal[d] * state[first_momentum_component + d] /
            state[density_component];
      nb += normal[d] * state[first_magnetic_component + d];
    }

  const double pressure = compute_pressure(state);
  const double a_s2     = adiabatic_index * pressure / state[density_component];
  const double c_a2     = b2 / state[density_component];
  const double d_n      = (a_s2 + c_a2) * (a_s2 + c_a2) -
                     4. * a_s2 * nb * nb / state[density_component];
  const double c_s2 = 0.5 * (a_s2 + c_a2 - std::sqrt(d_n));
  const double c_f2 = 0.5 * (a_s2 + c_a2 + std::sqrt(d_n));

  eigenvalues[0] = nu - std::sqrt(c_f2);
  eigenvalues[1] = nu - std::sqrt(c_a2);
  eigenvalues[2] = nu - std::sqrt(c_s2);
  eigenvalues[3] = nu;
  eigenvalues[4] = nu;
  eigenvalues[5] = nu + std::sqrt(c_s2);
  eigenvalues[6] = nu + std::sqrt(c_a2);
  eigenvalues[7] = nu + std::sqrt(c_f2);
}



void
sapphirepp::MHD::MHDEquations::convert_primitive_to_conserved(
  const state_type &primitive_state,
  state_type       &conserved_state) const
{
  AssertDimension(primitive_state.size(), n_components);
  AssertDimension(conserved_state.size(), n_components);

  double u2 = 0.;
  double b2 = 0.;
  for (unsigned int d = 0; d < spacedim; ++d)
    {
      u2 += primitive_state[first_velocity_component + d] *
            primitive_state[first_velocity_component + d];
      b2 += primitive_state[first_magnetic_component + d] *
            primitive_state[first_magnetic_component + d];
    }

  const double energy =
    0.5 * primitive_state[density_component] * u2 +
    primitive_state[pressure_component] / (adiabatic_index - 1.) + 0.5 * b2;

  conserved_state[density_component] = primitive_state[density_component];
  conserved_state[energy_component]  = energy;
  for (unsigned int d = 0; d < spacedim; ++d)
    {
      conserved_state[first_momentum_component + d] =
        primitive_state[density_component] *
        primitive_state[first_velocity_component + d];
      conserved_state[first_magnetic_component + d] =
        primitive_state[first_magnetic_component + d];
    }
}



void
sapphirepp::MHD::MHDEquations::convert_conserved_to_primitive(
  const state_type &conserved_state,
  state_type       &primitive_state) const
{
  AssertDimension(conserved_state.size(), n_components);
  AssertDimension(primitive_state.size(), n_components);

  primitive_state[density_component]  = conserved_state[density_component];
  primitive_state[pressure_component] = compute_pressure(conserved_state);
  for (unsigned int d = 0; d < spacedim; ++d)
    {
      primitive_state[first_velocity_component + d] =
        conserved_state[first_momentum_component + d] /
        conserved_state[density_component];
      primitive_state[first_magnetic_component + d] =
        conserved_state[first_magnetic_component + d];
    }
}
