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
 * @file mhd-postprocessor.cpp
 * @author Florian Schulze (florian.schulze@mpi-hd.mpg.de)
 * @brief Implement @ref sapphirepp::MHD::MHDPostprocessor
 */

#include "mhd-postprocessor.h"

#include <deal.II/base/exceptions.h>



template <unsigned int dim, bool divergence_cleaning>
sapphirepp::MHD::MHDPostprocessor<dim, divergence_cleaning>::MHDPostprocessor(
  const MHDEqs      &mhd_equations,
  const std::string &prefix)
  : mhd_equations{mhd_equations}
  , prefix(prefix)
{}



template <unsigned int dim, bool divergence_cleaning>
std::vector<std::string>
sapphirepp::MHD::MHDPostprocessor<dim, divergence_cleaning>::get_names() const
{
  std::vector<std::string> component_names(n_components_out);
  const char vec_component_name[MHDEqs::n_vec_components] = {'x', 'y', 'z'};

  component_names[pressure_component_out] = prefix + "P";
  for (unsigned int d = 0; d < dim; ++d)
    {
      component_names[first_velocity_component_out + d] = prefix + "u";
    }
  for (unsigned int d = dim; d < MHDEqs::n_vec_components; ++d)
    {
      component_names[first_velocity_component_out + d] =
        prefix + "u_" + vec_component_name[d];
    }

  return component_names;
}


template <unsigned int dim, bool divergence_cleaning>
std::vector<dealii::DataComponentInterpretation::DataComponentInterpretation>
sapphirepp::MHD::MHDPostprocessor<dim, divergence_cleaning>::
  get_data_component_interpretation() const
{
  std::vector<dealii::DataComponentInterpretation::DataComponentInterpretation>
    data_component_interpretation(n_components_out);

  data_component_interpretation[pressure_component_out] =
    dealii::DataComponentInterpretation::component_is_scalar;
  // Interpret components in the plain as vector components
  for (unsigned int d = 0; d < dim; ++d)
    {
      data_component_interpretation[first_velocity_component_out + d] =
        dealii::DataComponentInterpretation::component_is_part_of_vector;
    }
  // Interpret components out of the plain as scalars
  for (unsigned int d = dim; d < MHDEqs::n_vec_components; ++d)
    {
      data_component_interpretation[first_velocity_component_out + d] =
        dealii::DataComponentInterpretation::component_is_scalar;
    }

  return data_component_interpretation;
}



template <unsigned int dim, bool divergence_cleaning>
dealii::UpdateFlags
sapphirepp::MHD::MHDPostprocessor<dim, divergence_cleaning>::
  get_needed_update_flags() const
{
  return update_values;
}



template <unsigned int dim, bool divergence_cleaning>
inline void
sapphirepp::MHD::MHDPostprocessor<dim, divergence_cleaning>::
  evaluate_vector_field(const DataPostprocessorInputs::Vector<dim> &inputs,
                        std::vector<Vector<double>> &computed_quantities) const
{
  AssertDimension(computed_quantities.size(), inputs.solution_values.size());

  for (unsigned int q_index = 0; q_index < computed_quantities.size();
       ++q_index)
    {
      AssertDimension(computed_quantities[q_index].size(), n_components_out);

      const state_type &state = inputs.solution_values[q_index];
      AssertDimension(state.size(), MHDEqs::n_components);

      computed_quantities[q_index][pressure_component_out] =
        mhd_equations.compute_pressure_unsafe(state); // No warnings for saving

      for (unsigned int d = 0; d < MHDEqs::n_vec_components; ++d)
        {
          computed_quantities[q_index][first_velocity_component_out + d] =
            state[MHDEqs::first_momentum_component + d] /
            state[MHDEqs::density_component];
        }
    }
}



// explicit instantiation
template class sapphirepp::MHD::MHDPostprocessor<1, false>;
template class sapphirepp::MHD::MHDPostprocessor<1, true>;
template class sapphirepp::MHD::MHDPostprocessor<2, false>;
template class sapphirepp::MHD::MHDPostprocessor<2, true>;
template class sapphirepp::MHD::MHDPostprocessor<3, false>;
template class sapphirepp::MHD::MHDPostprocessor<3, true>;
