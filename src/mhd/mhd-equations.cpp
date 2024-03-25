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

template <unsigned int dim>
sapphirepp::MHD::MHDEquations<dim>::MHDEquations() = default;



template <unsigned int dim>
std::vector<std::string>
sapphirepp::MHD::MHDEquations<dim>::create_component_name_list(
  const std::string &prefix)
{
  std::vector<std::string> component_names(MHDEquations<dim>::n_components);

  component_names[MHDEquations<dim>::density_component] = prefix + "rho";
  component_names[MHDEquations<dim>::energy_component]  = prefix + "E";
  for (unsigned int d = 0; d < dim; ++d)
    component_names[MHDEquations<dim>::momentum_component + d] =
      prefix + "u_" + std::to_string(d + 1);
  for (unsigned int d = 0; d < 3; ++d)
    component_names[MHDEquations<dim>::magnetic_component + d] =
      prefix + "B_" + std::to_string(d + 1);

  return component_names;
}



// Explicit instantiations
template class sapphirepp::MHD::MHDEquations<1>;
template class sapphirepp::MHD::MHDEquations<2>;
template class sapphirepp::MHD::MHDEquations<3>;