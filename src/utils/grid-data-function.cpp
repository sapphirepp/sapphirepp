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
 * @file grid-data-function.cpp
 * @author Florian Schulze (florian.schulze@mpi-hd.mpg.de)
 * @brief Implement @ref sapphirepp::Utils::GridDataFunction
 */

#include "grid-data-function.h"

// #include <deal.II/grid/grid_out.h>

#include <fstream>
#include <iostream>

#include "sapphirepp-logstream.h"



template <unsigned int dim, unsigned int spacedim>
sapphirepp::Utils::GridDataFunction<dim, spacedim>::GridDataFunction(
  const unsigned int n_components,
  const double       inital_time)
  : Function<spacedim>(n_components, inital_time)
{
  saplog << "Constructor" << std::endl;
}



template <unsigned int dim, unsigned int spacedim>
double
sapphirepp::Utils::GridDataFunction<dim, spacedim>::value(
  const Point<spacedim> &p,
  const unsigned int     component) const
{
  Point<dim> point;
  for (unsigned int d = 0; d < dim; ++d)
    point[d] = p[d];

  static_cast<void>(point);
  static_cast<void>(component);
  return 0.0;
}



template <unsigned int dim, unsigned int spacedim>
dealii::Tensor<1, spacedim>
sapphirepp::Utils::GridDataFunction<dim, spacedim>::gradient(
  const Point<spacedim> &p,
  const unsigned int     component) const
{
  Point<dim> point;
  for (unsigned int d = 0; d < dim; ++d)
    point[d] = p[d];

  static_cast<void>(point);
  static_cast<void>(component);

  Tensor<1, dim> gradient;

  Tensor<1, spacedim> return_gradient;
  for (unsigned int d = 0; d < dim; ++d)
    return_gradient[d] = gradient[d];

  return return_gradient;
}



template <unsigned int dim, unsigned int spacedim>
void
sapphirepp::Utils::GridDataFunction<dim, spacedim>::set_time(
  const double new_time)
{
  Function<spacedim>::set_time(new_time);
}



template class sapphirepp::Utils::GridDataFunction<1, 1>;
template class sapphirepp::Utils::GridDataFunction<1, 2>;
template class sapphirepp::Utils::GridDataFunction<1, 3>;
template class sapphirepp::Utils::GridDataFunction<2, 2>;
template class sapphirepp::Utils::GridDataFunction<2, 3>;
template class sapphirepp::Utils::GridDataFunction<3, 3>;
