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



namespace
{
  using namespace dealii;

  // interpolate a data value from a table where ix denotes
  // the (lower) left endpoint of the interval to interpolate
  // in, and p_unit denotes the point in unit coordinates to do so.
  double
  interpolate(const Table<1, double> &data_values,
              const TableIndices<1>  &ix,
              const Point<1>         &xi)
  {
    return ((1 - xi[0]) * data_values[ix[0]] + xi[0] * data_values[ix[0] + 1]);
  }

  double
  interpolate(const Table<2, double> &data_values,
              const TableIndices<2>  &ix,
              const Point<2>         &p_unit)
  {
    return (((1 - p_unit[0]) * data_values[ix[0]][ix[1]] +
             p_unit[0] * data_values[ix[0] + 1][ix[1]]) *
              (1 - p_unit[1]) +
            ((1 - p_unit[0]) * data_values[ix[0]][ix[1] + 1] +
             p_unit[0] * data_values[ix[0] + 1][ix[1] + 1]) *
              p_unit[1]);
  }

  double
  interpolate(const Table<3, double> &data_values,
              const TableIndices<3>  &ix,
              const Point<3>         &p_unit)
  {
    return ((((1 - p_unit[0]) * data_values[ix[0]][ix[1]][ix[2]] +
              p_unit[0] * data_values[ix[0] + 1][ix[1]][ix[2]]) *
               (1 - p_unit[1]) +
             ((1 - p_unit[0]) * data_values[ix[0]][ix[1] + 1][ix[2]] +
              p_unit[0] * data_values[ix[0] + 1][ix[1] + 1][ix[2]]) *
               p_unit[1]) *
              (1 - p_unit[2]) +
            (((1 - p_unit[0]) * data_values[ix[0]][ix[1]][ix[2] + 1] +
              p_unit[0] * data_values[ix[0] + 1][ix[1]][ix[2] + 1]) *
               (1 - p_unit[1]) +
             ((1 - p_unit[0]) * data_values[ix[0]][ix[1] + 1][ix[2] + 1] +
              p_unit[0] * data_values[ix[0] + 1][ix[1] + 1][ix[2] + 1]) *
               p_unit[1]) *
              p_unit[2]);
  }


  // Interpolate the gradient of a data value from a table where ix
  // denotes the lower left endpoint of the interval to interpolate
  // in, p_unit denotes the point in unit coordinates, and dx
  // denotes the width of the interval in each dimension.
  Tensor<1, 1>
  gradient_interpolate(const Table<1, double> &data_values,
                       const TableIndices<1>  &ix,
                       const Point<1>         &p_unit,
                       const Point<1>         &dx)
  {
    (void)p_unit;
    Tensor<1, 1> grad;
    grad[0] = (data_values[ix[0] + 1] - data_values[ix[0]]) / dx[0];
    return grad;
  }


  Tensor<1, 2>
  gradient_interpolate(const Table<2, double> &data_values,
                       const TableIndices<2>  &ix,
                       const Point<2>         &p_unit,
                       const Point<2>         &dx)
  {
    Tensor<1, 2> grad;
    double u00 = data_values[ix[0]][ix[1]], u01 = data_values[ix[0] + 1][ix[1]],
           u10 = data_values[ix[0]][ix[1] + 1],
           u11 = data_values[ix[0] + 1][ix[1] + 1];

    grad[0] = ((1 - p_unit[1]) * (u01 - u00) + p_unit[1] * (u11 - u10)) / dx[0];
    grad[1] = ((1 - p_unit[0]) * (u10 - u00) + p_unit[0] * (u11 - u01)) / dx[1];
    return grad;
  }


  Tensor<1, 3>
  gradient_interpolate(const Table<3, double> &data_values,
                       const TableIndices<3>  &ix,
                       const Point<3>         &p_unit,
                       const Point<3>         &dx)
  {
    Tensor<1, 3> grad;
    double       u000 = data_values[ix[0]][ix[1]][ix[2]],
           u001       = data_values[ix[0] + 1][ix[1]][ix[2]],
           u010       = data_values[ix[0]][ix[1] + 1][ix[2]],
           u100       = data_values[ix[0]][ix[1]][ix[2] + 1],
           u011       = data_values[ix[0] + 1][ix[1] + 1][ix[2]],
           u101       = data_values[ix[0] + 1][ix[1]][ix[2] + 1],
           u110       = data_values[ix[0]][ix[1] + 1][ix[2] + 1],
           u111       = data_values[ix[0] + 1][ix[1] + 1][ix[2] + 1];

    grad[0] = ((1 - p_unit[2]) *
                 ((1 - p_unit[1]) * (u001 - u000) + p_unit[1] * (u011 - u010)) +
               p_unit[2] * ((1 - p_unit[1]) * (u101 - u100) +
                            p_unit[1] * (u111 - u110))) /
              dx[0];
    grad[1] = ((1 - p_unit[2]) *
                 ((1 - p_unit[0]) * (u010 - u000) + p_unit[0] * (u011 - u001)) +
               p_unit[2] * ((1 - p_unit[0]) * (u110 - u100) +
                            p_unit[0] * (u111 - u101))) /
              dx[1];
    grad[2] = ((1 - p_unit[1]) *
                 ((1 - p_unit[0]) * (u100 - u000) + p_unit[0] * (u101 - u001)) +
               p_unit[1] * ((1 - p_unit[0]) * (u110 - u010) +
                            p_unit[0] * (u111 - u011))) /
              dx[2];

    return grad;
  }
} // namespace



template <int dim>
sapphirepp::Utils::InterpolatedUniformGridData2<dim>::
  InterpolatedUniformGridData2(
    const std::array<std::pair<double, double>, dim> &interval_endpoints,
    const std::array<unsigned int, dim>              &n_subintervals,
    const dealii::Table<dim, double>                 &data_values)
  : interval_endpoints(interval_endpoints)
  , n_subintervals(n_subintervals)
  , data_values(data_values)
{
  for (unsigned int d = 0; d < dim; ++d)
    {
      Assert(n_subintervals[d] >= 1,
             ExcMessage("There needs to be at least one subinterval in each "
                        "coordinate direction."));
      Assert(interval_endpoints[d].first < interval_endpoints[d].second,
             ExcMessage("The interval in each coordinate direction needs "
                        "to have positive size"));
      Assert(data_values.size()[d] == n_subintervals[d] + 1,
             ExcMessage("The data table does not have the correct size."));
    }
}



template <int dim>
sapphirepp::Utils::InterpolatedUniformGridData2<dim>::
  InterpolatedUniformGridData2(
    std::array<std::pair<double, double>, dim> &&interval_endpoints,
    std::array<unsigned int, dim>              &&n_subintervals,
    dealii::Table<dim, double>                 &&data_values)
  : interval_endpoints(std::move(interval_endpoints))
  , n_subintervals(std::move(n_subintervals))
  , data_values(std::move(data_values))
{
  for (unsigned int d = 0; d < dim; ++d)
    {
      Assert(this->n_subintervals[d] >= 1,
             ExcMessage("There needs to be at least one subinterval in each "
                        "coordinate direction."));
      Assert(this->interval_endpoints[d].first <
               this->interval_endpoints[d].second,
             ExcMessage("The interval in each coordinate direction needs "
                        "to have positive size"));
      Assert(this->data_values.size()[d] == this->n_subintervals[d] + 1,
             ExcMessage("The data table does not have the correct size."));
    }
}



template <int dim>
double
sapphirepp::Utils::InterpolatedUniformGridData2<dim>::value(
  const Point<dim>  &p,
  const unsigned int component) const
{
  (void)component;
  Assert(
    component == 0,
    ExcMessage(
      "This is a scalar function object, the component can only be zero."));

  // find out where this data point lies, relative to the given
  // subdivision points
  TableIndices<dim> ix;
  for (unsigned int d = 0; d < dim; ++d)
    {
      const double delta_x =
        ((interval_endpoints[d].second - interval_endpoints[d].first) /
         n_subintervals[d]);
      if (p[d] <= interval_endpoints[d].first)
        ix[d] = 0;
      else if (p[d] >= interval_endpoints[d].second - delta_x)
        ix[d] = n_subintervals[d] - 1;
      else
        ix[d] = static_cast<unsigned int>((p[d] - interval_endpoints[d].first) /
                                          delta_x);
    }

  // now compute the relative point within the interval/rectangle/box
  // defined by the point coordinates found above. truncate below and
  // above to accommodate points that may lie outside the range
  Point<dim> p_unit;
  for (unsigned int d = 0; d < dim; ++d)
    {
      const double delta_x =
        ((interval_endpoints[d].second - interval_endpoints[d].first) /
         n_subintervals[d]);
      p_unit[d] = std::max(std::min((p[d] - interval_endpoints[d].first -
                                     static_cast<double>(ix[d]) * delta_x) /
                                      delta_x,
                                    1.),
                           0.);
    }

  return interpolate(data_values, ix, p_unit);
}



template <int dim>
dealii::Tensor<1, dim>
sapphirepp::Utils::InterpolatedUniformGridData2<dim>::gradient(
  const Point<dim>  &p,
  const unsigned int component) const
{
  (void)component;
  Assert(
    component == 0,
    ExcMessage(
      "This is a scalar function object, the component can only be zero."));

  // find out where this data point lies, relative to the given
  // subdivision points
  TableIndices<dim> ix;
  for (unsigned int d = 0; d < dim; ++d)
    {
      const double delta_x = ((this->interval_endpoints[d].second -
                               this->interval_endpoints[d].first) /
                              this->n_subintervals[d]);
      if (p[d] <= this->interval_endpoints[d].first)
        ix[d] = 0;
      else if (p[d] >= this->interval_endpoints[d].second - delta_x)
        ix[d] = this->n_subintervals[d] - 1;
      else
        ix[d] = static_cast<unsigned int>(
          (p[d] - this->interval_endpoints[d].first) / delta_x);
    }

  // now compute the relative point within the interval/rectangle/box
  // defined by the point coordinates found above. truncate below and
  // above to accommodate points that may lie outside the range
  Point<dim> p_unit;
  Point<dim> delta_x;
  for (unsigned int d = 0; d < dim; ++d)
    {
      delta_x[d] = ((this->interval_endpoints[d].second -
                     this->interval_endpoints[d].first) /
                    this->n_subintervals[d]);
      p_unit[d]  = std::max(std::min((p[d] - this->interval_endpoints[d].first -
                                     static_cast<double>(ix[d]) * delta_x[d]) /
                                      delta_x[d],
                                    1.),
                           0.);
    }

  return gradient_interpolate(this->data_values, ix, p_unit, delta_x);
}



template <int dim>
std::size_t
sapphirepp::Utils::InterpolatedUniformGridData2<dim>::memory_consumption() const
{
  return sizeof(*this) + data_values.memory_consumption() - sizeof(data_values);
}



template <int dim>
const dealii::Table<dim, double> &
sapphirepp::Utils::InterpolatedUniformGridData2<dim>::get_data() const
{
  return data_values;
}



template class sapphirepp::Utils::InterpolatedUniformGridData2<1>;
template class sapphirepp::Utils::InterpolatedUniformGridData2<2>;
template class sapphirepp::Utils::InterpolatedUniformGridData2<3>;



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
