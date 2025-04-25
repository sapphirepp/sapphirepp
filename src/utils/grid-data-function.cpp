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

#include <deal.II/base/exceptions.h>
#include <deal.II/base/table_indices.h>
#include <deal.II/base/utilities.h>

#include <fstream>
#include <iostream>
#include <sstream>

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



template <int dim>
void
sapphirepp::Utils::InterpolatedUniformGridData2<dim>::set_data(
  const std::array<std::pair<double, double>, dim> &new_interval_endpoints,
  const std::array<unsigned int, dim>              &new_n_subintervals,
  const dealii::Table<dim, double>                 &new_data_values)
{
  interval_endpoints = new_interval_endpoints;
  n_subintervals     = new_n_subintervals;
  data_values        = new_data_values;

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
void
sapphirepp::Utils::InterpolatedUniformGridData2<dim>::set_data(
  std::array<std::pair<double, double>, dim> &&new_interval_endpoints,
  std::array<unsigned int, dim>              &&new_n_subintervals,
  dealii::Table<dim, double>                 &&new_data_values)

{
  interval_endpoints = std::move(new_interval_endpoints);
  n_subintervals     = std::move(new_n_subintervals);
  data_values        = std::move(new_data_values);

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



template class sapphirepp::Utils::InterpolatedUniformGridData2<1>;
template class sapphirepp::Utils::InterpolatedUniformGridData2<2>;
template class sapphirepp::Utils::InterpolatedUniformGridData2<3>;



template <unsigned int dim>
sapphirepp::Utils::GridDataFunction<dim>::GridDataFunction(
  const std::filesystem::path &input_path,
  const std::string           &base_filename,
  const unsigned int           n_components,
  const double                 inital_time,
  const unsigned int           n_components_input)
  : Function<dim>(n_components, inital_time)
  , n_components_data(n_components_input == 0 ? n_components :
                                                n_components_input)
  , input_path(input_path)
  , base_filename(base_filename)
  , time_index(0)
  , grid_functions(n_components_data)
{
  {
    LogStream::Prefix p("GridDataFunction", saplog);
    read_data_hst(input_path, base_filename + ".hst", time_series);
  }

  const std::string filename = base_filename + ".block0.out1." +
                               Utilities::int_to_string(time_index, 5) + ".tab";
  load_data_from_file(filename);
  set_time(inital_time);
}



template <unsigned int dim>
double
sapphirepp::Utils::GridDataFunction<dim>::value(
  const Point<dim>  &p,
  const unsigned int component) const
{
  if (component >= n_components_data)
    return 0.0;

  return grid_functions[component].value(p);
}



template <unsigned int dim>
dealii::Tensor<1, dim>
sapphirepp::Utils::GridDataFunction<dim>::gradient(
  const Point<dim>  &p,
  const unsigned int component) const
{
  if (component >= n_components_data)
    return Tensor<1, dim>();

  return grid_functions[component].gradient(p);
}



template <unsigned int dim>
void
sapphirepp::Utils::GridDataFunction<dim>::set_time(const double new_time)
{
  Function<dim>::set_time(new_time);

  // Skip if in same time interval
  if (time_index == time_series.size() - 1)
    {
      if (new_time >= time_series[time_index])
        return;
    }
  else
    {
      AssertIndexRange(time_index + 1, time_series.size());
      if ((new_time >= time_series[time_index]) &&
          (new_time < time_series[time_index + 1]))
        return;
    }

  // Find new index
  // TODO: Improve
  for (time_index = 0; time_index < time_series.size(); ++time_index)
    {
      if (time_series[time_index] > new_time)
        {
          time_index--;
          break;
        }
    }
  if (time_index >= time_series.size())
    time_index = static_cast<unsigned int>(time_series.size() - 1);
  AssertIndexRange(time_index, time_series.size());

  // Load new data
  const std::string filename = base_filename + ".block0.out1." +
                               Utilities::int_to_string(time_index, 5) + ".tab";
  load_data_from_file(filename);
}



template <unsigned int dim>
void
sapphirepp::Utils::GridDataFunction<dim>::load_data_from_file(
  const std::string &filename)
{
  LogStream::Prefix p("GridDataFunction", saplog);
  saplog << "Load data from file: " << filename << std::endl;

  std::array<std::pair<double, double>, dim> interval_endpoints;
  std::array<unsigned int, dim>              n_subintervals;
  std::vector<Table<dim, double>>            data_values(n_components_data);
  read_data_tab(input_path,
                filename,
                interval_endpoints,
                n_subintervals,
                data_values,
                n_components_data);


  if (n_components_data == 5)
    {
      // Density
      grid_functions[0].set_data(std::move(interval_endpoints),
                                 std::move(n_subintervals),
                                 std::move(data_values[0]));
      // Energy
      grid_functions[4].set_data(std::move(interval_endpoints),
                                 std::move(n_subintervals),
                                 std::move(data_values[1]));
      // Momentum
      for (unsigned int d = 0; d < 3; ++d)
        grid_functions[1 + d].set_data(std::move(interval_endpoints),
                                       std::move(n_subintervals),
                                       std::move(data_values[2 + d]));
    }
  if (n_components_data == 8)
    {
      // Density
      grid_functions[0].set_data(std::move(interval_endpoints),
                                 std::move(n_subintervals),
                                 std::move(data_values[0]));
      // Energy
      grid_functions[4].set_data(std::move(interval_endpoints),
                                 std::move(n_subintervals),
                                 std::move(data_values[1]));
      for (unsigned int d = 0; d < 3; ++d)
        {
          // Momentum
          grid_functions[1 + d].set_data(std::move(interval_endpoints),
                                         std::move(n_subintervals),
                                         std::move(data_values[2 + d]));
          // Magnetic field
          grid_functions[5 + d].set_data(std::move(interval_endpoints),
                                         std::move(n_subintervals),
                                         std::move(data_values[5 + d]));
        }
    }
  else
    {
      for (unsigned int c = 0; c < n_components_data; c++)
        grid_functions[c].set_data(std::move(interval_endpoints),
                                   std::move(n_subintervals),
                                   std::move(data_values[c]));
    }
}



template <unsigned int dim>
void
sapphirepp::Utils::GridDataFunction<dim>::read_data_tab(
  const std::filesystem::path                &input_path,
  const std::string                          &filename,
  std::array<std::pair<double, double>, dim> &interval_endpoints,
  std::array<unsigned int, dim>              &n_subintervals,
  std::vector<Table<dim, double>>            &data_values,
  unsigned int                                n_components)
{
  AssertDimension(data_values.size(), n_components);
  Assert(dim == 1, ExcNotImplemented());
  saplog << "Read data from file: " << input_path / filename << std::endl;
  std::ifstream input_file(input_path / filename);
  AssertThrow(input_file.is_open(), ExcFileNotOpen(input_path / filename));

  // Table<2, double> data_table;
  // data_table.resize(N, M);

  std::vector<std::vector<double>> values(n_components);
  unsigned int                     index = 0;
  std::vector<std::string>         string_list;
  std::vector<double>              double_list;
  double                           x_min  = 0.;
  double                           x_max  = 0.;
  double                           x_last = 0.;
  unsigned int                     N;

  std::string line;
  while (std::getline(input_file, line))
    {
      // Skip lines that start with '#'
      if (line.empty() || line[0] == '#')
        continue;

      // string_list = Utilities::split_string_list(line, " ");
      string_list.clear();
      std::istringstream iss(line);
      std::string        token;
      while (iss >> token)
        string_list.push_back(token);

      double_list = Utilities::string_to_double(string_list);
      AssertDimension(double_list.size(), 1 + dim + n_components);

      for (unsigned int c = 0; c < n_components; ++c)
        values[c].push_back(double_list[2 + c]);

      const double x_coord = double_list[1];
      if (index == 0)
        {
          x_min  = x_coord;
          x_max  = x_coord;
          x_last = x_coord;
        }
      else
        {
          x_max                  = std::max(x_max, x_coord);
          N                      = index;
          const double dx        = x_coord - x_last;
          const double global_dx = (x_max - x_min) / static_cast<double>(N);
          const double eps       = 1e-6;
          AssertThrow(std::abs(dx - global_dx) < eps,
                      ExcMessage(
                        "Non-uniform grid detected in " +
                        std::string(input_path / filename) +
                        ": x_n - x_{n-1} = " + std::to_string(x_coord) + " - " +
                        std::to_string(x_last) + " = " + std::to_string(dx) +
                        " != dx = " + std::to_string(global_dx)));
        }

      x_last = x_coord;
      index++;
    }
  N = index;


  TableIndices<dim> table_indices;
  for (unsigned int d = 0; d < dim; ++d)
    table_indices[d] = 0;
  table_indices[0]             = N;
  n_subintervals[0]            = N - 1;
  interval_endpoints[0].first  = x_min;
  interval_endpoints[0].second = x_max;

  for (unsigned int c = 0; c < n_components; ++c)
    {
      AssertDimension(values[c].size(), N);
      data_values[c].reinit(table_indices);
      data_values[c].fill(values[c].begin());
    }
}



template <unsigned int dim>
void
sapphirepp::Utils::GridDataFunction<dim>::read_data_hst(
  const std::filesystem::path &input_path,
  const std::string           &filename,
  std::vector<double>         &time_series)
{
  saplog << "Read time series from file: " << input_path / filename
         << std::endl;
  time_series.clear();
  std::ifstream input_file(input_path / filename);
  AssertThrow(input_file.is_open(), ExcFileNotOpen(input_path / filename));

  unsigned int index  = 0;
  double       t_last = 0.;

  std::string line;
  while (std::getline(input_file, line))
    {
      // Skip lines that start with '#'
      if (line.empty() || line[0] == '#')
        continue;

      saplog << line << std::endl;

      std::istringstream iss(line);
      std::string        token;
      iss >> token;

      const double t_now = Utilities::string_to_double(token);

      if (index > 0)
        AssertThrow(t_now > t_last,
                    ExcMessage("Non-increasing time series detected in " +
                               std::string(input_path / filename) +
                               ": t_n = " + std::to_string(t_now) +
                               ",  t_{n-1} =  " + std::to_string(t_last)));

      time_series.push_back(t_now);
      t_last = t_now;
      index++;
    }
}



template class sapphirepp::Utils::GridDataFunction<1>;
template class sapphirepp::Utils::GridDataFunction<2>;
template class sapphirepp::Utils::GridDataFunction<3>;
