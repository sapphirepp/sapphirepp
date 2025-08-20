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

#include <deal.II/base/exceptions.h>
#include <deal.II/base/table_indices.h>
#include <deal.II/base/utilities.h>

#include <fstream>
#include <iostream>
#include <sstream>

#include "sapphirepp-logstream.h"
#include "tools.h"



template <unsigned int dim>
sapphirepp::Utils::GridDataFunction<dim>::GridDataFunction(
  const std::filesystem::path &input_path,
  const std::string           &base_filename,
  const unsigned int           n_components,
  const double                 inital_time,
  const unsigned int           n_components_data,
  const std::string           &delimiter,
  const unsigned int           col_start_coordinates,
  const unsigned int           col_start_data,
  const bool                   last_coordinate_runs_fastest,
  const bool                   athena_ordering)
  : Function<dim>(n_components, inital_time)
  , input_path{input_path}
  , base_filename{base_filename}
  , delimiter{delimiter}
  , col_start_coordinates{col_start_coordinates}
  , col_start_data{col_start_data}
  , last_coordinate_runs_fastest{last_coordinate_runs_fastest}
  , athena_ordering{athena_ordering}
  , time_series{read_hst_to_time_series(input_path / (base_filename + ".hst"))}
  , time_index{0}
  , grid_functions(n_components_data == 0 ? n_components : n_components_data)
{
  LogStream::Prefix pre0("Setup", saplog);
  const std::string filename = base_filename + ".block0.out1." +
                               Utilities::int_to_string(time_index, 5) + ".tab";
  load_data_from_file(input_path / filename);
  set_time(inital_time);
}



template <unsigned int dim>
sapphirepp::Utils::GridDataFunction<dim>::GridDataFunction(
  const std::filesystem::path &filename,
  const unsigned int           n_components,
  const unsigned int           n_components_data,
  const std::string           &delimiter,
  const unsigned int           col_start_coordinates,
  const unsigned int           col_start_data,
  const bool                   last_coordinate_runs_fastest,
  const bool                   athena_ordering)
  : Function<dim>(n_components)
  , input_path{filename}
  , base_filename{""}
  , delimiter{delimiter}
  , col_start_coordinates{col_start_coordinates}
  , col_start_data{col_start_data}
  , last_coordinate_runs_fastest{last_coordinate_runs_fastest}
  , athena_ordering{athena_ordering}
  , time_series(1, 0.)
  , time_index{0}
  , grid_functions(n_components_data == 0 ? n_components : n_components_data)
{
  LogStream::Prefix pre0("Setup", saplog);
  load_data_from_file(filename);
}



template <unsigned int dim>
double
sapphirepp::Utils::GridDataFunction<dim>::value(
  const Point<dim>  &p,
  const unsigned int component) const
{
  if (component >= grid_functions.size())
    return 0.0;

  return grid_functions[component]->value(p);
}



template <unsigned int dim>
dealii::Tensor<1, dim>
sapphirepp::Utils::GridDataFunction<dim>::gradient(
  const Point<dim>  &p,
  const unsigned int component) const
{
  if (component >= grid_functions.size())
    return Tensor<1, dim>();

  return grid_functions[component]->gradient(p);
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
  load_data_from_file(input_path / filename);
}



template <unsigned int dim>
void
sapphirepp::Utils::GridDataFunction<dim>::load_data_from_file(
  const std::filesystem::path &filename)
{
  LogStream::Prefix p("GridDataFunction", saplog);
  saplog << "Load data from file: " << filename << std::endl;

  std::array<std::vector<double>, dim> coordinate_values;
  std::vector<Table<dim, double>>      data_values(grid_functions.size());
  const unsigned int                   n_columns =
    col_start_data + static_cast<unsigned int>(grid_functions.size());
  Tools::read_dat_to_tensor_product_grid_data<dim>(
    filename,
    n_columns,
    coordinate_values,
    data_values,
    delimiter,
    col_start_coordinates,
    col_start_data,
    last_coordinate_runs_fastest);

  std::array<std::pair<double, double>, dim> interval_endpoints;
  std::array<unsigned int, dim>              n_subintervals;
  for (unsigned int d = 0; d < dim; ++d)
    {
      const unsigned int coords_size =
        static_cast<unsigned int>(coordinate_values[d].size());
      n_subintervals[d]            = coords_size - 1;
      interval_endpoints[d].first  = coordinate_values[d][0];
      interval_endpoints[d].second = coordinate_values[d][coords_size - 1];
    }

  if (athena_ordering)
    {
      switch (grid_functions.size())
        {
          case 8:
            {
              for (unsigned int d = 0; d < 3; ++d)
                {
                  // Magnetic field
                  std::unique_ptr<
                    dealii::Functions::InterpolatedUniformGridData<dim>>
                    tmp_4(
                      new dealii::Functions::InterpolatedUniformGridData<dim>(
                        std::move(interval_endpoints),
                        std::move(n_subintervals),
                        std::move(data_values[5 + d])));
                  grid_functions[5 + d].swap(tmp_4);
                }
              [[fallthrough]];
            }
          case 5:
            {
              // Density
              std::unique_ptr<
                dealii::Functions::InterpolatedUniformGridData<dim>>
                tmp_1(new dealii::Functions::InterpolatedUniformGridData<dim>(
                  std::move(interval_endpoints),
                  std::move(n_subintervals),
                  std::move(data_values[0])));
              grid_functions[0].swap(tmp_1);
              // Energy
              std::unique_ptr<
                dealii::Functions::InterpolatedUniformGridData<dim>>
                tmp_2(new dealii::Functions::InterpolatedUniformGridData<dim>(
                  std::move(interval_endpoints),
                  std::move(n_subintervals),
                  std::move(data_values[1])));
              grid_functions[4].swap(tmp_2);
              for (unsigned int d = 0; d < 3; ++d)
                {
                  // Momentum
                  std::unique_ptr<
                    dealii::Functions::InterpolatedUniformGridData<dim>>
                    tmp_3(
                      new dealii::Functions::InterpolatedUniformGridData<dim>(
                        std::move(interval_endpoints),
                        std::move(n_subintervals),
                        std::move(data_values[2 + d])));
                  grid_functions[1 + d].swap(tmp_3);
                }
              break;
            }
          default:
            Assert(false,
                   ExcMessage(
                     "Expect either 5 or 8 components for Athena++ ordering."));
        }
    }
  else
    for (unsigned int c = 0; c < grid_functions.size(); c++)
      {
        std::unique_ptr<dealii::Functions::InterpolatedUniformGridData<dim>>
          tmp(new dealii::Functions::InterpolatedUniformGridData<dim>(
            std::move(interval_endpoints),
            std::move(n_subintervals),
            std::move(data_values[c])));
        grid_functions[c].swap(tmp);
      }
}



template <unsigned int dim>
std::vector<double>
sapphirepp::Utils::GridDataFunction<dim>::read_hst_to_time_series(
  const std::filesystem::path &filename)
{
  LogStream::Prefix pre0("Setup", saplog);
  LogStream::Prefix pew1("GridDataFunction", saplog);
  saplog << "Load time series from file" << filename << std::endl;
  const unsigned int               n_columns = 13;
  std::vector<std::vector<double>> data_vector(n_columns);

  Tools::read_dat_to_vector(filename, n_columns, data_vector, " ");
  const std::vector<double> &time_series = data_vector[0];

  for (unsigned int n = 1; n < time_series.size(); ++n)
    AssertThrow(time_series[n] > time_series[n - 1],
                ExcMessage("Non-increasing time series detected in " +
                           std::string(filename) + ": t_n = " +
                           std::to_string(time_series[n]) + ", t_{n-1} =  " +
                           std::to_string(time_series[n - 1])));

  return time_series;
}



template class sapphirepp::Utils::GridDataFunction<1>;
template class sapphirepp::Utils::GridDataFunction<2>;
template class sapphirepp::Utils::GridDataFunction<3>;
