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
 * @file tools.cpp
 * @author Florian Schulze (florian.schulze@mpi-hd.mpg.de)
 * @brief Implement functions in namespace @ref sapphirepp::Utils::Tools.
 */


#include "tools.h"

#include <deal.II/base/exceptions.h>
#include <deal.II/base/mpi.h>
#include <deal.II/base/utilities.h>

#include <algorithm>

#include "sapphirepp-logstream.h"



std::vector<double>
sapphirepp::Utils::Tools::create_linear_range(const double       start,
                                              const double       stop,
                                              const unsigned int num)
{
  std::vector<double> values(num);
  if (num == 0)
    return values;

  for (unsigned int i = 0; i < num; ++i)
    values[i] = start + i * (stop - start) / (num - 1);
  // Sanitize
  values[0]       = start;
  values[num - 1] = stop;
  return values;
}



template <unsigned int dim>
std::string
sapphirepp::Utils::Tools::tensor_to_string(
  const dealii::Tensor<1, dim, double> &value)
{
  std::string s;
  for (unsigned int d = 0; d < dim; ++d)
    s += dealii::Utilities::to_string(value[d]) + " ";
  if (!s.empty())
    s.erase(s.size() - 1); // Remove trailing separator
  return s;
}


/** @cond sapinternal */
// Explicit instantiations of to_sting functions
template std::string
sapphirepp::Utils::Tools::tensor_to_string<1>(
  const dealii::Tensor<1, 1, double> &);
template std::string
sapphirepp::Utils::Tools::tensor_to_string<2>(
  const dealii::Tensor<1, 2, double> &);
template std::string
sapphirepp::Utils::Tools::tensor_to_string<3>(
  const dealii::Tensor<1, 3, double> &);

template std::string
sapphirepp::Utils::Tools::tensor_list_to_string<
  std::vector<dealii::Tensor<1, 1, double>>,
  1>(const std::vector<dealii::Tensor<1, 1, double>> &);
template std::string
sapphirepp::Utils::Tools::tensor_list_to_string<
  std::vector<dealii::Tensor<1, 2, double>>,
  2>(const std::vector<dealii::Tensor<1, 2, double>> &);
template std::string
sapphirepp::Utils::Tools::tensor_list_to_string<
  std::vector<dealii::Tensor<1, 3, double>>,
  3>(const std::vector<dealii::Tensor<1, 3, double>> &);
/** @endcond */



dealii::Table<2, double>
sapphirepp::Utils::Tools::read_csv_to_table(
  const std::filesystem::path &filename,
  const unsigned int           n_rows,
  const unsigned int           n_cols,
  const std::string           &delimiter,
  const bool                   transpose,
  const bool                   reverse_rows,
  const bool                   reverse_cols)
{
  using namespace dealii;
  const MPI_Comm     mpi_communicator = MPI_COMM_WORLD;
  const unsigned int root_rank        = 0;

  Table<2, double> data_table;

  if (Utilities::MPI::this_mpi_process(mpi_communicator) == root_rank)
    {
      saplog << "Read file " << filename << std::endl;
      std::ifstream input_file(filename);
      AssertThrow(input_file.is_open(), dealii::ExcFileNotOpen(filename));

      if (transpose)
        data_table.reinit(n_cols, n_rows);
      else
        data_table.reinit(n_rows, n_cols);

      std::string  line;
      unsigned int i = 0;
      while (std::getline(input_file, line))
        {
          // Skip lines that start with '#'
          line = Utilities::trim(line);
          if (line.empty() || Utilities::match_at_string_start(line, "#"))
            continue;

          // Clean line (needed e.g. for multiple spaces with ' ' delimiter)
          line = Utilities::replace_in_string(line, "    ", " ");
          line = Utilities::replace_in_string(line, "   ", " ");
          line = Utilities::replace_in_string(line, "  ", " ");

          const std::vector<double> values = Utilities::string_to_double(
            Utilities::split_string_list(line, delimiter));
          AssertThrow(values.size() == n_cols,
                      dealii::ExcDimensionMismatch(values.size(), n_cols));
          AssertThrow(i < n_rows,
                      dealii::ExcIndexRange(values.size(), 0, n_rows));

          for (unsigned int j = 0; j < data_table.size(1); ++j)
            if (transpose)
              data_table(reverse_cols ? n_cols - (j + 1) : j,
                         reverse_rows ? n_rows - (i + 1) : i) = values[j];
            else
              data_table(reverse_rows ? n_rows - (i + 1) : i,
                         reverse_cols ? n_cols - (j + 1) : j) = values[j];

          ++i;
        }
      AssertThrow(i == n_rows, dealii::ExcDimensionMismatch(i, n_rows));

      input_file.close();
    }

  // Now distribute to all processes
  data_table.replicate_across_communicator(mpi_communicator, root_rank);

  return data_table;
}



void
sapphirepp::Utils::Tools::read_dat_to_vector(
  const std::filesystem::path      &filename,
  const unsigned int                n_columns,
  std::vector<std::vector<double>> &data_vector,
  const std::string                &delimiter)
{
  using namespace dealii;
  AssertDimension(data_vector.size(), n_columns);

  saplog << "Read file " << filename << std::endl;
  std::ifstream input_file(filename);
  AssertThrow(input_file.is_open(), dealii::ExcFileNotOpen(filename));

  for (unsigned int j = 0; j < n_columns; ++j)
    data_vector[j].clear();

  std::string line;
  while (std::getline(input_file, line))
    {
      // Skip lines that start with '#'
      line = Utilities::trim(line);
      if (line.empty() || Utilities::match_at_string_start(line, "#"))
        continue;

      // Clean line (needed e.g. for multiple spaces with ' ' delimiter)
      line = Utilities::replace_in_string(line, "    ", " ");
      line = Utilities::replace_in_string(line, "   ", " ");
      line = Utilities::replace_in_string(line, "  ", " ");

      const std::vector<double> values = Utilities::string_to_double(
        Utilities::split_string_list(line, delimiter));

      AssertThrow(values.size() == n_columns,
                  dealii::ExcDimensionMismatch(values.size(), n_columns));

      for (unsigned int j = 0; j < n_columns; ++j)
        data_vector[j].push_back(values[j]);
    }

  input_file.close();
}



template <unsigned int dim>
void
sapphirepp::Utils::Tools::read_dat_to_tensor_product_grid_data(
  const std::filesystem::path             &filename,
  const unsigned int                       n_columns,
  std::array<std::vector<double>, dim>    &coordinate_values,
  std::vector<dealii::Table<dim, double>> &data_values,
  const std::string                       &delimiter,
  const unsigned int                       col_start_coordinates,
  const unsigned int                       col_start_data,
  const bool                               last_coordinate_runs_fastest)
{
  using namespace dealii;
  const MPI_Comm     mpi_communicator = MPI_COMM_WORLD;
  const unsigned int root_rank        = 0;
  LogStream::Prefix  pre("ReadTensorGrid", saplog);

  Assert(col_start_coordinates + dim <= col_start_data,
         ExcMessage("Coordinates must be given in columns before data."));
  AssertDimension(data_values.size(), n_columns - col_start_data);

  if (Utilities::MPI::this_mpi_process(mpi_communicator) == root_rank)
    {
      std::vector<std::vector<double>> data_vector(n_columns);
      for (auto &vector_comp : data_vector)
        vector_comp.reserve(data_values[0].n_elements());

      read_dat_to_vector(filename, n_columns, data_vector, delimiter);

      AssertThrow(data_vector[0].size() > 0,
                  ExcMessage("Empty file: " + std::string(filename)));

      std::array<dealii::ArrayView<double>, dim> coordinates;
      for (unsigned int d = 0; d < dim; ++d)
        coordinates[d] =
          make_array_view(data_vector[col_start_coordinates + d]);
      coordinates_to_tensor_grid<dim>(coordinates,
                                      coordinate_values,
                                      last_coordinate_runs_fastest);

      TableIndices<dim> table_indices;
      for (unsigned int d = 0; d < dim; ++d)
        table_indices[d] = coordinate_values[d].size();

      saplog << "Grid of size " << table_indices << " with "
             << data_vector[0].size() << " points" << std::endl;
      saplog << "Grid starts at ";
      for (unsigned int d = 0; d < dim; ++d)
        saplog << coordinate_values[d][0] << " ";
      saplog << "and ends at ";
      for (unsigned int d = 0; d < dim; ++d)
        saplog << coordinate_values[d][coordinate_values[d].size() - 1] << " ";
      saplog << std::endl;


      for (unsigned int c = 0; c < data_values.size(); ++c)
        {
          data_values[c].reinit(table_indices);
          data_values[c].fill(data_vector[col_start_data + c].begin(),
                              last_coordinate_runs_fastest);
        }
    }

  // Now distribute to all processes
  for (unsigned int d = 0; d < dim; ++d)
    {
      const std::size_t count =
        Utilities::MPI::broadcast(mpi_communicator,
                                  coordinate_values[d].size(),
                                  root_rank);
      coordinate_values[d].resize(count);
      Utilities::MPI::broadcast(coordinate_values[d].data(),
                                count,
                                root_rank,
                                mpi_communicator);
    }

  for (unsigned int c = 0; c < data_values.size(); ++c)
    data_values[c].replicate_across_communicator(mpi_communicator, root_rank);
}


/** @cond sapinternal */
// Explicit instantiation
template void
sapphirepp::Utils::Tools::read_dat_to_tensor_product_grid_data<1>(
  const std::filesystem::path &filename,
  const unsigned int,
  std::array<std::vector<double>, 1> &,
  std::vector<dealii::Table<1, double>> &,
  const std::string &,
  const unsigned int,
  const unsigned int,
  const bool);
template void
sapphirepp::Utils::Tools::read_dat_to_tensor_product_grid_data<2>(
  const std::filesystem::path &filename,
  const unsigned int,
  std::array<std::vector<double>, 2> &,
  std::vector<dealii::Table<2, double>> &,
  const std::string &,
  const unsigned int,
  const unsigned int,
  const bool);
template void
sapphirepp::Utils::Tools::read_dat_to_tensor_product_grid_data<3>(
  const std::filesystem::path &filename,
  const unsigned int,
  std::array<std::vector<double>, 3> &,
  std::vector<dealii::Table<3, double>> &,
  const std::string &,
  const unsigned int,
  const unsigned int,
  const bool);
/** @endcond */



template <unsigned int dim>
void
sapphirepp::Utils::Tools::coordinates_to_tensor_grid(
  const std::array<dealii::ArrayView<double>, dim> &coordinates_in,
  std::array<std::vector<double>, dim>             &tensor_grid,
  const bool last_coordinate_runs_fastest)
{
  const size_t n_rows = coordinates_in[0].size();

  Assert(n_rows > 0, dealii::ExcMessage("Empty coordinate vector"));
  if constexpr (dim > 1)
    AssertDimension(coordinates_in[1].size(), coordinates_in[0].size());
  if constexpr (dim > 2)
    AssertDimension(coordinates_in[2].size(), coordinates_in[0].size());

  std::array<dealii::ArrayView<double>, dim> coordinates = coordinates_in;
  if (last_coordinate_runs_fastest)
    std::reverse(coordinates.begin(), coordinates.end());

  std::array<unsigned int, dim> i_coord;
  for (unsigned int d = 0; d < dim; ++d)
    {
      i_coord[d] = 0;
      tensor_grid[d].clear();
      tensor_grid[d].push_back(coordinates[d][0]);
    }

  for (unsigned int i = 1; i < n_rows; ++i)
    {
      // Find wich component changed
      unsigned int d;
      for (d = 0; d < dim; ++d)
        if (coordinates[d][i] > tensor_grid[d][i_coord[d]])
          break;
      AssertThrow(d < dim,
                  dealii::ExcMessage("Wrong coordinate ordering at i=" +
                                     std::to_string(i)));

      // Faster running components rest to 0
      for (unsigned int k = 0; k < d; ++k)
        {
          i_coord[k] = 0;
          AssertThrow(std::abs(coordinates[k][i] - tensor_grid[k][i_coord[k]]) <
                        epsilon_d,
                      dealii::ExcMessage("Wrong coordinate ordering at i=" +
                                         std::to_string(i)));
        }

      // Slower running components stay the same
      bool new_value = true;
      for (unsigned int k = d + 1; k < dim; ++k)
        {
          AssertThrow(std::abs(coordinates[k][i] - tensor_grid[k][i_coord[k]]) <
                        epsilon_d,
                      dealii::ExcMessage("Wrong coordinate ordering at i=" +
                                         std::to_string(i)));
          new_value = new_value && (i_coord[k] == 0);
        }

      // Only add if new value (i.e. slower components still at 0)
      if (new_value)
        tensor_grid[d].push_back(coordinates[d][i]);
      ++i_coord[d];
    }

  if (last_coordinate_runs_fastest)
    std::reverse(tensor_grid.begin(), tensor_grid.end());
}


/** @cond sapinternal */
// Explicit instantiation
template void
sapphirepp::Utils::Tools::coordinates_to_tensor_grid<1>(
  const std::array<dealii::ArrayView<double>, 1> &,
  std::array<std::vector<double>, 1> &,
  const bool);
template void
sapphirepp::Utils::Tools::coordinates_to_tensor_grid<2>(
  const std::array<dealii::ArrayView<double>, 2> &,
  std::array<std::vector<double>, 2> &,
  const bool);
template void
sapphirepp::Utils::Tools::coordinates_to_tensor_grid<3>(
  const std::array<dealii::ArrayView<double>, 3> &,
  std::array<std::vector<double>, 3> &,
  const bool);
/** @endcond */
