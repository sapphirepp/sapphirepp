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
          if (line.empty() || line[0] == '#')
            continue;

          std::vector<double> values = Utilities::string_to_double(
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
