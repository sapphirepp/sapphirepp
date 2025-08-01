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
 * @file tools.h
 * @author Florian Schulze (florian.schulze@mpi-hd.mpg.de)
 * @brief Define functions in namespace @ref sapphirepp::Utils::Tools.
 */

#ifndef UTILS_TOOLS_H
#define UTILS_TOOLS_H


#include <deal.II/base/point.h>
#include <deal.II/base/table.h>

#include <array>
#include <filesystem>
#include <string>
#include <vector>



namespace sapphirepp
{
  namespace Utils
  {

    /**
     * @namespace sapphirepp::Utils::Tools
     * @brief Namespace for general utility functions
     */
    namespace Tools
    {
      /**
       * @addtogroup numerical-parameters Numerical parameters
       * Parameters for the numerical algorithms.
       * @warning !DO NOT CHANGE unless you are aware what the parameters do!
       * @{
       */
      /** Precision for double / zero comparision. */
      static constexpr double epsilon_d = 1e-8;
      /** @} */



      /**
       * @brief Create a vector with `num` linear spaced values.
       *
       * @param start Starting value
       * @param stop End value
       * @param num Number of values
       * @return std::vector<double>
       */
      std::vector<double>
      create_linear_range(const double       start,
                          const double       stop,
                          const unsigned int num);



      /**
       * @brief Convert a point of dim1 into a point of dim2.
       *
       * For dim1 > dim2 it neglects the additional components,
       * for dim2 > dim1 it sets the additional components to 0.
       *
       * @tparam dim1 Dimension of input point
       * @tparam dim2 Dimension of output point
       * @param point @dealref{Point} of dimension dim1
       * @return dealii::Point<dim2> @dealref{Point} of dimension dim2
       */
      template <unsigned int dim1, unsigned int dim2>
      inline dealii::Point<dim2>
      convert_point_dimension(const dealii::Point<dim1> &point)
      {
        if constexpr (dim1 == dim2)
          return point;

        dealii::Point<dim2> point_out;
        for (unsigned int d = 0; d < dim1 && d < dim2; ++d)
          point_out[d] = point[d];
        return point_out;
      };



      /**
       * @brief Read a csv file into a @dealref{Table}.
       *
       * Skips empty lines and lines starting with '#'.
       *
       * @note The table is only filled in one processor,
       *       and then replicated across all processors.
       *
       * @param filename Path to the csv file
       * @param n_rows Number of rows
       * @param n_columns Number of columns
       * @param delimiter Delimiter for the csv file
       * @param transpose Transpose the table?
       * @param reverse_rows Reverse the ordering of rows?
       * @param reverse_columns Reverse the ordering of columns?
       * @return @dealref{dealii::Table<2\, double>,classTable_3_012_00_01T_01_4}
       *         of size n_rows x n_columns.
       *         (For transposed table: n_columns x n_rows.)
       */
      dealii::Table<2, double>
      read_csv_to_table(const std::filesystem::path &filename,
                        const unsigned int           n_rows,
                        const unsigned int           n_columns,
                        const std::string           &delimiter       = ",",
                        const bool                   transpose       = false,
                        const bool                   reverse_rows    = false,
                        const bool                   reverse_columns = false);

    } // namespace Tools
  } // namespace Utils
} // namespace sapphirepp
#endif
