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


#include <deal.II/base/array_view.h>
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
       * @brief Convert a `dealii::Tensor` to a string.
       *
       * @tparam dim Dimension of the tensor.
       * @param value The tensor to convert.
       * @return std::string String representation of the tensor.
       */
      template <unsigned int dim>
      std::string
      tensor_to_string(const dealii::Tensor<1, dim, double> &value);



      /**
       * @brief Convert a list of tensors to a string.
       *
       * @tparam Container Some container of `dealii::Tensor`,
       *         e.g. `std::vector<dealii::Tensor<1,dim,double>>`.
       * @tparam dim Dimension of the tensor
       * @param values The list of tensors.
       * @return std::string String representation of the list of tensors.
       */
      template <typename Container, unsigned int dim>
      std::enable_if_t<std::is_same_v<typename Container::value_type,
                                      dealii::Tensor<1, dim, double>>,
                       std::string>
      tensor_list_to_string(const Container &values)
      {
        std::string s;
        for (const dealii::Tensor<1, dim, double> &v : values)
          s += tensor_to_string<dim>(v) + ";  ";
        if (!s.empty() && s.size() >= 3)
          s.erase(s.size() - 3); // Remove trailing separator
        return s;
      };



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
       * @brief Read a `.csv` file into a @dealref{Table}.
       *
       * Skips empty lines and lines starting with '#'.
       *
       * @note The table is only filled in one processor,
       *       and then replicated across all processors.
       *
       * @param filename Path to the `.csv` file.
       * @param n_rows Number of rows.
       * @param n_columns Number of columns.
       * @param delimiter Delimiter for the `.csv` file.
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



      /**
       * @brief Read a `.dat` file into a vector.
       *
       * The file should be structured the following way:
       *
       * ```data
       * # Column0  Column1    ... ColumnJ    ... ColumnN-1
       * data[0][0] data[1][0] ... data[j][0] ... data[N-1][0]
       * data[0][1] data[1][1] ... data[j][1] ... data[N-1][1]
       * ...
       * ```
       *
       * The number of columns must be know,
       * the number if rows is determined automatically.
       *
       * Skips empty lines and lines starting with '#'.
       *
       * @param filename Path to the `.dat` file.
       * @param n_columns Number of columns.
       * @param data_vector The return vector with the data,
       *        `data_vector[j][i]`
       *         with `j` the column index and `i` the row index.
       * @param delimiter Delimiter for the `.dat` file.
       */
      void
      read_dat_to_vector(const std::filesystem::path      &filename,
                         const unsigned int                n_columns,
                         std::vector<std::vector<double>> &data_vector,
                         const std::string                &delimiter = " ");



      /**
       * @brief Read a `.dat` file into as a
       *        multicomponent tensor product grid data,
       *        in the format as defined in the
       *        @dealref{dealii::InterpolatedTensorProductGridData,classFunctions_1_1InterpolatedTensorProductGridData}.
       *
       * The file should be a table of the following format:
       *
       * | `x_1` | `x_2` | `x_3` | `comp_1` | `comp_2` |
       * | ----- | ----- | ----- | -------- | -------- |
       * | x0    | y0    | z0    | v1(p0)   | v2(p0)   |
       * | x1    | y0    | z0    | v1(p1)   | v2(p1)   |
       * | x2    | y0    | z0    | v1(p2)   | v2(p2)   |
       * | x0    | y1    | z0    | v1(p3)   | v2(p3)   |
       *
       * @tparam dim Dimension of the data.
       * @param filename Path to the `.dat` file.
       * @param n_columns Number of total columns in the `.dat` file.
       * @param coordinate_values An array of dim arrays
       *        to be filled with the coordinate values.
       *        Each of the inner arrays contains the coordinate values
       *        \f$ x_0 , \dots , x_{K−1} \f$
       *        and similarly for the other coordinate directions.
       *        These arrays need not have the same size.
       * @param data_values An vector of dim-dimensional tables
       *        for each of the data components
       *        to be filled with data at each of the mesh points
       *        defined by the coordinate arrays above.
       * @param delimiter Delimiter for the `.dat` file.
       * @param col_start_coordinates Column index with the first coordinate.
       * @param col_start_data Column index with the first data component.
       */
      template <unsigned int dim>
      void
      read_dat_to_tensor_product_grid_data(
        const std::filesystem::path             &filename,
        const unsigned int                       n_columns,
        std::array<std::vector<double>, dim>    &coordinate_values,
        std::vector<dealii::Table<dim, double>> &data_values,
        const std::string                       &delimiter             = " ",
        const unsigned int                       col_start_coordinates = 0,
        const unsigned int                       col_start_data        = dim);



      /**
       * @brief Convert coordinate values to a tensor grid.
       *
       * The `coordinates` are the coordinates of the grid points,
       * and must be given with the first component running fastest:
       *
       * | `coordinates[0]` | `coordinates[1]` | `coordinates[2]` |
       * | ---------------- | ---------------- | ---------------- |
       * | \f$ x_0 \f$      | \f$ y_0 \f$      |  \f$ z_0 \f$     |
       * | \f$ x_1 \f$      | \f$ y_0 \f$      |  \f$ z_0 \f$     |
       * | \f$ x_2 \f$      | \f$ y_0 \f$      |  \f$ z_0 \f$     |
       * | \f$ x_0 \f$      | \f$ y_1 \f$      |  \f$ z_0 \f$     |
       * | \f$ x_1 \f$      | \f$ y_1 \f$      |  \f$ z_0 \f$     |
       * | \f$ x_2 \f$      | \f$ y_1 \f$      |  \f$ z_0 \f$     |
       * | \f$ \dots \f$    | \f$ \dots \f$    |  \f$ \dots \f$   |
       * | \f$ x_0 \f$      | \f$ y_0 \f$      |  \f$ z_1 \f$     |
       * | \f$ x_1 \f$      | \f$ y_0 \f$      |  \f$ z_1 \f$     |
       * | \f$ \dots \f$    | \f$ \dots \f$    |  \f$ \dots \f$   |
       *
       * The `tensor_grid` of the coordinates is extracted as:
       *
       * ```cpp
       *  tensor_grid[0] = [x_0, x_1, x_2]
       *  tensor_grid[1] = [y_0, y_1, ...]
       *  tensor_grid[2] = [z_0, z_1, ...]
       *  ```
       *
       * @tparam dim Dimension of the grid.
       * @param coordinates The coordinate values of the grid points.
       * @param tensor_grid The tensor grid of the coordinates.
       * @param last_coordinate_runs_fastest If set to true,
       *        the last coordinate runs the fastest
       *        and the first coordinate runs slowest.
       */
      template <unsigned int dim>
      void
      coordinates_to_tensor_grid(
        const std::array<dealii::ArrayView<double>, dim> &coordinates,
        std::array<std::vector<double>, dim>             &tensor_grid,
        const bool last_coordinate_runs_fastest = false);

    } // namespace Tools
  } // namespace Utils
} // namespace sapphirepp
#endif
