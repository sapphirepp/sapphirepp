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
 * @file grid-data-function.h
 * @author Florian Schulze (florian.schulze@mpi-hd.mpg.de)
 * @brief Define @ref sapphirepp::Utils::GridDataFunction
 */

#ifndef UTILS_GRIDDATAFUNCTION_H
#define UTILS_GRIDDATAFUNCTION_H

#include <deal.II/base/function.h>
#include <deal.II/base/function_lib.h>
#include <deal.II/base/table.h>
#include <deal.II/base/tensor.h>

#include <filesystem>
#include <string>
#include <vector>

namespace sapphirepp
{
  namespace Utils
  {
    using namespace dealii;

    /**
     * @brief Multicomponent @dealref{Function} created from data in a file.
     *
     * This function reads data from a file and uses the
     * @dealref{InterpolatedTensorProductGridData,classFunctions_1_1InterpolatedTensorProductGridData}
     * class
     * and @ref Tools::read_dat_to_tensor_product_grid_data()
     * to create a function with multiple components.
     *
     * @tparam dim Dimension of the function.
     */
    template <unsigned int dim>
    class GridDataFunction : public Function<dim>
    {
    public:
      /**
       * @brief Creates a time dependent function
       *        from a time series saved in a `.hst` file
       *        and the data saved in `.tab` files.
       *
       * Assumes that the data is saved as specified in the Athena++ output,
       * with the time series file `input_path/base_filename.hst`
       * and the data files `input_path/base_filename.block0.out1.XXXXX.tab`.
       *
       * @param input_path Folder where the files are located.
       * @param base_filename Base filename of the files.
       * @param n_components Number of components of function.
       * @param inital_time Initial time.
       * @param n_components_data Number of components saved in the file.
       *        Can be used if the number of components in the file
       *        is different from the number of components for the function.
       *        `n_components_data=0` means that it is equal to `n_components`.
       * @param delimiter Delimiter for the data files.
       * @param col_start_coordinates Column index with the first coordinate
       *        in the data files.
       * @param col_start_data Column index with the first data component
       *        in the data files.
       * @param last_coordinate_runs_fastest Last coordinate in the data files
       *        runs fastest?
       * @param athena_ordering Data files use Athena++ ordering of components?
       */
      GridDataFunction(const std::filesystem::path &input_path,
                       const std::string           &base_filename,
                       const unsigned int           n_components          = 1,
                       const double                 inital_time           = 0.0,
                       const unsigned int           n_components_data     = 0,
                       const std::string           &delimiter             = " ",
                       const unsigned int           col_start_coordinates = 0,
                       const unsigned int           col_start_data        = dim,
                       const bool last_coordinate_runs_fastest = false,
                       const bool athena_ordering              = false);



      /**
       * @brief Creates a time independent function loaded from one data file.
       *
       * @param filename Path to the data file.
       * @param n_components Number of components of function.
       * @param n_components_data Number of components saved in the file.
       *        Can be used if the number of components in the file
       *        is different from the number of components for the function.
       *        `n_components_data=0` means that it is equal to `n_components`.
       * @param delimiter Delimiter for the data files.
       * @param col_start_coordinates Column index with the first coordinate
       *        in the data files.
       * @param col_start_data Column index with the first data component
       *        in the data files.
       * @param last_coordinate_runs_fastest Last coordinate in the data files
       *        runs fastest?
       * @param athena_ordering Data files use Athena++ ordering of components?
       */
      GridDataFunction(const std::filesystem::path &filename,
                       const unsigned int           n_components          = 1,
                       const unsigned int           n_components_data     = 0,
                       const std::string           &delimiter             = " ",
                       const unsigned int           col_start_coordinates = 0,
                       const unsigned int           col_start_data        = dim,
                       const bool last_coordinate_runs_fastest = false,
                       const bool athena_ordering              = false);



      /**
       * @brief Value of the interpolated function at a point.
       *
       * @param p @dealref{Point}
       * @param component Component
       * @return Interpolated value at the point.
       */
      virtual double
      value(const Point<dim>  &p,
            const unsigned int component = 0) const override;



      /**
       * @brief Gradient of the interpolated function at a point.
       *
       * @param p @dealref{Point}
       * @param component Component
       * @return Tensor<1, dim> Gradient of interpolated function at the point/
       */
      virtual Tensor<1, dim>
      gradient(const Point<dim>  &p,
               const unsigned int component = 0) const override;



      /**
       * @brief Set time.
       *
       * Read in new data for the corresponding timestep if needed.
       *
       * @param new_time New time.
       */
      virtual void
      set_time(const double new_time) override;



      /**
       * @brief Load data from the file.
       *
       * Uses @ref Tools::read_dat_to_tensor_product_grid_data()
       * to read the data from a file.
       * The dataformat must be as specified there.
       *
       * @param filename Path to the data file.
       */
      void
      load_data_from_file(const std::filesystem::path &filename);



      /**
       * @brief Load time series information from a Athena++ history file.
       *
       * @param filename Path to the `.hst` file.
       * @return time_series Returns the time series.
       */
      static std::vector<double>
      read_hst_to_time_series(const std::filesystem::path &filename);



    private:
      /** Folder where the files are located. */
      const std::filesystem::path input_path;
      /** Base filename of the files. */
      const std::string base_filename;

      /** Delimiter for the data files. */
      const std::string delimiter;
      /** Column index with the first coordinate in the data files. */
      const unsigned int col_start_coordinates;
      /** Column index with the first coordinate in the data files. */
      const unsigned int col_start_data;
      /** Last coordinate in the data files runs fastest? */
      const bool last_coordinate_runs_fastest;
      /** Files use Athena++ ordering of components? */
      const bool athena_ordering;

      /** Vector of time stamps for the data. */
      const std::vector<double> time_series;
      /** Time index of currently loaded data. */
      unsigned int time_index;

      /**
       * @dealref{InterpolatedTensorProductGridData,classFunctions_1_1InterpolatedTensorProductGridData}
       * for each component.
       */
      std::vector<
        std::unique_ptr<dealii::Functions::InterpolatedUniformGridData<dim>>>
        grid_functions;
    };

  } // namespace Utils
} // namespace sapphirepp
#endif
