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
 * @file output-parameters.h
 * @author Nils Schween (nils.schween@mpi-hd.mpg.de)
 * @author Florian Schulze (florian.schulze@mpi-hd.mpg.de)
 * @brief Define @ref sapphirepp::Utils::OutputParameters
 */

#ifndef UTILS_OUTPUTPARAMETERS_H
#define UTILS_OUTPUTPARAMETERS_H

#include <deal.II/base/parameter_handler.h>

#include <deal.II/grid/tria.h>

#include <deal.II/numerics/data_out.h>

#include <mpi.h>

#include <filesystem>
#include <string>
#include <vector>

namespace sapphirepp
{
  namespace Utils
  {
    using namespace dealii;

    /** @brief Enum for output formats */
    enum class OutputFormat
    {
      /** VTU output, recommended for short runs */
      vtu,
      /** VTU output with `.pvtu` record for parallel runs */
      pvtu,
      /** HDF5 output with XDMF record, recommended for large simulations */
      hdf5
    };

    /**
     * @brief Parameters class to handle output parameters
     *
     * This function handles the output of @sapphire. It defines all the output
     * related parameters and provides functions to write the results.
     */
    class OutputParameters
    {
    public:
      /** Only put out every n-th step */
      unsigned int output_frequency = 1;
      /** Output folder = results_path + simulation_id */
      std::filesystem::path output_path;
      /** Path to the results directory */
      std::filesystem::path results_path;
      /** Name of the simulation run, i.e. subfolder for the simulation */
      std::string simulation_id = "01";
      /** File name for the solution */
      std::string base_file_name = "solution";
      /** Number of digits for time steps */
      unsigned int n_digits_for_counter = 4;
      /** Output format */
      OutputFormat format;
      /** Append user defined input functions to output as debug information? */
      bool debug_input_functions = false;

      /** @brief Constructor */
      OutputParameters();

      /**
       * @brief Delcare parameters in parameter file
       *
       * @param prm @dealref{ParameterHandler}
       */
      void
      declare_parameters(ParameterHandler &prm);

      /**
       * @brief Parse parameters from parameter file
       *
       * @param prm @dealref{ParameterHandler}
       */
      void
      parse_parameters(ParameterHandler &prm);

      /**
       * @brief Write results to file
       *
       * @tparam dim Dimension of the solution
       * @param data_out @dealref{DataOut} class with the solution.
       *        It must already contain all the data that should be written.
       *        This function only handles the last step of writing to the file.
       * @param time_step_number Number of the time step
       * @param cur_time Simulation time associated with this time step
       * @param filename Overwrite for the @ref base_file_name
       * @note This function should be a const member. The problem is, that in
       *       case of HDF5 output, the xdmf_entries are changed.
       */
      template <unsigned int dim>
      void
      write_results(DataOut<dim>      &data_out,
                    const unsigned int time_step_number = 0,
                    const double cur_time = std::numeric_limits<double>::min(),
                    const std::string &filename = "");

      /**
       * @brief Write a grid to a file in ucd format
       *
       * @tparam dim Dimension of the grid
       * @param triangulation Grid
       * @param filename Name of the file
       */
      template <unsigned int dim>
      void
      write_grid(const Triangulation<dim> &triangulation,
                 const std::string        &filename = "grid.ucd") const;


    private:
      /**
       * @brief Callback function for parsing parameters
       *
       * Creates output folder and writes parameter log.
       *
       * @param prm @dealref{ParameterHandler}
       */
      void
      parse_parameters_callback(ParameterHandler &prm) const;

      /** mpi_communicator */
      const MPI_Comm mpi_communicator;
      /** XDMF entries for the HDF5 record */
      std::vector<XDMFEntry> xdmf_entries;
    };

  } // namespace Utils
} // namespace sapphirepp
#endif
