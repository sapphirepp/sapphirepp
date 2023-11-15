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
 * @file output_parameters.h
 * @author Florian Schulze (florian.schulze@mpi-hd.mpg.de)
 * @brief Define OutputParameters class
 * @version 0.1
 * @date 2023-07-17
 */

#ifndef UTILS_OUTPUTPARAMETERS_H
#define UTILS_OUTPUTPARAMETERS_H

#include <deal.II/base/mpi.h>
#include <deal.II/base/parameter_handler.h>

#include <deal.II/numerics/data_out.h>

#include <mpi.h>

#include <filesystem>

namespace Sapphire
{
  namespace Utils
  {
    using namespace dealii;

    enum class OutputFormat
    {
      vtu,
      pvtu,
      hdf5
    };

    template <int dim>
    class OutputParameters
    {
    public:
      OutputParameters();

      void
      declare_parameters(ParameterHandler &prm);

      void
      parse_parameters(ParameterHandler &prm);

      void
      write_results(DataOut<dim>      &data_out,
                    const unsigned int time_step_number);

      void
      write_grid(const Triangulation<dim> &triangulation,
                 const std::string        &filename = "grid.vtu") const;

      unsigned int          output_frequency;
      std::filesystem::path results_path;
      std::string           simulation_id;
      std::filesystem::path output_path;
      std::string           base_file_name;
      unsigned int          n_digits_for_counter;
      OutputFormat          format;

    private:
      void
      parse_parameters_callback(ParameterHandler &prm) const;

      const MPI_Comm         mpi_communicator;
      std::vector<XDMFEntry> xdmf_entries;
    };

  } // namespace Utils
} // namespace Sapphire
#endif
