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
 * @file output-parameters.cpp
 * @author Nils Schween (nils.schween@mpi-hd.mpg.de)
 * @author Florian Schulze (florian.schulze@mpi-hd.mpg.de)
 * @brief Implement @ref sapphirepp::Utils::OutputParameters
 */

#include "output-parameters.h"

#include <deal.II/grid/grid_out.h>

#include <fstream>
#include <iostream>

#include "sapphirepp-logstream.h"



sapphirepp::Utils::OutputParameters::OutputParameters()
  : mpi_communicator(MPI_COMM_WORLD)
{}



void
sapphirepp::Utils::OutputParameters::declare_parameters(ParameterHandler &prm)
{
  LogStream::Prefix pre1("Startup", saplog);
  LogStream::Prefix pre2("OutputParameters", saplog);
  saplog << "Declaring parameters" << std::endl;

  prm.enter_subsection("Output");

  prm.declare_entry(
    "Results folder",
    "./results",
    Patterns::DirectoryName(),
    "Path to the folder in which the simulation results will be stored. "
    "Without a trailing slash.");
  prm.add_parameter("Simulation identifier",
                    simulation_id,
                    "Name of the simulation run. It will be used to create a "
                    "subdirectory in the results folder.");
  prm.add_parameter("Base file name",
                    base_file_name,
                    "The base file name for the output files.");
  prm.add_parameter("Number of digits for counter",
                    n_digits_for_counter,
                    "The number of digits used for the counter in the "
                    "output file names.");
  prm.declare_entry("Format",
                    "pvtu",
                    Patterns::Selection("vtu|pvtu|hdf5"),
                    "The format in which the simulation "
                    "output will be stored.");
  prm.add_parameter("Output frequency",
                    output_frequency,
                    "The frequency at which output files will be written. "
                    "(In units of time steps)");
  prm.add_parameter("Debug input functions",
                    debug_input_functions,
                    "Append the user defined input_functions "
                    "like magnetic field and scattering frequency "
                    "to the output as debug information.");

  prm.leave_subsection();
}



void
sapphirepp::Utils::OutputParameters::parse_parameters(ParameterHandler &prm)
{
  LogStream::Prefix pre1("Startup", saplog);
  LogStream::Prefix pre2("OutputParameters", saplog);
  saplog << "Parsing parameters" << std::endl;
  std::string s;
  prm.enter_subsection("Output");

  results_path = prm.get("Results folder");
  output_path  = this->results_path / this->simulation_id;

  s = prm.get("Format");
  if (s == "vtu")
    format = sapphirepp::Utils::OutputFormat::vtu;
  else if (s == "pvtu")
    format = sapphirepp::Utils::OutputFormat::pvtu;
  else if (s == "hdf5")
    format = sapphirepp::Utils::OutputFormat::hdf5;
  else
    Assert(false, ExcNotImplemented());

  prm.leave_subsection();

  parse_parameters_callback(prm);
}



void
sapphirepp::Utils::OutputParameters::parse_parameters_callback(
  ParameterHandler &prm) const
{
  // create output directory
  if (Utilities::MPI::this_mpi_process(mpi_communicator) == 0)
    {
      saplog << "Create results folder " << output_path << std::endl;
      std::filesystem::create_directories(output_path);

      saplog << "Log parameters" << std::endl;
      prm.log_parameters(saplog);
      prm.print_parameters(output_path / static_cast<std::string>("log.prm"),
                           ParameterHandler::ShortPRM);
    }
}



template <unsigned int dim>
void
sapphirepp::Utils::OutputParameters::write_results(
  DataOut<dim>      &data_out,
  const unsigned int time_step_number,
  const double       cur_time,
  const std::string &filename)
{
  LogStream::Prefix p("OutputParameters", saplog);
  saplog << "Writing results at time_step " << time_step_number << std::endl;

  const std::string tmp_base_file_name =
    filename.empty() ? base_file_name : filename;

  switch (format)
    {
      case OutputFormat::vtu:
        {
          DataOutBase::VtkFlags vtk_flags(cur_time, time_step_number);
          data_out.set_flags(vtk_flags);
          const std::string filename_vtk =
            tmp_base_file_name + "_" +
            Utilities::int_to_string(time_step_number, n_digits_for_counter) +
            ".vtu";
          data_out.write_vtu_in_parallel(output_path / filename_vtk,
                                         mpi_communicator);
          break;
        }
      case OutputFormat::pvtu:
        {
          DataOutBase::VtkFlags vtk_flags(cur_time, time_step_number);
          data_out.set_flags(vtk_flags);
          data_out.write_vtu_with_pvtu_record(output_path / "",
                                              tmp_base_file_name,
                                              time_step_number,
                                              mpi_communicator,
                                              n_digits_for_counter);
          break;
        }
      case OutputFormat::hdf5:
        {
          Assert(
            dim > 1,
            ExcMessage(
              "HDF5 only supports datasets that live in 2 or 3 dimensions."));

          // I follow this pull request:
          // https://github.com/dealii/dealii/pull/14958
          const std::string filename_h5 =
            tmp_base_file_name +
            Utilities::int_to_string(time_step_number, n_digits_for_counter) +
            ".h5";
          const std::string filename_mesh = "mesh.h5";
          const std::string xdmf_file     = tmp_base_file_name + ".xdmf";
          // https://dealii.org/developer/doxygen/deal.II/structDataOutBase_1_1DataOutFilterFlags.html
          // Whether or not to filter out duplicate vertices and associated
          // values. Setting this value to true will drastically reduce the
          // output data size but will result in an output file that does
          // not faithfully represent the actual data if the data
          // corresponds to discontinuous fields. In particular, along
          // subdomain boundaries the data will still be discontinuous,
          // while it will look like a continuous field inside of the
          // subdomain. NOTE: I suspect that discontinuous elements produce
          // discontinuous fields, but I do not know.
          const bool filter_duplicate_vertices{false};
          const bool xdmf_hdf5_output{true};

          DataOutBase::DataOutFilterFlags flags(filter_duplicate_vertices,
                                                xdmf_hdf5_output);
          DataOutBase::DataOutFilter      data_filter(flags);
          // / Filter the data and store it in data_filter
          data_out.write_filtered_data(data_filter);
          // set the HDF5 compression level Future versions of dealii will
          // allow to set the compression level: compare pull request.
          // dealii::DataOutBase::Hdf5Flags hdf5Flags;
          // hdf5Flags.compression_level =
          // dealii::DataOutBase::CompressionLevel::best_compression;
          // dataOut.set_flags(dealii::DataOutBase::Hdf5Flags(hdf5Flags));
          const bool write_mesh_hdf5 = time_step_number == 0 ? true : false;

          data_out.write_hdf5_parallel(data_filter,
                                       write_mesh_hdf5,
                                       output_path / filename_mesh,
                                       output_path / filename_h5,
                                       mpi_communicator);

          XDMFEntry new_xdmf_entry = data_out.create_xdmf_entry(
            data_filter,
            filename_mesh,
            tmp_base_file_name +
              Utilities::int_to_string(time_step_number, n_digits_for_counter) +
              ".h5",
            cur_time,
            mpi_communicator);
          /** @todo Append lines to xdmf file instead of rewriting the file */
          xdmf_entries.push_back(new_xdmf_entry);
          // NOTE: For now I a writing the xdmf file in every time step.
          // That is unnecessary. There is missing a function add entry to
          // xdmf_file
          data_out.write_xdmf_file(xdmf_entries,
                                   output_path / xdmf_file,
                                   mpi_communicator);
          break;
        }
      default:
        Assert(false, ExcNotImplemented("Output format not implemented."));
    }
}

// Explicit instantiation
template void
sapphirepp::Utils::OutputParameters::write_results<1>(DataOut<1> &,
                                                      const unsigned int,
                                                      const double,
                                                      const std::string &);
template void
sapphirepp::Utils::OutputParameters::write_results<2>(DataOut<2> &,
                                                      const unsigned int,
                                                      const double,
                                                      const std::string &);
template void
sapphirepp::Utils::OutputParameters::write_results<3>(DataOut<3> &,
                                                      const unsigned int,
                                                      const double,
                                                      const std::string &);



template <unsigned int dim>
void
sapphirepp::Utils::OutputParameters::write_grid(
  const Triangulation<dim> &triangulation,
  const std::string        &filename) const
{
  LogStream::Prefix p("OutputParameters", saplog);
  saplog << "Write grid " << filename << std::endl;

  GridOutFlags::Ucd ucd_flags;
  ucd_flags.write_faces = true;
  ucd_flags.write_lines = true;

  GridOut grid_out;
  grid_out.set_flags(ucd_flags);
  std::ofstream output(output_path / filename);
  grid_out.write_ucd(triangulation, output);
}

// Explicit
template void
sapphirepp::Utils::OutputParameters::write_grid<1>(const Triangulation<1> &,
                                                   const std::string &) const;
template void
sapphirepp::Utils::OutputParameters::write_grid<2>(const Triangulation<2> &,
                                                   const std::string &) const;
template void
sapphirepp::Utils::OutputParameters::write_grid<3>(const Triangulation<3> &,
                                                   const std::string &) const;
