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

#include "output-module.h"

#include <deal.II/grid/grid_out.h>

#include <fstream>
#include <iostream>

#include "sapphirepp-logstream.h"

template <int dim>
Sapphire::Utils::OutputModule<dim>::OutputModule()
  : mpi_communicator(MPI_COMM_WORLD){};

template <int dim>
void
Sapphire::Utils::OutputModule<dim>::declare_parameters(ParameterHandler &prm)
{
  LogStream::Prefix p("OutputModule", saplog);
  saplog << "Declaring parameters" << std::endl;

  prm.enter_subsection("Output");

  prm.declare_entry(
    "Results folder",
    "./results",
    Patterns::Anything(),
    "Path to the folder in which the simulation results will be stored. "
    "Without a trailing slash.");
  prm.declare_entry("Simulation identifier",
                    "",
                    Patterns::Anything(),
                    "Name of the simulation run. It will be "
                    "used to create a subdirectory "
                    "in the results folder.");
  prm.declare_entry("Base file name",
                    "solution",
                    Patterns::Anything(),
                    "The base file name for the output files.");
  prm.declare_entry("Number of digits for counter",
                    "4",
                    Patterns::Integer(0),
                    "The number of digits used for the counter in the "
                    "output file names.");
  prm.declare_entry(
    "Format",
    "pvtu",
    Patterns::Selection("vtu|pvtu|hdf5"),
    "The format in which the simulation output will be stored.");
  prm.declare_entry("Output frequency",
                    "1",
                    Patterns::Integer(0),
                    "The frequency at which output files will "
                    "be written. (In units of time steps)");

  prm.leave_subsection();
}

template <int dim>
void
Sapphire::Utils::OutputModule<dim>::parse_parameters(ParameterHandler &prm)
{
  LogStream::Prefix p("OutputModule", saplog);
  saplog << "Parsing parameters" << std::endl;
  std::string s;
  prm.enter_subsection("Output");

  results_path         = prm.get("Results folder");
  simulation_id        = prm.get("Simulation identifier");
  output_path          = this->results_path / this->simulation_id;
  base_file_name       = prm.get("Base file name");
  n_digits_for_counter = prm.get_integer("Number of digits for counter");

  s = prm.get("Format");
  if (s == "vtu")
    format = Sapphire::Utils::OutputFormat::vtu;
  else if (s == "pvtu")
    format = Sapphire::Utils::OutputFormat::pvtu;
  else if (s == "hdf5")
    format = Sapphire::Utils::OutputFormat::hdf5;
  else
    AssertThrow(false, ExcNotImplemented());

  output_frequency = prm.get_integer("Output frequency");

  prm.leave_subsection();

  parse_parameters_callback(prm);
}

template <int dim>
void
Sapphire::Utils::OutputModule<dim>::parse_parameters_callback(
  ParameterHandler &prm) const
{
  // create output directory
  if (Utilities::MPI::this_mpi_process(mpi_communicator) == 0)
    {
      saplog << "Create results folder " << output_path << std::endl;
      std::filesystem::create_directories(output_path);

      saplog << "Log parameters" << std::endl;
      prm.log_parameters(saplog);
      prm.print_parameters(output_path / "log.prm", ParameterHandler::ShortPRM);
    }


  if (format == OutputFormat::vtu)
    {
      int comm_size;
      MPI_Comm_size(mpi_communicator, &comm_size);
      AssertThrow(comm_size == 1,
                  ExcMessage("'vtu' format only works with a single "
                             "processor. Use 'pvtu' instead."));
    }
  else if (format == OutputFormat::hdf5)
    {
      AssertThrow(
        dim > 1,
        ExcMessage(
          "HDF5 only supports datasets that live in 2 or 3 dimensions. "
          "Use 'vtu' or 'pvtu' instead."));
    }
}

template <int dim>
void
Sapphire::Utils::OutputModule<dim>::write_results(
  DataOut<dim>      &data_out,
  const unsigned int time_step_number)
{
  if (time_step_number % output_frequency == 0)
    {
      LogStream::Prefix p("OutputModule", saplog);
      saplog << "Writing results at time_step " << time_step_number
             << std::endl;
      const unsigned int counter = time_step_number / output_frequency;
      switch (format)
        {
          case OutputFormat::vtu:
            {
              DataOutBase::VtkFlags vtk_flags;
              vtk_flags.compression_level =
                DataOutBase::VtkFlags::ZlibCompressionLevel::best_speed;
              data_out.set_flags(vtk_flags);
              const std::string filename_vtk =
                base_file_name + "_" +
                Utilities::int_to_string(counter, n_digits_for_counter) +
                ".vtu";
              std::ofstream output(output_path / filename_vtk);
              data_out.write_vtu(output);
              break;
            }
          case OutputFormat::pvtu:
            {
              data_out.write_vtu_with_pvtu_record(output_path / "",
                                                  base_file_name,
                                                  counter,
                                                  mpi_communicator,
                                                  n_digits_for_counter);
              break;
            }
          case OutputFormat::hdf5:
            {
              // I follow this pull request:
              // https://github.com/dealii/dealii/pull/14958
              const std::string filename_h5 =
                base_file_name +
                Utilities::int_to_string(counter, n_digits_for_counter) + ".h5";
              const std::string filename_mesh = "mesh.h5";
              const std::string xdmf_file     = base_file_name + ".xdmf";
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
              const bool write_mesh_hdf5 = counter == 0 ? true : false;

              data_out.write_hdf5_parallel(data_filter,
                                           write_mesh_hdf5,
                                           output_path / filename_mesh,
                                           output_path / filename_h5,
                                           mpi_communicator);

              XDMFEntry new_xdmf_entry = data_out.create_xdmf_entry(
                data_filter,
                filename_mesh,
                base_file_name +
                  Utilities::int_to_string(counter, n_digits_for_counter) +
                  ".h5",
                counter,
                mpi_communicator);
              // TODO: Change this, append lines to xdmf file instead
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
            {
              Assert(false,
                     ExcNotImplemented("Output format not implemented."));
              break;
            }
        }
    }
}

template <int dim>
void
Sapphire::Utils::OutputModule<dim>::write_grid(
  const Triangulation<dim> &triangulation,
  const std::string        &filename) const
{
  LogStream::Prefix p("OutputModule", saplog);
  saplog << "Write grid " << filename << std::endl;

  GridOutFlags::Ucd ucd_flags;
  ucd_flags.write_faces = true;
  ucd_flags.write_lines = true;

  GridOut grid_out;
  grid_out.set_flags(ucd_flags);
  std::ofstream output(output_path / filename);
  grid_out.write_ucd(triangulation, output);
}

// explicit instantiation
template class Sapphire::Utils::OutputModule<1>;
template class Sapphire::Utils::OutputModule<2>;
template class Sapphire::Utils::OutputModule<3>;
