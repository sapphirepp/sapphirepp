#include "output-module.h"

#include <fstream>
#include <iostream>

#include "parameter-flags.h"
#include "sapphire-logstream.h"

template <int dim>
Sapphire::Utils::OutputModule<dim>::OutputModule(
  const Sapphire::Utils::ParameterParser &prm)
  : output_frequency(prm.out_output_frequency)
  , results_path(prm.out_results_path)
  , simulation_id(prm.out_simulation_id)
  , output_path(this->results_path / this->simulation_id)
  , base_file_name(prm.out_base_file_name)
  , n_digits_for_counter(prm.out_n_digits_for_counter)
  , format(prm.out_format)
  , mpi_communicator(MPI_COMM_WORLD)
{
  init(prm);
};

template <int dim>
void
Sapphire::Utils::OutputModule<dim>::init(
  const Sapphire::Utils::ParameterParser &prm) const
{
  LogStream::Prefix p("OutputModule", saplog);
  // create output directory
  if (Utilities::MPI::this_mpi_process(mpi_communicator) == 0)
    {
      saplog << "Create results folder " << output_path << std::endl;
      std::filesystem::create_directory(output_path);

      prm.write_parameters(output_path / "log");
    }


  if (format == OutputFormat::vtu)
    {
      int comm_size;
      MPI_Comm_size(mpi_communicator, &comm_size);
      Assert(comm_size == 1,
             ExcMessage("'vtu' format only works with a single "
                        "processor. Use 'pvtu' instead."));
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
              // TODO: Change this, append lines to xmdf file instead
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

// explicit instantiation
template class Sapphire::Utils::OutputModule<1>;
template class Sapphire::Utils::OutputModule<2>;
template class Sapphire::Utils::OutputModule<3>;
