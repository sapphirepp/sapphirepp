#include "output_module.h"

template <int dim> void Sapphire::Utils::OutputModule<dim>::init() const {
  // create output directory
  if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0) {
    std::filesystem::create_directory(output_path);
  }
}

template <int dim>
std::filesystem::path Sapphire::Utils::OutputModule<dim>::get_filename(
    const unsigned int time_step_number) const {

  return output_path / (base_file_name +
                        Utilities::int_to_string(time_step_number, 3) + ".vtu");
}

template <int dim>
dealii::DataOutBase::VtkFlags
Sapphire::Utils::OutputModule<dim>::get_vtk_flags() const {
  DataOutBase::VtkFlags vtk_flags;
  vtk_flags.compression_level =
      DataOutBase::VtkFlags::ZlibCompressionLevel::best_speed;
  return vtk_flags;
}

// explicit instantiation
template class Sapphire::Utils::OutputModule<1>;
template class Sapphire::Utils::OutputModule<2>;
template class Sapphire::Utils::OutputModule<3>;
