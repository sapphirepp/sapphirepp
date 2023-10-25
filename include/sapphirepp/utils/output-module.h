/**
 * @file output_module.h
 * @author Florian Schulze (florian.schulze@mpi-hd.mpg.de)
 * @brief Define OutputModule class
 * @version 0.1
 * @date 2023-07-17
 */

#ifndef UTILS_OUTPUTMODULE_H
#define UTILS_OUTPUTMODULE_H

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
    class OutputModule
    {
    public:
      OutputModule();

      void
      declare_parameters(ParameterHandler &prm);

      void
      parse_parameters(ParameterHandler &prm);

      void
      init(ParameterHandler &prm) const;

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
      const MPI_Comm         mpi_communicator;
      std::vector<XDMFEntry> xdmf_entries;
    };

  } // namespace Utils
} // namespace Sapphire
#endif
