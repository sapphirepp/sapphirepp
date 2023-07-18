/**
 * @file output_module.h
 * @author Florian Schulze (florian.schulze@mpi-hd.mpg.de)
 * @brief Define OutputModule class
 * @version 0.1
 * @date 2023-07-17
 */

#ifndef SAPPHIREUTISL_OUTPUTMODULE_H
#define SAPPHIREUTISL_OUTPUTMODULE_H

#include <deal.II/base/mpi.h>
#include <deal.II/base/parameter_handler.h>

#include <deal.II/numerics/data_out.h>

#include <mpi.h>

#include <filesystem>

#include "parameter-flags.h"
#include "parameter-parser.h"

namespace Sapphire
{
  namespace Utils
  {
    using namespace dealii;

    template <int dim>
    class OutputModule
    {
    public:
      OutputModule(const Sapphire::Utils::ParameterParser &prm);

      void
      init(const ParameterHandler &prm) const;
      void
      write_results(DataOut<dim>      &data_out,
                    const unsigned int time_step_number);

      const unsigned int output_frequency;

    private:
      const std::filesystem::path results_path;
      const std::string           simulation_id;
      const std::filesystem::path output_path;
      const std::string           base_file_name;
      const unsigned int          n_digits_for_counter;
      const OutputFormat          format;

      const MPI_Comm         mpi_communicator;
      std::vector<XDMFEntry> xdmf_entries;
    };

  } // namespace Utils
} // namespace Sapphire
#endif
