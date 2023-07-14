/**
 * @file output_module.h
 * @author Florian Schulze (florian.schulze@mpi-hd.mpg.de)
 * @brief Define OutputModule class
 * @version 0.1
 * @date 2023-07-12
 */

#ifndef SAPPHIREUTISL_OUTPUTMODULE_H
#define SAPPHIREUTISL_OUTPUTMODULE_H

#include <deal.II/base/mpi.h>
#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/smartpointer.h>
#include <deal.II/lac/vector.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/vector_tools.h>

#include <mpi.h>

#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>

#include <iostream>

#define DEBUG_LEVEL 1

#define DEBUG_PRINT(ostream, level, message)                                   \
  if (level <= DEBUG_LEVEL) {                                                  \
    for (int i = 0; i < level; ++i) {                                          \
      ostream << "  ";                                                         \
    }                                                                          \
    ostream << message << std::endl;                                           \
  }

namespace Sapphire {
namespace Utils {
using namespace dealii;

enum class OutputFormat { vtu, pvtu, hdf5 };

template <int dim> class OutputModule {
public:
  OutputModule(const std::string &results_path,
               const std::string &simulation_id, const OutputFormat &format,
               const unsigned int &output_frequency)
      : results_path(results_path), simulation_id(simulation_id),
        output_path(this->results_path / this->simulation_id), format(format),
        output_frequency(output_frequency){};
  OutputModule(const OutputModule<dim> &output_module) = default;
  ~OutputModule() = default;

  static void declare_parameters(ParameterHandler &prm) {
    prm.enter_subsection("Output");
    prm.declare_entry(
        "Results folder", "./results", dealii::Patterns::Anything(),
        "Path to the folder in which the simulation results will be stored. "
        "Without a trailing slash.");
    prm.declare_entry("Simulation identifier", "001",
                      dealii::Patterns::Anything(),
                      "Name of the simulation run. It will be "
                      "used to create a subdirectory "
                      "in the results folder.");
    prm.declare_entry(
        "Format", "vtu", dealii::Patterns::Selection("vtu|pvtu|hdf5"),
        "The format in which the simulation output will be stored.");
    prm.declare_entry("Output frequency", "1", dealii::Patterns::Integer(0),
                      "The frequence at which output files will "
                      "be written. (In units of time steps)");
    prm.leave_subsection(); // subsection Output
  };

  static OutputModule parse_parameters(ParameterHandler &prm) {
    std::string s;

    prm.enter_subsection("Output");
    // std::string results_path = "";
    std::string results_path = prm.get("Results folder");
    std::string simulation_id = prm.get("Simulation identifier");
    OutputFormat format;
    s = prm.get("Format");
    if (s == "vtu")
      format = OutputFormat::vtu;
    else if (s == "pvtu")
      format = OutputFormat::pvtu;
    else if (s == "hdf5")
      format = OutputFormat::hdf5;
    else
      AssertThrow(false, ExcNotImplemented());
    unsigned int output_frequency = prm.get_integer("Output frequency");
    prm.leave_subsection(); // subsection Output

    OutputModule output_module(results_path, simulation_id, format,
                               output_frequency);
    output_module.init(prm);
    return output_module;
  };

  void init(const ParameterHandler &prm) const;
  void write_results(DataOut<dim> &data_out,
                     const unsigned int time_step_number,
                     const MPI_Comm &mpi_communicator,
                     std::vector<XDMFEntry> &xdmf_entries) const;
  void write_results(DataOut<dim> &data_out,
                     const unsigned int time_step_number,
                     const MPI_Comm &mpi_communicator = MPI_COMM_WORLD) const {
    Assert(format != OutputFormat::hdf5, ExcNotImplemented());
    std::vector<XDMFEntry> xdmf_entries(0);
    write_results(data_out, time_step_number, mpi_communicator, xdmf_entries);
  };

private:
  const std::filesystem::path results_path;
  const std::string simulation_id;
  const std::string base_file_name = "solution";
  const unsigned int n_digits_for_counter = 4;
  const std::filesystem::path output_path;
  const OutputFormat format;

public:
  const unsigned int output_frequency;
};

} // namespace Utils
} // namespace Sapphire
#endif