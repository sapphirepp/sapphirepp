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

template <int dim> class OutputModule {
public:
  OutputModule(const std::string &output_path,
               const std::string &base_file_name)
      : output_path(output_path), base_file_name(base_file_name){};
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
        "Format", "vtu", dealii::Patterns::Selection("vtu|hdf5"),
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
    s = prm.get("Format");
    std::string format = s;
    unsigned int output_frequency = prm.get_integer("Output frequency");
    prm.leave_subsection(); // subsection Output

    (void)output_frequency;

    OutputModule output_module(results_path, simulation_id);
    output_module.init(prm);
    return output_module;
  };

  static void save_template_parameter(
      const ParameterHandler &prm,
      const std::string &template_filename = "parameter-template.prm") {
    prm.print_parameters(template_filename, ParameterHandler::PRM);
  };

  // TODO_BE: print log paramfile
  void init(const ParameterHandler &prm) const;
  std::filesystem::path get_filename(const unsigned int time_step_number) const;
  DataOutBase::VtkFlags get_vtk_flags() const;

private:
  // TODO_BE: change
  const std::filesystem::path output_path;
  const std::string base_file_name;
};

} // namespace Utils
} // namespace Sapphire
#endif