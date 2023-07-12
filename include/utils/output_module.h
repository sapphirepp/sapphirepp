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
  OutputModule(const std::string &output_path = "results",
               const std::string &base_file_name = "solution-")
      : output_path(output_path), base_file_name(base_file_name){};
  OutputModule(const OutputModule<dim> &output_module) = default;
  ~OutputModule() = default;

  void init() const;
  std::filesystem::path get_filename(const unsigned int time_step_number) const;
  DataOutBase::VtkFlags get_vtk_flags() const;

private:
  const std::filesystem::path output_path;
  const std::string base_file_name;
};

} // namespace Utils
} // namespace Sapphire
#endif