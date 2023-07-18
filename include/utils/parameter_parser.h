/**
 * @file parameter_parser.h
 * @author Florian Schulze (florian.schulze@mpi-hd.mpg.de)
 * @brief Implement ParameterParser class
 * @version 0.1
 * @date 2023-07-17
 */

#ifndef SAPPHIREUTISL_PARAMETERPARSER_H
#define SAPPHIREUTISL_PARAMETERPARSER_H

#include <deal.II/base/parameter_handler.h>

#include <string>

#include "parameter-flags.h"

namespace Sapphire
{
  namespace Utils
  {
    using namespace dealii;

    class ParameterParser
    {
    public:
      ParameterParser(const std::string &prm_file_name);

      const ParameterHandler &
      get_parameter_handler() const
      {
        return prm;
      }

      /**Output parameter*/
      std::string                   out_results_path;
      std::string                   out_simulation_id;
      std::string                   out_output_path;
      std::string                   out_base_file_name;
      unsigned int                  out_n_digits_for_counter;
      Sapphire::Utils::OutputFormat out_format;
      unsigned int                  out_output_frequency;

    private:
      void
      declare_parameters();
      void
      parse_parameters();

      ParameterHandler prm;
    };

  } // namespace Utils
} // namespace Sapphire
#endif