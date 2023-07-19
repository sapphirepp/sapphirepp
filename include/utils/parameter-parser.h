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

      void
      write_parameters(const std::string &filename) const;
      void
      write_template_parameters(const std::string &filename);

      /**Output parameter*/
      std::string                   out_results_path;
      std::string                   out_simulation_id;
      std::string                   out_output_path;
      std::string                   out_base_file_name;
      unsigned int                  out_n_digits_for_counter;
      Sapphire::Utils::OutputFormat out_format;
      unsigned int                  out_output_frequency;

      // VFP Control Parameter
      // Spherical harmonic expansion
      int vfp_expansion_order;
      // Mesh
      std::string vfp_p1;
      std::string vfp_p2;
      std::string vfp_n_cells;
      std::string vfp_periodicity;
      // Finite element
      unsigned int vfp_polynomial_degree;
      // Time stepping
      std::string vfp_time_stepping_method;
      double      vfp_theta;
      double      vfp_time_step;
      double      vfp_final_time;

    private:
      void
      declare_parameters();
      void
      parse_parameters();

      ParameterHandler  prm;
      const std::string prm_file_name;
    };

  } // namespace Utils
} // namespace Sapphire
#endif
