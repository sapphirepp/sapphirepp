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

      /**HDSolverControl parameter*/
      Sapphire::Hydro::TimeSteppingScheme    hdsolver_scheme;
      Sapphire::Hydro::FluxType              hdsolver_flux_type;
      Sapphire::Hydro::SlopeLimiter          hdsolver_limiter;
      Sapphire::Hydro::SlopeLimiterCriterion hdsolver_limiter_criterion;

      unsigned int hdsolver_fe_degree;
      double       hdsolver_time_step;
      double       hdsolver_end_time;
      unsigned int hdsolver_refinement_level;

      unsigned int hdsolver_max_iterations;
      double       hdsolver_tolerance;

      // VFP Control Parameter
      // Spherical harmonic expansion
      int expansion_order;
      // Mesh
      std::string p1;
      std::string p2;
      std::string n_cells;
      std::string periodicity;
      // Finite element
      unsigned int polynomial_degree;
      // Time stepping
      std::string time_stepping_method;
      double      theta;
      double      time_step;
      double      final_time;

    private:
      void
      declare_parameters();
      void
      parse_parameters();

      ParameterHandler prm;
      const std::string prm_file_name;
    };

  } // namespace Utils
} // namespace Sapphire
#endif
