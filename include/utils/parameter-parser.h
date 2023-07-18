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

namespace Sapphire
{

  namespace Utils
  {
    enum class OutputFormat
    {
      vtu,
      pvtu,
      hdf5
    };
  } // namespace Utils

  namespace Hydro
  { // TODO: Should this be in Hydro namespace?
    enum class TimeSteppingScheme
    {
      ForwardEuler,
      ExplicitRK
    };
    enum class FluxType
    {
      Central,
      Upwind,
      LaxFriedrichs
    };
    enum class SlopeLimiter
    {
      NoLimiter,
      LinearReconstruction,
      MinMod,
      MUSCL
    };
    enum class SlopeLimiterCriterion
    {
      Never,
      Always,
      GerneralizedSlopeLimiter
    };
  } // namespace Hydro

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