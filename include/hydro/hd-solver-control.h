/**
 * @file hd-solver-control.h
 * @author Florian Schulze (florian.schulze@mpi-hd.mpg.de)
 * @brief Implement numerical functions to solve the hydrodynamics equations
 * @version 0.1
 * @date 2023-07-12
 */

#ifndef HYDROSOLVER_HDSOLVERCONTROL_H
#define HYDROSOLVER_HDSOLVERCONTROL_H

#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_control.h>
#include <deal.II/lac/vector.h>

#include "parameter-parser.h"

namespace Sapphire
{
  namespace Hydro
  {
    using namespace dealii;

    class HDSolverControl
    {
    public:
      HDSolverControl(const Sapphire::Utils::ParameterParser &prm)
        : scheme(prm.hdsolver_scheme)
        , fe_degree(prm.hdsolver_fe_degree)
        , time_step(prm.hdsolver_time_step)
        , end_time(prm.hdsolver_end_time)
        , refinement_level(prm.hdsolver_refinement_level)
        , max_iterations(prm.hdsolver_max_iterations)
        , tolerance(prm.hdsolver_tolerance){};

      const TimeSteppingScheme scheme;

      const unsigned int fe_degree;
      const double       time_step;
      const double       end_time;
      const unsigned int refinement_level;

      const unsigned int max_iterations;
      const double       tolerance;
    };

  } // namespace Hydro
} // namespace Sapphire
#endif