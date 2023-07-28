#include "numerics.h"

#include <iostream>

Sapphire::Hydro::HDSolverControl::HDSolverControl(
  const Sapphire::Utils::ParameterParser &prm)
  : scheme(prm.hdsolver_scheme)
  , limiter_criterion(prm.hdsolver_limiter_criterion)
  , fe_degree(prm.hdsolver_fe_degree)
  , time_step(prm.hdsolver_time_step)
  , end_time(prm.hdsolver_end_time)
  , refinement_level(prm.hdsolver_refinement_level)
  , max_iterations(prm.hdsolver_max_iterations)
  , tolerance(prm.hdsolver_tolerance)
{}
