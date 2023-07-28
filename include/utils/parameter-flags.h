#ifndef UTILS_PARAMETERFLAGS_H
#define UTILS_PARAMETERFLAGS_H

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
    enum class TimeSteppingMethod
    {
      crank_nicolson,
      forward_euler,
      backward_euler,
      erk4,
      lserk
    };
  } // namespace Utils


  namespace Hydro
  {
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
    enum class SlopeLimiterType
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
} // namespace Sapphire
#endif
