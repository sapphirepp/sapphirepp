#ifndef VFPEQUATION_VFPSOLVERCONTROL_H
#define VFPEQUATION_VFPSOLVERCONTROL_H

#include <deal.II/base/parameter_handler.h>

#include <ostream>
#include <string>

#include "compile-time-flags.h"

namespace VFPEquation {
class VFPSolverControl {
 public:
  // NOTE: Member variable needs be constexpr to be used as template arguments.
  // But it can only be constexpr if it is static ,i.e. if it is the same for
  // all class instances. If it was not static, it would be determined when
  // constructing an instance, which happens at runtime.

  // compile time settings
  static constexpr TermFlags terms =
      TermFlags::spatial_advection | TermFlags::magnetic;
  // Deactivating the spatial advection term is equivalent to assuming a
  // homogeneous distribution function (i.e. a distribution function which does
  // not depend on x,y z). In this program this equivalent to set dimension of
  // the configuration to zero.
  static constexpr int dim_configuration_space =
      ((terms & TermFlags::spatial_advection) == TermFlags::none ? 0 : 2);
  // Compute the total dimension and set momentum to true
  static constexpr bool momentum =
      ((terms & TermFlags::momentum) != TermFlags::none ? true : false);
  static constexpr int dim = dim_configuration_space + momentum;
  static_assert(dim <= 3, "The total dimension must be smaller or equal three");

  // Runtime settings
  // These settings are read from a parameter file
  int expansion_order;
  unsigned int num_refinements;
  unsigned int polynomial_degree;
  double time_step;
  double final_time;

  VFPSolverControl(const std::string& file_path);
  void print_settings(std::ostream& os) const;

 private:
  void declare_parameters();
  void parse_parameters();
  void get_parameters();
  const std::string parameter_file;
  dealii::ParameterHandler parameter_handler;
};

}  // namespace VFPEquation
#endif
