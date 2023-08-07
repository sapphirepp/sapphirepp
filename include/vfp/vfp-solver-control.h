#ifndef VFP_VFPSOLVERCONTROL_H
#define VFP_VFPSOLVERCONTROL_H

#include <deal.II/base/parameter_handler.h>

#include <ostream>
#include <string>
#include <vector>

#include "config.h"
#include "parameter-flags.h"
#include "parameter-parser.h"

namespace Sapphire
{
  namespace VFP
  {
    class VFPSolverControl
    {
    public:
      VFPSolverControl(const Sapphire::Utils::ParameterParser &prm);

      void
      print_settings(std::ostream &os) const;

      // NOTE: A member variable needs be constexpr to be used as template
      // arguments. But it can only be constexpr if it is static ,i.e. if it is
      // the same for all class instances. If it was not static, it would be
      // determined when constructing an instance, which happens at runtime.

      // compile time settings
      static constexpr TermFlags terms =
        TermFlags::spatial_advection | TermFlags::source;

      // This variabale controls if p is linear or logarithmic
      static constexpr bool logarithmic_p = true;

      // Deactivating the spatial advection term is equivalent to assuming a
      // homogeneous distribution function (i.e. a distribution function which
      // does not depend on x,y z). In this program this is equivalent to set
      // dimension of the configuration to zero.
      static constexpr int dim_configuration_space = 2;

      // If the background velocity field and the the magnetic field do not
      // depend on time, the time stepping methods can be accelerated a lot: In
      // this case it is not necessary to reassamble the spatial discretisation
      // matrix in every stage of the Runge-Kutta method. Actually it only has
      // to be assembled once at time zero.
      static constexpr bool time_dependent_fields = false;

      // If the source term depends on time, it needs to be projected onto the
      // FEM space in the time stepping methods.
      static constexpr bool time_dependent_source = true;

      // The following static_assert uses an exlusive or (xor),
      // represented in C/C++ as "!=" for expressions of boolean type: Either
      // the spatial advection term is included in the equation or the dimension
      // of the configuration space is set to zero.
      static_assert(
        (((terms & TermFlags::spatial_advection) != TermFlags::none) !=
         (dim_configuration_space == 0)),
        "If the spatial advection term is deactivated, the "
        "distribtuion function is assumed to be homogeneous, i.e. the "
        "dimension of the configuratin space needs to be set to zero.");
      // Compute the total dimension and set momentum to true
      static constexpr bool momentum =
        ((terms & TermFlags::momentum) != TermFlags::none ? true : false);
      static constexpr int dim = dim_configuration_space + momentum;
      static_assert(
        dim <= 3,
        "The total dimension must be greater than or equal to one and "
        "smaller or equal to three.");

      // Runtime settings
      // These settings are read from a parameter file
      int expansion_order;
      // Mesh
      std::vector<bool>         periodicity;
      dealii::Point<dim>        p1;
      dealii::Point<dim>        p2;
      std::vector<unsigned int> n_cells;
      // Finite element
      unsigned int polynomial_degree;
      // Time stepping
      Utils::TimeSteppingMethod time_stepping_method;
      double                    theta;
      double                    time_step;
      double                    final_time;

    private:
      // void
      // declare_parameters();
      // void
      // parse_parameters();
      // void
      //                          get_parameters();
      // const std::string        parameter_file;
      // dealii::ParameterHandler parameter_handler;
    };
  } // namespace VFP
} // namespace Sapphire
#endif
