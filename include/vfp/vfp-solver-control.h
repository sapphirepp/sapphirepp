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

      // The following static_assert uses an exlusive or (xor),
      // represented in C/C++ as "!=" for expressions of boolean type: Either
      // the spatial advection term is included in the equation or the dimension
      // of the configuration space is set to zero.
      static_assert(
        (((vfp_terms & TermFlags::spatial_advection) != TermFlags::none) !=
         (dim_configuration_space == 0)),
        "If the spatial advection term is deactivated, the "
        "distribtuion function is assumed to be homogeneous, i.e. the "
        "dimension of the configuratin space needs to be set to zero.");
      // Compute the total dimension and set momentum to true
      static constexpr bool momentum =
        ((vfp_terms & TermFlags::momentum) != TermFlags::none ? true : false);
      static constexpr int dim = dim_configuration_space + momentum;
      static_assert(
        dim <= 3,
        "The total dimension must be greater than or equal to one and "
        "smaller or equal to three.");

      // Runtime settings
      // These settings are read from a parameter file
      int expansion_order;
      // Mesh
      Utils::GridType           grid_type;
      std::string               grid_file;
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
