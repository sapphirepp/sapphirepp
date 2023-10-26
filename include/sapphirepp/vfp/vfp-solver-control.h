#ifndef VFP_VFPSOLVERCONTROL_H
#define VFP_VFPSOLVERCONTROL_H

#include <deal.II/base/parameter_handler.h>

#include <ostream>
#include <string>
#include <vector>

#include "config.h"
#include "vfp-flags.h"

namespace Sapphire
{
  namespace VFP
  {
    using namespace dealii;

    template <int dim>
    class VFPSolverControl
    {
    public:
      VFPSolverControl();

      void
      declare_parameters(ParameterHandler &prm);
      void
      parse_parameters(ParameterHandler &prm);

      // Runtime settings
      // These settings are read from a parameter file
      int expansion_order;
      // Mesh
      GridType                        grid_type;
      std::string                     grid_file;
      dealii::Point<dim>              p1;
      dealii::Point<dim>              p2;
      std::vector<unsigned int>       n_cells;
      std::vector<BoundaryConditions> boundary_conditions;
      // Finite element
      unsigned int polynomial_degree;
      // Time stepping
      TimeSteppingMethod time_stepping_method;
      double             theta;
      double             time_step;
      double             final_time;

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
