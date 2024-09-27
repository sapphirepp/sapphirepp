// -----------------------------------------------------------------------------
//
// Copyright (C) 2023 by the Sapphire++ authors
//
// This file is part of Sapphire++.
//
// Sapphire++ is free software: you can redistribute it and/or modify it under
// the terms of the GNU Lesser General Public License as published by the Free
// Software Foundation, either version 3 of the License, or (at your option) any
// later version.
//
// Sapphire++ is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
// A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more
// details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with Sapphire++. If not, see <https://www.gnu.org/licenses/>.
//
// -----------------------------------------------------------------------------

/**
 * @file examples/mhd/shock-tube/config.h
 * @author Florian Schulze (florian.schulze@mpi-hd.mpg.de)
 * @brief Implement physical setup for shock-tube example
 */

#ifndef CONFIG_H
#define CONFIG_H

#include <deal.II/base/exceptions.h>
#include <deal.II/base/function.h>
#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/point.h>
#include <deal.II/base/tensor.h>

#include <deal.II/lac/vector.h>

#include <cmath>
#include <sstream>
#include <vector>

#include "mhd-equations.h"
#include "mhd-flags.h"
#include "sapphirepp-logstream.h"



namespace sapphirepp
{
  class PhysicalParameters
  {
  public:
    /** [Define runtime parameter] */
    unsigned int test_case;
    // Copy of MHD parameters for InitialValueFunction
    double adiabatic_index;
    /** [Define runtime parameter] */

    PhysicalParameters() = default;



    void
    declare_parameters(dealii::ParameterHandler &prm)
    {
      dealii::LogStream::Prefix pre1("Startup", saplog);
      dealii::LogStream::Prefix pre2("Physical parameters", saplog);
      saplog << "Declaring parameters" << std::endl;
      prm.enter_subsection("Physical parameters");

      /** [Declare runtime parameter] */
      prm.declare_entry("Test case",
                        "0",
                        dealii::Patterns::Integer(0, 1),
                        "Test case to run: 0 - Sod Shock Tube, "
                        "1 - Brio & Wu Shock Tube");
      /** [Declare runtime parameter] */

      prm.leave_subsection();
    }



    void
    parse_parameters(dealii::ParameterHandler &prm)
    {
      dealii::LogStream::Prefix pre1("Startup", saplog);
      dealii::LogStream::Prefix pre2("PhysicalParameters", saplog);
      saplog << "Parsing parameters" << std::endl;
      std::string s;
      prm.enter_subsection("Physical parameters");

      /** [Parse runtime parameter]  */
      test_case = static_cast<unsigned int>(prm.get_integer("Test case"));
      /** [Parse runtime parameter]  */

      prm.leave_subsection();
    }
  };



  namespace MHD
  {
    /** [MHD Dimension] */
    /** Specify mhd configuration space dimension \f$ (\mathbf{x}) \f$ */
    static constexpr unsigned int dim_mhd = 1;
    /** [MHD Dimension] */



    /** [MHD Flags] */
    /** Specify which MHD flags should be active */
    static constexpr MHDFlags mhd_flags = MHDFlags::none;
    /** [MHD Flags] */



    template <unsigned int spacedim>
    class InitialConditionMHD : public dealii::Function<spacedim>
    {
    public:
      InitialConditionMHD(const PhysicalParameters &physical_parameters)
        : dealii::Function<spacedim>(MHDEquations::n_components)
        , prm{physical_parameters}
        , mhd_equations(prm.adiabatic_index)
        , primitive_left_state(MHDEquations::n_components)
        , primitive_right_state(MHDEquations::n_components)
      {
        const double rho_l = 1.0;
        const double u_l   = 0.0;
        const double P_l   = 1.0;
        const double rho_r = 0.125;
        const double u_r   = 0.0;
        const double P_r   = 0.1;

        primitive_left_state                                         = 0.;
        primitive_left_state[MHDEquations::density_component]        = rho_l;
        primitive_left_state[MHDEquations::first_velocity_component] = u_l;
        primitive_left_state[MHDEquations::pressure_component]       = P_l;
        saplog << "Primitive left state: " << std::endl;
        saplog << primitive_left_state << std::endl;

        primitive_right_state                                         = 0.;
        primitive_right_state[MHDEquations::density_component]        = rho_r;
        primitive_right_state[MHDEquations::first_velocity_component] = u_r;
        primitive_right_state[MHDEquations::pressure_component]       = P_r;
        saplog << "Primitive right state: " << std::endl;
        saplog << primitive_right_state << std::endl;
      }



      void
      vector_value(const dealii::Point<spacedim> &point,
                   dealii::Vector<double>        &f) const override
      {
        AssertDimension(f.size(), this->n_components);
        static_cast<void>(point); // suppress compiler warning

        /** [MHD Initial condition] */
        if (point[0] < 0.5)
          mhd_equations.convert_primitive_to_conserved(primitive_left_state, f);
        else
          mhd_equations.convert_primitive_to_conserved(primitive_right_state,
                                                       f);
        /** [MHD Initial condition] */
      }



    private:
      const PhysicalParameters          prm;
      const MHDEquations                mhd_equations;
      typename MHDEquations::state_type primitive_left_state;
      typename MHDEquations::state_type primitive_right_state;
    };

  } // namespace MHD
} // namespace sapphirepp
#endif
