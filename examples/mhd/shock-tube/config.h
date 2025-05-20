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
    unsigned int test_case = 0;
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
      prm.add_parameter("Test case",
                        test_case,
                        "Test case to run: "
                        "0 - Sod Shock Tube, "
                        "1 - Brio & Wu Shock Tube",
                        dealii::Patterns::Integer(0, 1));
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

      /** [Parse runtime parameter] */
      // Parameters are automatically parsed by add_parameter()
      /** [Parse runtime parameter] */

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



    template <unsigned int dim>
    class InitialConditionMHD : public dealii::Function<dim>
    {
    public:
      static constexpr bool hdc =
        (mhd_flags & MHDFlags::hyperbolic_divergence_cleaning) !=
        MHDFlags::none;

      InitialConditionMHD(const PhysicalParameters &physical_parameters,
                          const double              adiabatic_index)
        : dealii::Function<dim>(MHDEquations<dim, hdc>::n_components)
        , prm{physical_parameters}
        , mhd_equations(adiabatic_index)
        , state_l(MHDEquations<dim, hdc>::n_components)
        , state_r(MHDEquations<dim, hdc>::n_components)
      {
        const unsigned int i_rho = MHDEquations<dim, hdc>::density_component;
        const unsigned int i_P   = MHDEquations<dim, hdc>::pressure_component;
        const unsigned int i_ux =
          MHDEquations<dim, hdc>::first_velocity_component;
        const unsigned int i_Bx =
          MHDEquations<dim, hdc>::first_magnetic_component;
        const unsigned int i_By =
          MHDEquations<dim, hdc>::first_magnetic_component + 1;

        state_l = 0.;
        state_r = 0.;

        switch (prm.test_case)
          {
            case 1:
              // Brio Wu
              state_l[i_Bx] = 0.75;
              state_l[i_By] = 1.0;

              state_r[i_Bx] = 0.75;
              state_r[i_By] = -1.0;
              [[fallthrough]]; // HD part is same as for Sod

            case 0:
              // Sod
              state_l[i_rho] = 1.0;
              state_l[i_ux]  = 0.0;
              state_l[i_P]   = 1.0;

              state_r[i_rho] = 0.125;
              state_r[i_ux]  = 0.0;
              state_r[i_P]   = 0.1;
              break;

            default:
              Assert(false, dealii::ExcNotImplemented());
              break;
          }

        saplog << "Primitive left state: " << std::endl;
        saplog << state_l << std::endl;
        saplog << "Primitive right state: " << std::endl;
        saplog << state_r << std::endl;
      }



      void
      vector_value(const dealii::Point<dim> &point,
                   dealii::Vector<double>   &f) const override
      {
        AssertDimension(f.size(), this->n_components);
        static_cast<void>(point); // suppress compiler warning

        /** [MHD Initial condition] */
        if (point[0] < 0.)
          mhd_equations.convert_primitive_to_conserved(state_l, f);
        else
          mhd_equations.convert_primitive_to_conserved(state_r, f);
        /** [MHD Initial condition] */
      }



    private:
      const PhysicalParameters                    prm;
      const MHDEquations<dim, hdc>                mhd_equations;
      typename MHDEquations<dim, hdc>::state_type state_l;
      typename MHDEquations<dim, hdc>::state_type state_r;
    };

  } // namespace MHD
} // namespace sapphirepp
#endif
