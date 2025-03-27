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
 * @file examples/mhd/orszag-tang/config.h
 * @author Florian Schulze (florian.schulze@mpi-hd.mpg.de)
 * @brief Implement physical setup for orszag-tang example
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
      /** [Parse runtime parameter]  */

      prm.leave_subsection();
    }
  };



  namespace MHD
  {
    /** [MHD Dimension] */
    /** Specify mhd configuration space dimension \f$ (\mathbf{x}) \f$ */
    static constexpr unsigned int dim_mhd = 2;
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
      {}



      void
      vector_value(const dealii::Point<dim> &point,
                   dealii::Vector<double>   &f) const override
      {
        AssertDimension(f.size(), this->n_components);
        static_cast<void>(point); // suppress compiler warning

        /** [MHD Initial condition] */
        typename MHDEquations<dim, hdc>::state_type primitive_state(
          MHDEquations<dim, hdc>::n_components);
        primitive_state = 0.;
        primitive_state[MHDEquations<dim, hdc>::density_component] =
          25. / (36. * M_PI);
        primitive_state[MHDEquations<dim, hdc>::pressure_component] =
          5. / (12. * M_PI);
        primitive_state[MHDEquations<dim, hdc>::first_velocity_component + 0] =
          -std::sin(2. * M_PI * point[1]);
        primitive_state[MHDEquations<dim, hdc>::first_velocity_component + 1] =
          std::sin(2. * M_PI * point[0]);

        const double b0 = 1. / std::sqrt(4. * M_PI);
        primitive_state[MHDEquations<dim, hdc>::first_magnetic_component + 0] =
          -b0 * std::sin(2. * M_PI * point[1]);
        primitive_state[MHDEquations<dim, hdc>::first_magnetic_component + 1] =
          b0 * std::sin(4. * M_PI * point[0]);


        mhd_equations.convert_primitive_to_conserved(primitive_state, f);
        /** [MHD Initial condition] */
      }



    private:
      const PhysicalParameters     prm;
      const MHDEquations<dim, hdc> mhd_equations;
    };

  } // namespace MHD
} // namespace sapphirepp
#endif
