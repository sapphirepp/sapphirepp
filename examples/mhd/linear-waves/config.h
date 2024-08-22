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
 * @file examples/mhd/liner-waves/config.h
 * @author Florian Schulze (florian.schulze@mpi-hd.mpg.de)
 * @brief Implement physical setup for liner-waves example
 */

#ifndef CONFIG_H
#define CONFIG_H

#include <deal.II/base/exceptions.h>
#include <deal.II/base/function.h>
#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/point.h>

#include <deal.II/lac/vector.h>

#include <cmath>
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
    double rho_0;
    double u_0;
    double P_0;
    double amplitude;
    // Copy of MHD parameters for InitialValueFunction
    double box_length_x;
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
      prm.declare_entry("rho_0",
                        "1.",
                        dealii::Patterns::Double(0),
                        "Background density");
      prm.declare_entry("u_0",
                        "1",
                        dealii::Patterns::Double(),
                        "Background velocity in x-direction");
      prm.declare_entry("P_0",
                        "0.6",
                        dealii::Patterns::Double(0),
                        "Background pressure");
      prm.declare_entry("Amplitude",
                        "1e-3",
                        dealii::Patterns::Double(0),
                        "Amplitude of the perturbation");
      /** [Declare runtime parameter] */

      prm.leave_subsection();
    }



    void
    parse_parameters(dealii::ParameterHandler &prm)
    {
      dealii::LogStream::Prefix pre1("Startup", saplog);
      dealii::LogStream::Prefix pre2("PhysicalParameters", saplog);
      saplog << "Parsing parameters" << std::endl;
      prm.enter_subsection("Physical parameters");

      /** [Parse runtime parameter]  */
      rho_0     = prm.get_double("rho_0");
      u_0       = prm.get_double("u_0");
      P_0       = prm.get_double("P_0");
      amplitude = prm.get_double("Amplitude");
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



    template <unsigned int dim>
    class InitialConditionMHD : public dealii::Function<dim>
    {
    public:
      InitialConditionMHD(const PhysicalParameters &physical_parameters)
        : dealii::Function<dim>(MHDEquations<dim>::n_components)
        , prm{physical_parameters}
        , mhd_equations(prm.adiabatic_index)
        , primitive_background_state(MHDEquations<dim>::n_components)
        , primitive_eigenvectors(MHDEquations<dim>::n_components)
      {
        primitive_background_state = 0;
        primitive_background_state[MHDEquations<dim>::density_component] =
          prm.rho_0;
        primitive_background_state
          [MHDEquations<dim>::first_velocity_component] = prm.u_0;
        primitive_background_state[MHDEquations<dim>::pressure_component] =
          prm.P_0;

        // Density entropy wave
        primitive_eigenvectors                                       = 0;
        primitive_eigenvectors[MHDEquations<dim>::density_component] = 1.;

        eigenvalues = prm.u_0;
      }



      void
      vector_value(const dealii::Point<dim> &point,
                   dealii::Vector<double>   &f) const override
      {
        AssertDimension(f.size(), this->n_components);
        static_cast<void>(point); // suppress compiler warning

        /** [MHD Initial condition] */
        const double t = this->get_time();
        const double x = point[0];

        typename MHDEquations<dim>::state_type primitive_state =
          primitive_background_state;

        for (unsigned int c = 0; c < MHDEquations<dim>::n_components; ++c)
          primitive_state[c] +=
            primitive_eigenvectors[c] * prm.amplitude *
            std::sin(2. * M_PI / prm.box_length_x * (x - eigenvalues * t));

        mhd_equations.convert_primitive_to_conserved(primitive_state, f);
        /** [MHD Initial condition] */
      }



    private:
      const PhysicalParameters               prm;
      const MHDEquations<dim>                mhd_equations;
      typename MHDEquations<dim>::state_type primitive_background_state;
      typename MHDEquations<dim>::state_type primitive_eigenvectors;
      double                                 eigenvalues;
    };

  } // namespace MHD
} // namespace sapphirepp
#endif
