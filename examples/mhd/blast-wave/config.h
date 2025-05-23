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
 * @file examples/mhd/blast-wave/config.h
 * @author Florian Schulze (florian.schulze@mpi-hd.mpg.de)
 * @brief Implement physical setup for blast-wave example
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
    double       radius    = 0.1;
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
                        "0 - HD, "
                        "1 - MHD",
                        dealii::Patterns::Integer(0, 1));
      prm.add_parameter("Radius",
                        radius,
                        "Initial radius of the blast wave",
                        dealii::Patterns::Double(0.));
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
      // Parameters are automatically parsed by add_parameter()
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
    static constexpr MHDFlags mhd_flags =
      MHDFlags::hyperbolic_divergence_cleaning;
    /** [MHD Flags] */



    template <unsigned int dim, bool divergence_cleaning>
    class InitialConditionMHD : public dealii::Function<dim>
    {
    public:
      /** Shorthand for @ref MHDEquations<dim, divergence_cleaning> */
      using MHDEqs = MHDEquations<dim, divergence_cleaning>;
      /** @ref MHDEquations::state_type */
      using state_type = typename MHDEqs::state_type;



      InitialConditionMHD(const PhysicalParameters &physical_parameters,
                          const MHDEqs             &mhd_equations)
        : dealii::Function<dim>(MHDEqs::n_components)
        , prm{physical_parameters}
        , mhd_equations{mhd_equations}
        , primitive_ambient_state(MHDEqs::n_components)
        , primitive_inner_state(MHDEqs::n_components)
      {
        primitive_ambient_state                             = 0.;
        primitive_ambient_state[MHDEqs::density_component]  = 1.0;
        primitive_ambient_state[MHDEqs::pressure_component] = 0.1;
        if (prm.test_case == 1)
          {
            primitive_ambient_state[MHDEqs::first_magnetic_component + 0] =
              M_SQRT1_2;
            primitive_ambient_state[MHDEqs::first_magnetic_component + 1] =
              M_SQRT1_2;
          }
        saplog << "Primitive ambient state: " << std::endl;
        saplog << primitive_ambient_state << std::endl;

        primitive_inner_state                             = 0.;
        primitive_inner_state[MHDEqs::density_component]  = 1.0;
        primitive_inner_state[MHDEqs::pressure_component] = 10.0;
        if (prm.test_case == 1)
          {
            primitive_inner_state[MHDEqs::first_magnetic_component + 0] =
              M_SQRT1_2;
            primitive_inner_state[MHDEqs::first_magnetic_component + 1] =
              M_SQRT1_2;
          }
        saplog << "Primitive inner state: " << std::endl;
        saplog << primitive_inner_state << std::endl;
      }



      void
      vector_value(const dealii::Point<dim> &point,
                   dealii::Vector<double>   &f) const override
      {
        AssertDimension(f.size(), this->n_components);
        static_cast<void>(point); // suppress compiler warning

        /** [MHD Initial condition] */
        if (point.norm() < prm.radius)
          mhd_equations.convert_primitive_to_conserved(primitive_inner_state,
                                                       f);
        else
          mhd_equations.convert_primitive_to_conserved(primitive_ambient_state,
                                                       f);
        /** [MHD Initial condition] */
      }



    private:
      const PhysicalParameters prm;
      const MHDEqs             mhd_equations;
      state_type               primitive_ambient_state;
      state_type               primitive_inner_state;
    };

  } // namespace MHD
} // namespace sapphirepp
#endif
