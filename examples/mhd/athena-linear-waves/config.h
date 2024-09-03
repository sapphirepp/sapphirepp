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
 * @file examples/mhd/athena-liner-waves/config.h
 * @author Florian Schulze (florian.schulze@mpi-hd.mpg.de)
 * @brief Implement physical setup for athena-liner-waves example
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
    double       amplitude;
    // Copy of MHD parameters for InitialValueFunction
    unsigned int        dimension;
    std::vector<double> box_length;
    double              adiabatic_index;
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
                        dealii::Patterns::Integer(0, 4),
                        "Test case to run: 0 - Sound Wave, 1 - Entropy Wave, "
                        "2 - Fast Wave, 3 - Alfven Wave, 4 - Slow Wave");
      prm.declare_entry("Amplitude",
                        "1e-4",
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
      std::string s;
      prm.enter_subsection("Physical parameters");

      /** [Parse runtime parameter]  */
      test_case = static_cast<unsigned int>(prm.get_integer("Test case"));
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



    template <unsigned int spacedim>
    class InitialConditionMHD : public dealii::Function<spacedim>
    {
    public:
      InitialConditionMHD(const PhysicalParameters &physical_parameters)
        : dealii::Function<spacedim>(MHDEquations::n_components)
        , prm{physical_parameters}
        , mhd_equations(prm.adiabatic_index)
        , background_state(MHDEquations::n_components)
        , eigenvector(MHDEquations::n_components)
      {
        const double density  = 1.;
        const double pressure = 1. / prm.adiabatic_index;
        double       u_x      = 0.;
        double       u_y      = 0.;
        double       u_z      = 0.;
        if (prm.test_case == 1)
          {
            u_x = 1.;
          }
        double b_x = 1.;
        double b_y = std::sqrt(2.);
        double b_z = 0.5;
        if (prm.test_case < 2)
          {
            b_x = 0.;
            b_y = 0.;
            b_z = 0.;
          }

        const double energy =
          0.5 * density * (u_x * u_x + u_y * u_y + u_z * u_z) +
          pressure / (prm.adiabatic_index - 1.) +
          0.5 * (b_x * b_x + b_y * b_y + b_z * b_z);

        // Background state in conserved variables
        background_state[MHDEquations::density_component] = density;
        background_state[MHDEquations::first_momentum_component + 0] =
          density * u_x;
        background_state[MHDEquations::first_momentum_component + 1] =
          density * u_y;
        background_state[MHDEquations::first_momentum_component + 2] =
          density * u_z;
        background_state[MHDEquations::energy_component]             = energy;
        background_state[MHDEquations::first_magnetic_component + 0] = b_x;
        background_state[MHDEquations::first_magnetic_component + 1] = b_y;
        background_state[MHDEquations::first_magnetic_component + 2] = b_z;

        saplog << "Conserved background state: " << std::endl;
        saplog << background_state << std::endl;

        // Eigenvector of linear wave
        eigenvector = 0;
        switch (prm.test_case)
          {
            case 0:
              {
                // Left going sound wave
                eigenvector[0] = 1.;
                eigenvector[1] = -1.;
                eigenvector[2] = 0.;
                eigenvector[3] = 0.;
                eigenvector[4] = 3. / 2.;
                eigenvalue     = -1.;
                saplog << "Left going sound wave, lambda = -a_s = "
                       << eigenvalue << std::endl;
                saplog << eigenvector << std::endl;
                break;
              }
            case 1:
              {
                // Density entropy wave
                eigenvector[0] = 1.;
                eigenvector[1] = 1.;
                eigenvector[2] = 0.;
                eigenvector[3] = 0.;
                eigenvector[4] = 0.5;
                eigenvalue     = 1.;
                saplog << "Density entropy wave, lambda = u = " << eigenvalue
                       << std::endl;
                saplog << eigenvector << std::endl;
                break;
              }
            case 2:
              {
                // Left going fast magneto-sonic wave
                eigenvector[0] = 4.472135954999580e-01;
                eigenvector[1] = -8.944271909999160e-01;
                eigenvector[2] = 4.216370213557840e-01;
                eigenvector[3] = 1.490711984999860e-01;
                eigenvector[4] = 2.012457825664615e+00;
                eigenvector[6] = 8.432740427115680e-01;
                eigenvector[7] = 2.981423969999720e-01;
                eigenvalue     = -2.;
                saplog << "Left going fast magneto-sonic wave, lambda = -c_f = "
                       << eigenvalue << std::endl;
                saplog << eigenvector << std::endl;
                break;
              }
            case 3:
              {
                // Left going Alfven wave
                eigenvector[0] = 0.000000000000000e+00;
                eigenvector[1] = 0.000000000000000e+00;
                eigenvector[2] = -3.333333333333333e-01;
                eigenvector[3] = 9.428090415820634e-01;
                eigenvector[4] = 0.000000000000000e+00;
                eigenvector[6] = -3.333333333333333e-01;
                eigenvector[7] = 9.428090415820634e-01;
                eigenvalue     = -1.;
                saplog << "Left going Alfven wave, lambda = -c_a = "
                       << eigenvalue << std::endl;
                saplog << eigenvector << std::endl;
                break;
              }
            case 4:
              {
                // Left going slow magneto-sonic wave
                eigenvector[0] = 8.944271909999159e-01;
                eigenvector[1] = -4.472135954999579e-01;
                eigenvector[2] = -8.432740427115680e-01;
                eigenvector[3] = -2.981423969999720e-01;
                eigenvector[4] = 6.708136850795449e-01;
                eigenvector[6] = -4.216370213557841e-01;
                eigenvector[7] = -1.490711984999860e-01;
                eigenvalue     = -0.5;
                saplog << "Left going slow magneto-sonic wave, lambda = -c_s = "
                       << eigenvalue << std::endl;
                saplog << eigenvector << std::endl;
                break;
              }
            default:
              AssertThrow(false, dealii::ExcNotImplemented());
          }

        wave_vector    = 0;
        wave_vector[0] = 2 * M_PI / prm.box_length[0];
        if (prm.dimension > 1)
          wave_vector[1] = 2 * M_PI / prm.box_length[1];

        wave_number = eigenvalue * std::sqrt(wave_vector * wave_vector);

        double theta = 0;
        if (prm.dimension == 2)
          theta = std::atan(prm.box_length[1] / prm.box_length[0]);

        saplog << prm.dimension << "D, k = " << wave_vector
               << ", omega = " << wave_number
               << ", T = " << 2 * M_PI / std::abs(wave_number)
               << ", theta = " << theta * 180 / M_PI << "deg" << std::endl;

        // Rotate state and eigenvector

        const double tmp_p0_x =
          background_state[MHDEquations::first_momentum_component + 0];
        const double tmp_p0_y =
          background_state[MHDEquations::first_momentum_component + 1];
        background_state[MHDEquations::first_momentum_component + 0] =
          std::cos(theta) * tmp_p0_x - std::sin(theta) * tmp_p0_y;
        background_state[MHDEquations::first_momentum_component + 1] =
          std::sin(theta) * tmp_p0_x + std::cos(theta) * tmp_p0_y;
        const double tmp_b0_x =
          background_state[MHDEquations::first_magnetic_component + 0];
        const double tmp_b0_y =
          background_state[MHDEquations::first_magnetic_component + 1];
        background_state[MHDEquations::first_magnetic_component + 0] =
          std::cos(theta) * tmp_b0_x - std::sin(theta) * tmp_b0_y;
        background_state[MHDEquations::first_magnetic_component + 1] =
          std::sin(theta) * tmp_b0_x + std::cos(theta) * tmp_b0_y;
        saplog << "Rotated background state: " << background_state << std::endl;

        const double tmp_p_x =
          eigenvector[MHDEquations::first_momentum_component + 0];
        const double tmp_p_y =
          eigenvector[MHDEquations::first_momentum_component + 1];
        eigenvector[MHDEquations::first_momentum_component + 0] =
          std::cos(theta) * tmp_p_x - std::sin(theta) * tmp_p_y;
        eigenvector[MHDEquations::first_momentum_component + 1] =
          std::sin(theta) * tmp_p_x + std::cos(theta) * tmp_p_y;
        const double tmp_b_x =
          eigenvector[MHDEquations::first_magnetic_component + 0];
        const double tmp_b_y =
          eigenvector[MHDEquations::first_magnetic_component + 1];
        eigenvector[MHDEquations::first_magnetic_component + 0] =
          std::cos(theta) * tmp_b_x - std::sin(theta) * tmp_b_y;
        eigenvector[MHDEquations::first_magnetic_component + 1] =
          std::sin(theta) * tmp_b_x + std::cos(theta) * tmp_b_y;
        saplog << "Rotated eigenvector: " << eigenvector << std::endl;
      }



      void
      vector_value(const dealii::Point<spacedim> &point,
                   dealii::Vector<double>        &f) const override
      {
        AssertDimension(f.size(), this->n_components);
        static_cast<void>(point); // suppress compiler warning

        /** [MHD Initial condition] */
        const double t = this->get_time();

        f = background_state;

        for (unsigned int c = 0; c < MHDEquations::n_components; ++c)
          {
            f[c] += prm.amplitude * eigenvector[c] *
                    std::sin(wave_vector * point - wave_number * t);
          }
        /** [MHD Initial condition] */
      }



    private:
      const PhysicalParameters            prm;
      const MHDEquations                  mhd_equations;
      typename MHDEquations::state_type   background_state;
      typename MHDEquations::state_type   eigenvector;
      double                              eigenvalue;
      dealii::Tensor<1, spacedim, double> wave_vector;
      double                              wave_number;
    };

  } // namespace MHD
} // namespace sapphirepp
#endif
