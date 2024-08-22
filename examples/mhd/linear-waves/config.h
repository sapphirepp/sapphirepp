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
    double                           rho_0;
    double                           P_0;
    std::vector<double>              u_0;
    std::vector<double>              B_0;
    double                           amplitude;
    std::vector<std::vector<double>> eigenmodes;
    // Copy of MHD parameters for InitialValueFunction
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
      prm.declare_entry("rho_0",
                        "1.",
                        dealii::Patterns::Double(0),
                        "Background density");
      prm.declare_entry("P_0",
                        "0.6",
                        dealii::Patterns::Double(0),
                        "Background pressure");
      prm.declare_entry("u_0",
                        "1., 0., 0.",
                        dealii::Patterns::Anything(),
                        "Background velocity in x, y, z direction");
      prm.declare_entry("B_0",
                        "0., 0., 0.",
                        dealii::Patterns::Anything(),
                        "Background magnetic field in x, y, z direction");
      prm.declare_entry("Amplitude",
                        "1e-3",
                        dealii::Patterns::Double(0),
                        "Amplitude of the perturbation");
      prm.declare_entry("Eigenmodes x",
                        "0, 0, 0, 0, 0, 0, 0, 0",
                        dealii::Patterns::Anything(),
                        "Select the eigenmodes in x direction");
      prm.declare_entry("Eigenmodes y",
                        "0, 0, 0, 0, 0, 0, 0, 0",
                        dealii::Patterns::Anything(),
                        "Select the eigenmodes in y direction");
      prm.declare_entry("Eigenmodes z",
                        "0, 0, 0, 0, 0, 0, 0, 0",
                        dealii::Patterns::Anything(),
                        "Select the eigenmodes in z direction");
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
      rho_0 = prm.get_double("rho_0");
      P_0   = prm.get_double("P_0");

      s = prm.get("u_0");
      std::stringstream u_0_string(s);
      for (std::string tmp; std::getline(u_0_string, tmp, ',');)
        u_0.push_back(std::stod(tmp));
      AssertThrow(u_0.size() == 3,
                  dealii::ExcMessage(
                    "Dimension of background velocity must be 3."));
      s = prm.get("B_0");
      std::stringstream B_0_string(s);
      for (std::string tmp; std::getline(B_0_string, tmp, ',');)
        B_0.push_back(std::stod(tmp));
      AssertThrow(B_0.size() == 3,
                  dealii::ExcMessage(
                    "Dimension of background magnetic field must be 3."));

      amplitude = prm.get_double("Amplitude");


      eigenmodes = std::vector<std::vector<double>>(3);

      s = prm.get("Eigenmodes x");
      std::stringstream em_x(s);
      for (std::string tmp; std::getline(em_x, tmp, ',');)
        eigenmodes[0].push_back(std::stod(tmp));
      AssertThrow(eigenmodes[0].size() == 8,
                  dealii::ExcMessage("Size of Eigenmodes x must be 8."));

      s = prm.get("Eigenmodes y");
      std::stringstream em_y(s);
      for (std::string tmp; std::getline(em_y, tmp, ',');)
        eigenmodes[1].push_back(std::stod(tmp));
      AssertThrow(eigenmodes[1].size() == 8,
                  dealii::ExcMessage("Size of Eigenmodes y must be 8."));

      s = prm.get("Eigenmodes z");
      std::stringstream em_z(s);
      for (std::string tmp; std::getline(em_z, tmp, ',');)
        eigenmodes[2].push_back(std::stod(tmp));
      AssertThrow(eigenmodes[2].size() == 8,
                  dealii::ExcMessage("Size of Eigenmodes z must be 8."));
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
        , conserved_background_state(MHDEquations<dim>::n_components)
        , primitive_eigenvectors(dim)
        , eigenvalues(dim)
      {
        // Create background state
        primitive_background_state[MHDEquations<dim>::density_component] =
          prm.rho_0;
        primitive_background_state[MHDEquations<dim>::pressure_component] =
          prm.P_0;
        for (unsigned int d = 0; d < MHDEquations<dim>::dim_uB; ++d)
          {
            primitive_background_state
              [MHDEquations<dim>::first_velocity_component + d] = prm.u_0[d];
            primitive_background_state
              [MHDEquations<dim>::first_magnetic_component + d] = prm.B_0[d];
          }

        saplog << "Primitive background state: " << std::endl;
        saplog << primitive_background_state << std::endl;

        mhd_equations.convert_primitive_to_conserved(
          primitive_background_state, conserved_background_state);
        saplog << "Conserved background state: " << std::endl;
        saplog << conserved_background_state << std::endl;

        // Define short hand definitions
        const double rho_0 =
          primitive_background_state[MHDEquations<dim>::density_component];
        const double P_0 =
          primitive_background_state[MHDEquations<dim>::pressure_component];
        dealii::Tensor<1, MHDEquations<dim>::dim_uB> u_0, B_0;
        for (unsigned int d = 0; d < MHDEquations<dim>::dim_uB; ++d)
          {
            u_0[d] = primitive_background_state
              [MHDEquations<dim>::first_velocity_component + d];
            B_0[d] = primitive_background_state
              [MHDEquations<dim>::first_magnetic_component + d];
          }
        (void)u_0;
        (void)B_0;

        const double a_s =
          std::sqrt(mhd_equations.adiabatic_index * P_0 / rho_0);
        (void)a_s;

        // Create normal vector
        dealii::Tensor<1, dim> normal;
        for (unsigned int d = 0; d < dim; ++d)
          {
            normal    = 0.;
            normal[d] = 1.;


            // Calculate eigenvalues
            eigenvalues[d] = dealii::Vector<double>(8);
            mhd_equations.compute_normale_eigenvalues(
              conserved_background_state, normal, eigenvalues[d]);
            saplog << "Eigenvalues in direction " << d << std::endl;
            saplog << eigenvalues[d] << std::endl;


            // Create eigenvectors
            primitive_eigenvectors[d] =
              std::vector<typename MHDEquations<dim>::state_type>(8);
            for (unsigned int i = 0; i < 8; ++i)
              {
                primitive_eigenvectors[d][i] =
                  typename MHDEquations<dim>::state_type(
                    MHDEquations<dim>::n_components);
                primitive_eigenvectors[d][i] = 0;
              }

            /** @todo Implement MHD waves */

            // Left going sound wave
            primitive_eigenvectors[d][0][MHDEquations<dim>::density_component] =
              rho_0;
            primitive_eigenvectors[d][0]
                                  [MHDEquations<dim>::first_velocity_component +
                                   d] = -a_s;
            primitive_eigenvectors[d][0]
                                  [MHDEquations<dim>::pressure_component] =
                                    rho_0 * a_s * a_s;
            saplog << "Left going sound wave in direction " << d
                   << ", eigenvalue=" << eigenvalues[d][0] << std::endl;
            saplog << primitive_eigenvectors[d][0] << std::endl;

            // Density entropy wave
            primitive_eigenvectors[d][3][MHDEquations<dim>::density_component] =
              1.;
            saplog << "Density entropy wave in direction " << d
                   << ", eigenvalue=" << eigenvalues[d][3] << std::endl;
            saplog << primitive_eigenvectors[d][3] << std::endl;

            // Right going sound wave
            primitive_eigenvectors[d][7][MHDEquations<dim>::density_component] =
              rho_0;
            primitive_eigenvectors[d][7]
                                  [MHDEquations<dim>::first_velocity_component +
                                   d] = a_s;
            primitive_eigenvectors[d][7]
                                  [MHDEquations<dim>::pressure_component] =
                                    rho_0 * a_s * a_s;
            saplog << "Right going sound wave in direction " << d
                   << ", eigenvalue=" << eigenvalues[d][7] << std::endl;
            saplog << primitive_eigenvectors[d][7] << std::endl;
          }
      }



      void
      vector_value(const dealii::Point<dim> &point,
                   dealii::Vector<double>   &f) const override
      {
        AssertDimension(f.size(), this->n_components);
        static_cast<void>(point); // suppress compiler warning

        /** [MHD Initial condition] */
        const double t = this->get_time();

        typename MHDEquations<dim>::state_type primitive_state =
          primitive_background_state;

        for (unsigned int d = 0; d < dim; ++d)
          for (unsigned int i = 0; i < 8; ++i)
            for (unsigned int c = 0; c < MHDEquations<dim>::n_components; ++c)
              primitive_state[c] +=
                primitive_eigenvectors[d][i][c] * prm.amplitude *
                prm.eigenmodes[d][i] *
                std::sin(2. * M_PI / prm.box_length[d] *
                         (point[d] - eigenvalues[d][i] * t));

        mhd_equations.convert_primitive_to_conserved(primitive_state, f);
        /** [MHD Initial condition] */
      }



    private:
      const PhysicalParameters               prm;
      const MHDEquations<dim>                mhd_equations;
      typename MHDEquations<dim>::state_type primitive_background_state;
      typename MHDEquations<dim>::state_type conserved_background_state;
      std::vector<std::vector<typename MHDEquations<dim>::state_type>>
                                          primitive_eigenvectors;
      std::vector<dealii::Vector<double>> eigenvalues;
    };

  } // namespace MHD
} // namespace sapphirepp
#endif
