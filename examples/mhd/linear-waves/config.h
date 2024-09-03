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
    double                           rho_0;
    double                           P_0;
    std::vector<double>              u_0;
    std::vector<double>              b_0;
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
      prm.declare_entry(
        "b_0",
        "0., 0., 0.",
        dealii::Patterns::Anything(),
        "Background magnetic field $b_0 = B_0 / sqrt(4 pi)$ in x, y, z direction");
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
      s = prm.get("b_0");
      std::stringstream b_0_string(s);
      for (std::string tmp; std::getline(b_0_string, tmp, ',');)
        b_0.push_back(std::stod(tmp));
      AssertThrow(b_0.size() == 3,
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



    template <unsigned int spacedim>
    class InitialConditionMHD : public dealii::Function<spacedim>
    {
    public:
      InitialConditionMHD(const PhysicalParameters &physical_parameters)
        : dealii::Function<spacedim>(MHDEquations::n_components)
        , prm{physical_parameters}
        , mhd_equations(prm.adiabatic_index)
        , primitive_background_state(MHDEquations::n_components)
        , conserved_background_state(MHDEquations::n_components)
        , primitive_eigenvectors(spacedim)
        , eigenvalues(spacedim)
      {
        // Helper variables
        const unsigned int rho = MHDEquations::density_component;
        const unsigned int P   = MHDEquations::pressure_component;
        const unsigned int u   = MHDEquations::first_velocity_component;
        const unsigned int b   = MHDEquations::first_magnetic_component;

        // Create background state
        primitive_background_state[rho] = prm.rho_0;
        primitive_background_state[P]   = prm.P_0;
        for (unsigned int d = 0; d < spacedim; ++d)
          {
            primitive_background_state[u + d] = prm.u_0[d];
            primitive_background_state[b + d] = prm.b_0[d];
          }

        saplog << "Primitive background state: " << std::endl;
        saplog << primitive_background_state << std::endl;

        mhd_equations.convert_primitive_to_conserved(
          primitive_background_state, conserved_background_state);
        saplog << "Conserved background state: " << std::endl;
        saplog << conserved_background_state << std::endl;

        // Define short hand definitions
        const double                rho_0 = primitive_background_state[rho];
        const double                P_0   = primitive_background_state[P];
        dealii::Tensor<1, spacedim> u_0, b_0;
        for (unsigned int d = 0; d < spacedim; ++d)
          {
            u_0[d] = primitive_background_state[u + d];
            b_0[d] = primitive_background_state[b + d];
          }



        // Create normal vector
        dealii::Tensor<1, spacedim> normal, eps_b;
        for (unsigned int d = 0; d < spacedim; ++d)
          {
            normal    = 0.;
            normal[d] = 1.;
            eps_b     = dealii::cross_product_3d<spacedim>(normal, b_0);

            const double b2 = b_0 * b_0;
            const double nb = normal * b_0;

            const double a_s2 = prm.adiabatic_index * P_0 / rho_0;
            const double c_a2 = b2 / rho_0;
            const double d_n =
              (a_s2 + c_a2) * (a_s2 + c_a2) - 4. * a_s2 * nb * nb / rho_0;
            const double c_s2 = 0.5 * (a_s2 + c_a2 - std::sqrt(d_n));
            const double c_f2 = 0.5 * (a_s2 + c_a2 + std::sqrt(d_n));

            saplog << "Normal vector: " << normal << std::endl;
            saplog << "a_s2=" << a_s2 << ", c_a2=" << c_a2 << ", d_n=" << d_n
                   << ", c_s2=" << c_s2 << ", c_f2=" << c_f2 << std::endl;

            // const double a_s   = std::sqrt(a_s2);
            const double c_s   = std::sqrt(c_s2);
            const double c_f   = std::sqrt(c_f2);
            const double del_s = rho_0 * c_s2 - nb * nb;
            const double del_f = rho_0 * c_f2 - nb * nb;

            saplog << "c_s=" << c_s << ", c_f=" << c_f << ", del_s=" << del_s
                   << ", del_f=" << del_f << std::endl;

            // Calculate eigenvalues
            eigenvalues[d] = dealii::Vector<double>(8);
            mhd_equations.compute_normale_eigenvalues(
              conserved_background_state, normal, eigenvalues[d]);
            saplog << std::endl
                   << "Eigenvalues in direction " << d << std::endl
                   << eigenvalues[d] << std::endl
                   << std::endl;


            // Create eigenvectors
            primitive_eigenvectors[d] =
              std::vector<typename MHDEquations::state_type>(8);
            for (unsigned int i = 0; i < 8; ++i)
              {
                primitive_eigenvectors[d][i] =
                  typename MHDEquations::state_type(MHDEquations::n_components);
                primitive_eigenvectors[d][i] = 0;
              }

            // Left going fast magneto-sonic wave
            primitive_eigenvectors[d][0][rho] = rho_0;
            primitive_eigenvectors[d][0][P]   = rho_0 * a_s2;
            for (unsigned int i = 0; i < spacedim; ++i)
              {
                primitive_eigenvectors[d][0][u + i] =
                  -c_f * normal[i] +
                  (1 - normal[i]) * c_f * b_0[d] * b_0[i] / del_f;
                primitive_eigenvectors[d][0][b + i] =
                  (1 - normal[i]) * rho_0 * c_f2 * b_0[i] / del_f;
              }
            saplog << "Left going fast magneto-sonic wave in direction " << d
                   << ", eigenvalue=" << eigenvalues[d][0] << std::endl;
            saplog << primitive_eigenvectors[d][0] << std::endl << std::endl;

            // Left going fast Alfven wave
            primitive_eigenvectors[d][1][rho] = 0;
            primitive_eigenvectors[d][1][P]   = 0;
            for (unsigned int i = 0; i < spacedim; ++i)
              {
                primitive_eigenvectors[d][1][u + i] =
                  -eps_b[i] / std::sqrt(rho_0);
                primitive_eigenvectors[d][1][b + i] = -eps_b[i];
              }
            saplog << "Left going Alfven wave in direction " << d
                   << ", eigenvalue=" << eigenvalues[d][1] << std::endl;
            saplog << primitive_eigenvectors[d][1] << std::endl << std::endl;


            // Left going slow magneto-sonic wave
            if (del_s > 0)
              {
                primitive_eigenvectors[d][2][rho] = rho_0;
                primitive_eigenvectors[d][2][P]   = rho_0 * a_s2;
                for (unsigned int i = 0; i < spacedim; ++i)
                  {
                    primitive_eigenvectors[d][2][u + i] =
                      -c_s * normal[i] +
                      (1 - normal[i]) * c_s * b_0[d] * b_0[i] / del_s;
                    primitive_eigenvectors[d][2][b + i] =
                      (1 - normal[i]) * rho_0 * c_s2 * b_0[i] / del_s;
                  }
                saplog << "Left going slow magneto-sonic wave in direction "
                       << d << ", eigenvalue=" << eigenvalues[d][2]
                       << std::endl;
                saplog << primitive_eigenvectors[d][2] << std::endl
                       << std::endl;
              }
            else
              {
                primitive_eigenvectors[d][2][rho] = 0;
                primitive_eigenvectors[d][2][P]   = 0;
                for (unsigned int i = 0; i < spacedim; ++i)
                  {
                    primitive_eigenvectors[d][2][u + i] = 0;
                    primitive_eigenvectors[d][2][b + i] = 0;
                  }
                primitive_eigenvectors[d][2][u + ((d - 1) % 3)] = 1;
                saplog
                  << "Left going transverse velocity entropy wave in direction "
                  << d << ", eigenvalue=" << eigenvalues[d][2] << std::endl;
                saplog << primitive_eigenvectors[d][2] << std::endl
                       << std::endl;
              }


            // Density entropy wave
            primitive_eigenvectors[d][3][rho] = 1.;
            saplog << "Density entropy wave in direction " << d
                   << ", eigenvalue=" << eigenvalues[d][3] << std::endl;
            saplog << primitive_eigenvectors[d][3] << std::endl << std::endl;

            // TODO: Remove this mode, it represents divB \neq 0
            // Magnetic entropy wave
            primitive_eigenvectors[d][4][b + d] = 1.;
            saplog << "Magnetic entropy wave in direction " << d
                   << ", eigenvalue=" << eigenvalues[d][4] << std::endl;
            saplog << primitive_eigenvectors[d][4] << std::endl << std::endl;


            // Right going slow magneto-sonic wave
            if (del_s > 0)
              {
                primitive_eigenvectors[d][5][rho] = rho_0;
                primitive_eigenvectors[d][5][P]   = rho_0 * a_s2;
                for (unsigned int i = 0; i < spacedim; ++i)
                  {
                    primitive_eigenvectors[d][5][u + i] =
                      c_s * normal[i] -
                      (1 - normal[i]) * c_s * b_0[d] * b_0[i] / del_s;
                    primitive_eigenvectors[d][5][b + i] =
                      (1 - normal[i]) * rho_0 * c_s2 * b_0[i] / del_s;
                  }
                saplog << "Right going slow magneto-sonic wave in direction "
                       << d << ", eigenvalue=" << eigenvalues[d][5]
                       << std::endl;
                saplog << primitive_eigenvectors[d][5] << std::endl
                       << std::endl;
              }
            else
              {
                primitive_eigenvectors[d][5][rho] = 0;
                primitive_eigenvectors[d][5][P]   = 0;
                for (unsigned int i = 0; i < spacedim; ++i)
                  {
                    primitive_eigenvectors[d][5][u + i] = 0;
                    primitive_eigenvectors[d][5][b + i] = 0;
                  }
                primitive_eigenvectors[d][5][u + ((d + 1) % 3)] = 1;
                saplog
                  << "Right going transverse velocity entropy wave in direction "
                  << d << ", eigenvalue=" << eigenvalues[d][5] << std::endl;
                saplog << primitive_eigenvectors[d][5] << std::endl
                       << std::endl;
              }

            // Right going fast Alfven wave
            primitive_eigenvectors[d][6][rho] = 0;
            primitive_eigenvectors[d][6][P]   = 0;
            for (unsigned int i = 0; i < spacedim; ++i)
              {
                primitive_eigenvectors[d][6][u + i] =
                  eps_b[i] / std::sqrt(rho_0);
                primitive_eigenvectors[d][6][b + i] = -eps_b[i];
              }
            saplog << "Right going Alfven wave in direction " << d
                   << ", eigenvalue=" << eigenvalues[d][6] << std::endl;
            saplog << primitive_eigenvectors[d][6] << std::endl << std::endl;

            // Right going fast magneto-sonic wave
            primitive_eigenvectors[d][7][rho] = rho_0;
            primitive_eigenvectors[d][7][P]   = rho_0 * a_s2;
            for (unsigned int i = 0; i < spacedim; ++i)
              {
                primitive_eigenvectors[d][7][u + i] =
                  c_f * normal[i] -
                  (1 - normal[i]) * c_f * b_0[d] * b_0[i] / del_f;
                primitive_eigenvectors[d][7][b + i] =
                  (1 - normal[i]) * rho_0 * c_f2 * b_0[i] / del_f;
              }
            saplog << "Right going fast magneto-sonic wave in direction " << d
                   << ", eigenvalue=" << eigenvalues[d][7] << std::endl;
            saplog << primitive_eigenvectors[d][7] << std::endl << std::endl;
          }
      }



      void
      vector_value(const dealii::Point<spacedim> &point,
                   dealii::Vector<double>        &f) const override
      {
        AssertDimension(f.size(), this->n_components);
        static_cast<void>(point); // suppress compiler warning

        /** [MHD Initial condition] */
        const double t = this->get_time();

        typename MHDEquations::state_type primitive_state =
          primitive_background_state;

        for (unsigned int d = 0; d < dim_mhd; ++d)
          for (unsigned int i = 0; i < 8; ++i)
            for (unsigned int c = 0; c < MHDEquations::n_components; ++c)
              primitive_state[c] +=
                primitive_eigenvectors[d][i][c] * prm.amplitude *
                prm.eigenmodes[d][i] *
                std::sin(2. * M_PI / prm.box_length[d] *
                         (point[d] - eigenvalues[d][i] * t));

        mhd_equations.convert_primitive_to_conserved(primitive_state, f);
        /** [MHD Initial condition] */
      }



    private:
      const PhysicalParameters          prm;
      const MHDEquations                mhd_equations;
      typename MHDEquations::state_type primitive_background_state;
      typename MHDEquations::state_type conserved_background_state;
      std::vector<std::vector<typename MHDEquations::state_type>>
                                          primitive_eigenvectors;
      std::vector<dealii::Vector<double>> eigenvalues;
    };

  } // namespace MHD
} // namespace sapphirepp
#endif
