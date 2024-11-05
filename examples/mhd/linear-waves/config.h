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
    double              rho_0;
    double              P_0;
    std::vector<double> u_0;
    std::vector<double> b_0;
    double              amplitude;
    std::vector<double> eigenmodes;
    std::vector<int>    direction;
    // Copy of MHD parameters for InitialValueFunction
    std::vector<double> box_length;
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
      prm.declare_entry("Eigenmodes",
                        "0, 0, 0, 0, 0, 0, 0, 0",
                        dealii::Patterns::Anything(),
                        "Select the eigenmodes");
      prm.declare_entry("Direction",
                        "1, 0, 0",
                        dealii::Patterns::Anything(),
                        "Integer multiple of wavelength in each direction: "
                        "n_x, n_y, n_z \n"
                        "The wave vector components are determined by"
                        "k_i = n_i 2 pi / L_i");
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

      s = prm.get("Eigenmodes");
      std::stringstream em(s);
      for (std::string tmp; std::getline(em, tmp, ',');)
        eigenmodes.push_back(std::stod(tmp));
      AssertThrow(eigenmodes.size() == 8,
                  dealii::ExcMessage("Size of Eigenmodes x must be 8."));

      s = prm.get("Direction");
      std::stringstream direction_string(s);
      for (std::string n; std::getline(direction_string, n, ',');)
        direction.push_back(std::stoi(n));
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



    template <unsigned int spacedim>
    class InitialConditionMHD : public dealii::Function<spacedim>
    {
    public:
      InitialConditionMHD(const PhysicalParameters &physical_parameters,
                          const double              adiabatic_index)

        : dealii::Function<spacedim>(MHDEquations::n_components)
        , prm{physical_parameters}
        , mhd_equations(adiabatic_index)
        , primitive_background_state(MHDEquations::n_components)
        , conserved_background_state(MHDEquations::n_components)
        , eigenvalues(MHDEquations::n_components)
        , eigenvectors(MHDEquations::n_components)
      {
        primitive_background_state[MHDEquations::density_component] = prm.rho_0;
        primitive_background_state[MHDEquations::pressure_component] = prm.P_0;
        for (unsigned int d = 0; d < spacedim; ++d)
          {
            primitive_background_state[MHDEquations::first_velocity_component +
                                       d] = prm.u_0[d];
            primitive_background_state[MHDEquations::first_magnetic_component +
                                       d] = prm.b_0[d];
          }

        saplog << "Primitive background state: " << std::endl;
        saplog << primitive_background_state << std::endl;

        mhd_equations.convert_primitive_to_conserved(
          primitive_background_state, conserved_background_state);
        saplog << "Conserved background state: " << std::endl;
        saplog << conserved_background_state << std::endl;

        AssertThrow(
          dim_mhd <= prm.direction.size(),
          dealii::ExcMessage(
            "Direction does not specify value in each dimension."
            "Please enter the multiple of wavelength in each direction: \n"
            "\tset Number of cells = n_x (, n_y) (, n_z)"));
        if (prm.direction.size() != dim_mhd)
          saplog << "WARNING: Direction specification does not match dimension!"
                 << std::endl;
        saplog << "Direction: ";
        for (unsigned int i = 0; i < prm.direction.size(); i++)
          saplog << prm.direction[i] << ", ";
        saplog << std::endl;
        AssertThrow(dim_mhd <= prm.box_length.size(),
                    dealii::ExcMessage(
                      "Box length does not specify value in each dimension."));
        if (prm.box_length.size() != dim_mhd)
          saplog
            << "WARNING: Box length specification does not match dimension!"
            << std::endl;
        saplog << "Box length: ";
        for (unsigned int i = 0; i < prm.box_length.size(); i++)
          saplog << prm.box_length[i] << ", ";
        saplog << std::endl;

        wave_vector = 0;
        for (unsigned int d = 0; d < dim_mhd; ++d)
          wave_vector[d] = prm.direction[d] * 2 * M_PI / prm.box_length[d];
        saplog << "Wave vector: " << wave_vector
               << ", norm = " << wave_vector.norm() << std::endl;

        wave_number = wave_vector.norm();
        saplog << "Wave number: " << wave_number << std::endl;

        const dealii::Tensor<1, spacedim> normal_vector =
          wave_vector / wave_number;
        saplog << "Normal vector: " << normal_vector
               << ", norm = " << normal_vector.norm() << std::endl;

        mhd_equations.compute_normale_eigenvalues(conserved_background_state,
                                                  normal_vector,
                                                  eigenvalues);
        saplog << "Eigenvalues: " << eigenvalues << std::endl;

        mhd_equations.compute_right_eigenvector_matrix(
          conserved_background_state, normal_vector, eigenvectors);

        saplog << "Left going fast magneto-sonic wave, eigenvalue=u-c_f="
               << eigenvalues[0] << std::endl;
        for (unsigned int c = 0; c < MHDEquations::n_components; ++c)
          saplog << eigenvectors[c][0] << " ";
        saplog << std::endl << std::endl;

        saplog << "Left going Alfven wave, eigenvalue=u-c_a=" << eigenvalues[1]
               << std::endl;
        for (unsigned int c = 0; c < MHDEquations::n_components; ++c)
          saplog << eigenvectors[c][1] << " ";
        saplog << std::endl << std::endl;

        saplog << "Left going slow magneto-sonic wave, eigenvalue=u-c_s="
               << eigenvalues[2] << std::endl;
        for (unsigned int c = 0; c < MHDEquations::n_components; ++c)
          saplog << eigenvectors[c][2] << " ";
        saplog << std::endl << std::endl;

        saplog << "Density entropy mode, eigenvalue=u=" << eigenvalues[3]
               << std::endl;
        for (unsigned int c = 0; c < MHDEquations::n_components; ++c)
          saplog << eigenvectors[c][3] << " ";
        saplog << std::endl << std::endl;

        saplog << "Magnetic entropy mode, eigenvalue=u=" << eigenvalues[4]
               << std::endl;
        for (unsigned int c = 0; c < MHDEquations::n_components; ++c)
          saplog << eigenvectors[c][4] << " ";
        saplog << std::endl << std::endl;

        saplog << "Right going slow magneto-sonic wave, eigenvalue=u+c_s="
               << eigenvalues[5] << std::endl;
        for (unsigned int c = 0; c < MHDEquations::n_components; ++c)
          saplog << eigenvectors[c][5] << " ";
        saplog << std::endl << std::endl;

        saplog << "Right going Alfven wave, eigenvalue=u+c_a=" << eigenvalues[6]
               << std::endl;
        for (unsigned int c = 0; c < MHDEquations::n_components; ++c)
          saplog << eigenvectors[c][6] << " ";
        saplog << std::endl << std::endl;

        saplog << "Right going fast magneto-sonic wave, eigenvalue=u+c_f="
               << eigenvalues[7] << std::endl;
        for (unsigned int c = 0; c < MHDEquations::n_components; ++c)
          saplog << eigenvectors[c][7] << " ";
        saplog << std::endl << std::endl;
      }



      void
      vector_value(const dealii::Point<spacedim> &point,
                   dealii::Vector<double>        &f) const override
      {
        AssertDimension(f.size(), this->n_components);
        static_cast<void>(point); // suppress compiler warning

        /** [MHD Initial condition] */
        const double t = this->get_time();

        f = conserved_background_state;

        for (unsigned int i = 0; i < MHDEquations::n_components; ++i)
          for (unsigned int c = 0; c < MHDEquations::n_components; ++c)
            f[c] +=
              prm.amplitude * prm.eigenmodes[i] * eigenvectors[c][i] *
              std::sin(wave_vector * point - eigenvalues[i] * wave_number * t);
        /** [MHD Initial condition] */
      }



    private:
      const PhysicalParameters          prm;
      const MHDEquations                mhd_equations;
      typename MHDEquations::state_type primitive_background_state;
      typename MHDEquations::state_type conserved_background_state;
      double                            wave_number;
      dealii::Tensor<1, spacedim>       wave_vector;
      dealii::Vector<double>            eigenvalues;
      dealii::FullMatrix<double>        eigenvectors;
    };

  } // namespace MHD
} // namespace sapphirepp
#endif
