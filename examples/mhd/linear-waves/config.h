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
    double                rho_0 = 1.;
    double                P_0   = 0.6;
    dealii::Tensor<1, 3>  u_0{{1., 0., 0.}};
    dealii::Tensor<1, 3>  b_0{{0., 0., 0.}};
    double                amplitude = 1e-4;
    std::array<double, 9> eigenmodes{{0, 0, 0, 0, 0, 0, 0, 0, 0}};
    std::vector<int>      direction{1, 0, 0};
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
      prm.add_parameter("rho_0",
                        rho_0,
                        "Background density",
                        dealii::Patterns::Double(0));
      prm.add_parameter("P_0",
                        P_0,
                        "Background pressure",
                        dealii::Patterns::Double(0));
      prm.add_parameter("u_0", u_0, "Background velocity in x, y, z direction");
      prm.add_parameter(
        "b_0",
        b_0,
        "Background magnetic field $b_0 = B_0 / sqrt(4 pi)$ in x, y, z direction");
      prm.add_parameter("Amplitude",
                        amplitude,
                        "Amplitude of the perturbation",
                        dealii::Patterns::Double(0));
      prm.add_parameter("Eigenmodes", eigenmodes, "Select the eigenmodes");
      prm.add_parameter("Direction",
                        direction,
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

      /** [Parse runtime parameter] */
      // Parameters are automatically parsed by add_parameter()

      // Copy MHD parameters
      std::vector<double> p1, p2;
      prm.enter_subsection("MHD");
      prm.enter_subsection("Mesh");
      {
        // Two diagonally opposite corner points of the grid
        dealii::Patterns::Tools::to_value(prm.get("Point 1"), p1);
        dealii::Patterns::Tools::to_value(prm.get("Point 2"), p2);
      } // Mesh
      prm.leave_subsection();
      prm.leave_subsection();

      const unsigned int dimension = static_cast<unsigned int>(p1.size());
      AssertDimension(p1.size(), dimension);
      AssertDimension(p2.size(), dimension);

      box_length = std::vector<double>(dimension, 1.);
      for (unsigned int d = 0; d < dimension; ++d)
        box_length[d] = std::abs(p1[d] - p2[d]);


      if (direction.size() > dimension)
        {
          saplog.print_warning(
            "Direction specification does not match dimension!");
          direction.resize(dimension);
        }
      AssertThrow(
        direction.size() == dimension,
        dealii::ExcMessage(
          "Direction does not specify value in each dimension."
          "Please enter the multiple of wavelength in each direction: \n"
          "\tset Direction = n_x (, n_y) (, n_z) \n"
          "You entered " +
          std::to_string(direction.size()) + " of " +
          std::to_string(dimension) + " dimensions."));
      /** [Parse runtime parameter] */
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
    static constexpr MHDFlags mhd_flags = MHDFlags::no_limiting;
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
        , primitive_background_state(MHDEqs::n_components)
        , conserved_background_state(MHDEqs::n_components)
        , eigenvalues(MHDEqs::n_components)
        , eigenvectors(MHDEqs::n_components)
      {
        dealii::LogStream::Prefix pre1("InitialConditionMHD", saplog);

        primitive_background_state[MHDEqs::density_component]  = prm.rho_0;
        primitive_background_state[MHDEqs::pressure_component] = prm.P_0;
        for (unsigned int d = 0; d < MHDEqs::n_vec_components; ++d)
          {
            primitive_background_state[MHDEqs::first_velocity_component + d] =
              prm.u_0[d];
            primitive_background_state[MHDEqs::first_magnetic_component + d] =
              prm.b_0[d];
          }

        saplog << "Primitive background state: " << std::endl;
        saplog << primitive_background_state << std::endl;

        mhd_equations.convert_primitive_to_conserved(
          primitive_background_state, conserved_background_state);
        saplog << "Conserved background state: " << std::endl;
        saplog << conserved_background_state << std::endl;

        saplog << "Direction: ";
        for (int d : prm.direction)
          saplog << d << ", ";
        saplog << std::endl;
        AssertDimension(prm.direction.size(), dim);

        saplog << "Box length: ";
        for (double l : prm.box_length)
          saplog << l << ", ";
        saplog << std::endl;
        AssertDimension(prm.box_length.size(), dim);

        wave_vector = 0;
        for (unsigned int d = 0; d < dim; ++d)
          wave_vector[d] = prm.direction[d] * 2 * M_PI / prm.box_length[d];
        saplog << "Wave vector: " << wave_vector
               << ", norm = " << wave_vector.norm() << std::endl;

        wave_number = wave_vector.norm();
        saplog << "Wave number: " << wave_number << std::endl;

        const dealii::Tensor<1, dim> normal_vector = wave_vector / wave_number;
        saplog << "Normal vector: " << normal_vector
               << ", norm = " << normal_vector.norm() << std::endl;

        mhd_equations.compute_normale_eigenvalues(conserved_background_state,
                                                  normal_vector,
                                                  eigenvalues);
        saplog << "Eigenvalues: " << eigenvalues << std::endl;

        mhd_equations.compute_right_eigenvector_matrix(
          conserved_background_state, normal_vector, eigenvectors);

        dealii::FullMatrix<double> left_eigenvectors(MHDEqs::n_components);
        mhd_equations.compute_left_eigenvector_matrix(
          conserved_background_state, normal_vector, left_eigenvectors);


        saplog << "Left going fast magneto-sonic wave, eigenvalue=u-c_f="
               << eigenvalues[0] << std::endl;
        for (unsigned int c = 0; c < MHDEqs::n_components; ++c)
          saplog << eigenvectors[c][0] << " ";
        saplog << std::endl;
        for (unsigned int c = 0; c < MHDEqs::n_components; ++c)
          saplog << left_eigenvectors[0][c] << " ";
        saplog << std::endl << std::endl;

        saplog << "Left going Alfven wave, eigenvalue=u-c_a=" << eigenvalues[1]
               << std::endl;
        for (unsigned int c = 0; c < MHDEqs::n_components; ++c)
          saplog << eigenvectors[c][1] << " ";
        saplog << std::endl;
        for (unsigned int c = 0; c < MHDEqs::n_components; ++c)
          saplog << left_eigenvectors[1][c] << " ";
        saplog << std::endl << std::endl;

        saplog << "Left going slow magneto-sonic wave, eigenvalue=u-c_s="
               << eigenvalues[2] << std::endl;
        for (unsigned int c = 0; c < MHDEqs::n_components; ++c)
          saplog << eigenvectors[c][2] << " ";
        saplog << std::endl;
        for (unsigned int c = 0; c < MHDEqs::n_components; ++c)
          saplog << left_eigenvectors[2][c] << " ";
        saplog << std::endl << std::endl;

        saplog << "Density entropy mode, eigenvalue=u=" << eigenvalues[3]
               << std::endl;
        for (unsigned int c = 0; c < MHDEqs::n_components; ++c)
          saplog << eigenvectors[c][3] << " ";
        saplog << std::endl;
        for (unsigned int c = 0; c < MHDEqs::n_components; ++c)
          saplog << left_eigenvectors[3][c] << " ";
        saplog << std::endl << std::endl;

        saplog << "Magnetic divergence mode, eigenvalue=";
        if constexpr (divergence_cleaning)
          saplog << "u-c_h=";
        else
          saplog << "u=";
        saplog << eigenvalues[4] << std::endl;
        for (unsigned int c = 0; c < MHDEqs::n_components; ++c)
          saplog << eigenvectors[c][4] << " ";
        saplog << std::endl;
        for (unsigned int c = 0; c < MHDEqs::n_components; ++c)
          saplog << left_eigenvectors[4][c] << " ";
        saplog << std::endl << std::endl;

        saplog << "Right going slow magneto-sonic wave, eigenvalue=u+c_s="
               << eigenvalues[5] << std::endl;
        for (unsigned int c = 0; c < MHDEqs::n_components; ++c)
          saplog << eigenvectors[c][5] << " ";
        saplog << std::endl;
        for (unsigned int c = 0; c < MHDEqs::n_components; ++c)
          saplog << left_eigenvectors[5][c] << " ";
        saplog << std::endl << std::endl;

        saplog << "Right going Alfven wave, eigenvalue=u+c_a=" << eigenvalues[6]
               << std::endl;
        for (unsigned int c = 0; c < MHDEqs::n_components; ++c)
          saplog << eigenvectors[c][6] << " ";
        saplog << std::endl;
        for (unsigned int c = 0; c < MHDEqs::n_components; ++c)
          saplog << left_eigenvectors[6][c] << " ";
        saplog << std::endl << std::endl;

        saplog << "Right going fast magneto-sonic wave, eigenvalue=u+c_f="
               << eigenvalues[7] << std::endl;
        for (unsigned int c = 0; c < MHDEqs::n_components; ++c)
          saplog << eigenvectors[c][7] << " ";
        saplog << std::endl;
        for (unsigned int c = 0; c < MHDEqs::n_components; ++c)
          saplog << left_eigenvectors[7][c] << " ";
        saplog << std::endl << std::endl;

        if constexpr (divergence_cleaning)
          {
            saplog << "Magnetic divergence mode 2, eigenvalue=u+ch="
                   << eigenvalues[8] << std::endl;
            for (unsigned int c = 0; c < MHDEqs::n_components; ++c)
              saplog << eigenvectors[c][8] << " ";
            saplog << std::endl;
            for (unsigned int c = 0; c < MHDEqs::n_components; ++c)
              saplog << left_eigenvectors[8][c] << " ";
            saplog << std::endl << std::endl;
          }


        saplog << "Test left and right eigenvector matrices:" << std::endl;
        saplog << "Right eigenvector matrix R:" << std::endl;
        eigenvectors.print(saplog, 12, 5);
        saplog << "Left eigenvector matrix L:" << std::endl;
        left_eigenvectors.print(saplog, 12, 5);

        dealii::FullMatrix<double> result(MHDEqs::n_components);
        dealii::IdentityMatrix     id_matrix(MHDEqs::n_components);
        dealii::FullMatrix<double> identity(id_matrix);
        const double               epsilon_d = 1e-6;

        left_eigenvectors.mmult(result, eigenvectors);
        saplog << "L*R:" << std::endl;
        result.print(saplog, 12, 5);
        for (unsigned int i = 0; i < MHDEqs::n_components; ++i)
          for (unsigned int j = 0; j < MHDEqs::n_components; ++j)
            AssertThrow(std::abs(result[i][j] - identity[i][j]) < epsilon_d,
                        dealii::ExcMessage("Problem with eigenvectors: l_" +
                                           std::to_string(i + 1) + " * r_" +
                                           std::to_string(j + 1) +
                                           " != delta_ij"));

        eigenvectors.mmult(result, left_eigenvectors);
        saplog << "R*L:" << std::endl;
        result.print(saplog, 12, 5);
        for (unsigned int i = 0; i < MHDEqs::n_components; ++i)
          for (unsigned int j = 0; j < MHDEqs::n_components; ++j)
            AssertThrow(std::abs(result[i][j] - identity[i][j]) < epsilon_d,
                        dealii::ExcMessage("Problem with eigenvectors: (R*L)_" +
                                           std::to_string(i + 1) +
                                           std::to_string(j + 1) +
                                           " != delta_ij"));
      }



      void
      vector_value(const dealii::Point<dim> &point,
                   dealii::Vector<double>   &f) const override
      {
        AssertDimension(f.size(), this->n_components);
        static_cast<void>(point); // suppress compiler warning

        /** [MHD Initial condition] */
        const double t = this->get_time();

        f = conserved_background_state;

        for (unsigned int i = 0; i < MHDEqs::n_components; ++i)
          for (unsigned int c = 0; c < MHDEqs::n_components; ++c)
            f[c] +=
              prm.amplitude * prm.eigenmodes[i] * eigenvectors[c][i] *
              std::sin(wave_vector * point - eigenvalues[i] * wave_number * t);
        /** [MHD Initial condition] */
      }



    private:
      const PhysicalParameters   prm;
      const MHDEqs               mhd_equations;
      state_type                 primitive_background_state;
      state_type                 conserved_background_state;
      double                     wave_number;
      dealii::Tensor<1, dim>     wave_vector;
      dealii::Vector<double>     eigenvalues;
      dealii::FullMatrix<double> eigenvectors;
    };



    template <unsigned int dim, bool divergence_cleaning>
    class BoundaryValueFunctionMHD : public dealii::Function<dim>
    {
    public:
      /** Shorthand for @ref MHDEquations */
      using MHDEqs = MHDEquations<dim, divergence_cleaning>;
      /** @ref MHDEquations::state_type */
      using state_type = typename MHDEqs::state_type;



      BoundaryValueFunctionMHD(const PhysicalParameters &physical_parameters,
                               const MHDEqs             &mhd_equations)
        : dealii::Function<dim>(MHDEqs::n_components)
        , prm{physical_parameters}
        , mhd_equations{mhd_equations}
      {}



      void
      bc_vector_value_list(const std::vector<dealii::Point<dim>> &points,
                           const unsigned int                     boundary_id,
                           std::vector<dealii::Vector<double>> &bc_values) const
      {
        AssertDimension(points.size(), bc_values.size());
        AssertDimension(bc_values[0].size(), this->n_components);
        static_cast<void>(points); // suppress compiler warning
        static_cast<void>(boundary_id);
        static_cast<void>(bc_values);

        for (unsigned int q_index = 0; q_index < points.size(); ++q_index)
          {
            /** [MHD Boundary value] */
            if (boundary_id == 0)
              {
                // lower x
              }
            else if (boundary_id == 1)
              {
                // upper x
              }
            else if (boundary_id == 2)
              {
                // lower y
              }
            else if (boundary_id == 3)
              {
                // upper y
              }
            else if (boundary_id == 4)
              {
                // lower z
              }
            else if (boundary_id == 5)
              {
                // upper z
              }
            /** [MHD Boundary value] */
          }
      }



    private:
      const PhysicalParameters prm;
      const MHDEqs             mhd_equations;
    };

  } // namespace MHD
} // namespace sapphirepp
#endif
