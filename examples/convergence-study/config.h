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
 * @file examples/convergence-study/config.h
 * @author Florian Schulze (florian.schulze@mpi-hd.mpg.de)
 * @brief Implement physical setup for convergence-study example
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

#include "pde-system.h"
#include "sapphirepp-logstream.h"
#include "vfp-flags.h"



namespace sapphirepp
{
  class PhysicalParameters
  {
  public:
    /** [Define runtime parameter] */
    double B0_2pi;
    double B0;
    /** [Define runtime parameter] */
    // Copy of VFP parameters for inital value function
    double velocity;
    double gamma;
    double charge;
    double mass;
    double box_length;

    PhysicalParameters() = default;



    void
    declare_parameters(dealii::ParameterHandler &prm)
    {
      dealii::LogStream::Prefix pre1("Startup", saplog);
      dealii::LogStream::Prefix pre2("Physical parameters", saplog);
      saplog << "Declaring parameters" << std::endl;
      prm.enter_subsection("Physical parameters");

      /** [Declare runtime parameter] */
      prm.declare_entry("B0/2pi",
                        "1.",
                        dealii::Patterns::Double(),
                        "Magnetic field strength");
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
      B0_2pi = prm.get_double("B0/2pi");
      B0     = B0_2pi * 2. * M_PI;
      /** [Parse runtime parameter]  */

      prm.leave_subsection();
    }
  };



  namespace VFP
  {
    /** [Dimension] */
    /** Specify reduced phase space dimension \f$ (\mathbf{x}, p) \f$ */
    static constexpr unsigned int dimension = 1;
    /** [Dimension] */



    /** [VFP Flags] */
    /** Specify which terms of the VFP equation should be active */
    static constexpr VFPFlags vfp_flags = VFPFlags::spatial_advection |
                                          VFPFlags::magnetic |
                                          VFPFlags::time_independent_fields;
    /** [VFP Flags] */



    template <unsigned int dim>
    class InitialValueFunction : public dealii::Function<dim>
    {
    public:
      InitialValueFunction(const PhysicalParameters &physical_parameters,
                           const unsigned int        exp_order)
        : dealii::Function<dim>((exp_order + 1) * (exp_order + 1))
        , prm{physical_parameters}
        , lms_indices{PDESystem::create_lms_indices(exp_order)}
        , V{prm.velocity}
        , L{prm.box_length}
        , omega{prm.B0 * prm.charge / (prm.gamma * prm.mass)}
      {}



      void
      vector_value(const dealii::Point<dim> &point,
                   dealii::Vector<double>   &f) const override
      {
        AssertDimension(f.size(), this->n_components);
        static_cast<void>(point); // suppress compiler warning

        /** [Initial value] */
        // 1D Fourier series analytic solution
        f = 0;

        const double t = this->get_time();
        const double x = point[0];

        for (unsigned int n = 1; n < 2; ++n)
          {
            const double A_n = 0.;
            const double B_n = 1.;
            const double k_n = n * 2. * M_PI / L;
            const double c_n = V * k_n / std::sqrt(3.);

            // f_000
            f[0] += c_n / std::sqrt(omega * omega + c_n * c_n) *
                    (std::cos(std::sqrt(omega * omega + c_n * c_n) * t) - 1. +
                     (omega * omega + c_n * c_n) / (c_n * c_n)) *
                    (A_n * std::cos(k_n * x) - B_n * std::sin(k_n * x));

            // f_110
            f[1] += omega / std::sqrt(omega * omega + c_n * c_n) *
                    (1. - std::cos(std::sqrt(omega * omega + c_n * c_n) * t)) *
                    (A_n * std::sin(k_n * x) + B_n * std::cos(k_n * x));

            // f_100
            f[2] += std::sin(std::sqrt(omega * omega + c_n * c_n) * t) *
                    (A_n * std::sin(k_n * x) + B_n * std::cos(k_n * x));

            // f_111
            f[3] += 0.;
          }
        /** [Initial value] */
      }



    private:
      const PhysicalParameters                       prm;
      const std::vector<std::array<unsigned int, 3>> lms_indices;
      const double                                   V;     /** Velocity */
      const double                                   L;     /** Length */
      const double                                   omega; /** Gyrofrequency */
    };



    template <unsigned int dim>
    class ScatteringFrequency : public dealii::Function<dim>
    {
    public:
      ScatteringFrequency(const PhysicalParameters &physical_parameters)
        : dealii::Function<dim>(1)
        , prm{physical_parameters}
      {}



      void
      value_list(const std::vector<dealii::Point<dim>> &points,
                 std::vector<double>                   &scattering_frequencies,
                 const unsigned int component = 0) const override
      {
        AssertDimension(scattering_frequencies.size(), points.size());
        static_cast<void>(component); // suppress compiler warning

        for (unsigned int i = 0; i < points.size(); ++i)
          {
            /** [Scattering frequency] */
            // No scattering frequency
            scattering_frequencies[i] = 0;
            /** [Scattering frequency] */
          }
      }



    private:
      const PhysicalParameters prm;
    };



    template <unsigned int dim>
    class Source : public dealii::Function<dim>
    {
    public:
      Source(const PhysicalParameters &physical_parameters,
             unsigned int              exp_order)
        : dealii::Function<dim>((exp_order + 1) * (exp_order + 1))
        , prm{physical_parameters}
        , lms_indices{PDESystem::create_lms_indices(exp_order)}
      {}



      void
      vector_value(const dealii::Point<dim> &point,
                   dealii::Vector<double>   &source_values) const override
      {
        AssertDimension(source_values.size(), this->n_components);
        static_cast<void>(point); // suppress compiler warning

        for (unsigned int i = 0; i < source_values.size(); ++i)
          {
            /** [Source] */
            // No Source
            source_values[i] = 0.;
            /** [Source] */
          }
      }



    private:
      const PhysicalParameters                       prm;
      const std::vector<std::array<unsigned int, 3>> lms_indices;
    };



    template <unsigned int dim>
    class MagneticField : public dealii::Function<dim>
    {
    public:
      MagneticField(const PhysicalParameters &physical_parameters)
        : dealii::Function<dim>(3)
        , prm{physical_parameters}
      {}



      void
      vector_value(const dealii::Point<dim> &point,
                   dealii::Vector<double>   &magnetic_field) const override
      {
        AssertDimension(magnetic_field.size(), this->n_components);
        static_cast<void>(point); // suppress compiler warning

        /** [Magnetic field] */
        magnetic_field[0] = 0.;     // B_x
        magnetic_field[1] = 0.;     // B_y
        magnetic_field[2] = prm.B0; // B_z
        /** [Magnetic field] */
      }



    private:
      const PhysicalParameters prm;
    };



    template <unsigned int dim>
    class BackgroundVelocityField : public dealii::Function<dim>
    {
    public:
      BackgroundVelocityField(const PhysicalParameters &physical_parameters)
        : dealii::Function<dim>(3)
        , prm{physical_parameters}
      {}



      void
      vector_value(const dealii::Point<dim> &point,
                   dealii::Vector<double>   &velocity) const override
      {
        AssertDimension(velocity.size(), this->n_components);
        static_cast<void>(point); // suppress compiler warning

        /** [Background velocity value] */
        // zero velocity field
        velocity[0] = 0.; // u_x
        velocity[1] = 0.; // u_y
        velocity[2] = 0.; // u_z
        /** [Background velocity value] */
      }



      void
      divergence_list(const std::vector<dealii::Point<dim>> &points,
                      std::vector<double>                   &divergence)
      {
        AssertDimension(divergence.size(), points.size());
        static_cast<void>(points); // suppress compiler warning

        for (unsigned int i = 0; i < divergence.size(); ++i)
          {
            /** [Background velocity divergence] */
            // div u
            divergence[i] = 0.;
            /** [Background velocity divergence] */
          }
      }



      void
      material_derivative_list(
        const std::vector<dealii::Point<dim>> &points,
        std::vector<dealii::Vector<double>>   &material_derivatives)
      {
        AssertDimension(material_derivatives.size(), points.size());
        AssertDimension(material_derivatives[0].size(), this->n_components);

        for (unsigned int i = 0; i < points.size(); ++i)
          {
            /** [Background velocity material derivative] */
            // zero velocity field
            material_derivatives[i][0] = 0.; // D/Dt u_x
            material_derivatives[i][1] = 0.; // D/Dt u_y
            material_derivatives[i][2] = 0.; // D/Dt u_z
            /** [Background velocity material derivative] */
          }
      }



      void
      jacobian_list(
        const std::vector<dealii::Point<dim>>            &points,
        std::vector<std::vector<dealii::Vector<double>>> &jacobians) const
      {
        AssertDimension(jacobians.size(), points.size());
        AssertDimension(jacobians[0].size(), this->n_components);
        AssertDimension(jacobians[0][0].size(), this->n_components);

        for (unsigned int i = 0; i < points.size(); ++i)
          {
            /** [Background velocity Jacobian] */
            // zero velocity field
            jacobians[i][0][0] = 0.; // \partial u_x / \partial x
            jacobians[i][0][1] = 0.; // \partial u_x / \partial y
            jacobians[i][0][2] = 0.; // \partial u_x / \partial z

            jacobians[i][1][0] = 0.; // \partial u_y / \partial x
            jacobians[i][1][1] = 0.; // \partial u_y / \partial y
            jacobians[i][1][2] = 0.; // \partial u_y / \partial z

            jacobians[i][2][0] = 0.; // \partial u_z / \partial x
            jacobians[i][2][1] = 0.; // \partial u_z / \partial y
            jacobians[i][2][2] = 0.; // \partial u_z / \partial z
            /** [Background velocity Jacobian] */
          }
      }



    private:
      const PhysicalParameters prm;
    };
  } // namespace VFP
} // namespace sapphirepp
#endif
