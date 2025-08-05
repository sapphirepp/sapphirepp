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
 * @file examples/vfp/convergence-study/config.h
 * @author Florian Schulze (florian.schulze@mpi-hd.mpg.de)
 * @brief Implement physical setup for convergence-study example
 */

#ifndef CONFIG_H
#define CONFIG_H

#include <deal.II/base/exceptions.h>
#include <deal.II/base/function.h>
#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/point.h>

#include <deal.II/lac/full_matrix.h>
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
    double B0;
    // Copy of VFP parameters for InitialValueFunction
    double box_length;
    double velocity;
    double gamma;
    double charge;
    double mass;
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

      /** [Parse runtime parameter] */
      B0 = prm.get_double("B0/2pi") * 2. * M_PI;
      /** [Parse runtime parameter] */

      prm.leave_subsection();
    }
  };



  namespace VFP
  {
    /** [Dimension] */
    /** Specify reduced phase space dimension \f$ (\mathbf{x}, p) \f$ */
    constexpr unsigned int dimension = 1;
    /** [Dimension] */



    /** [VFP Flags] */
    /** Specify which terms of the VFP equation should be active */
    constexpr VFPFlags vfp_flags =
      VFPFlags::time_evolution | VFPFlags::spatial_advection |
      VFPFlags::rotation | VFPFlags::time_independent_fields;
    /** [VFP Flags] */



    template <unsigned int dim>
    class InitialValueFunction : public dealii::Function<dim>
    {
    public:
      InitialValueFunction(const PhysicalParameters &physical_parameters,
                           const unsigned int        system_size)
        : dealii::Function<dim>(system_size)
        , prm{physical_parameters}
        , lms_indices{PDESystem::create_lms_indices(system_size)}
      {}



      void
      vector_value([[maybe_unused]] const dealii::Point<dim> &point,
                   dealii::Vector<double>                    &f) const override
      {
        AssertDimension(f.size(), this->n_components);

        /** [Initial value constants] */
        const double t     = this->get_time();
        const double x     = point[0];
        const double omega = prm.B0 * prm.charge / (prm.gamma * prm.mass);

        // Coefficients
        const double A[2] = {2., 0.};
        const double B[2] = {0., 1.};
        /** [Initial value constants] */

        /** [Initial value] */
        f = 0; // Initialize to zero
        for (unsigned int n = 0; n < 2; ++n)
          {
            const double k_n = n * 2. * M_PI / prm.box_length;
            const double c_n = prm.velocity / std::sqrt(3.) * k_n;

            // f_000
            f[0] +=
              (c_n * c_n / std::sqrt(omega * omega + c_n * c_n) *
                 (std::cos(std::sqrt(omega * omega + c_n * c_n) * t) - 1.) +
               std::sqrt(omega * omega + c_n * c_n)) *
              (A[n] * std::cos(k_n * x) - B[n] * std::sin(k_n * x));

            // f_110
            f[1] += omega * c_n / std::sqrt(omega * omega + c_n * c_n) *
                    (1. - std::cos(std::sqrt(omega * omega + c_n * c_n) * t)) *
                    (A[n] * std::sin(k_n * x) + B[n] * std::cos(k_n * x));

            // f_100
            f[2] += c_n * std::sin(std::sqrt(omega * omega + c_n * c_n) * t) *
                    (A[n] * std::sin(k_n * x) + B[n] * std::cos(k_n * x));

            // f_111
            f[3] += 0.;
          }
        /** [Initial value] */
      }



    private:
      const PhysicalParameters                       prm;
      const std::vector<std::array<unsigned int, 3>> lms_indices;
    };



    template <unsigned int dim>
    class BoundaryValueFunction : public dealii::Function<dim>
    {
    public:
      BoundaryValueFunction(const PhysicalParameters &physical_parameters,
                            const unsigned int        system_size)
        : dealii::Function<dim>(system_size)
        , prm{physical_parameters}
        , lms_indices{PDESystem::create_lms_indices(system_size)}
      {}



      void
      bc_vector_value_list(
        [[maybe_unused]] const std::vector<dealii::Point<dim>> &points,
        [[maybe_unused]] const unsigned int                     boundary_id,
        [[maybe_unused]] std::vector<dealii::Vector<double>>   &bc_values) const
      {
        AssertDimension(points.size(), bc_values.size());
        AssertDimension(bc_values[0].size(), this->n_components);

        for (unsigned int q_index = 0; q_index < points.size(); ++q_index)
          {
            /** [Boundary value] */
            // No inflow boundary
            /** [Boundary value] */
          }
      }



    private:
      const PhysicalParameters                       prm;
      const std::vector<std::array<unsigned int, 3>> lms_indices;
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
      value_list(
        const std::vector<dealii::Point<dim>> &points,
        std::vector<double>                   &scattering_frequencies,
        [[maybe_unused]] const unsigned int    component = 0) const override
      {
        AssertDimension(scattering_frequencies.size(), points.size());

        for (unsigned int q_index = 0; q_index < points.size(); ++q_index)
          {
            /** [Scattering frequency] */
            // No scattering frequency
            scattering_frequencies[q_index] = 0;
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
             unsigned int              system_size)
        : dealii::Function<dim>(system_size)
        , prm{physical_parameters}
        , lms_indices{PDESystem::create_lms_indices(system_size)}
      {}



      void
      vector_value([[maybe_unused]] const dealii::Point<dim> &point,
                   dealii::Vector<double> &source_values) const override
      {
        AssertDimension(source_values.size(), this->n_components);

        // NOLINTNEXTLINE(modernize-loop-convert)
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
      vector_value([[maybe_unused]] const dealii::Point<dim> &point,
                   dealii::Vector<double> &magnetic_field) const override
      {
        AssertDimension(magnetic_field.size(), this->n_components);

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
      vector_value([[maybe_unused]] const dealii::Point<dim> &point,
                   dealii::Vector<double> &velocity) const override
      {
        AssertDimension(velocity.size(), this->n_components);

        /** [Background velocity value] */
        // zero velocity field
        velocity[0] = 0.; // u_x
        velocity[1] = 0.; // u_y
        velocity[2] = 0.; // u_z
        /** [Background velocity value] */
      }



      void
      divergence_list(
        [[maybe_unused]] const std::vector<dealii::Point<dim>> &points,
        std::vector<double> &divergence) const
      {
        AssertDimension(divergence.size(), points.size());

        for (unsigned int q_index = 0; q_index < points.size(); ++q_index)
          {
            /** [Background velocity divergence] */
            // div u
            divergence[q_index] = 0.;
            /** [Background velocity divergence] */
          }
      }



      void
      material_derivative_list(
        const std::vector<dealii::Point<dim>> &points,
        std::vector<dealii::Vector<double>>   &material_derivatives) const
      {
        AssertDimension(material_derivatives.size(), points.size());
        AssertDimension(material_derivatives[0].size(), this->n_components);

        for (unsigned int q_index = 0; q_index < points.size(); ++q_index)
          {
            /** [Background velocity material derivative] */
            // zero velocity field
            material_derivatives[q_index][0] = 0.; // D/Dt u_x
            material_derivatives[q_index][1] = 0.; // D/Dt u_y
            material_derivatives[q_index][2] = 0.; // D/Dt u_z
            /** [Background velocity material derivative] */
          }
      }



      void
      jacobian_list(const std::vector<dealii::Point<dim>>   &points,
                    std::vector<dealii::FullMatrix<double>> &jacobians) const
      {
        AssertDimension(jacobians.size(), points.size());
        AssertDimension(jacobians[0].m(), this->n_components);
        AssertDimension(jacobians[0].n(), this->n_components);

        for (unsigned int q_index = 0; q_index < points.size(); ++q_index)
          {
            /** [Background velocity Jacobian] */
            // zero velocity field
            jacobians[q_index][0][0] = 0.; // \partial u_x / \partial x
            jacobians[q_index][0][1] = 0.; // \partial u_x / \partial y
            jacobians[q_index][0][2] = 0.; // \partial u_x / \partial z

            jacobians[q_index][1][0] = 0.; // \partial u_y / \partial x
            jacobians[q_index][1][1] = 0.; // \partial u_y / \partial y
            jacobians[q_index][1][2] = 0.; // \partial u_y / \partial z

            jacobians[q_index][2][0] = 0.; // \partial u_z / \partial x
            jacobians[q_index][2][1] = 0.; // \partial u_z / \partial y
            jacobians[q_index][2][2] = 0.; // \partial u_z / \partial z
            /** [Background velocity Jacobian] */
          }
      }



    private:
      const PhysicalParameters prm;
    };
  } // namespace VFP
} // namespace sapphirepp
#endif
