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
/// @file examples/gyro-motion-f0/config.h
/// @author Florian Schulze (florian.schulze@mpi-hd.mpg.de)
/// @brief Implement physical setup for gyro-motion-f0 example

/// [Includes]
#ifndef CONFIG_H
#  define CONFIG_H

#  include <deal.II/base/conditional_ostream.h>
#  include <deal.II/base/exceptions.h>
#  include <deal.II/base/function.h>
#  include <deal.II/base/parameter_handler.h>
#  include <deal.II/base/point.h>

#  include <deal.II/lac/vector.h>

#  include <cmath>
#  include <ostream>
#  include <vector>

#  include "sapphirepp-logstream.h"
#  include "vfp-flags.h"

namespace Sapphire
{
  /// [Includes]
  /// [PhysicalParameters]
  class PhysicalParameters
  {
  public:
    PhysicalParameters() = default;

    void
    declare_parameters(dealii::ParameterHandler &prm)
    {
      dealii::LogStream::Prefix pre1("Startup", saplog);
      dealii::LogStream::Prefix pre2("Physical parameters", saplog);
      saplog << "Declaring parameters" << std::endl;
      prm.enter_subsection("Physical parameters");

      prm.declare_entry("B0/2pi",
                        "1.",
                        dealii::Patterns::Double(),
                        "Magnetic field strength");
      prm.declare_entry("sigma",
                        "1.",
                        dealii::Patterns::Double(0),
                        "Spread of the initial distribution");

      prm.leave_subsection();
    }

    void
    parse_parameters(dealii::ParameterHandler &prm)
    {
      dealii::LogStream::Prefix pre1("Startup", saplog);
      dealii::LogStream::Prefix pre2("PhysicalParameters", saplog);
      saplog << "Parsing parameters" << std::endl;
      prm.enter_subsection("Physical parameters");
      B0    = prm.get_double("B0/2pi") * 2. * M_PI;
      sigma = prm.get_double("sigma");

      prm.leave_subsection();
    }

    double B0;
    double sigma;
  };
  /// [PhysicalParameters]
  /// [Namespace VFP]
  namespace VFP
  {
    /// [Namespace VFP]
    /// [Dimension]
    static constexpr int dimension = 2;
    /// [Dimension]

    /// [VFP flags]
    static constexpr VFPFlags vfp_flags = VFPFlags::spatial_advection |
                                          VFPFlags::magnetic |
                                          VFPFlags::time_independent_fields;
    /// [VFP flags]

    /// [InitialValueFunction constructor]
    template <unsigned int dim>
    class InitialValueFunction : public dealii::Function<dim>
    {
    public:
      InitialValueFunction(const PhysicalParameters &physical_parameters,
                           unsigned int              exp_order)
        : dealii::Function<dim>((exp_order + 1) * (exp_order + 1))
        , sigma(physical_parameters.sigma)
      {}
      /// [InitialValueFunction constructor]

      /// [InitialValueFunction value]
      void
      vector_value(const dealii::Point<dim> &point,
                   dealii::Vector<double>   &f) const override
      {
        AssertDimension(f.size(), this->n_components);

        std::fill(f.begin(), f.end(), 0.);
        f = 0.;

        f[0] = std::exp(-point.norm_square() / (2. * sigma * sigma));
      }

      const double sigma;
    };
    /// [InitialValueFunction value]

    /// [MagneticField constructor]
    template <unsigned int dim>
    class MagneticField : public dealii::Function<dim>
    {
    public:
      MagneticField(const PhysicalParameters &physical_parameters)
        : dealii::Function<dim>(3)
        , B0(physical_parameters.B0)
      {}
      /// [MagneticField constructor]

      /// [MagneticField value]
      void
      vector_value(const dealii::Point<dim> &point,
                   dealii::Vector<double>   &magnetic_field) const override
      {
        AssertDimension(magnetic_field.size(), this->n_components);
        static_cast<void>(point); // suppress compiler warning

        // constant magnetic field
        magnetic_field[0] = 0.;
        magnetic_field[1] = 0.;
        magnetic_field[2] = B0;
      }

      const double B0;
    };
    /// [MagneticField value]

    /// [BackgroundVelocityField]
    template <unsigned int dim>
    class BackgroundVelocityField : public dealii::Function<dim>
    {
    public:
      BackgroundVelocityField(const PhysicalParameters &physical_parameters)
        : dealii::Function<dim>(3)
      {
        static_cast<void>(physical_parameters); // suppress compiler warning
      }

      void
      vector_value(const dealii::Point<dim> &point,
                   dealii::Vector<double>   &velocity) const override
      {
        AssertDimension(velocity.size(), this->n_components);
        static_cast<void>(point); // suppress compiler warning

        // zero velocity field
        velocity[0] = 0.; // u_x
        velocity[1] = 0.; // u_y
        velocity[2] = 0.; // u_z
      }

      // Divergence
      void
      divergence_list(const std::vector<dealii::Point<dim>> &points,
                      std::vector<double>                   &divergence)
      {
        AssertDimension(divergence.size(), points.size());
        static_cast<void>(points); // suppress compiler warning

        // zero velocity
        std::fill(divergence.begin(), divergence.end(), 0.);
      }

      // Material derivative
      void
      material_derivative_list(
        const std::vector<dealii::Point<dim>> &points,
        std::vector<dealii::Vector<double>>   &material_derivatives)
      {
        AssertDimension(material_derivatives[0].size(), this->n_components);
        AssertDimension(material_derivatives.size(), points.size());
        static_cast<void>(points); // suppress compiler warning

        for (unsigned int i = 0; i < points.size(); ++i)
          {
            // zero velocity
            material_derivatives[i][0] = 0.; // D/Dt u_x
            material_derivatives[i][1] = 0.; // D/Dt u_y
            material_derivatives[i][2] = 0.; // D/Dt u_z
          }
      }

      // Jacobian matrix
      void
      jacobian_list(
        const std::vector<dealii::Point<dim>>            &points,
        std::vector<std::vector<dealii::Vector<double>>> &jacobians) const
      {
        AssertDimension(jacobians[0].size(), this->n_components);
        AssertDimension(jacobians.size(), points.size());
        static_cast<void>(points); // suppress compiler warning

        for (unsigned int i = 0; i < points.size(); ++i)
          {
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
          }
      }
    };
    /// [BackgroundVelocityField]

    /// [ScatteringFrequency]
    template <unsigned int dim>
    class ScatteringFrequency : public dealii::Function<dim>
    {
    public:
      ScatteringFrequency(const PhysicalParameters &physical_parameters)
        : dealii::Function<dim>(1)
      {
        static_cast<void>(physical_parameters); // suppress compiler warning
      }

      void
      value_list(const std::vector<dealii::Point<dim>> &points,
                 std::vector<double>                   &scattering_frequencies,
                 const unsigned int component = 0) const override
      {
        // suppress compiler warning
        static_cast<void>(points);
        static_cast<void>(scattering_frequencies);
        static_cast<void>(component);
      }
    };
    /// [ScatteringFrequency]

    /// [Source]
    template <unsigned int dim>
    class Source : public dealii::Function<dim>
    {
    public:
      Source(const PhysicalParameters &physical_parameters,
             unsigned int              exp_order)
        : dealii::Function<dim>((exp_order + 1) * (exp_order + 1))
      {
        static_cast<void>(physical_parameters); // suppress compiler warning
      }

      void
      vector_value(const dealii::Point<dim> &point,
                   dealii::Vector<double>   &values) const override
      {
        // suppress compiler warning
        static_cast<void>(point);
        static_cast<void>(values);
      }
    };
    /// [Source]
    /// [Close namespaces]
  } // namespace VFP
} // namespace Sapphire
#endif
/// [Close namespaces]
