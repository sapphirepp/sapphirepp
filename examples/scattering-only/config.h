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
 * @file examples/scattering-only/config.h
 * @author Florian Schulze (florian.schulze@mpi-hd.mpg.de)
 * @brief Implement physical setup for scattering-only example
 */

/** [Includes] */
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

#  include "pde-system.h"
#  include "sapphirepp-logstream.h"
#  include "vfp-flags.h"
/** [Includes] */

/** [Namespace Sapphire] */
namespace sapphirepp
{
  /** [Namespace sapphirepp] */

  /** [PhysicalParameters] */
  class PhysicalParameters
  {
  public:
    double nu;
    double f0;

    PhysicalParameters() = default;
    /** [PhysicalParameters] */
    /** [Declare parameters] */
    void
    declare_parameters(dealii::ParameterHandler &prm)
    {
      prm.enter_subsection("Physical parameters");

      prm.declare_entry("nu",
                        "0.1",
                        dealii::Patterns::Double(),
                        "Scattering frequency");
      prm.declare_entry("f0",
                        "1.",
                        dealii::Patterns::Double(),
                        "Initial value of the expansion coefficients");

      prm.leave_subsection();
    }
    /** [Declare parameters] */

    /** [Parse parameters] */
    void
    parse_parameters(dealii::ParameterHandler &prm)
    {
      prm.enter_subsection("Physical parameters");

      nu = prm.get_double("nu");
      f0 = prm.get_double("f0");

      prm.leave_subsection();
    }
  };
  /** [Parse parameters] */

  /** [Namespace VFP] */
  namespace VFP
  {
    /** [Namespace VFP] */

    /** [Dimension] */
    static constexpr int dimension = 1;
    /** [Dimension] */

    /** [VFP flags] */
    static constexpr VFPFlags vfp_flags = VFPFlags::collision;
    /** [VFP flags] */

    /** [InitialValueFunction constructor] */
    template <unsigned int dim>
    class InitialValueFunction : public dealii::Function<dim>
    {
    public:
      InitialValueFunction(const PhysicalParameters &physical_parameters,
                           unsigned int              exp_order)
        : dealii::Function<dim>((exp_order + 1) * (exp_order + 1))
        , f0(physical_parameters.f0)
        , nu(physical_parameters.nu)
      {
        PDESystem::create_lms_indices(exp_order, lms_indices);
      }
      /** [InitialValueFunction constructor] */

      /** [InitialValueFunction value] */
      void
      vector_value(const dealii::Point<dim> &p,
                   dealii::Vector<double>   &f) const override
      {
        static_cast<void>(p); // suppress compiler warning

        for (unsigned int i = 0; i < InitialValueFunction<dim>::n_components;
             i++)
          {
            const unsigned int l = lms_indices[i][0];
            const double       t = this->get_time();
            f[i]                 = f0 * std::exp(-nu * l * (l + 1) / 2. * t);
          }
      }

    private:
      const double                             f0;
      const double                             nu;
      std::vector<std::array<unsigned int, 3>> lms_indices;
    };
    /** [InitialValueFunction value] */

    /** [ScatteringFrequency] */
    template <unsigned int dim>
    class ScatteringFrequency : public dealii::Function<dim>
    {
    public:
      ScatteringFrequency(const PhysicalParameters &physical_parameters)
        : dealii::Function<dim>(1)
        , nu(physical_parameters.nu)
      {}

      void
      value_list(const std::vector<dealii::Point<dim>> &points,
                 std::vector<double>                   &scattering_frequencies,
                 const unsigned int component = 0) const override
      {
        Assert(scattering_frequencies.size() == points.size(),
               dealii::ExcDimensionMismatch(scattering_frequencies.size(),
                                            points.size()));
        // suppress compiler warning
        static_cast<void>(component);
        static_cast<void>(points);

        std::fill(scattering_frequencies.begin(),
                  scattering_frequencies.end(),
                  nu);
      }

    private:
      const double nu;
    };
    /** [ScatteringFrequency] */

    /** [Source] */
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
      vector_value(const dealii::Point<dim> &p,
                   dealii::Vector<double>   &values) const override
      {
        static_cast<void>(p); // suppress compiler warning

        values[0] = 0.; // unused
      }
    };
    /** [Source] */

    /** [MagneticField] */
    template <unsigned int dim>
    class MagneticField : public dealii::Function<dim>
    {
    public:
      MagneticField(const PhysicalParameters &physical_parameters)
        : dealii::Function<dim>(3)
      {
        static_cast<void>(physical_parameters); // suppress compiler warning
      }

      void
      vector_value(const dealii::Point<dim> &point,
                   dealii::Vector<double>   &magnetic_field) const override
      {
        static_cast<void>(point); // suppress compiler warning

        // zero magnetic field
        magnetic_field[0] = 0.;
        magnetic_field[1] = 0.;
        magnetic_field[2] = 0.;
      }
    };
    /** [MagneticField] */

    /** [BackgroundVelocityField] */
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
        Assert(velocity.size() == BackgroundVelocityField<dim>::n_components,
               dealii::ExcDimensionMismatch(
                 velocity.size(), BackgroundVelocityField<dim>::n_components));
        static_cast<void>(point); // suppress compiler warning

        // zero velocity field
        velocity[0] = .0; // u_x
        velocity[1] = .0; // u_y
        velocity[2] = .0; // u_z
      }

      // Divergence
      void
      divergence_list(const std::vector<dealii::Point<dim>> &points,
                      std::vector<double>                   &divergence)
      {
        Assert(divergence.size() == points.size(),
               dealii::ExcDimensionMismatch(divergence.size(), points.size()));
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
        Assert(material_derivatives[0].size() == 3,
               dealii::ExcDimensionMismatch(material_derivatives[0].size(), 3));
        Assert(material_derivatives.size() == points.size(),
               dealii::ExcDimensionMismatch(material_derivatives.size(),
                                            points.size()));
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
        Assert(jacobians.size() == points.size(),
               dealii::ExcDimensionMismatch(jacobians.size(), points.size()));
        Assert(jacobians[0].size() == 3,
               dealii::ExcDimensionMismatch(jacobians[0].size(), 3));
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
    /** [BackgroundVelocityField] */
    /** [Close namespaces] */
  } // namespace VFP
} // namespace sapphirepp
#endif
/** [Close namespaces] */
