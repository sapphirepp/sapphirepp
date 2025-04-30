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
 * @file examples/vfp/steady-state-parallel-shock/config.h
 * @author Nils Schween (nils.schween@mpi-hd.mpg.de)
 * @brief Implement physical setup for steady-state-parallel-shock example
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
    double u_sh;
    double B0;
    double compression_ratio;
    double shock_width;
    double nu0;

    // Source
    double Q;
    double p_inj;
    double x_inj;
    double sig_p;
    double sig_x;
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
      prm.declare_entry("u_sh",
                        "0.1",
                        dealii::Patterns::Double(),
                        "The shock velocity.");
      prm.declare_entry("B0",
                        "1.",
                        dealii::Patterns::Double(),
                        "The magnetic field strength upstream.");
      prm.declare_entry("compression ratio",
                        "4.",
                        dealii::Patterns::Double(),
                        "The compression ratio of the shock.");
      prm.declare_entry("shock width",
                        "0.04",
                        dealii::Patterns::Double(),
                        "The width of the shock.");
      prm.declare_entry("nu0",
                        "0.1",
                        dealii::Patterns::Double(),
                        "The scattering frequency.");

      prm.enter_subsection("Source");
      prm.declare_entry("Q",
                        "0.1",
                        dealii::Patterns::Double(),
                        "The injection rate.");
      prm.declare_entry("p_inj",
                        "2.",
                        dealii::Patterns::Double(),
                        "The injection momentum.");
      prm.declare_entry("x_inj",
                        "0.0",
                        dealii::Patterns::Double(),
                        "The injection position.");
      prm.declare_entry("sig_p",
                        "0.125",
                        dealii::Patterns::Double(),
                        "The width of the source in momentum space.");
      prm.declare_entry("sig_x",
                        "0.125",
                        dealii::Patterns::Double(),
                        "The width of the source in configuration space.");
      prm.leave_subsection();
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
      u_sh              = prm.get_double("u_sh");
      B0                = prm.get_double("B0");
      compression_ratio = prm.get_double("compression ratio");
      shock_width       = prm.get_double("shock width");
      nu0               = prm.get_double("nu0");

      prm.enter_subsection("Source");
      p_inj = prm.get_double("p_inj");
      x_inj = prm.get_double("x_inj");
      sig_x = prm.get_double("sig_x");
      sig_p = prm.get_double("sig_p");
      Q     = prm.get_double("Q");
      prm.leave_subsection();
      /** [Parse runtime parameter] */

      prm.leave_subsection();
    }
  };



  namespace VFP
  {
    /** [Dimension] */
    /** Specify reduced phase space dimension \f$ (\mathbf{x}, p) \f$ */
    constexpr unsigned int dimension = 2;
    /** [Dimension] */



    /** [VFP Flags] */
    /** Specify which terms of the VFP equation should be active */
    constexpr VFPFlags vfp_flags = VFPFlags::spatial_advection |
                                   VFPFlags::momentum | VFPFlags::collision |
                                   VFPFlags::rotation | VFPFlags::source;
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
      vector_value(const dealii::Point<dim> &point,
                   dealii::Vector<double>   &f) const override
      {
        AssertDimension(f.size(), this->n_components);
        static_cast<void>(point); // suppress compiler warning

        for (unsigned int i = 0; i < f.size(); ++i)
          {
            /** [Initial value] */
            // No initial value
            f[i] = 0.;
            /** [Initial value] */
          }
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
      bc_vector_value_list(const std::vector<dealii::Point<dim>> &points,
                           const unsigned int                     boundary_id,
                           std::vector<dealii::Vector<double>> &bc_values) const
      {
        AssertDimension(points.size(), bc_values.size());
        AssertDimension(bc_values[0].size(), this->n_components);

        for (unsigned int i = 0; i < points.size(); ++i)
          {
            // !!!EDIT HERE!!!
            if (boundary_id == 0)
              bc_values[i][0] =
                std::sqrt(4 * M_PI) * std::exp(-points[i][0] * points[i][0]);
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
      value_list(const std::vector<dealii::Point<dim>> &points,
                 std::vector<double>                   &scattering_frequencies,
                 const unsigned int component = 0) const override
      {
        AssertDimension(scattering_frequencies.size(), points.size());
        static_cast<void>(component); // suppress compiler warning

        for (unsigned int q_index = 0; q_index < points.size(); ++q_index)
          {
            /** [Scattering frequency] */
            // Bohm limit
            scattering_frequencies[q_index] =
              prm.nu0 * prm.B0 * std::exp(-points[q_index][1]);
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
      vector_value(const dealii::Point<dim> &point,
                   dealii::Vector<double>   &source_values) const override
      {
        AssertDimension(source_values.size(), this->n_components);
        static_cast<void>(point); // suppress compiler warning

        for (unsigned int i = 0; i < source_values.size(); ++i)
          {
            /** [Source] */
            if (i == 0)
              {
                const double p = std::exp(point[1]);
                const double x = point[0];

                // S_000 = sqrt(4 pi) * S
                source_values[0] =
                  prm.Q /
                  (4 * std::pow(M_PI, 1.5) * prm.sig_p * prm.sig_x * p * p) *
                  std::exp(-(p - prm.p_inj) * (p - prm.p_inj) /
                           (2. * prm.sig_p * prm.sig_p)) *
                  std::exp(-(x - prm.x_inj) * (x - prm.x_inj) /
                           (2. * prm.sig_x * prm.sig_x));
              }
            else
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
        magnetic_field[0] = prm.B0; // B_x
        magnetic_field[1] = 0.;     // B_y
        magnetic_field[2] = 0.;     // B_z
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
        // u(x) = u_sh/2r * ((1-r)*tanh(x/x_s) + (1+r))

        // u_x
        velocity[0] =
          prm.u_sh / (2 * prm.compression_ratio) *
          ((1 - prm.compression_ratio) * std::tanh(point[0] / prm.shock_width) +
           (1 + prm.compression_ratio));
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

        for (unsigned int q_index = 0; q_index < points.size(); ++q_index)
          {
            /** [Background velocity divergence] */
            // u(x) = u_sh/2r * ((1-r)*tanh(x/x_s) + (1+r))
            // => d/dx u(x) = u_sh/2r 1/x_s (1-r) (1-tanh(x/x_s)^2)

            // div u
            divergence[q_index] =
              prm.u_sh / (2 * prm.compression_ratio) *
              (1 - prm.compression_ratio) / prm.shock_width *
              (1 - std::tanh(points[q_index][0] / prm.shock_width) *
                     std::tanh(points[q_index][0] / prm.shock_width));
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

        for (unsigned int q_index = 0; q_index < points.size(); ++q_index)
          {
            /** [Background velocity material derivative] */
            // u(x) = u_sh/2r * ((1-r)*tanh(x/x_s) + (1+r))
            // => D/Dt u(x) = d/dt u(x) + u d/dx u(x)
            //     = (u_sh/2r)^2 * ((1-r)*tanh(x/x_s) + (1+r)) *
            //        1/x_s (1-r) (1-tanh(x/x_s)^2)

            // D/Dt u_x
            material_derivatives[q_index][0] =
              prm.u_sh * prm.u_sh /
              (4 * prm.compression_ratio * prm.compression_ratio) /
              prm.shock_width * (1 - prm.compression_ratio) *
              ((1 - prm.compression_ratio) *
                 std::tanh(points[q_index][0] / prm.shock_width) +
               (1 + prm.compression_ratio)) *
              (1 - std::tanh(points[q_index][0] / prm.shock_width) *
                     std::tanh(points[q_index][0] / prm.shock_width));
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
            //  u(x) = u_sh/2r * ((1-r)*tanh(x/x_s) + (1+r))
            //  => u_00 = du/dx = u_sh/2r 1/x_s (1-r) (1-tanh(x/x_s)^2)

            // \partial u_x / \partial x
            jacobians[q_index][0][0] =
              prm.u_sh / (2 * prm.compression_ratio) *
              (1 - prm.compression_ratio) / prm.shock_width *
              (1 - std::tanh(points[q_index][0] / prm.shock_width) *
                     std::tanh(points[q_index][0] / prm.shock_width));
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
