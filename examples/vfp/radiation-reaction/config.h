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
 * @file examples/vfp/radiation-reaction/config.h
 * @author Nils Schween (nils.schween@mpi-hd.mpg.de)
 * @author Harsh Goyal (harsh.goyal@mpi-hd.mpg.de)
 * @brief Implement physical setup for radiation-reaction example
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
#include <numbers>
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
    double u_sh              = 0.1;
    double B0                = 1.;
    double compression_ratio = 4.;
    double shock_width       = 0.04;
    double nu0               = 0.1;

    // Source
    double Q     = 0.1;
    double p_inj = 2.;
    double x_inj = 0.;
    double sig_p = 0.125;
    double sig_x = 0.125;
    /** [Define runtime parameter] */



    PhysicalParameters() = default;



    void
    declare_parameters(dealii::ParameterHandler &prm)
    {
      dealii::LogStream::Prefix prefix_startup("Startup", saplog);
      dealii::LogStream::Prefix prefix("PhysicalParameters", saplog);
      saplog << "Declaring parameters" << std::endl;
      prm.enter_subsection("Physical parameters");

      /** [Declare runtime parameter] */
      prm.add_parameter("u_sh", u_sh, "The shock velocity.");
      prm.add_parameter("B0", B0, "The magnetic field strength upstream.");
      prm.add_parameter("compression ratio",
                        compression_ratio,
                        "The compression ratio of the shock.");
      prm.add_parameter("shock width", shock_width, "The width of the shock.");
      prm.add_parameter("nu0", nu0, "The scattering frequency.");

      // Source
      prm.add_parameter("Q", Q, "The injection rate.");
      prm.add_parameter("p_inj", p_inj, "The injection momentum.");
      prm.add_parameter("x_inj", x_inj, "The injection position.");
      prm.add_parameter("sig_p",
                        sig_p,
                        "The width of the source in momentum space.");
      prm.add_parameter("sig_x",
                        sig_x,
                        "The width of the source in configuration space.");
      /** [Declare runtime parameter] */

      prm.leave_subsection();
    }



    void
    parse_parameters(dealii::ParameterHandler &prm)
    {
      dealii::LogStream::Prefix prefix_startup("Startup", saplog);
      dealii::LogStream::Prefix prefix("PhysicalParameters", saplog);
      saplog << "Parsing parameters" << std::endl;
      prm.enter_subsection("Physical parameters");

      /** [Parse runtime parameter] */
      // Parameters are automatically parsed by add_parameter()
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
    constexpr VFPFlags vfp_flags = VFPFlags::time_evolution |               //
                                   VFPFlags::momentum |                     //
                                   VFPFlags::radiation_reaction |           //
                                   VFPFlags::spatial_advection |            //
                                   VFPFlags::scaled_distribution_function | //
                                   VFPFlags::rotation |                     //
                                   VFPFlags::time_independent_fields;
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

        /** [Initial value] */
        double p     = std::exp(point[0]) / 1e7;
        double p_max = (0.5 * 1e8) / 1e7;
        f[0]         = 0; // f_000
        f[2]         = 0; // f_100
        f[1] = 0.5 * 1 / std::sqrt(6.0) * std::pow(p, 3) * std::pow(p, -4) *
               std::exp(-p / p_max); // f_110
        f[3] = 0.5 * 1 / std::sqrt(6.0) * std::pow(p, 3) * std::pow(p, -4) *
               std::exp(-p / p_max); // f_111
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
            scattering_frequencies[q_index] = 100.;
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
            // Define parameters
            const double p0    = std::log(1e7);
            const double Q     = 0.1;
            const double sig_p = 0.1;
            const double sig_x = 0.001;

            // Get x and p value
            const double x = point[0];
            const double p = std::exp(point[1]); // point[1] = log(p)

            // Convert i -> l, m, s
            const unsigned int l = lms_indices[i][0];
            const unsigned int m = lms_indices[i][1];
            const unsigned int s = lms_indices[i][2];

            // isotropic part
            if (l == 0 && m == 0 && s == 0)
              {
                // shifted Gaussian in x and p
                // s_000 = sqrt(4 pi) * s
                source_values[i] =
                  Q / (std::sqrt(std::numbers::pi) * sig_p * sig_x) *
                  std::exp(-(p - p0) * (p - p0) / (2. * sig_p * sig_p)) *
                  std::exp(-x * x / (2. * sig_x * sig_x));
                // const double mu    = p0;     //
                // const double sigma = sig_p;            // std dev in ln p

                // const double z      = (std::log(p) - mu) / sigma;

                // source_values[i] = Q * (1.0 / (p * sigma * std::sqrt(2.0 *
                // std::numbers::pi))) * std::exp(-0.5 * z * z);
              }
            // vanishing anisotropic part
            else
              {
                source_values[i] = 0.;
              }
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
        magnetic_field[0] = 1e9; // B_x
        magnetic_field[1] = 0;   // B_y
        magnetic_field[2] = 0;   // B_z = 1 µG
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
        // u_x
        velocity[0] = 0.;
        // u_sh / (2 * r) * ((1 - r) * std::tanh(point[0] / d_sh) + (1 + r));
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
            // u_sh / (2 * r) * (1 - r) / d_sh *
            // (1 - std::tanh(x / d_sh) * std::tanh(x / d_sh));
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
            // D/Dt u_x
            material_derivatives[q_index][0] = 0.;
            // // u_sh * u_sh / (4 * r * r) / d_sh * (1 - r) *
            // ((1 - r) * std::tanh(x / d_sh) + (1 + r)) *
            // (1 - std::tanh(x / d_sh) * std::tanh(x / d_sh));
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
            // \partial u_x / \partial x
            jacobians[q_index][0][0] = 0.;
            // u_sh / (2 * r) * (1 - r) / d_sh *
            // (1 - std::tanh(x / d_sh) * std::tanh(x / d_sh));
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
