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
 * @file examples/vfp/parallel-shock/config.h
 * @author Florian Schulze (florian.schulze@mpi-hd.mpg.de)
 * @brief Implement physical setup for parallel-shock example
 */

#ifndef CONFIG_H
#define CONFIG_H

#include <deal.II/base/exceptions.h>
#include <deal.II/base/function.h>
#include <deal.II/base/function_lib.h>
#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/point.h>
#include <deal.II/base/table.h>

#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/vector.h>

#include <cmath>
#include <filesystem>
#include <fstream>
#include <vector>

#include "pde-system.h"
#include "sapphirepp-logstream.h"
#include "tools.h"
#include "vfp-flags.h"



namespace sapphirepp
{
  class PhysicalParameters
  {
  public:
    static constexpr unsigned int dim_cs = 2;
    /** [Define runtime parameter] */
    double B0;
    double nu0;
    /** [Define runtime parameter] */
    const unsigned int                            N = 500; // N entries
    std::array<unsigned int, dim_cs>              n_subintervals;
    std::array<std::pair<double, double>, dim_cs> interval_endpoints;
    std::filesystem::path                         input_path;



    PhysicalParameters() = default;



    void
    declare_parameters(dealii::ParameterHandler &prm)
    {
      dealii::LogStream::Prefix pre1("Startup", saplog);
      dealii::LogStream::Prefix pre2("Physical parameters", saplog);
      saplog << "Declaring parameters" << std::endl;
      prm.enter_subsection("Physical parameters");

      /** [Declare runtime parameter] */
      prm.declare_entry("B0",
                        "1.",
                        dealii::Patterns::Double(),
                        "The magnetic field strength upstream.");
      prm.declare_entry("nu0",
                        "0.1",
                        dealii::Patterns::Double(),
                        "The scattering frequency.");
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
      B0  = prm.get_double("B0");
      nu0 = prm.get_double("nu0");
      /** [Parse runtime parameter] */
      prm.leave_subsection();

      input_path = "examples/vfp/magnetic-mirroring";
      prm.enter_subsection("VFP");
      prm.enter_subsection("Mesh");
      {
        dealii::Point<dim_cs> p1;
        dealii::Point<dim_cs> p2;

        // Two diagonally opposite corner points of the grid
        dealii::Patterns::Tools::to_value(prm.get("Point 1"), p1);
        dealii::Patterns::Tools::to_value(prm.get("Point 2"), p2);

        saplog << p1 << p2 << std::endl;

        for (unsigned int d = 0; d < dim_cs; ++d)
          {
            interval_endpoints[d] = {p1[d], p2[d]};
            n_subintervals[d]     = N - 1;
          }
      } // Mesh
      prm.leave_subsection();
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
    constexpr VFPFlags vfp_flags =
      VFPFlags::spatial_advection | VFPFlags::collision | VFPFlags::rotation |
      VFPFlags::time_independent_fields | VFPFlags::source |
      VFPFlags::time_independent_source;
    /** [VFP Flags] */



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
      vector_value(const dealii::Point<dim> &point,
                   dealii::Vector<double>   &bc) const override
      {
        AssertDimension(bc.size(), this->n_components);
        static_cast<void>(point); // suppress compiler warning

        // constant isotropic part only
        bc[0] = 0.1;
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
        /** [Scattering frequency] */
        for (unsigned int q_index = 0; q_index < points.size(); ++q_index)
          {
            const double B_0   = 50;
            const double Gamma = 1070376;
            const double nu_b  = B_0 / (2 * M_PI * Gamma);

            // Constant scattering frequency
            scattering_frequencies[q_index] = prm.nu0 * nu_b;
          }
        /** [Scattering frequency] */
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
                source_values[i] = 0.;
                // shifted Gaussian in x and p
                const double x_inj = 1600000 / 2.;
                const double sig_x = 1600000 / 250. * 3.;
                const double Q     = 0.1;
                const double x     = point[0] - x_inj;

                // s_000 = sqrt(4 pi) * s
                source_values[0] = Q / (std::sqrt(M_PI) * sig_x) *
                                   std::exp(-x * x / (2. * sig_x * sig_x));
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
      static constexpr unsigned int dim_cs = PhysicalParameters::dim_cs;

      MagneticField(const PhysicalParameters &physical_parameters)
        : dealii::Function<dim>(3)
        , prm{physical_parameters}
        , table_Bx{Utils::Tools::read_csv(prm.input_path / "AcsvBx.csv",
                                          prm.N,
                                          prm.N,
                                          true,
                                          false,
                                          true)}
        , table_By{Utils::Tools::read_csv(prm.input_path / "AcsvBy.csv",
                                          prm.N,
                                          prm.N,
                                          true,
                                          false,
                                          true)}
        // , Bx(prm.interval_endpoints, prm.n_subintervals, std::move(table_Bx))
        , Bx(prm.interval_endpoints,
             prm.n_subintervals,
             std::move(Utils::Tools::read_csv(prm.input_path / "AcsvBx.csv",
                                              prm.N,
                                              prm.N,
                                              true,
                                              false,
                                              true)))
        // , Bx(std::move(prm.interval_endpoints),
        // std::move(prm.n_subintervals), std::move(table_Bx))
        // TODO: Debug to see if this is using the correct move constructor
        , By(prm.interval_endpoints, prm.n_subintervals, std::move(table_By))
        , Az(std::move(prm.interval_endpoints),
             std::move(prm.n_subintervals),
             std::move(Utils::Tools::read_csv(prm.input_path / "AcsvAz.csv",
                                              prm.N,
                                              prm.N,
                                              true,
                                              false,
                                              true)))
      {}



      void
      vector_value(const dealii::Point<dim> &point,
                   dealii::Vector<double>   &magnetic_field) const override
      {
        AssertDimension(magnetic_field.size(), this->n_components);
        static_cast<void>(point); // suppress compiler warning


        const double B_0 = 50;
        // const double rescale = 1e-10;

        /** [Magnetic field] */
        // magnetic_field[0] = prm.B0; // B_x
        magnetic_field[0] = B_0; // B_x
        magnetic_field[1] = 0;   // B_y
        magnetic_field[2] = 0.;  // B_z
        // magnetic_field[2] = Az.value(p); // B_z
        /** [Magnetic field] */
      }



    private:
      const PhysicalParameters                                     prm;
      const dealii::Table<dim_cs, double>                          table_Bx;
      const dealii::Table<dim_cs, double>                          table_By;
      const dealii::Functions::InterpolatedUniformGridData<dim_cs> Bx;
      const dealii::Functions::InterpolatedUniformGridData<dim_cs> By;
      const dealii::Functions::InterpolatedUniformGridData<dim_cs> Az;
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
        velocity[0] = 0.1; // u_x
        velocity[1] = 0.;  // u_y
        velocity[2] = 0.;  // u_z
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
            // div u
            divergence[q_index] = 0.;
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



    template <unsigned int dim>
    class InitialValueFunction : public dealii::Function<dim>
    {
    public:
      InitialValueFunction(const PhysicalParameters &physical_parameters,
                           const unsigned int        system_size)
        : dealii::Function<dim>(system_size)
        , prm{physical_parameters}
        , lms_indices{PDESystem::create_lms_indices(system_size)}
        , magnetic_field(prm)
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

        // // DEBUG: Output magnetic field for visualization
        // dealii::Vector<double> magnetic_field_value(3);
        // magnetic_field.vector_value(point, magnetic_field_value);
        // for (unsigned int i = 0; i < 3; ++i)
        //   f[i] = magnetic_field_value[i];
      }



    private:
      const PhysicalParameters                       prm;
      const std::vector<std::array<unsigned int, 3>> lms_indices;
      const MagneticField<dim>                       magnetic_field;
    };
  } // namespace VFP
} // namespace sapphirepp
#endif
