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

/** [Includes] */
#ifndef CONFIG_H
#  define CONFIG_H

#  include <deal.II/base/exceptions.h>
#  include <deal.II/base/function.h>
#  include <deal.II/base/parameter_handler.h>
#  include <deal.II/base/point.h>

#  include <deal.II/lac/full_matrix.h>
#  include <deal.II/lac/vector.h>

#  include <cmath>
#  include <vector>

#  include "pde-system.h"
#  include "sapphirepp-logstream.h"
#  include "vfp-flags.h"
/** [Includes] */



/** [Namespace sapphirepp] */
namespace sapphirepp
{
  /** [Namespace sapphirepp] */
  /** [PhysicalParameters] */
  class PhysicalParameters
  {
  public:
    // !!!EDIT HERE!!!
    double nu;
    double f0;

    PhysicalParameters() = default;
    /** [PhysicalParameters] */
    /** [Declare parameters] */
    void
    declare_parameters(dealii::ParameterHandler &prm)
    {
      dealii::LogStream::Prefix pre1("Startup", saplog);
      dealii::LogStream::Prefix pre2("Physical parameters", saplog);
      saplog << "Declaring parameters" << std::endl;
      prm.enter_subsection("Physical parameters");

      // !!!EDIT HERE!!!
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
      dealii::LogStream::Prefix pre1("Startup", saplog);
      dealii::LogStream::Prefix pre2("PhysicalParameters", saplog);
      saplog << "Parsing parameters" << std::endl;
      prm.enter_subsection("Physical parameters");

      // !!!EDIT HERE!!!
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
    // !!!EDIT HERE!!!
    /** Specify reduced phase space dimension \f$ (\mathbf{x}, p) \f$ */
    constexpr unsigned int dimension = 1;
    /** [Dimension] */



    /** [VFP Flags] */
    // !!!EDIT HERE!!!
    /** Specify which terms of the VFP equation should be active */
    constexpr VFPFlags vfp_flags =
      VFPFlags::time_evolution | VFPFlags::collision;
    /** [VFP Flags] */



    /** [InitialValueFunction constructor] */
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
      /** [InitialValueFunction constructor] */



      /** [InitialValueFunction value] */
      void
      vector_value(const dealii::Point<dim> &point,
                   dealii::Vector<double>   &f) const override
      {
        AssertDimension(f.size(), this->n_components);
        static_cast<void>(point); // suppress compiler warning

        for (unsigned int i = 0; i < f.size(); ++i)
          {
            // !!!EDIT HERE!!!
            f[i] = prm.f0;
          }
      }



    private:
      const PhysicalParameters                       prm;
      const std::vector<std::array<unsigned int, 3>> lms_indices;
    };
    /** [InitialValueFunction value] */



    /** [BoundaryValueFunction constructor] */
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
      /** [BoundaryValueFunction constructor] */



      /** [BoundaryValueFunction value] */
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
            // !!!EDIT HERE!!!
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
          }
      }



    private:
      const PhysicalParameters                       prm;
      const std::vector<std::array<unsigned int, 3>> lms_indices;
    };
    /** [BoundaryValueFunction value] */



    /** [ScatteringFrequency constructor] */
    template <unsigned int dim>
    class ScatteringFrequency : public dealii::Function<dim>
    {
    public:
      ScatteringFrequency(const PhysicalParameters &physical_parameters)
        : dealii::Function<dim>(1)
        , prm{physical_parameters}
      {}
      /** [ScatteringFrequency constructor] */



      /** [ScatteringFrequency value] */
      void
      value_list(const std::vector<dealii::Point<dim>> &points,
                 std::vector<double>                   &scattering_frequencies,
                 const unsigned int component = 0) const override
      {
        AssertDimension(scattering_frequencies.size(), points.size());
        static_cast<void>(component); // suppress compiler warning

        for (unsigned int q_index = 0; q_index < points.size(); ++q_index)
          {
            // !!!EDIT HERE!!!
            scattering_frequencies[q_index] = prm.nu;
          }
      }



    private:
      const PhysicalParameters prm;
    };
    /** [ScatteringFrequency value] */



    /** [Source] */
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
            // !!!EDIT HERE!!!
            source_values[i] = 0.;
          }
      }



    private:
      const PhysicalParameters                       prm;
      const std::vector<std::array<unsigned int, 3>> lms_indices;
    };
    /** [Source] */



    /** [MagneticField] */
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

        // !!!EDIT HERE!!!
        magnetic_field[0] = 0.; // B_x
        magnetic_field[1] = 0.; // B_y
        magnetic_field[2] = 0.; // B_z
      }



    private:
      const PhysicalParameters prm;
    };
    /** [MagneticField] */



    /** [BackgroundVelocityField] */
    template <unsigned int dim>
    class BackgroundVelocityField : public dealii::Function<dim>
    {
    public:
      BackgroundVelocityField(const PhysicalParameters &physical_parameters)
        : dealii::Function<dim>(3)
        , prm{physical_parameters}
      {}
      /** [BackgroundVelocityField] */



      /** [BackgroundVelocityField value] */
      void
      vector_value(const dealii::Point<dim> &point,
                   dealii::Vector<double>   &velocity) const override
      {
        AssertDimension(velocity.size(), this->n_components);
        static_cast<void>(point); // suppress compiler warning

        // !!!EDIT HERE!!!
        velocity[0] = 0.; // u_x
        velocity[1] = 0.; // u_y
        velocity[2] = 0.; // u_z
      }
      /** [BackgroundVelocityField value] */



      /** [BackgroundVelocityField divergence] */
      void
      divergence_list(const std::vector<dealii::Point<dim>> &points,
                      std::vector<double>                   &divergence) const
      {
        AssertDimension(divergence.size(), points.size());
        static_cast<void>(points); // suppress compiler warning

        for (unsigned int q_index = 0; q_index < points.size(); ++q_index)
          {
            // !!!EDIT HERE!!!
            divergence[q_index] = 0.; // div u
          }
      }
      /** [BackgroundVelocityField divergence] */



      /** [BackgroundVelocityField material derivative] */
      void
      material_derivative_list(
        const std::vector<dealii::Point<dim>> &points,
        std::vector<dealii::Vector<double>>   &material_derivatives) const
      {
        AssertDimension(material_derivatives.size(), points.size());
        AssertDimension(material_derivatives[0].size(), this->n_components);

        for (unsigned int q_index = 0; q_index < points.size(); ++q_index)
          {
            // !!!EDIT HERE!!!
            material_derivatives[q_index][0] = 0.; // D/Dt u_x
            material_derivatives[q_index][1] = 0.; // D/Dt u_y
            material_derivatives[q_index][2] = 0.; // D/Dt u_z
          }
      }
      /** [BackgroundVelocityField material derivative] */



      /** [BackgroundVelocityField Jacobian] */
      void
      jacobian_list(const std::vector<dealii::Point<dim>>   &points,
                    std::vector<dealii::FullMatrix<double>> &jacobians) const
      {
        AssertDimension(jacobians.size(), points.size());
        AssertDimension(jacobians[0].m(), this->n_components);
        AssertDimension(jacobians[0].n(), this->n_components);

        for (unsigned int q_index = 0; q_index < points.size(); ++q_index)
          {
            // !!!EDIT HERE!!!
            jacobians[q_index][0][0] = 0.; // \partial u_x / \partial x
            jacobians[q_index][0][1] = 0.; // \partial u_x / \partial y
            jacobians[q_index][0][2] = 0.; // \partial u_x / \partial z

            jacobians[q_index][1][0] = 0.; // \partial u_y / \partial x
            jacobians[q_index][1][1] = 0.; // \partial u_y / \partial y
            jacobians[q_index][1][2] = 0.; // \partial u_y / \partial z

            jacobians[q_index][2][0] = 0.; // \partial u_z / \partial x
            jacobians[q_index][2][1] = 0.; // \partial u_z / \partial y
            jacobians[q_index][2][2] = 0.; // \partial u_z / \partial z
          }
      }



    private:
      const PhysicalParameters prm;
    };
    /** [BackgroundVelocityField Jacobian] */
    /** [Close namespaces] */
  } // namespace VFP
} // namespace sapphirepp
#endif
/** [Close namespaces] */
