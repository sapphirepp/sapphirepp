#ifndef CONFIG_H
#define CONFIG_H

#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/base/function.h>
#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/point.h>

#include <deal.II/lac/vector.h>

#include <cmath>
#include <ostream>
#include <vector>

#include "sapphire-logstream.h"
#include "vfp-flags.h"

namespace Sapphire
{
  class PhysicalProperties
  {
  public:
    PhysicalProperties() = default;

    void
    declare_parameters(dealii::ParameterHandler &prm)
    {
      dealii::LogStream::Prefix p("Physical properties", saplog);
      saplog << "Declaring parameters" << std::endl;
      prm.enter_subsection("Physical properties");

      prm.leave_subsection();
    };

    void
    parse_parameters(dealii::ParameterHandler &prm)
    {
      dealii::LogStream::Prefix p("PhysicalProperties", saplog);
      saplog << "Parsing parameters" << std::endl;
      prm.enter_subsection("Physical properties");

      prm.leave_subsection();
    };

    const double B0 = 2 * M_PI;
    const double R0 = 1.;
  };


  namespace VFP
  {
    static constexpr VFPFlags vfp_flags = VFPFlags::spatial_advection |
                                          VFPFlags::magnetic |
                                          VFPFlags::time_independent_fields;

    static constexpr int dim_configuration_space = 2;

    struct ParticleProperties
    {
      const double mass   = 1.;
      const double charge = 1.;
    };

    struct TransportOnly
    {
      const double gamma    = 2.0;
      const double velocity = std::sqrt(1 - 1 / std::pow(gamma, 2));
    };

    template <int dim>
    class InitialValueFunction : public dealii::Function<dim>
    {
    public:
      InitialValueFunction(const PhysicalProperties &physical_properties,
                           int                       exp_order)
        : dealii::Function<dim>((exp_order + 1) * (exp_order + 1))
        , R0(physical_properties.R0)
      {}

      void
      vector_value(const dealii::Point<dim> &p,
                   dealii::Vector<double>   &f) const override
      {
        Assert(dim <= 3, dealii::ExcNotImplemented());
        Assert(f.size() == InitialValueFunction<dim>::n_components,
               dealii::ExcDimensionMismatch(
                 f.size(), InitialValueFunction<dim>::n_components));

        // Constant disc
        const double radius      = p.norm();
        const double shock_width = 0.2;
        f[0] = 0.5 - 0.5 * std::tanh((radius - R0) / shock_width);
      }

      const double R0;
    };

    template <int dim>
    class ScatteringFrequency : public dealii::Function<dim>
    {
    public:
      ScatteringFrequency(const PhysicalProperties &physical_properties)
        : dealii::Function<dim>(1)
      {
        (void)physical_properties;
      }

      void
      value_list(const std::vector<dealii::Point<dim>> &points,
                 std::vector<double>                   &scattering_frequencies,
                 const unsigned int component = 0) const override
      {
        Assert(scattering_frequencies.size() == points.size(),
               dealii::ExcDimensionMismatch(scattering_frequencies.size(),
                                            points.size()));
        static_cast<void>(component); // suppress compiler warning

        std::fill(scattering_frequencies.begin(),
                  scattering_frequencies.end(),
                  0.0);
      }
    };

    template <int dim>
    class Source : public dealii::Function<dim>
    {
    public:
      Source(const PhysicalProperties &physical_properties,
             unsigned int              exp_order)
        : dealii::Function<dim>((exp_order + 1) * (exp_order + 1))
      {
        (void)physical_properties;
      }

      void
      vector_value(const dealii::Point<dim> &p,
                   dealii::Vector<double>   &values) const override
      {
        Assert(values.size() == Source<dim>::n_components,
               dealii::ExcDimensionMismatch(values.size(),
                                            Source<dim>::n_components));
        static_cast<void>(p); // suppress compiler warning

        values[0] = 0.; // unused
      }
    };

    template <int dim>
    class MagneticField : public dealii::Function<dim>
    {
    public:
      // set n_components
      MagneticField(const PhysicalProperties &physical_properties)
        : dealii::Function<dim>(3)
        , B0(physical_properties.B0)
      {}

      void
      vector_value(const dealii::Point<dim> &point,
                   dealii::Vector<double>   &magnetic_field) const override
      {
        Assert(magnetic_field.size() == MagneticField<dim>::n_components,
               dealii::ExcDimensionMismatch(magnetic_field.size(),
                                            MagneticField<dim>::n_components));
        static_cast<void>(point); // suppress compiler warning

        // constant magnetic field
        magnetic_field[0] = 0.;
        magnetic_field[1] = 0.;
        magnetic_field[2] = B0;
      }

      const double B0;
    };

    template <int dim>
    class BackgroundVelocityField : public dealii::Function<dim>
    {
    public:
      BackgroundVelocityField(const PhysicalProperties &physical_properties)
        : dealii::Function<dim>(3)
      {
        (void)physical_properties;
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
  } // namespace VFP
} // namespace Sapphire
#endif
