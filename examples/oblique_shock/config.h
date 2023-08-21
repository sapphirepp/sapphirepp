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

#include "parameter-flags.h"
#include "sapphire-logstream.h"

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
  };


  namespace VFP
  {
    static constexpr TermFlags vfp_terms =
      TermFlags::spatial_advection | TermFlags::momentum |
      TermFlags::collision | TermFlags::magnetic | TermFlags::source;

    static constexpr bool logarithmic_p = true;

    static constexpr int dim_configuration_space = 1;

    static constexpr bool time_dependent_fields = false;

    static constexpr bool time_dependent_source = false;

    struct ParticleProperties
    {
      const double mass   = 1.;
      const double charge = 1.;
    };

    struct TransportOnly
    {
      const double gamma    = 3.;
      const double velocity = std::sqrt(1 - 1 / std::pow((gamma * gamma), 2));
    };

    template <int dim>
    class InitialValueFunction : public dealii::Function<dim>
    {
    public:
      InitialValueFunction(const PhysicalProperties &physical_properties,
                           int                       exp_order)
        : dealii::Function<dim>((exp_order + 1) * (exp_order + 1))
      {
        (void)physical_properties;
      }

      void
      vector_value(const dealii::Point<dim> &p,
                   dealii::Vector<double>   &f) const override
      {
        Assert(dim <= 3, dealii::ExcNotImplemented());
        Assert(f.size() == InitialValueFunction<dim>::n_components,
               dealii::ExcDimensionMismatch(
                 f.size(), InitialValueFunction<dim>::n_components));
        (void)p;
        for (unsigned int i = 0; i < f.size(); ++i)
          f[i] = 0;
      }
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
        static_cast<void>(component);
        static_cast<void>(points);

        // Constant scattering frequency
        std::fill(scattering_frequencies.begin(),
                  scattering_frequencies.end(),
                  alpha);

        // TODO: Scattering frequency propotional to omega_g
        // for (unsigned int i = 0; i < points.size(); ++i)
        //   {
        //     // w_g  = c * e * B / p
        //     const double B_x = B_0;
        //     const double B_z =
        //       std::sin(theta) * B_0 * (2 * compression_ratio) /
        //       ((1 - compression_ratio) * std::tanh(points[i][0] /
        //       shock_width) +
        //        (1 + compression_ratio));
        //     const double B = std::sqrt(B_x * B_x + B_z * B_z);
        //     // const double B =
        //     //   B_0 * std::sqrt(
        //     //           1 + std::pow(std::tan(theta), 2) *
        //     //                 std::pow(2 * compression_ratio, 2) /
        //     //                 std::pow((1 - compression_ratio) *
        //     //                              std::tanh(points[i][0] /
        //     //                              shock_width) +
        //     //                            (1 + compression_ratio),
        //     //                          2));
        //     const double w_g = B / std::exp(points[i][1]);

        //     scattering_frequencies[i] = alpha * w_g;
        //   }
      }

    private:
      const double alpha = 0.1;
      const double B_0   = 1.;
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
      vector_value(const dealii::Point<dim> &point,
                   dealii::Vector<double>   &values) const override
      {
        // Gausian in x and p
        // TODO: Maybe shift x injection a bit into upstream: x -> x-0.3
        const double x = point[0] - x_inj;
        const double p = std::exp(point[1]) - p_inj;

        values[0] = normalization * std::exp(-x * x / (2. * sig_x)) *
                    std::exp(-p * p / (2. * sig_p));
      }

    private:
      const double p_inj         = 1.;
      const double x_inj         = 0.5;
      const double sig_x         = (1. / 4.) * (1. / 4.);
      const double sig_p         = (1. / 4.) * (1. / 4.);
      const double normalization = 0.1 / (2 * M_PI * std::sqrt(sig_p * sig_x));
    };

    template <int dim>
    class MagneticField : public dealii::Function<dim>
    {
    public:
      MagneticField(const PhysicalProperties &physical_properties)
        : dealii::Function<dim>(3)
      {
        (void)physical_properties;
      }

      void
      vector_value(const dealii::Point<dim> &point,
                   dealii::Vector<double>   &magnetic_field) const override
      {
        Assert(magnetic_field.size() == MagneticField<dim>::n_components,
               dealii::ExcDimensionMismatch(magnetic_field.size(),
                                            MagneticField<dim>::n_components));

        magnetic_field[0] = std::cos(theta) * B_0;
        magnetic_field[1] = 0.;
        // B_z(x) = B0 * sin(theta) * 2r / ((1-r)*tanh(x/x_s) + (1+r))
        magnetic_field[2] =
          std::sin(theta) * B_0 * (2 * compression_ratio) /
          ((1 - compression_ratio) * std::tanh(point[0] / shock_width) +
           (1 + compression_ratio));
      }

    private:
      const double B_0               = 1.;
      const double theta             = 0.;
      const double compression_ratio = 4.;
      const double shock_width       = 0.04;
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

      // Velocity field
      void
      vector_value(const dealii::Point<dim> &point,
                   dealii::Vector<double>   &velocity) const override
      {
        Assert(velocity.size() == BackgroundVelocityField<dim>::n_components,
               dealii::ExcDimensionMismatch(
                 velocity.size(), BackgroundVelocityField<dim>::n_components));

        // u(x) = u_sh/2r * ((1-r)*tanh(x/x_s) + (1+r))
        velocity[0] =
          u_sh / (2 * compression_ratio) *
          ((1 - compression_ratio) * std::tanh(point[0] / shock_width) +
           (1 + compression_ratio));
        velocity[1] = 0.;
        velocity[2] = 0.;
      }

      // Divergence
      void
      divergence_list(const std::vector<dealii::Point<dim>> &points,
                      std::vector<double>                   &divergence)
      {
        Assert(divergence.size() == points.size(),
               dealii::ExcDimensionMismatch(divergence.size(), points.size()));
        // u(x) = u_sh/2r * ((1-r)*tanh(x/x_s) + (1+r))
        // => d/dx u(x) = u_sh/2r 1/x_s (1-r) (1-tanh(x/x_s)^2)
        for (unsigned int i = 0; i < points.size(); ++i)
          divergence[i] = u_sh / (2 * compression_ratio) *
                          (1 - compression_ratio) / shock_width *
                          (1 - std::tanh(points[i][0] / shock_width) *
                                 std::tanh(points[i][0] / shock_width));
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
        // u(x) = u_sh/2r * ((1-r)*tanh(x/x_s) + (1+r))
        // => D/Dt u(x) = d/dt u(x) + u d/dx u(x)
        //     = (u_sh/2r)^2 * ((1-r)*tanh(x/x_s) + (1+r)) *
        //        1/x_s (1-r) (1-tanh(x/x_s)^2)
        for (unsigned int i = 0; i < points.size(); ++i)
          {
            material_derivatives[i][0] =
              u_sh * u_sh / (4 * compression_ratio * compression_ratio) /
              shock_width * (1 - compression_ratio) *
              ((1 - compression_ratio) * std::tanh(points[i][0] / shock_width) +
               (1 + compression_ratio)) *
              (1 - std::tanh(points[i][0] / shock_width) *
                     std::tanh(points[i][0] / shock_width));
            material_derivatives[i][1] = 0.;
            material_derivatives[i][2] = 0.;
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
        //  u(x) = u_sh/2r * ((1-r)*tanh(x/x_s) + (1+r))
        //  => u_00 = du/dx = u_sh/2r 1/x_s (1-r) (1-tanh(x/x_s)^2)
        for (unsigned int i = 0; i < points.size(); ++i)
          {
            jacobians[i][0][0] = u_sh / (2 * compression_ratio) *
                                 (1 - compression_ratio) / shock_width *
                                 (1 - std::tanh(points[i][0] / shock_width) *
                                        std::tanh(points[i][0] / shock_width));
            jacobians[i][0][1] = 0.;
            jacobians[i][0][2] = 0.;

            jacobians[i][1][0] = 0.;
            jacobians[i][1][1] = 0.;
            jacobians[i][1][2] = 0.;

            jacobians[i][2][0] = 0.;
            jacobians[i][2][1] = 0.;
            jacobians[i][2][2] = 0.;
          }
      }

    private:
      const double u_sh              = 0.1;
      const double compression_ratio = 4.;
      const double shock_width       = 0.04;
    };
  } // namespace VFP
} // namespace Sapphire
#endif
