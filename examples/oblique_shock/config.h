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

      prm.declare_entry("u_sh",
                        "0.1",
                        dealii::Patterns::Double(),
                        "The shock velocity.");
      prm.declare_entry("B_0",
                        "1.",
                        dealii::Patterns::Double(),
                        "The magnetic field strength upstream.");
      prm.declare_entry(
        "theta",
        "0.",
        dealii::Patterns::Double(),
        "The angle between the magnetic field and the shock normal.");
      prm.declare_entry("compression ratio",
                        "4.",
                        dealii::Patterns::Double(),
                        "The compression ratio of the shock.");
      prm.declare_entry("shock width",
                        "0.04",
                        dealii::Patterns::Double(),
                        "The width of the shock.");
      prm.declare_entry(
        "alpha",
        "0.1",
        dealii::Patterns::Double(),
        "The ratio of the gyrofrequency to the scattering frequency.");

      prm.enter_subsection("Mesh");
      prm.declare_entry("p_min",
                        "0.1",
                        dealii::Patterns::Double(),
                        "The minimum momentum.");
      prm.declare_entry("p_max",
                        "100",
                        dealii::Patterns::Double(),
                        "The maximum momentum.");
      prm.declare_entry("n_cells_p",
                        "64",
                        dealii::Patterns::Integer(1),
                        "The number of cells in momentum direction.");
      prm.declare_entry("n_cells_downstream",
                        "64",
                        dealii::Patterns::Integer(1),
                        "The number of cells downstream of the shock.");
      prm.declare_entry("n_cells_upstream",
                        "64",
                        dealii::Patterns::Integer(1),
                        "The number of cells upstream of the shock.");
      prm.declare_entry("n_cells_shock",
                        "8",
                        dealii::Patterns::Integer(1),
                        "The number of cells in the shock.");
      prm.declare_entry("scaling delta_x",
                        "1.5",
                        dealii::Patterns::Double(0),
                        "The scaling factor for the cell sizes in the shock.");
      prm.leave_subsection();

      prm.enter_subsection("Source");
      prm.declare_entry("p_inj",
                        "1.",
                        dealii::Patterns::Double(),
                        "The injection momentum.");
      prm.declare_entry("x_inj",
                        "0.0",
                        dealii::Patterns::Double(),
                        "The injection position.");
      prm.declare_entry("sig_x",
                        "0.1",
                        dealii::Patterns::Double(),
                        "The width of the source in configuration space.");
      prm.declare_entry("sig_p",
                        "0.1",
                        dealii::Patterns::Double(),
                        "The width of the source in momentum space.");
      prm.declare_entry("normalization",
                        "1.",
                        dealii::Patterns::Double(),
                        "The normalization of the source.");
      prm.leave_subsection();

      prm.leave_subsection();
    };

    void
    parse_parameters(dealii::ParameterHandler &prm)
    {
      dealii::LogStream::Prefix p("PhysicalProperties", saplog);
      saplog << "Parsing parameters" << std::endl;
      prm.enter_subsection("Physical properties");

      u_sh              = prm.get_double("u_sh");
      B_0               = prm.get_double("B_0");
      theta             = prm.get_double("theta");
      compression_ratio = prm.get_double("compression ratio");
      shock_width       = prm.get_double("shock width");
      alpha             = prm.get_double("alpha");

      prm.enter_subsection("Mesh");
      p_min              = prm.get_double("p_min");
      p_max              = prm.get_double("p_max");
      n_cells_p          = prm.get_integer("n_cells_p");
      n_cells_downstream = prm.get_integer("n_cells_downstream");
      n_cells_upstream   = prm.get_integer("n_cells_upstream");
      n_cells_shock      = prm.get_integer("n_cells_shock");
      scaling_delta_x    = prm.get_double("scaling delta_x");
      prm.leave_subsection();

      prm.enter_subsection("Source");
      p_inj         = prm.get_double("p_inj");
      x_inj         = prm.get_double("x_inj");
      sig_x         = prm.get_double("sig_x");
      sig_p         = prm.get_double("sig_p");
      normalization = prm.get_double("normalization");
      prm.leave_subsection();

      prm.leave_subsection();
    };

    double u_sh;
    double B_0;
    double theta;
    double compression_ratio;
    double shock_width;
    double alpha;

    // Grid
    double       p_min;
    double       p_max;
    unsigned int n_cells_p;
    unsigned int n_cells_downstream;
    unsigned int n_cells_upstream;
    unsigned int n_cells_shock;
    double       scaling_delta_x;

    // Source
    double p_inj;
    double x_inj;
    double sig_x;
    double sig_p;
    double normalization;
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

    struct TransportOnly // unused
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
        , B_0(physical_properties.B_0)
        , theta(physical_properties.theta)
        , compression_ratio(physical_properties.compression_ratio)
        , shock_width(physical_properties.shock_width)
        , alpha(physical_properties.alpha)
      {}

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
      const double B_0;
      const double theta;
      const double compression_ratio;
      const double shock_width;
      const double alpha;
    };

    template <int dim>
    class Source : public dealii::Function<dim>
    {
    public:
      Source(const PhysicalProperties &physical_properties,
             unsigned int              exp_order)
        : dealii::Function<dim>((exp_order + 1) * (exp_order + 1))
        , p_inj(physical_properties.p_inj)
        , x_inj(physical_properties.x_inj)
        , sig_x(physical_properties.sig_x)
        , sig_p(physical_properties.sig_p)
        , normalization(physical_properties.normalization)
      {}

      void
      vector_value(const dealii::Point<dim> &point,
                   dealii::Vector<double>   &values) const override
      {
        // shifted Gausian in x and p
        const double x = point[0] - x_inj;
        const double p = std::exp(point[1]) - p_inj;

        values[0] = normalization * std::exp(-x * x / (2. * sig_x)) *
                    std::exp(-p * p / (2. * sig_p));
      }

    private:
      const double p_inj;
      const double x_inj;
      const double sig_x;
      const double sig_p;
      const double normalization;
      // const double normalization = 0.1 / (2 * M_PI * std::sqrt(sig_p *
      // sig_x)); const double normalization = 1.; const double p_inj = 1.;
      // const double sig_x         = 1.;
      // const double sig_p         = 1.;
    };

    template <int dim>
    class MagneticField : public dealii::Function<dim>
    {
    public:
      MagneticField(const PhysicalProperties &physical_properties)
        : dealii::Function<dim>(3)
        , B_0(physical_properties.B_0)
        , theta(physical_properties.theta)
        , compression_ratio(physical_properties.compression_ratio)
        , shock_width(physical_properties.shock_width)
      {}

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
      const double B_0;
      const double theta;
      const double compression_ratio;
      const double shock_width;
    };

    template <int dim>
    class BackgroundVelocityField : public dealii::Function<dim>
    {
    public:
      BackgroundVelocityField(const PhysicalProperties &physical_properties)
        : dealii::Function<dim>(3)
        , u_sh(physical_properties.u_sh)
        , compression_ratio(physical_properties.compression_ratio)
        , shock_width(physical_properties.shock_width)
      {}

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
      const double u_sh;
      const double compression_ratio;
      const double shock_width;
    };
  } // namespace VFP
} // namespace Sapphire
#endif
