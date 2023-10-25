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
  };


  namespace VFP
  {
    // Physical setup

    // NOTE: A member variable needs be constexpr to be used as template
    // arguments. But it can only be constexpr if it is static ,i.e. if it is
    // the same for all class instances. If it was not static, it would be
    // determined when constructing an instance, which happens at runtime.

    // compile time settings
    static constexpr TermFlags vfp_terms =
      TermFlags::spatial_advection | TermFlags::source;

    // This variabale controls if p is linear or logarithmic
    static constexpr bool logarithmic_p = true;

    // Deactivating the spatial advection term is equivalent to assuming a
    // homogeneous distribution function (i.e. a distribution function which
    // does not depend on x,y z). In this program this is equivalent to set
    // dimension of the configuration to zero.
    static constexpr int dim_configuration_space = 2;

    // If the background velocity field and the the magnetic field do not
    // depend on time, the time stepping methods can be accelerated a lot: In
    // this case it is not necessary to reassamble the spatial discretisation
    // matrix in every stage of the Runge-Kutta method. Actually it only has
    // to be assembled once at time zero.
    static constexpr bool time_dependent_fields = false;

    // If the source term depends on time, it needs to be projected onto the
    // FEM space in the time stepping methods.
    static constexpr bool time_dependent_source = true;

    // NOTE: All physical quantities are dimensionless. The reference values are
    // defined in the reference-values.h header.
    struct ParticleProperties
    {
      const double mass   = 1.;
      const double charge = 1.;
    };

    struct TransportOnly
    {
      // In the transport-only case (i.e. no dependence on p) the energy of the
      // particles has to be given.
      const double gamma = 3.;
      // Compute the partilces Lorentz factor gamma and its velocity
      // TODO: here is a mistake: it must be gamma^2, not gama^4 ?!
      const double velocity = std::sqrt(1 - 1 / std::pow((gamma * gamma), 2));
    };

    // Initial values
    template <int dim>
    class InitialValueFunction : public dealii::Function<dim>
    {
    public:
      InitialValueFunction(const PhysicalProperties &physical_properties,
                           int                       exp_order)
        : // set the number of components with the constructor of the base class
        dealii::Function<dim>((exp_order + 1) * (exp_order + 1))
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
        // NOTE: The zeroth component of values corresponds to f_000, the first
        // component to f_110 etc.
        //
        // NOTE: If the momentum term is activated the last component
        // of point, i.e. p[dim] is the p coordinate.
        //
        // EXAMPLES
        // DIM = 1
        // constant function
        static_cast<void>(p);
        f[0] = 0.;
        // Gaussian
        // f[0] = 1. * std::exp(-(std::pow(p[0], 2)));

        // Constant disc
        // if (std::abs(p[0]) < 1.) f[0] = 1.;

        // space/momentum dependent distribution
        // f[0] = 1. + 0.2 * p[0];

        // sinusoidal distribution
        // f[0] = std::sin((1. * 3.14159265359) / 2 * p[0]) + 1.;

        // DIM = 2
        // constant
        // static_cast<void>(p);
        // f[0] = 0.;
        // Gaussian
        // f[0] = 1. * std::exp(-((std::pow(p[0], 2) + std::pow(p[1], 2))));

        // Constant disc
        // if (p.norm() <= 1.) f[0] = 1.;

        // Constant cross (e.g. for a rigid rotator)
        // if (std::abs(p[0]) <= 3 && std::abs(p[1]) <= 0.5) f[0] = 1.;
        // if (std::abs(p[0]) <= 0.5 && std::abs(p[1]) <= 3) f[0] = 1.;

        // Initialise f_000, f_110, f_100, f111, ... with a Gaussian
        // std::fill(
        //     f.begin(), f.end(),
        //     1. * std::exp(-((std::pow(p[0] - 0.5, 2) + std::pow(p[1] - 0.5,
        //     2)) /
        //                     0.01)));

        // DIM 3
        // f[0] = 1. * std::exp(-((std::pow(p[0], 2) + std::pow(p[1], 2) +
        //                                std::pow(p[2], 2))));
      }
    };

    // Scattering frequency
    template <int dim>
    class ScatteringFrequency : public dealii::Function<dim>
    {
    public:
      // Set the variable  n_components of the base class
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
        // EXAMPLES:
        // Constant scattering frequency
        std::fill(scattering_frequencies.begin(),
                  scattering_frequencies.end(),
                  0.5);
      }
    };

    // Source term
    template <int dim>
    class Source : public dealii::Function<dim>
    {
    public:
      Source(const PhysicalProperties &physical_properties,
             unsigned int              exp_order)
        : // set n_components
        dealii::Function<dim>((exp_order + 1) * (exp_order + 1))
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
        // The zeroth component of values corresponds to isotropic part of the
        // distribution

        // EXAMPLES:
        // Transport only
        // 1D Gaussian (isotropic)
        // values[0] = 0.01 * std::exp(-(std::pow(p[0], 2)));
        // 2D Gaussian (isotropic)
        // values[0] = .01 * std::exp(-(std::pow(p[0], 2) + std::pow(p[1], 2)));
        // 2D pulsating Gaussian (isotropic, time dependent)
        values[0] = 0.01 * (std::sin(this->get_time()) + 1.) *
                    std::exp(-(std::pow(p[0], 2) + std::pow(p[1], 2)));
      }
    };

    // Magnetic field
    template <int dim>
    class MagneticField : public dealii::Function<dim>
    {
    public:
      // set n_components
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

        // EXAMPLES:
        // constant magnetic field
        static_cast<void>(point); // suppress compiler warning
        magnetic_field[0] = 0.;
        magnetic_field[1] = 0.;
        magnetic_field[2] = 3.;
      }
    };

    // Background velocity field
    template <int dim>
    class BackgroundVelocityField : public dealii::Function<dim>
    {
    public:
      // set n_components
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
        // EXAMPLES:
        // constant velocity field
        static_cast<void>(point); // suppress compiler warning
        velocity[0] = .0;         // u_x
        velocity[1] = .0;         // u_y
        velocity[2] = .0;         // u_z

        // space dependent velocity
        // u_x = 1/5 * x
        // velocity[0] = 1. / 5 * point[0];
        // velocity[1] = .0;
        // velocity[2] = .0;

        // u_x = 0.1 * cos(pi/2 * x)
        // velocity[0] = 0.1 * std::cos(pi / 2 * point[0]);
        // velocity[1] = .0;
        // velocity[2] = .0;

        // time-dependent velocity-field
        // u_x = 1/10 * t
        // velocity[0] = 1. / 10 * this->get_time();
        // velocity[1] = .0;
        // velocity[2] = .0;

        // u_x = 1/2 * sin(t)
        // velocity[0] = 1. / 2 * std::sin(this->get_time());
        // velocity[1] = .0;
        // velocity[2] = .0;

        // time- and space dependent velocity
        // u_x = -1/10 * sin(t) * cos(pi/2 * x)
        // velocity[0] = -0.1 * std::sin(this->get_time()) * std::cos(pi / 2 *
        // point[0]); velocity[1] = .0; velocity[2] = .0;

        // u_x = 1/10 * x * (1 - e^(-t))
        // velocity[0] = 0.1 * point[0] * (1 - std::exp(-this->get_time()));
        // velocity[1] = .0;
        // velocity[2] = .0;

        // rigid rotator
        // velocity[0] = -point[1];
        // velocity[1] = point[0];
        // velocity[2] = 0.;
      }

      // Divergence
      void
      divergence_list(const std::vector<dealii::Point<dim>> &points,
                      std::vector<double>                   &divergence)
      {
        Assert(divergence.size() == points.size(),
               dealii::ExcDimensionMismatch(divergence.size(), points.size()));
        // EXAMPLES:
        // constant velocity
        std::fill(divergence.begin(), divergence.end(), 0.);

        // space-dependent velocity field
        // u_x = 1/5 * x => div u = 1/5
        // std::fill(divergence.begin(), divergence.end(), 1. / 5);

        // u_x = 0.1 * cos(pi/2 * x) => div u = -0.1 * pi/2 * sin(pi/2 * x)
        // for (unsigned int i = 0; i < points.size(); ++i)
        //   divergence[i] = -0.1 * pi / 2 * std::sin(pi / 2 * points[i][0]);

        // time-dependent velocity field
        // u_x = 1/10 * t => div u = 0
        // std::fill(divergence.begin(), divergence.end(), 0.);

        // u_x = 1/2 * sin(t) => div u = 0
        // std::fill(divergence.begin(), divergence.end(), 0.);

        // time- and space-dependent velocity field
        // u_x = -1/10 * sin(t) * cos(pi/2 * x) =>
        // div u = 0.1 * pi/2 * sin(t) * sin(pi/2 * x)
        // for (unsigned int i = 0; i < points.size(); ++i)
        //   divergence[i] = 0.1 * pi / 2 * std::sin(this->get_time()) *
        //                   std::sin(pi / 2 * points[i][0]);

        // u_x = 1/10 * x * (1 - e^(-t))  => div u = 0.1 (1  - e^(-t))
        // std::fill(divergence.begin(), divergence.end(),
        //           0.1 * (1 - std::exp(-this->get_time())));
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
            // EXAMPLES:
            // constant velocity
            material_derivatives[i][0] = 0.; // D/Dt u_x
            material_derivatives[i][1] = 0.; // D/Dt u_y
            material_derivatives[i][2] = 0.; // D/Dt u_z

            // space dependent velocity field
            // u_x = 1/5 * x => D/Dt u_x = 1/25 * x
            // material_derivatives[i][0] = 1. / 25 * points[i][0];
            // material_derivatives[i][1] = 0.;
            // material_derivatives[i][2] = 0.;

            // u_x = 0.1 * cos(pi/2 * x)
            // => D/Dt u_x = 1/100 * pi/2 * cos(pi/2*x) * sin(pi/2 * x)
            // material_derivatives[i][0] = 1. / 100 * pi / 2 *
            //                              std::sin(pi / 2 * points[i][0]) *
            //                              std::cos(pi / 2 * points[i][0]);
            // material_derivatives[i][1] = 0.;
            // material_derivatives[i][2] = 0.;

            // time-dependent velocity field
            // u_x = 1/10 * t => D/Dt u_x = 1/10
            // material_derivatives[i][0] = 1. / 10;
            // material_derivatives[i][1] = 0.;
            // material_derivatives[i][2] = 0.;

            // u_x = 1/2 * sin(t) => D/Dt u_x = 1/2 * cos(t)
            // material_derivatives[i][0] = 1. / 2 * std::cos(this->get_time());
            // material_derivatives[i][1] = 0.;
            // material_derivatives[i][2] = 0.;

            // time- and space-dependent velocity field
            // u_x = -1/10 * sin(t) * cos(pi/2 * x)
            // => D/Dt u_x = -0.1 * cos(pi/2 * x) * (cos(t) + 0.1 * pi/2 *
            // sin(t) * sin(pi/2 * x)) material_derivatives[i][0] =
            //     -0.1 * std::cos((pi / 2 * points[i][0])) *
            //     (std::cos(this->get_time()) +
            //      0.1 * pi / 2 * std::sin(this->get_time()) *
            //          std::sin(this->get_time()) * std::sin(pi / 2 *
            //          points[i][0]));
            // material_derivatives[i][1] = 0.;
            // material_derivatives[i][2] = 0.;

            // u_x = 1/10 * x * (1 - e^(-t))  => D/Dt u_x = 0.1 * x * (e^(-t)  +
            // (1 - e^(-t))^2) material_derivatives[i][0] = 0.1 * points[i][0] *
            //                              (std::exp(-this->get_time()) +
            //                               0.1 * (1 -
            //                               std::exp(-this->get_time())) *
            //                                   (1 -
            //                                   std::exp(-this->get_time())));
            // material_derivatives[i][1] = 0.;
            // material_derivatives[i][2] = 0.;
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
            // EXAMPLES:
            // constant velocity field
            jacobians[i][0][0] = 0.; // \partial u_x / \partial x
            jacobians[i][0][1] = 0.; // \partial u_x / \partial y
            jacobians[i][0][2] = 0.; // \partial u_x / \partial z

            jacobians[i][1][0] = 0.; // \partial u_y / \partial x
            jacobians[i][1][1] = 0.; // \partial u_y / \partial y
            jacobians[i][1][2] = 0.; // \partial u_y / \partial z

            jacobians[i][2][0] = 0.; // \partial u_z / \partial x
            jacobians[i][2][1] = 0.; // \partial u_z / \partial y
            jacobians[i][2][2] = 0.; // \partial u_z / \partial z

            // space dependent velocity field
            // u_x = 1/5 * x
            // jacobians[i][0][0] = 1. / 5;
            // jacobians[i][0][1] = 0.;
            // jacobians[i][0][2] = 0.;

            // jacobians[i][1][0] = 0.;
            // jacobians[i][1][1] = 0.;
            // jacobians[i][1][2] = 0.;

            // jacobians[i][2][0] = 0.;
            // jacobians[i][2][1] = 0.;
            // jacobians[i][2][2] = 0.;

            // u_x = 0.1 * cos(pi/2 * x)
            // jacobians[i][0][0] = 0.1 * pi/2 * std::sin(pi/2 * points[i][0]);
            // jacobians[i][0][1] = 0.;
            // jacobians[i][0][2] = 0.;

            // jacobians[i][1][0] = 0.;
            // jacobians[i][1][1] = 0.;
            // jacobians[i][1][2] = 0.;

            // jacobians[i][2][0] = 0.;
            // jacobians[i][2][1] = 0.;
            // jacobians[i][2][2] = 0.;

            // time-dependent velocity field
            // u_x = 1/10 * t and u_x = 1/2 * sin(t)
            // jacobians[i][0][0] = 0.;
            // jacobians[i][0][1] = 0.;
            // jacobians[i][0][2] = 0.;

            // jacobians[i][1][0] = 0.;
            // jacobians[i][1][1] = 0.;
            // jacobians[i][1][2] = 0.;

            // jacobians[i][2][0] = 0.;
            // jacobians[i][2][1] = 0.;
            // jacobians[i][2][2] = 0.;

            // time- and space dependent velocity field
            // u_x = -1/10 * sin(t) * cos(pi/2 * x)
            // jacobians[i][0][0] = 0.1 * pi / 2 * std::sin(this->get_time()) *
            //                      std::sin(pi / 2 * points[i][0]);
            // jacobians[i][0][1] = 0.;
            // jacobians[i][0][2] = 0.;

            // jacobians[i][1][0] = 0.;
            // jacobians[i][1][1] = 0.;
            // jacobians[i][1][2] = 0.;

            // jacobians[i][2][0] = 0.;
            // jacobians[i][2][1] = 0.;
            // jacobians[i][2][2] = 0.;

            // u_x = 1/10 * x * (1 - e^(-t))
            // jacobians[i][0][0] = 0.1 * (1 - std::exp(-this->get_time()));x
            // jacobians[i][0][1] = 0.;
            // jacobians[i][0][2] = 0.;

            // jacobians[i][1][0] = 0.;
            // jacobians[i][1][1] = 0.;
            // jacobians[i][1][2] = 0.;

            // jacobians[i][2][0] = 0.;
            // jacobians[i][2][1] = 0.;
            // jacobians[i][2][2] = 0.;
          }
      }

    private:
      // Numerical constants
      double pi = 2 * std::acos(0.);
    };
  } // namespace VFP
} // namespace Sapphire
#endif
