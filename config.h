#ifndef CONFIG_H
#define CONFIG_H

#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/base/function.h>
#include <deal.II/base/point.h>
#include <deal.II/base/tensor_function.h>

#include <deal.II/lac/vector.h>

#include <cmath>
#include <ostream>
#include <vector>

namespace Sapphire
{
  namespace VFP
  {
    enum class TermFlags
    {
      none              = 0,
      spatial_advection = 1 << 0,
      collision         = 1 << 1,
      magnetic          = 1 << 2,
      momentum          = 1 << 3,
      source            = 1 << 4
    };

    constexpr TermFlags
    operator|(TermFlags f1, TermFlags f2)
    {
      return static_cast<TermFlags>(static_cast<int>(f1) |
                                    static_cast<int>(f2));
    }

    constexpr TermFlags
    operator&(TermFlags f1, TermFlags f2)
    {
      return static_cast<TermFlags>(static_cast<int>(f1) &
                                    static_cast<int>(f2));
    }

    template <typename StreamType>
    inline StreamType &
    operator<<(StreamType &os, TermFlags f)
    {
      os << "Term flags: \n";
      if ((f & TermFlags::spatial_advection) != TermFlags::none)
        os << "	 - Spatial Advection\n";
      if ((f & TermFlags::collision) != TermFlags::none)
        os << "	 - Collision\n";
      if ((f & TermFlags::magnetic) != TermFlags::none)
        os << "	 - Magnetic\n";
      if ((f & TermFlags::momentum) != TermFlags::none)
        os << "	 - Momentum\n";
      if ((f & TermFlags::source) != TermFlags::none)
        os << "	 - Source\n";
      return os;
    }

    // explicit instantiation
    template std::ostream &
    operator<<(std::ostream &os, TermFlags f);
    template dealii::ConditionalOStream &
    operator<<(dealii::ConditionalOStream &os, TermFlags f);


    // Physical setup

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
      const double velocity = std::sqrt(1 - 1 / std::pow((gamma * gamma), 2));
    };

    // Initial values
    template <int dim>
    class InitialValueFunction : public dealii::Function<dim>
    {
    public:
      InitialValueFunction(unsigned int exp_order)
        : // set the number of components with the constructor of the base class
        dealii::Function<dim>((exp_order + 1) * (exp_order + 1))
      {}

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
      ScatteringFrequency()
        : dealii::Function<dim>(1)
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
      Source(unsigned int exp_order)
        : // set n_components
        dealii::Function<dim>((exp_order + 1) * (exp_order + 1))
      {}

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
      MagneticField()
        : dealii::Function<dim>(3)
      {}

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
      BackgroundVelocityField()
        : dealii::Function<dim>(3)
      {}

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

  namespace Hydro
  {
    using namespace dealii;

    /**
     * @brief Velocity field for a rigid rotator.
     *
     * \f$ \mathbf{\beta}(\v{x}, t) = \omega (y \mathbf{e}_x - x \mathbf{e}_y)
     * \f$
     *
     * \tparam dim: space dimension
     */
    template <int dim>
    class VelocityFieldRigidRotator : public TensorFunction<1, dim, double>
    {
    public:
      VelocityFieldRigidRotator(const double &omega, const double time = 0.0)
        : TensorFunction<1, dim, double>(time)
        , omega(omega)
      {}

      Tensor<1, dim>
      value(const Point<dim> &p) const override
      {
        AssertDimension(dim, 2);

        Tensor<1, dim> values;
        values[0] = omega * p[1];
        values[1] = -omega * p[0];
        return values;
      }

      Tensor<2, dim>
      gradient(const Point<dim> &p) const override
      {
        (void)p; // suppress unused parameter warning
        AssertDimension(dim, 2);

        Tensor<2, dim> values;
        values[0][0] = 0.0;
        values[0][1] = omega;
        values[1][0] = -omega;
        values[1][1] = 0.0;
        return values;
      }

    private:
      const double omega;
    };

    /**
     * @brief Exact analytical solution of rigid rotator.
     *
     * \tparam dim: space dimension
     */
    template <int dim>
    class ExactSolutionRigidRotator : public Function<dim>
    {
    public:
      ExactSolutionRigidRotator(const double &omega, const double time = 0.0)
        : Function<dim>(1, time)
        , omega(omega)
      {}

      double
      value(const Point<dim>  &p,
            const unsigned int component = 0) const override
      {
        AssertIndexRange(component, 1);
        AssertDimension(dim, 2);

        Point<dim> x;
        x[0] = p[0] * std::cos(omega * this->get_time()) -
               p[1] * std::sin(omega * this->get_time());
        x[1] = p[0] * std::sin(omega * this->get_time()) +
               p[1] * std::cos(omega * this->get_time());

        const double width  = 0.2;
        const double length = 0.7;

        if ((std::abs(x[0]) < width) && (std::abs(x[1]) < length))
          return 1.0;
        if ((std::abs(x[1]) < width) && (std::abs(x[0]) < length))
          return 1.0;

        return 0.0;
      }

    private:
      const double omega;
    };

    /**
     * @brief Constant velocity field.
     *
     * \tparam dim: space dimension
     */
    template <int dim>
    class ConstantVelocityField : public TensorFunction<1, dim, double>
    {
    public:
      ConstantVelocityField(const Tensor<1, dim> &beta)
        : TensorFunction<1, dim, double>()
        , beta(beta)
      {}

      Tensor<1, dim>
      value(const Point<dim> &p) const override
      {
        (void)p; // suppress unused parameter warning

        return beta;
      }

      Tensor<2, dim>
      gradient(const Point<dim> &p) const override
      {
        (void)p; // suppress unused parameter warning

        Tensor<2, dim> values;
        values[0][0] = 0.0;
        values[0][1] = 0.0;
        values[1][0] = 0.0;
        values[1][1] = 0.0;
        return values;
      }

    private:
      const Tensor<1, dim> beta;
    };

    /**
     * @brief Exact analytical solution of constant linear advection equation.
     *
     * \f$ u(x, t) = u_0(x - \beta t) \f$
     * with
     * <!-- \f$ u_0(\mathbf{x}) =  1 \f$ -->
     * <!-- \f$ u_0(\mathbf{x}) =  \sin(\pi * \hat{n} \cdot \mathbf{x}) \f$
     * --> \f$ u_0(\mathbf{x}) =  \exp(-\mathbf{x} \cdot \mathbf{x} / (2
     * \sigma^2)) \f$
     *
     * \tparam dim: space dimension
     */
    template <int dim>
    class ExactSolutionConstantAdvection : public Function<dim>
    {
    public:
      ExactSolutionConstantAdvection(const Tensor<1, dim> &beta,
                                     const double          time = 0.0)
        : Function<dim>(1, time)
        , beta(beta)
      {}

      double
      value(const Point<dim>  &p,
            const unsigned int component = 0) const override
      {
        AssertIndexRange(component, 1);

        const Point<dim> x = p - beta * this->get_time();

        // return 1.0;
        // Tensor<1, dim> normal;
        // normal[0] = 1.0;
        // normal[1] = 1.0;
        // normal /= normal.norm();
        // return std::sin(numbers::PI * normal * x);
        const double sigma = 0.1;
        return std::exp(-(x * x) / (2.0 * sigma * sigma));
      }

    private:
      const Tensor<1, dim> beta;
    };

    /**
     * @brief Exact solution of Burgers' equation.
     *
     * \f$ u(x, t) = \Theta(x - t/2) \f$
     * with \f$ \Theta(x) \f$ being the Heaviside step function.
     *
     * \tparam dim: space dimension (must be 1)
     */
    template <int dim>
    class ExactSolutionBurgersEq : public Function<dim>
    {
    public:
      ExactSolutionBurgersEq(const double time = 0.0)
        : Function<dim>(1, time)
      {
        AssertDimension(dim, 1);
      }

      double
      value(const Point<dim>  &p,
            const unsigned int component = 0) const override
      {
        (void)component; // suppress unused parameter warning
        Point<dim> x;

        // return -std::sin(numbers::PI * p[0]);

        // x[0] = p[0];
        // return 0.25 + 0.5 * std::sin(numbers::PI * (2 * x[0] - 1));

        // x[0] = p[0] - this->get_time() / 2.0;
        // if (x[0] > 0.0)
        //   return 1.0;
        // else
        //   return 0.0;

        x[0] = p[0] - this->get_time() * 3.0;
        if (x[0] > -0.5)
          return 1.0;
        else
          return 2.0;
      }
    };

    /**
     * @brief Initial condition extracted from an exact solution.
     *
     * \f$ u_0(\mathbf{x}) = u(\mathbf{x}, t = 0) \f$
     *
     * \tparam dim: space dimension
     */
    template <int dim>
    class InitialConditionBurgersEq : public Function<dim>
    {
    public:
      InitialConditionBurgersEq()
        : Function<dim>(1)
        , exact_solution()
      {
        exact_solution.set_time(0.0);
      }

      double
      value(const Point<dim>  &p,
            const unsigned int component = 0) const override
      {
        return exact_solution.value(p, component);
      }

    private:
      ExactSolutionBurgersEq<dim> exact_solution;
    };

    /**
     * @brief Boundary values extracted from an exact solution.
     *
     * \f$ u_b(\mathbf{x}, t) = u(\mathbf{x}, t) \f$
     *
     * @tparam dim: space dimension
     */
    template <int dim>
    class BoundaryValuesBurgersEq : public Function<dim>
    {
    public:
      BoundaryValuesBurgersEq(const double time = 0.0)
        : Function<dim>(1, time)
        , exact_solution()
      {}

      double
      value(const Point<dim>  &p,
            const unsigned int component = 0) const override
      {
        return exact_solution.value(p, component);
      }

      void
      set_time(const double new_time) override
      {
        exact_solution.set_time(new_time);
      }



    private:
      ExactSolutionBurgersEq<dim> exact_solution;
    };

    // /**
    //  * @brief Initial condition extracted from an exact solution.
    //  *
    //  * \f$ u_0(\mathbf{x}) = u(\mathbf{x}, t = 0) \f$
    //  *
    //  * \tparam dim: space dimension
    //  */
    // template <int dim>
    // class InitialCondition : public Function<dim>
    // {
    // public:
    //   InitialCondition(Function<dim> *exact_solution)
    //     : Function<dim>(1)
    //     , exact_solution(exact_solution)
    //   {}

    //   double
    //   value(const Point<dim>  &p,
    //         const unsigned int component = 0) const override
    //   {
    //     exact_solution->set_time(0.0);
    //     return exact_solution->value(p, component);
    //   }

    // private:
    //   const SmartPointer<Function<dim>> exact_solution;
    // };

    // /**
    //  * @brief Boundary values extracted from an exact solution.
    //  *
    //  * \f$ u_b(\mathbf{x}, t) = u(\mathbf{x}, t) \f$
    //  *
    //  * @tparam dim: space dimension
    //  */
    // template <int dim>
    // class BoundaryValues : public Function<dim>
    // {
    // public:
    //   BoundaryValues(Function<dim> *exact_solution, const double time =
    //   0.0)
    //     : Function<dim>(1, time)
    //     , exact_solution(exact_solution)
    //   {}

    //   double
    //   value(const Point<dim>  &p,
    //         const unsigned int component = 0) const override
    //   {
    //     exact_solution->set_time(this->get_time());
    //     return exact_solution->value(p, component);
    //   }

    // private:
    //   const SmartPointer<Function<dim>> exact_solution;
    // };
  } // namespace Hydro
} // namespace Sapphire
#endif
