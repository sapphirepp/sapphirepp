/// \file examples/scattering/config.h
/// \author Florian Schulze (florian.schulze@mpi-hd.mpg.de)
/// \brief config.h file for the scattering example
/// \date 2023-10-27
/// @copyright Copyright (c) 2023

/// \page scattering
/// \subsection config config.h
/// \dontinclude scattering/config.h
/// We start by going line by line trough the `config.h` file.
/// First, we have to make sure that the file is only included onces, and then
/// import some dependencies:
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

#  include "sapphire-logstream.h"
#  include "vfp-flags.h"

/// \page scattering
/// \skip ifndef
/// \until CODEBLOCK_DELIMITER
/// Everything implemented in `Sapphire++` is part of the namespace Sapphire.
/// \todo How should we name Sapphire in the docs? Sapphire, Sapphire++,
/// `Sapphire++`, ...?
namespace Sapphire
{
  /// \page scattering
  /// \until CODEBLOCK_DELIMITER
  /// Often we parametrize the physical setup with some runtime parameter.
  /// Since it is setup dependent what these parameters are, they have to be
  /// specified by the user. The PhysicalProperties class allows for this. It
  /// uses the deal.II  concept of a ParameterHandler, for more details see
  /// ... \todo Link to deal.II, and how to name deal.II ?
  ///
  /// The PhysicalProperties class consists of __public__ variables for the
  /// user defined runtime parameter, a default constructor and tho functions
  /// to __delcare__ and __parse__ the parameter from the parameter file.
  /// In this example, we have two parameter that we want to specify, the
  /// scattering frequency \f$ \nu \f$ and the initial value of the expansion
  /// coefficients, \f$ f_{lms, 0} \f$. We will call these parameters `nu` and
  /// `f0` respectively, assuming all expansion coefficients have the same
  /// initial value.
  class PhysicalProperties
  {
  public:
    PhysicalProperties() = default;

    /// \page scattering
    /// \until CODEBLOCK_DELIMITER
    /// The declare_parameters function is using the dealii::ParameterHandler
    /// class to delcare parameter in the parameter file. We sort all parameter
    /// in a subsection "Physical properties".
    /// In addition, we write a message to the custom Sapphire logstream
    /// `saplog`. The LogStream::Prefix ensures, that the message is prefixed
    /// and only shown, if detailed output is requested.
    /// The declaration of the parameter is straight forward, using the
    /// `declare_entry` function of the ParameterHandler. It takes the name of
    /// the parameter, a default value and its type/pattern. Additionally, can
    /// give a description of the parameter.
    void
    declare_parameters(dealii::ParameterHandler &prm)
    {
      dealii::LogStream::Prefix p("Physical properties", saplog);
      saplog << "Declaring parameters" << std::endl;
      prm.enter_subsection("Physical properties");

      prm.declare_entry("nu",
                        "1",
                        dealii::Patterns::Double(),
                        "Scattering frequency");
      prm.declare_entry("f0",
                        "1.",
                        dealii::Patterns::Double(),
                        "Initial value of the expansion coefficients");

      prm.leave_subsection();
    };

    /// \page scattering
    /// \until CODEBLOCK_DELIMITER
    /// The parsing of the parameter is equally simple. We use the
    /// `get_double()` function of the ParameterParser to get the value for a
    /// previously declared parameter.
    void
    parse_parameters(dealii::ParameterHandler &prm)
    {
      dealii::LogStream::Prefix p("PhysicalProperties", saplog);
      saplog << "Parsing parameters" << std::endl;
      prm.enter_subsection("Physical properties");

      nu = prm.get_double("nu");
      f0 = prm.get_double("f0");

      prm.leave_subsection();
    };


    /// \page scattering
    /// \until CODEBLOCK_DELIMITER
    /// At the end, we define the runtime parameter as __public__ variables of
    /// the class, so that subsequent functions have easy access to them.
    double nu;
    double f0;
  };


  /// \page scattering
  /// \until CODEBLOCK_DELIMITER
  /// Next, we define static variables and functions related to the VFP
  /// equation. We therefore use the namespace `VFP` to collect them in one
  /// place.
  namespace VFP
  {
    /// \page scattering
    /// \until CODEBLOCK_DELIMITER
    /// First, we define the dimensionality of the problem. Since the solution
    /// does not depend on either \f$ \mathbf{x} \f$ nor \f$ p \f$, we just use
    /// one space dimension.
    /// \todo Use dimension = 0 if possible
    static constexpr int dimension = 1;

    /// \page scattering
    /// \until CODEBLOCK_DELIMITER
    /// Next, we define which terms of the VFP equation we use. To this end, we
    /// define a `static constexpr` which the __compiler__ can use to determine
    /// which terms are activated. Selecting only the relevant terms here
    /// results in big performance gains, even so setting the respective terms
    /// to 0 at runtime would produce the same results.
    /// In this example, we only have a scattering term, and no \f$ p \f$
    /// dependence. We therefore only activate the `collision` term, while all
    /// other terms (and the \f$ p \f$ dependence are turned off by default.)
    static constexpr VFPFlags vfp_flags = VFPFlags::collision;


    /// \page scattering
    /// \until CODEBLOCK_DELIMITER
    /// To define the initial values, we use the class `Function` provided by
    /// `deal.II`. We define our own class `InitialValueFunction` functions that
    /// inherits all properties of the parent class `Function`. Keeping the
    /// style of `deal.II`, we keep the dimension as a template parameter (even
    /// so it is at the moment defined by the compile time variable
    /// `dimension`).
    ///
    /// As we have seen before, the inital condition depends on only one
    /// parameter, `f0`. But in this example, we will extend the scope of this
    /// function, by using it as the analytic solution we can compare to.
    /// Therefore, the function will also depend on `nu` as a function of time.
    /// The function has to provide a value for each component \f$ f_{lms} \f$.
    /// Therefore it is a __vector-valued__ fuction with \f$ (l_{\rm max} +1)^2
    /// \f$ components.
    template <int dim>
    class InitialValueFunction : public dealii::Function<dim>
    {
    public:
      InitialValueFunction(const PhysicalProperties &physical_properties,
                           int                       exp_order)
        : dealii::Function<dim>((exp_order + 1) * (exp_order + 1))
        , f0(physical_properties.f0)
        , nu(physical_properties.nu)
      {}

      /// \page scattering
      /// \until CODEBLOCK_DELIMITER
      /// After defining the constructor of the function, we have to define its
      /// value. For this we override the `vector_value` function of the parent
      /// class. At a *point* `p` it has to provide a *value* `f` for each
      /// component.
      /// In this example, value is given by
      /// \f[ f_{lms}(t) = f_{lms, 0} \exp\left(-\nu \frac{l(l + 1)}{2} t
      /// \right)
      /// \,. \f]
      void
      vector_value(const dealii::Point<dim> &p,
                   dealii::Vector<double>   &f) const override
      {
        static_cast<void>(p); // suppress compiler warning

        for (unsigned int i = 0; i < InitialValueFunction<dim>::n_components;
             i++)
          {
            const int    l = 0; // TODO
            const double t = this->get_time();
            f[i]           = f0 * std::exp(-nu * l * (l + 1) / 2. * t);
          }
      }

      const double f0;
      const double nu;
    };

    /// \page scattering
    /// \until CODEBLOCK_DELIMITER
    /// The definition of the scattering frequency works similar. The only
    /// difference is, that the scattering frequency is a scalar function,
    /// therefore we use the function `value_list` to get the scattering
    /// frequency at multiple points in one function call.
    template <int dim>
    class ScatteringFrequency : public dealii::Function<dim>
    {
    public:
      ScatteringFrequency(const PhysicalProperties &physical_properties)
        : dealii::Function<dim>(1)
        , nu(physical_properties.nu)
      {}

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
                  nu);
      }

      const double nu;
    };

    /// \page scattering
    /// \until CODEBLOCK_DELIMITER
    /// Since `Sapphire++` a definition of the function `Source`,
    /// `MagneticField` and `BackgroundVelocityField` in the `config.h`, we have
    /// to implement them here. We can however leave the implementation empty,
    /// since the functions won't be used in this example.
    template <int dim>
    class Source : public dealii::Function<dim>
    {
    public:
      Source(const PhysicalProperties &physical_properties,
             unsigned int              exp_order)
        : dealii::Function<dim>((exp_order + 1) * (exp_order + 1))
      {
        static_cast<void>(physical_properties); // suppress compiler warning
      }

      void
      vector_value(const dealii::Point<dim> &p,
                   dealii::Vector<double>   &values) const override
      {
        static_cast<void>(p); // suppress compiler warning

        values[0] = 0.; // unused
      }
    };

    /// \page scattering
    /// \until CODEBLOCK_DELIMITER
    /// Empty implementation of the magnetic field:
    template <int dim>
    class MagneticField : public dealii::Function<dim>
    {
    public:
      MagneticField(const PhysicalProperties &physical_properties)
        : dealii::Function<dim>(3)
      {
        static_cast<void>(physical_properties); // suppress compiler warning
      }

      void
      vector_value(const dealii::Point<dim> &point,
                   dealii::Vector<double>   &magnetic_field) const override
      {
        static_cast<void>(point); // suppress compiler warning

        // constant magnetic field
        magnetic_field[0] = 0.;
        magnetic_field[1] = 0.;
        magnetic_field[2] = 0.;
      }
    };

    /// \page scattering
    /// \until CODEBLOCK_DELIMITER
    /// Empty implementation of the background velocity field:
    template <int dim>
    class BackgroundVelocityField : public dealii::Function<dim>
    {
    public:
      BackgroundVelocityField(const PhysicalProperties &physical_properties)
        : dealii::Function<dim>(3)
      {
        static_cast<void>(physical_properties); // suppress compiler warning
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

    /// \page scattering
    /// \until CODEBLOCK_DELIMITER
    /// Last, we have to close the namespaces again and end the include guard.
  } // namespace VFP
} // namespace Sapphire
#endif

/// \page scattering
/// \until END_OF_FILE
/// This concludes the `config.h` file. Next, we have to implement the `main` in
/// the `scattering.cpp` file.
