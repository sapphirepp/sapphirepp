#ifndef VFPEQUATION_PHYSICALSETUP_H
#define VFPEQUATION_PHYSICALSETUP_H

#include <deal.II/base/exceptions.h>
#include <deal.II/base/function.h>
#include <deal.II/base/point.h>
#include <deal.II/lac/vector.h>

#include <cmath>
#include <vector>
// NOTE: All physical quantities are dimensionless. The reference values are
// defined in the reference-values.h header.
#include "reference-values.h"

namespace Sapphire {

struct ParticleProperties {
  const double mass = 1.;
  const double charge = 1.;
};

struct TransportOnly {
  // In the transport-only case (i.e. no dependence on p) the energy of the
  // particles has to be given.
  const double energy = 1.;
  // Compute the partilces Lorentz factor gamma and its velocity
  ParticleProperties particle_properties;
  const double gamma = energy / particle_properties.mass;
  ReferenceValues reference_values;
  const double velocity =
      std::sqrt(1 - 1 / std::pow((reference_values.gamma * gamma), 2));
};
// Initial values
template <int dim>
class InitialValueFunction : public dealii::Function<dim> {
 public:
  InitialValueFunction(unsigned int exp_order)
      :  // set the number of components with the constructor of the base class
        dealii::Function<dim>((exp_order + 1) * (exp_order + 1)) {}
  void vector_value(const dealii::Point<dim> &p,
                    dealii::Vector<double> &values) const override;
};

// Scattering frequency
template <int dim>
class ScatteringFrequency : public dealii::Function<dim> {
 public:
  // Set the variable  n_components of the base class
  ScatteringFrequency() : dealii::Function<dim>(1) {}
  void value_list(const std::vector<dealii::Point<dim>> &p,
                  std::vector<double> &scattering_frequencies,
                  const unsigned int component = 0) const override;
};

// Source term
template <int dim>
class Source : public dealii::Function<dim> {
 public:
  Source(unsigned int exp_order)
      :  // set n_components
        dealii::Function<dim>((exp_order + 1) * (exp_order + 1)) {}
  void vector_value(const dealii::Point<dim> &p,
                    dealii::Vector<double> &values) const override;
};

// Magnetic field
template <int dim>
class MagneticField : public dealii::Function<dim> {
 public:
  // set n_components
  MagneticField() : dealii::Function<dim>(3) {}
  void vector_value(const dealii::Point<dim> &point,
                    dealii::Vector<double> &magnetic_field) const override;
};

// Background velocity field
template <int dim>
class BackgroundVelocityField : public dealii::Function<dim> {
 public:
  // set n_components
  BackgroundVelocityField() : dealii::Function<dim>(3) {}
  // Velocity field
  void vector_value(const dealii::Point<dim> &point,
                    dealii::Vector<double> &velocity) const override;
  // Divergence
  void divergence_list(const std::vector<dealii::Point<dim>> &points,
                       std::vector<double> &divergence);

  // Material derivative
  void material_derivative_list(
      const std::vector<dealii::Point<dim>> &points,
      std::vector<dealii::Vector<double>> &material_derivatives);

  // Jacobian matrix
  void jacobian_list(
      const std::vector<dealii::Point<dim>> &points,
      std::vector<std::vector<dealii::Vector<double>>> &jacobians) const;

 private:
  // Numerical constants
  double pi = 2 * std::acos(0.);
};
}  // namespace Sapphire
#endif
