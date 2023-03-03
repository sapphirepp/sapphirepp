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

namespace VFPEquation {
const double scattering_frequency = 1.;

struct ParticleProperties {
  const double mass = 1.;
  const double charge = 1.;

  // In the transport-only case (i.e. no dependence on p) the energy of the
  // particles has to be given.
  const double energy = 1.;
  const double gamma = energy / mass;
  // Compute the particle velocity
  ReferenceValues reference_values;
  const double velocity =
      std::sqrt(1 - 1 / std::pow((reference_values.gamma * gamma), 2));
};

// Source term
template <int dim>
class Source : public dealii::Function<dim> {
 public:
  Source(unsigned int exp_order) : expansion_order{exp_order} {}
  void vector_value(const dealii::Point<dim> &p,
                    dealii::Vector<double> &values) const override;

 private:
  const unsigned int expansion_order;
};

// Magnetic field
template <int dim>
class MagneticField : public dealii::Function<dim> {
 public:
  void vector_value(const dealii::Point<dim> &point,
                    dealii::Vector<double> &magnetic_field) const override;
};

// Background velocity field
template <int dim>
class BackgroundVelocityField : public dealii::Function<dim> {
 public:
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
}  // namespace VFPEquation
#endif
