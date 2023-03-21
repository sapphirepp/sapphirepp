#include "particle-functions.h"

#include <cmath>

template <int dim>
void Sapphire::ParticleVelocity<dim>::value_list(
    const std::vector<dealii::Point<dim>> &points,
    std::vector<double> &velocities, unsigned int component) const {
  Assert(velocities.size() == points.size(),
         dealii::ExcDimensionMismatch(velocities.size(), points.size()));
  static_cast<void>(component);
  for (unsigned int i = 0; i < points.size(); ++i) {
    // NOTE: It is assumed that magnitude of p is always the last component of a
    // point (i.e the coordinates of the phase space are x,(y,z), p)
    double p =
        (logarithmic_p ? std::exp(points[i][dim - 1]) : points[i][dim - 1]);
    velocities[i] = p / std::sqrt(p * p + 1 / (reference_values.gamma *
                                               reference_values.gamma));
  }
}

// explicit instantiation
template class Sapphire::ParticleVelocity<1>;
template class Sapphire::ParticleVelocity<2>;
template class Sapphire::ParticleVelocity<3>;

template <int dim>
void Sapphire::ParticleGamma<dim>::value_list(
    const std::vector<dealii::Point<dim>> &points, std::vector<double> &gammas,
    unsigned int component) const {
  Assert(gammas.size() == points.size(),
         dealii::ExcDimensionMismatch(gammas.size(), points.size()));
  static_cast<void>(component);

  for (unsigned int i = 0; i < points.size(); ++i) {
    double p =
        (logarithmic_p ? std::exp(points[i][dim - 1]) : points[i][dim - 1]);
    gammas[i] = std::sqrt(
        p * p + 1 / (reference_values.gamma * reference_values.gamma));
  }
}

// explicit instantiation
template class Sapphire::ParticleGamma<1>;
template class Sapphire::ParticleGamma<2>;
template class Sapphire::ParticleGamma<3>;
