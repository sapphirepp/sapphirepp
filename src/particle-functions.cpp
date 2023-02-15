#include "particle-functions.h"

#include <cmath>

template <int dim>
void VFPEquation::ParticleVelocity<dim>::value_list(
    const std::vector<dealii::Point<dim>> &points,
    std::vector<double> &velocities, unsigned int component) const {
  Assert(velocities.size() == points.size(),
         dealii::ExcDimensionMismatch(velocities.size(), points.size()));
  static_cast<void>(component);
  for (unsigned int i = 0; i < points.size(); ++i) {
    // NOTE: It is assumed that magnitude of p is always the last component of a
    // point (i.e the coordinates of the phase space are x,(y,z), p)
    velocities[i] =
        points[i][dim - 1] /
        std::sqrt(points[i][dim - 1] * points[i][dim - 1] +
                  1 / (reference_values.gamma * reference_values.gamma));
  }
}

// explicit instantiation
template class VFPEquation::ParticleVelocity<1>;
template class VFPEquation::ParticleVelocity<2>;
template class VFPEquation::ParticleVelocity<3>;

template <int dim>
void VFPEquation::ParticleGamma<dim>::value_list(
    const std::vector<dealii::Point<dim>> &points, std::vector<double> &gammas,
    unsigned int component) const {
  Assert(gammas.size() == points.size(),
         dealii::ExcDimensionMismatch(gammas.size(), points.size()));
  static_cast<void>(component);

  for (unsigned int i = 0; i < points.size(); ++i) {
    gammas[i] =
        std::sqrt(points[i][dim - 1] * points[i][dim - 1] +
                  1 / (reference_values.gamma * reference_values.gamma));
  }
}

// explicit instantiation
template class VFPEquation::ParticleGamma<1>;
template class VFPEquation::ParticleGamma<2>;
template class VFPEquation::ParticleGamma<3>;
