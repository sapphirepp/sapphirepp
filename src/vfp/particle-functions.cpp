#include "particle-functions.h"

#include <cmath>

#include "config.h"

template <int dim, bool logarithmic_p>
void
Sapphire::VFP::ParticleVelocity<dim, logarithmic_p>::value_list(
  const std::vector<dealii::Point<dim>> &points,
  std::vector<double>                   &velocities,
  unsigned int                           component) const
{
  Assert(velocities.size() == points.size(),
         dealii::ExcDimensionMismatch(velocities.size(), points.size()));
  static_cast<void>(component);
  ParticleProperties particle_properties;

  for (unsigned int i = 0; i < points.size(); ++i)
    {
      // NOTE: It is assumed that magnitude of p is always the last component of
      // a point (i.e the coordinates of the phase space are x,(y,z), p)
      double p =
        (logarithmic_p ? std::exp(points[i][dim - 1]) : points[i][dim - 1]);
      velocities[i] = 1. / std::sqrt(particle_properties.mass *
                                       particle_properties.mass / (p * p) +
                                     1);
      // add mass tomorrow
    }
}

// explicit instantiation
template class Sapphire::VFP::ParticleVelocity<1, true>;
template class Sapphire::VFP::ParticleVelocity<1, false>;
template class Sapphire::VFP::ParticleVelocity<2, true>;
template class Sapphire::VFP::ParticleVelocity<2, false>;
template class Sapphire::VFP::ParticleVelocity<3, true>;
template class Sapphire::VFP::ParticleVelocity<3, false>;


template <int dim, bool logarithmic_p>
void
Sapphire::VFP::ParticleGamma<dim, logarithmic_p>::value_list(
  const std::vector<dealii::Point<dim>> &points,
  std::vector<double>                   &gammas,
  unsigned int                           component) const
{
  Assert(gammas.size() == points.size(),
         dealii::ExcDimensionMismatch(gammas.size(), points.size()));
  static_cast<void>(component);
  ParticleProperties particle_properties;
  for (unsigned int i = 0; i < points.size(); ++i)
    {
      double p =
        (logarithmic_p ? std::exp(points[i][dim - 1]) : points[i][dim - 1]);
      gammas[i] = std::sqrt(
        p * p / (particle_properties.mass * particle_properties.mass) + 1.);
      // add mass tomorrow
    }
}

// explicit instantiation
template class Sapphire::VFP::ParticleGamma<1, true>;
template class Sapphire::VFP::ParticleGamma<1, false>;
template class Sapphire::VFP::ParticleGamma<2, true>;
template class Sapphire::VFP::ParticleGamma<2, false>;
template class Sapphire::VFP::ParticleGamma<3, true>;
template class Sapphire::VFP::ParticleGamma<3, false>;
