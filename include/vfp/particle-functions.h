#ifndef VFPEQUATION_PARTICLEFUNCTIONS_H
#define VFPEQUATION_PARTICLEFUNCTIONS_H

#include <deal.II/base/exceptions.h>
#include <deal.II/base/function.h>
#include <deal.II/base/point.h>

#include <vector>

#include "reference-values.h"

// Functions to compute the velocity and the gamma of a particle at the point
// x,(y,z) of p in phase space
namespace Sapphire {
template <int dim>
class ParticleVelocity : public dealii::Function<dim> {
 public:
  ParticleVelocity(bool log_p) : logarithmic_p {log_p} {}
  void value_list(const std::vector<dealii::Point<dim>> &points,
                  std::vector<double> &velocities,
                  unsigned int component = 0) const override;

 private:
  bool logarithmic_p;
  ReferenceValues reference_values;
};

template <int dim>
class ParticleGamma : public dealii::Function<dim> {
 public:
  ParticleGamma(bool log_p) : logarithmic_p {log_p} {}
  void value_list(const std::vector<dealii::Point<dim>> &points,
                  std::vector<double> &gammas,
                  unsigned int component = 0) const override;

 private:
  ReferenceValues reference_values;
  bool logarithmic_p;
};
}  // namespace Sapphire
#endif
