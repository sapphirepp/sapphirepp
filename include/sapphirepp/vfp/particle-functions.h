#ifndef VFP_PARTICLEFUNCTIONS_H
#define VFP_PARTICLEFUNCTIONS_H

#include <deal.II/base/exceptions.h>
#include <deal.II/base/function.h>
#include <deal.II/base/point.h>

#include <vector>

// Functions to compute the velocity and the gamma of a particle at the point
// x,(y,z) of p in phase space
namespace Sapphire
{
  namespace VFP
  {
    template <int dim, bool logarithmic_p>
    class ParticleVelocity : public dealii::Function<dim>
    {
    public:
      ParticleVelocity()
      {}
      void
      value_list(const std::vector<dealii::Point<dim>> &points,
                 std::vector<double>                   &velocities,
                 unsigned int component = 0) const override;
    };

    template <int dim, bool logarithmic_p>
    class ParticleGamma : public dealii::Function<dim>
    {
    public:
      ParticleGamma()
      {}
      void
      value_list(const std::vector<dealii::Point<dim>> &points,
                 std::vector<double>                   &gammas,
                 unsigned int component = 0) const override;
    };
  } // namespace VFP
} // namespace Sapphire
#endif
