// -----------------------------------------------------------------------------
//
// Copyright (C) 2023 by the Sapphire++ authors
//
// This file is part of Sapphire++.
//
// Sapphire++ is free software: you can redistribute it and/or modify it under
// the terms of the GNU Lesser General Public License as published by the Free
// Software Foundation, either version 3 of the License, or (at your option) any
// later version.
//
// Sapphire++ is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
// A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more
// details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with Sapphire++. If not, see <https://www.gnu.org/licenses/>.
//
// -----------------------------------------------------------------------------

/**
 * @file particle-functions.h
 * @author Nils Schween (nils.schween@mpi-hd.mpg.de)
 * @author Florian Schulze (florian.schulze@mpi-hd.mpg.de)
 * @brief Define @ref sapphirepp::VFP::ParticleVelocity and
 *        @ref sapphirepp::VFP::ParticleGamma
 */

#ifndef VFP_PARTICLEFUNCTIONS_H
#define VFP_PARTICLEFUNCTIONS_H

#include <deal.II/base/exceptions.h>
#include <deal.II/base/function.h>
#include <deal.II/base/point.h>

#include <vector>

// Functions to compute the velocity and the gamma of a particle at the point
// x,(y,z) of p in phase space
namespace sapphirepp
{
  namespace VFP
  {
    template <unsigned int dim, bool logarithmic_p>
    class ParticleVelocity : public dealii::Function<dim>
    {
    public:
      ParticleVelocity(const double &mass);
      void
      value_list(const std::vector<dealii::Point<dim>> &points,
                 std::vector<double>                   &velocities,
                 unsigned int component = 0) const override;

    private:
      const double mass;
    };

    template <unsigned int dim, bool logarithmic_p>
    class ParticleGamma : public dealii::Function<dim>
    {
    public:
      ParticleGamma(const double &mass);
      void
      value_list(const std::vector<dealii::Point<dim>> &points,
                 std::vector<double>                   &gammas,
                 unsigned int component = 0) const override;

    private:
      const double mass;
    };
  } // namespace VFP
} // namespace sapphirepp
#endif
