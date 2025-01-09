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
 * @file particle-functions.cpp
 * @author Nils Schween (nils.schween@mpi-hd.mpg.de)
 * @author Florian Schulze (florian.schulze@mpi-hd.mpg.de)
 * @brief Implement @ref sapphirepp::VFP::ParticleVelocity and
 *        @ref sapphirepp::VFP::ParticleGamma
 */

#include "particle-functions.h"

#include <cmath>



template <unsigned int dim, bool logarithmic_p>
sapphirepp::VFP::ParticleVelocity<dim, logarithmic_p>::ParticleVelocity(
  const double &mass)
  : dealii::Function<dim>()
  , mass(mass)
{}



template <unsigned int dim, bool logarithmic_p>
void
sapphirepp::VFP::ParticleVelocity<dim, logarithmic_p>::value_list(
  const std::vector<dealii::Point<dim>> &points,
  std::vector<double>                   &velocities,
  unsigned int                           component) const
{
  Assert(velocities.size() == points.size(),
         dealii::ExcDimensionMismatch(velocities.size(), points.size()));
  static_cast<void>(component);

  for (unsigned int i = 0; i < points.size(); ++i)
    {
      // NOTE: It is assumed that magnitude of p is always the last component of
      // a point (i.e the coordinates of the phase space are x,(y,z), p)
      double p =
        (logarithmic_p ? std::exp(points[i][dim - 1]) : points[i][dim - 1]);
      velocities[i] = 1. / std::sqrt(mass * mass / (p * p) + 1);
    }
}

// explicit instantiation
template class sapphirepp::VFP::ParticleVelocity<1, true>;
template class sapphirepp::VFP::ParticleVelocity<1, false>;
template class sapphirepp::VFP::ParticleVelocity<2, true>;
template class sapphirepp::VFP::ParticleVelocity<2, false>;
template class sapphirepp::VFP::ParticleVelocity<3, true>;
template class sapphirepp::VFP::ParticleVelocity<3, false>;



template <unsigned int dim, bool logarithmic_p>
sapphirepp::VFP::ParticleGamma<dim, logarithmic_p>::ParticleGamma(
  const double &mass)
  : dealii::Function<dim>()
  , mass(mass)
{}



template <unsigned int dim, bool logarithmic_p>
void
sapphirepp::VFP::ParticleGamma<dim, logarithmic_p>::value_list(
  const std::vector<dealii::Point<dim>> &points,
  std::vector<double>                   &gammas,
  unsigned int                           component) const
{
  Assert(gammas.size() == points.size(),
         dealii::ExcDimensionMismatch(gammas.size(), points.size()));
  static_cast<void>(component);
  for (unsigned int i = 0; i < points.size(); ++i)
    {
      double p =
        (logarithmic_p ? std::exp(points[i][dim - 1]) : points[i][dim - 1]);
      gammas[i] = std::sqrt(p * p / (mass * mass) + 1.);
    }
}

// explicit instantiation
template class sapphirepp::VFP::ParticleGamma<1, true>;
template class sapphirepp::VFP::ParticleGamma<1, false>;
template class sapphirepp::VFP::ParticleGamma<2, true>;
template class sapphirepp::VFP::ParticleGamma<2, false>;
template class sapphirepp::VFP::ParticleGamma<3, true>;
template class sapphirepp::VFP::ParticleGamma<3, false>;
