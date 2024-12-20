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



namespace sapphirepp
{
  namespace VFP
  {

    /**
     * @brief A function to compute the particle velocity at a point in
     *        reduced phase space \f$ (\mathbf{x}, p) \f$.
     *
     * @tparam dim Dimension of the reduced phase space \f$ (\mathbf{x}, p) \f$
     * @tparam logarithmic_p Do we use a logarithmic momentum variable?
     */
    template <unsigned int dim, bool logarithmic_p>
    class ParticleVelocity : public dealii::Function<dim>
    {
    public:
      /**
       * @brief Constructor
       *
       * @param mass Mass of the particles
       */
      ParticleVelocity(const double &mass);

      /**
       * @brief Compute the particle velocity at a list of points.
       *
       * @param points List of points in reduced phase space
       *               \f$ (\mathbf{x}, p) \f$
       * @param velocities List of (absolute) velocities
       * @param component unused
       */
      void
      value_list(const std::vector<dealii::Point<dim>> &points,
                 std::vector<double>                   &velocities,
                 unsigned int component = 0) const override;


    private:
      /** Mass of the particles */
      const double mass;
    };



    /**
     * @brief Function to compute the particle gamma factor \f$ \gamma \f$
     *        at a point in reduced phase space \f$ (\mathbf{x}, p) \f$.
     *
     * @tparam dim Dimension of the reduced phase space \f$ (\mathbf{x}, p) \f$
     * @tparam logarithmic_p Do we use a logarithmic momentum variable?
     */
    template <unsigned int dim, bool logarithmic_p>
    class ParticleGamma : public dealii::Function<dim>
    {
    public:
      /**
       * @brief Constructor
       *
       * @param mass Mass of the particles
       */
      ParticleGamma(const double &mass);


      /**
       * @brief Compute the particle gamma factor at a list of points.
       *
       * @param points List of points in reduced phase space
       *               \f$ (\mathbf{x}, p) \f$
       * @param gammas List of gamma factors
       * @param component unused
       */
      void
      value_list(const std::vector<dealii::Point<dim>> &points,
                 std::vector<double>                   &gammas,
                 unsigned int component = 0) const override;

    private:
      /** Mass of the particles */
      const double mass;
    };
  } // namespace VFP
} // namespace sapphirepp
#endif
