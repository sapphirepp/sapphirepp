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
 * @file numerical-flux.h
 * @author Florian Schulze (florian.schulze@mpi-hd.mpg.de)
 * @brief Define @ref sapphirepp::MHD::NumericalFlux
 */

#ifndef MHD_NUMERICALFLUX_H
#define MHD_NUMERICALFLUX_H

#include <deal.II/base/tensor.h>

#include <deal.II/lac/vector.h>

#include <vector>

#include "mhd-equations.h"
#include "mhd-flags.h"

namespace sapphirepp
{
  namespace MHD
  {
    using namespace dealii;



    /**
     * @brief Class for computing numerical flux in MHD simulations.
     *
     * The flux is computed based on the given normal vector and the states of
     * the system.
     *
     * @tparam dim Dimension of the configuration space \f$ (\mathbf{x}) \f$,
     *         `dim`
     */
    template <unsigned int dim>
    class NumericalFlux
    {
    public:
      /**
       * @brief Constructor
       *
       * @param mhd_equations Instance of the underlying @ref MHDEquations.
       */
      NumericalFlux(const MHDEquations<dim> &mhd_equations);



      /**
       * @brief Compute the numerical flux given the normal vector and states.
       *
       * Uses a local Lax-Friedrichs flux to compute the numerical flux,
       * \f[
       *   \mathbf{F}^{\text{num}} =
       *   \frac{1}{2} \hat{\mathbf{n}} \cdot
       *   \left( \mathbf{F}(\mathbf{w}_-) + \mathbf{F}(\mathbf{w}_+) \right)
       *   - \frac{\max(\lambda_{-, \hat{\mathbf{n}}},
       *   \lambda_{+, \hat{\mathbf{n}}})}{2}
       *   \left( \mathbf{w}_+ - \mathbf{w}_- \right) \,.
       * \f]
       *
       * @param normal Normal vector of the face \f$ \hat{\mathbf{n}} \f$.
       * @param state_1 Left state \f$ \mathbf{w}_- \f$.
       * @param state_2 Right state \f$ \mathbf{w}_+ \f$.
       * @param numerical_normal_flux Resulting numerical flux
       *        \f$ \mathbf{F}^{\text{num}} \f$.
       */
      void
      compute_numerical_normal_flux(
        const dealii::Tensor<1, dim>                 &normal,
        const typename MHDEquations<dim>::state_type &state_1,
        const typename MHDEquations<dim>::state_type &state_2,
        typename MHDEquations<dim>::state_type &numerical_normal_flux) const;



    private:
      /** @ref MHDEquations */
      const MHDEquations<dim> mhd_equations;
    };
  } // namespace MHD
} // namespace sapphirepp
#endif
