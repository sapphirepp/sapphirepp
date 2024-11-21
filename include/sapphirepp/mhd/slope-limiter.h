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
 * @file slope-limiter.h
 * @author Florian Schulze (florian.schulze@mpi-hd.mpg.de)
 * @brief Define utility functions related to slope limiting
 */

#ifndef MHD_SLOPELIMITER_H
#define MHD_SLOPELIMITER_H

#include <deal.II/base/tensor.h>

#include <deal.II/lac/vector.h>

#include <array>
#include <vector>

#include "mhd-equations.h"
#include "mhd-flags.h"

namespace sapphirepp
{
  namespace MHD
  {

    /**
     * @namespace sapphirepp::MHD::SlopeLimiter
     * @brief Namespace for functions related to slope limiting
     */
    namespace SlopeLimiter
    {

      /**
       * @brief Calculate the minmod function of a list of values.
       *
       * The minmod function is defined as,
       * \f[
       *   {\rm minmod}(v_1, \dots, v_k) =
       *   \begin{cases}
       *     s \min(|v_1|, \dots, |v_k|)
       *      & {\rm if} \quad s = {\rm sgn}(v_1) = \dots {\rm sgn}(v_k) \\
       *     0
       *      & {\rm otherwise}
       *   \end{cases} \,.
       * \f]
       *
       * @param values Vector of values \f$ (v_1, \dots, v_k) \f$.
       * @return double The minmod function of the input values,
       *                \f$ {\rm minmod}(v_1, \dots, v_k) \f$.
       */
      double
      minmod(const std::vector<double> &values);



      /**
       * @brief Calculate the component and direction wise minmod function for a
       *        list of gradients.
       *
       * This functions ignores `nan`/`inf` values, and only takes the finite
       * entires into account.
       *
       * @param cell_gradient Average gradient on the cell.
       * @param neighbor_gradients Gradients to neighbor cells.
       * @param limited_gradient Component and direction wise minmod function,
       *                         `minmod(cell_gradient, neighbor_gradients)`.
       * @return double Average difference between cell_gradient and limited_gradient
       *                per component and direction.
       */
      double
      minmod_gradients(
        const MHDEquations::flux_type              &cell_gradient,
        const std::vector<MHDEquations::flux_type> &neighbor_gradients,
        MHDEquations::flux_type                    &limited_gradient);
    } // namespace SlopeLimiter
  }   // namespace MHD
} // namespace sapphirepp
#endif
