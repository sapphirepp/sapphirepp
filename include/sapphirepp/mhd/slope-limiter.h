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

#include <deal.II/base/exceptions.h>
#include <deal.II/base/tensor.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/lac/vector.h>

#include <array>
#include <vector>

#include "mhd-equations.h"
#include "mhd-flags.h"

namespace sapphirepp
{
  namespace MHD
  {
    using namespace dealii;



    /**
     * @brief Collect functions related to slope limiting
     *
     * @tparam dim Dimension of the configuration space \f$ (\mathbf{x}) \f$,
     *         `dim`
     * @tparam divergence_cleaning Use Lagrange multiplier \f$ \psi \f$
     *         for hyperbolic divergence cleaning.
     */
    template <unsigned int dim, bool divergence_cleaning>
    class SlopeLimiter
    {
    public:
      /** Constructor */
      SlopeLimiter() = default;



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
      static double
      minmod(const std::vector<double> &values);



      /**
       * @brief Calculate the component and direction wise modified minmod
       *        function for a list of gradients.
       *
       * This functions ignores `nan`/`inf` values, and only takes the finite
       * entires into account. It furthermore uses only the modified minmod
       * function.
       *
       * @param cell_gradient Average gradient on the cell.
       * @param neighbor_gradients Gradients to neighbor cells.
       * @param limited_gradient Component and direction wise minmod function,
       *                         `minmod(cell_gradient, neighbor_gradients)`.
       * @param dx Cell size \f$ \Delta x \f$ for modified minmod.
       * @return double Average difference between cell_gradient and limited_gradient
       *                per component and direction.
       */
      static double
      minmod_gradients(
        const typename MHDEquations<dim, divergence_cleaning>::flux_type
          &cell_gradient,
        const std::vector<
          typename MHDEquations<dim, divergence_cleaning>::flux_type>
          &neighbor_gradients,
        typename MHDEquations<dim, divergence_cleaning>::flux_type
                    &limited_gradient,
        const double dx);



      /**
       * @brief Computes the distance form the neighbor cell to this cell in
       *        direction of face `face_no`, accounting for periodic boundaries.
       *
       * @param cell Current cell.
       * @param face_no Face number for direction of neighbor.
       * @return Tensor<1, dim> Distance between this cell and the
       *         neighbor cell, accounting for periodic boundaries.
       */
      static inline Tensor<1, dim>
      compute_periodic_distance_cell_neighbor(
        const typename DoFHandler<dim>::active_cell_iterator &cell,
        const unsigned int                                    face_no)
      {
        if (cell->has_periodic_neighbor(face_no))
          return (cell->face(face_no)->center() - cell->center()) * 2.;

        Assert(!cell->at_boundary(face_no),
               ExcMessage("Can not compute distance for non-periodic "
                          "boundary cells."));
        return cell->neighbor(face_no)->center() - cell->center();
      }



      /**
       * @brief Enforce that the limited gradients
       *        have a divergence free magnetic field
       *
       * @param limited_gradient limited gradient
       */
      static void
      enforce_divergence_free_limited_gradient(
        typename MHDEquations<dim, divergence_cleaning>::flux_type
          &limited_gradient);
    };
  } // namespace MHD
} // namespace sapphirepp
#endif
