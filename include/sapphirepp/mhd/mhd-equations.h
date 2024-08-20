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
 * @file mhd-equations.h
 * @author Florian Schulze (florian.schulze@mpi-hd.mpg.de)
 * @brief Define @ref sapphirepp::MHD::MHDEquations
 */

#ifndef MHD_MHDEQUATIONS_H
#define MHD_MHDEQUATIONS_H

#include <deal.II/base/tensor.h>

#include <deal.II/lac/vector.h>

#include <array>
#include <string>
#include <vector>

namespace sapphirepp
{
  namespace MHD
  {
    /**
     * @brief Define functions related to MHD equations.
     *
     * This class calculates everything related to the PDE system resulting from
     * the MHD equations, including the mapping between the components and the
     * corresponding physical quantities.
     *
     * @tparam dim Dimension of the configuration space \f$ (\mathbf{x}) \f$,
     *         `dim`
     */
    template <unsigned int dim>
    class MHDEquations
    {
    public:
      /** @{ */
      /** Dimension of the velocity and magnetic field (always fixed to 3). */
      static constexpr unsigned int dim_uB = 3;
      /** Number of components `c`. */
      static constexpr unsigned int n_components = 2 + dim_uB + dim_uB;
      /** Index of the density component \f$ \rho \f$. */
      static constexpr unsigned int density_component = 0;
      /**
       * Starting index of the momentum components \f$ \mathbf{p} \f$
       * (@ref dim_uB components).
       */
      static constexpr unsigned int first_momentum_component = 1;
      /** Index of the energy component \f$ \mathcal{E} \f$. */
      static constexpr unsigned int energy_component = dim + 1;
      /**
       * Starting index of the magnetic components \f$ \mathbf{B} \f$
       * (@ref dim_uB components).
       */
      static constexpr unsigned int first_magnetic_component = dim + 2;
      /** @} */

      /**
       * @brief Type definition for MHD states \f$ \mathbf{w} \f$.
       *
       * The states in conservative form are defined as
       * \f[
       *   \mathbf{w} =
       *   \begin{pmatrix}
       *     \rho        \\
       *     \mathbf{p}  \\
       *     \mathcal{E} \\
       *     \mathbf{B}
       *   \end{pmatrix} \,.
       * \f]
       * They are indexed as `state[c]` with the index `c` corresponding to the
       * component.
       */
      using state_type = dealii::Vector<double>;
      /**
       * @brief Type definition for MHD flux matrices \f$ \mathbf{F} \f$.
       *
       * They are indexed as `flux_matrix[c][d]` with the first index `c`
       * corresponding to the component and the second index `d` corresponding
       * to the spatial dimension.
       */
      using flux_type =
        std::array<dealii::Tensor<1, dim, double>, n_components>;


      /** Constructor */
      MHDEquations();


      /** @{ */
      /**
       * @brief Create a list of component names: `rho`, `p_i`, `E`, `B_i`.
       *
       * @param prefix Prefix for the component names.
       * @return std::vector<std::string> component_names.
       */
      static std::vector<std::string>
      create_component_name_list(const std::string &prefix = "");
      /** @} */


      /** @{ */
      /**
       * @brief Compute the MHD flux matrix.
       *
       * @param state The MHD state \f$ \mathbf{w} \f$.
       * @param flux_matrix The MHD flux \f$ \mathbf{F}(\mathbf{w}) \f$.
       */
      void
      compute_flux_matrix(const state_type &state,
                          flux_type        &flux_matrix) const;



      /**
       * @brief Computes the maximal eigenvalue in normal direction.
       *
       * @param state The MHD state \f$ \mathbf{w} \f$.
       * @param normal The normal vector \f$ \hat{\mathbf{n}} \f$.
       * @return double The maximal eigenvalue in normal direction.
       */
      double
      compute_maximal_eigenvalue_normal(
        const state_type             &state,
        const dealii::Tensor<1, dim> &normal) const;
      /** @} */
    };
  } // namespace MHD
} // namespace sapphirepp
#endif
