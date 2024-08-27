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
     */
    class MHDEquations
    {
    public:
      /** @{ */
      /**
       * Dimension of the space in which the equations operate, i.e. the
       * dimension of the velocity and magnetic field.
       */
      static constexpr unsigned int spacedim = 3;
      /** Number of components `c`. */
      static constexpr unsigned int n_components = 2 + spacedim + spacedim;
      /** Index of the density component \f$ \rho \f$. */
      static constexpr unsigned int density_component = 0;
      /**
       * Starting index of the momentum components \f$ \mathbf{p} \f$
       * (@ref spacedim components).
       */
      static constexpr unsigned int first_momentum_component = 1;
      /** Index of the energy component \f$ \mathcal{E} \f$. */
      static constexpr unsigned int energy_component = spacedim + 1;
      /**
       * Starting index of the magnetic components \f$ \mathbf{B} \f$
       * (@ref spacedim components).
       */
      static constexpr unsigned int first_magnetic_component = spacedim + 2;
      /**
       * Only in primitive states: Starting index of the velocity components
       * \f$ \mathbf{u} \f$ (@ref spacedim components).
       */
      static constexpr unsigned int first_velocity_component =
        first_momentum_component;
      /** Only in primitive states: Index of the energy component \f$ P \f$. */
      static constexpr unsigned int pressure_component = energy_component;
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
        std::array<dealii::Tensor<1, spacedim, double>, n_components>;


      /** @{ */
      /** Adiabatic index \f$ \gamma \f$ */
      const double adiabatic_index;
      /** @} */


      /** Constructor */
      MHDEquations(const double adiabatic_index);


      /** @{ */
      /**
       * @brief Create a list of component names corresponding to the conserved
       *        MHD state.
       *
       * Returns a list: `rho`, `p_i`, `E`, `B_i`.
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
       * @param state The MHD state in conservative form \f$ \mathbf{w} \f$.
       * @param flux_matrix The MHD flux \f$ \mathbf{F}(\mathbf{w}) \f$.
       */
      void
      compute_flux_matrix(const state_type &state,
                          flux_type        &flux_matrix) const;



      /**
       * @brief Computes the maximum eigenvalue in normal direction.
       *
       * @param state The MHD state in conservative form \f$ \mathbf{w} \f$.
       * @param normal The normal vector \f$ \hat{\mathbf{n}} \f$.
       * @return double The maximum eigenvalue in normal direction.
       */
      double
      compute_maximum_normal_eigenvalue(
        const state_type                  &state,
        const dealii::Tensor<1, spacedim> &normal) const;
      /** @} */


      /** @{ */
      /**
       * @brief Compute the pressure \f$ P \f$ from the MHD state
       *        \f$ \mathbf{w} \f$.
       *
       * @param state The MHD state in conservative form \f$ \mathbf{w} \f$.
       * @return double The pressure \f$ P \f$.
       */
      double
      compute_pressure(const state_type &state) const;


      /**
       * @brief Compute all eigenvalues corresponding to the MHD state
       *        \f$ \mathbf{w} \f$.
       *
       * Returns a list of all eigenvalues in normal direction, ordered in the
       * following way,
       * \f[
       *   \lambda_n =
       *   \begin{pmatrix}
       *     u - c_f \\
       *     u - c_a \\
       *     u - c_s \\
       *     u       \\
       *     u       \\
       *     u + c_s \\
       *     u + c_a \\
       *     u + c_f
       *   \end{pmatrix} \,,
       * \f]
       * where \f$ c_f \f$ and \f$ c_s \f$ are the fast and slow magnetosonic
       * speeds respectively, and \f$ c_a \f$ is the Alfven speed.
       *
       * @param state The MHD state in conservative form \f$ \mathbf{w} \f$.
       * @param normal The normal vector \f$ \hat{\mathbf{n}} \f$.
       * @param eigenvalues List of the eigenvalues.
       *
       */
      void
      compute_normale_eigenvalues(const state_type                  &state,
                                  const dealii::Tensor<1, spacedim> &normal,
                                  dealii::Vector<double> &eigenvalues) const;
      /** @} */


      /** @{ */
      /**
       * @brief Convert primitive to conserved state.
       *
       * @param primitive_state Primitive state
       *        \f[
       *          \tilde{\mathbf{w}} =
       *          \begin{pmatrix}
       *            \rho        \\
       *            \mathbf{u}  \\
       *            P \\
       *            \mathbf{B}
       *          \end{pmatrix} \,.
       *        \f]
       * @param conserved_state Returns the conserved state
       *        \f[
       *          \mathbf{w} =
       *          \begin{pmatrix}
       *            \rho        \\
       *            \mathbf{p}  \\
       *            \mathcal{E} \\
       *            \mathbf{B}
       *          \end{pmatrix} \,.
       *        \f]
       */
      void
      convert_primitive_to_conserved(const state_type &primitive_state,
                                     state_type       &conserved_state) const;


      /**
       * @brief Convert conserved to primitive state.
       *
       * @param conserved_state Conserved state
       *        \f[
       *          \mathbf{w} =
       *          \begin{pmatrix}
       *            \rho        \\
       *            \mathbf{p}  \\
       *            \mathcal{E} \\
       *            \mathbf{B}
       *          \end{pmatrix} \,.
       *        \f]
       * @param primitive_state Returns the primitive state
       *        \f[
       *          \tilde{\mathbf{w}} =
       *          \begin{pmatrix}
       *            \rho        \\
       *            \mathbf{u}  \\
       *            P \\
       *            \mathbf{B}
       *          \end{pmatrix} \,.
       *        \f]
       */
      void
      convert_conserved_to_primitive(const state_type &conserved_state,
                                     state_type       &primitive_state) const;
      /** @} */
    };
  } // namespace MHD
} // namespace sapphirepp
#endif
