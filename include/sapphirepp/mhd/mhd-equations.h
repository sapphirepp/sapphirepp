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

#include <deal.II/base/exceptions.h>
#include <deal.II/base/point.h>
#include <deal.II/base/tensor.h>

#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/vector.h>

#include <deal.II/numerics/data_component_interpretation.h>

#include <array>
#include <exception>
#include <ostream>
#include <string>
#include <vector>

namespace sapphirepp
{
  namespace MHD
  {
    /**
     * @brief This exception is raised if a non-admissible MHD state is
     *        encountered, e.g. due to negative energy or pressure.
     *
     * @tparam dim Dimension of the configuration space \f$ (\mathbf{x}) \f$,
     *         `dim`
     */
    template <unsigned int dim = 1>
    class ExcNonAdmissibleState : public dealii::ExceptionBase
    {
    public:
      /**
       * @brief Constructor
       *
       * @param state The @ref MHDEquations::state_type "MHD state"
       * @param msg Error message
       * @param point Point on the grid where the error occurred
       */
      ExcNonAdmissibleState(
        const dealii::Vector<double>     &state,
        const std::string                &msg   = "",
        const dealii::Point<dim, double> &point = dealii::Point<dim, double>())
        : state(state)
        , msg(msg)
        , point(point)
      {}

      /** Destructor */
      virtual ~ExcNonAdmissibleState() noexcept
      {}



      /**
       * @brief Print error message
       *
       * @param out Output stream
       */
      virtual void
      print_info(std::ostream &out) const override
      {
        if (!msg.empty())
          out << "    " << msg << "\n";
        out << "    "
            << "Non-admissible MHD state encountered";
        if (point.norm_square() > 0.)
          out << " at (" << point << ")";
        out << ": \n"
            << "    " << state << std::endl;
      }



    private:
      const dealii::Vector<double>     state;
      const std::string                msg;
      const dealii::Point<dim, double> point;
    };



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
      /**
       * Number of components in a vector,
       * e.g. for momentum or magnetic field.
       */
      static constexpr unsigned int n_vec_components = 3;
      /** Number of components `c`. */
      static constexpr unsigned int n_components =
        2 + n_vec_components + n_vec_components;
      /** Index of the density component \f$ \rho \f$. */
      static constexpr unsigned int density_component = 0;
      /**
       * Starting index of the momentum components \f$ \mathbf{p} \f$
       * (@ref n_vec_components components).
       */
      static constexpr unsigned int first_momentum_component = 1;
      /** Index of the energy component \f$ \mathcal{E} \f$. */
      static constexpr unsigned int energy_component = n_vec_components + 1;
      /**
       * Starting index of the magnetic components
       * \f$ \mathbf{b} = \frac{\mathbf{B}}{\sqrt{4\pi}} \f$
       * (@ref n_vec_components components).
       */
      static constexpr unsigned int first_magnetic_component =
        n_vec_components + 2;
      /**
       * Only in primitive states: Starting index of the velocity components
       * \f$ \mathbf{u} \f$ (@ref n_vec_components components).
       */
      static constexpr unsigned int first_velocity_component =
        first_momentum_component;
      /**
       * Only in primitive states:
       * Index of the pressure component \f$ P \f$.
       */
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
       *     \mathbf{b}
       *   \end{pmatrix} \,,
       * \f]
       * with \f$ \mathbf{b} = \frac{\mathbf{B}}{\sqrt{4\pi}} \f$.
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
       * Returns a list: `rho`, `p`, `E`, `b`.
       *
       * @param prefix Prefix for the component names.
       * @return std::vector<std::string> component_names.
       */
      static std::vector<std::string>
      create_component_name_list(const std::string &prefix = "");

      /**
       * @brief Create a list of component interpretation corresponding to the
       *        MHD state.
       *
       * Returns a list of
       * @dealref{DataComponentInterpretation,namespaceDataComponentInterpretation},
       * so the velocity and magnetic field components are interpreted as
       * vectors in the output.
       *
       * @return std::vector<DataComponentInterpretation>
       *         data_component_interpretation.
       */
      static std::vector<
        dealii::DataComponentInterpretation::DataComponentInterpretation>
      create_component_interpretation_list();
      /** @} */


      /** @{ */
      /**
       * @brief Compute the MHD flux matrix.
       *
       * @param state The @ref state_type "MHD state" in conservative form
       *              \f$ \mathbf{w} \f$.
       * @param flux_matrix The MHD flux \f$ \mathbf{F}(\mathbf{w}) \f$.
       */
      void
      compute_flux_matrix(const state_type &state,
                          flux_type        &flux_matrix) const;



      /**
       * @brief Computes the absolute value of maximum eigenvalue in normal
       *        direction.
       *
       * @param state The @ref state_type "MHD state" in conservative form
       *              \f$ \mathbf{w} \f$.
       * @param normal The normal vector \f$ \hat{\mathbf{n}} \f$.
       * @return double The absolute value of maximum eigenvalue in normal
       *         direction.
       */
      double
      compute_maximum_normal_eigenvalue(
        const state_type             &state,
        const dealii::Tensor<1, dim> &normal) const;



      /**
       * @brief Computes an upper limit for the maximum eigenvalue of the state.
       *
       * @param state The @ref state_type "MHD state" in conservative form
       *              \f$ \mathbf{w} \f$.
       * @return double Upper limit for the maximum eigenvalue.
       */
      double
      compute_maximum_eigenvalue(const state_type &state) const;
      /** @} */


      /** @{ */
      /**
       * @brief Compute the pressure \f$ P \f$ from the MHD state
       *        \f$ \mathbf{w} \f$.
       *
       * @param state The @ref state_type "MHD state" in conservative form
       *              \f$ \mathbf{w} \f$.
       * @return double The pressure \f$ P \f$.
       */
      double
      compute_pressure(const state_type &state) const;

      /**
       * @brief Same as @ref compute_pressure() but without pre- and
       *        post-conditions.
       *
       * @note Use this function with caution. It is intended for use in
       *       scenarios where negative energy or pressure warnings are handled
       *       elsewhere in the code.
       *
       * @param state The @ref state_type "MHD state" in conservative form
       *              \f$ \mathbf{w} \f$.
       * @return double The pressure \f$ P \f$.
       */
      double
      compute_pressure_unsafe(const state_type &state) const;


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
       * @param state The @ref state_type "MHD state" in conservative form
       *              \f$ \mathbf{w} \f$.
       * @param normal The normal vector \f$ \hat{\mathbf{n}} \f$.
       * @param eigenvalues List of the eigenvalues.
       *
       */
      void
      compute_normale_eigenvalues(const state_type             &state,
                                  const dealii::Tensor<1, dim> &normal,
                                  dealii::Vector<double> &eigenvalues) const;


      /**
       * @brief Compute the right eigenvector matrix \f$ \mathbf{R} \f$ for the
       *        MHD equations in conserved variables.
       *
       * This function computes the right eigenvector matrix \f$ \mathbf{R} \f$
       * in the normal direction. The matrix \f$ \mathbf{R} \f$ is defined as:
       * \f[
       *   \mathbf{R} =
       *   \begin{pmatrix}
       *      \mathbf{r}_{1} \,, &
       *      \mathbf{r}_{2} \,, &
       *      \dots          \,, &
       *      \mathbf{r}_{8}
       *   \end{pmatrix} \,,
       * \f]
       * where \f$ \mathbf{r}_{n} \f$ is the right eigenvector corresponding to
       * the eigenvalue \f$ \lambda_{n} \f$.
       *
       * The eigenvectors correspond to the following modes:
       * - \f$ \mathbf{r}_{1/8} \f$: Left- and right-going fast magneto-sonic
       *   modes.
       * - \f$ \mathbf{r}_{2/7} \f$: Alfven modes.
       * - \f$ \mathbf{r}_{3/6} \f$: Slow magneto-sonic modes.
       * - \f$ \mathbf{r}_{4} \f$: Density entropy mode.
       * - \f$ \mathbf{r}_{5} \f$: Unphysical \f$ \nabla \cdot \mathbf{B} \f$
       *   mode.
       *
       * @param state The @ref state_type "MHD state" in conservative form
       *              \f$ \mathbf{w} \f$.
       * @param normal The normal vector \f$ \hat{\mathbf{n}} \f$.
       * @param eigenvectors Returns a @dealref{FullMatrix} of the right
       *                     eigenvectors, \f$ \mathbf{R} \f$.
       *
       * @see compute_normale_eigenvalues()
       */
      void
      compute_right_eigenvector_matrix(
        const state_type             &state,
        const dealii::Tensor<1, dim> &normal,
        dealii::FullMatrix<double>   &eigenvectors) const;

      /**
       * @brief Compute the left eigenvector matrix \f$ \mathbf{L} \f$ for the
       *        MHD equations in conserved variables.
       *
       * This function computes the left eigenvector matrix \f$ \mathbf{L} \f$
       * in the normal direction. The matrix \f$ \mathbf{L} \f$ is defined as:
       * \f[
       *   \mathbf{L} =
       *   \begin{pmatrix}
       *      \mathbf{l}_{1} \\
       *      \mathbf{l}_{2} \\
       *      \vdots         \\
       *      \mathbf{l}_{8}
       *   \end{pmatrix} \,,
       * \f]
       * where \f$ \mathbf{l}_{n} \f$ is the left eigenvector (row-vector)
       * corresponding to the eigenvalue \f$ \lambda_{n} \f$.
       *
       * The vectors are normalized such that
       * \f[
       *   \mathbf{L} \mathbf{R} = \mathbf{1} \,.
       * \f]
       *
       * @param state The @ref state_type "MHD state" in conservative form
       *              \f$ \mathbf{w} \f$.
       * @param normal The normal vector \f$ \hat{\mathbf{n}} \f$.
       * @param eigenvectors Returns a @dealref{FullMatrix} of the left
       *                     eigenvectors, \f$ \mathbf{L} \f$.
       *
       * @see compute_normale_eigenvalues()
       * @see compute_right_eigenvector_matrix()
       */
      void
      compute_left_eigenvector_matrix(
        const state_type             &state,
        const dealii::Tensor<1, dim> &normal,
        dealii::FullMatrix<double>   &eigenvectors) const;

      /**
       * @brief Computes the transformation matrices to convert between
       *        conserved variables \f$ \mathbf{w} \f$ and characteristic
       *        variables \f$ \omega \f$ in each dimension \f$ j \f$.
       *
       * @param state Conserved state \f$ \mathbf{w} \f$ to use for the
       *        linearized computation of the eigenstates.
       * @param left_matrices Array of left eigenvector matrices
       *        \f$ \mathbf{L}_j \f$ in each dimension \f$ j \f$.
       * @param right_matrices Array of right eigenvector matrices
       *        \f$ \mathbf{R}_j \f$ in each dimension \f$ j \f$.
       *
       * @see compute_right_eigenvector_matrix()
       * @see compute_left_eigenvector_matrix()
       */
      void
      compute_transformation_matrices(
        const state_type                            &state,
        std::array<dealii::FullMatrix<double>, dim> &left_matrices,
        std::array<dealii::FullMatrix<double>, dim> &right_matrices) const;
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
       *            \mathbf{b}
       *          \end{pmatrix} \,.
       *        \f]
       * @param conserved_state Returns the conserved state
       *        \f[
       *          \mathbf{w} =
       *          \begin{pmatrix}
       *            \rho        \\
       *            \mathbf{p}  \\
       *            \mathcal{E} \\
       *            \mathbf{b}
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
       *            \mathbf{b}
       *          \end{pmatrix} \,.
       *        \f]
       * @param primitive_state Returns the primitive state
       *        \f[
       *          \tilde{\mathbf{w}} =
       *          \begin{pmatrix}
       *            \rho        \\
       *            \mathbf{u}  \\
       *            P \\
       *            \mathbf{b}
       *          \end{pmatrix} \,.
       *        \f]
       */
      void
      convert_conserved_to_primitive(const state_type &conserved_state,
                                     state_type       &primitive_state) const;



      /**
       * @brief Convert characteristic to conserved state.
       *
       * @param characteristic_state Characteristic state
       *        \f[
       *          \omega_{\hat{\mathbf{n}}} =
       *          \begin{pmatrix}
       *            \omega_{1, \hat{\mathbf{n}}} \\
       *            \omega_{2, \hat{\mathbf{n}}} \\
       *            \vdots                       \\
       *            \omega_{8, \hat{\mathbf{n}}} \\
       *          \end{pmatrix} \,.
       *        \f]
       * @param normal Normal vector \f$ \hat{\mathbf{n}} \f$
       * @param conserved_state Returns the conserved state
       *        \f[
       *          \mathbf{w} =
       *          \begin{pmatrix}
       *            \rho        \\
       *            \mathbf{p}  \\
       *            \mathcal{E} \\
       *            \mathbf{b}
       *          \end{pmatrix} \,.
       *        \f]
       *       This state is also used as input to calculate the transformation
       *       matrix.
       *
       * @see compute_right_eigenvector_matrix()
       */
      void
      convert_characteristic_to_conserved(
        const state_type             &characteristic_state,
        const dealii::Tensor<1, dim> &normal,
        state_type                   &conserved_state) const;


      /**
       * @brief Convert conserved to characteristic state.
       *
       * @param conserved_state Conserved state
       *        \f[
       *          \mathbf{w} =
       *          \begin{pmatrix}
       *            \rho        \\
       *            \mathbf{p}  \\
       *            \mathcal{E} \\
       *            \mathbf{b}
       *          \end{pmatrix} \,.
       *        \f]
       * @param normal Normal vector \f$ \hat{\mathbf{n}} \f$
       * @param characteristic_state Returns the characteristic state
       *        \f[
       *          \omega_{\hat{\mathbf{n}}} =
       *          \begin{pmatrix}
       *            \omega_{1, \hat{\mathbf{n}}} \\
       *            \omega_{2, \hat{\mathbf{n}}} \\
       *            \vdots                       \\
       *            \omega_{8, \hat{\mathbf{n}}} \\
       *          \end{pmatrix} \,.
       *        \f]
       *
       * @see compute_left_eigenvector_matrix()
       */
      void
      convert_conserved_to_characteristic(
        const state_type             &conserved_state,
        const dealii::Tensor<1, dim> &normal,
        state_type                   &characteristic_state) const;


      /**
       * @brief Convert gradient in characteristic to conserved variables.
       *
       * @param characteristic_gradient Gradient in characteristic variables,
       *        \f$ \partial_j \omega \f$.
       * @param right_matrices Array of right eigenvector matrices
       *        \f$ \mathbf{R}_j \f$ in each space direction \f$ j \f$.
       * @param conserved_gradient Returns the gradient in conserved variables,
       *        \f$ \partial_j \mathbf{w} \f$.
       *
       * @see compute_transformation_matrices()
       */
      void
      convert_gradient_characteristic_to_conserved(
        const flux_type                             &characteristic_gradient,
        std::array<dealii::FullMatrix<double>, dim> &right_matrices,
        flux_type                                   &conserved_gradient) const;


      /**
       * @brief Convert gradient in conserved to characteristic variables.
       *
       * @param conserved_gradient Gradient in conserved variables,
       *        \f$ \partial_j \mathbf{w} \f$.
       * @param left_matrices Array of left eigenvector matrices
       *        \f$ \mathbf{L}_j \f$ in each space direction \f$ j \f$.
       * @param characteristic_gradient Returns the gradient in characteristic
       *        variables, \f$ \partial_j \omega \f$.
       *
       * @see compute_transformation_matrices()
       */
      void
      convert_gradient_conserved_to_characteristic(
        const flux_type                             &conserved_gradient,
        std::array<dealii::FullMatrix<double>, dim> &left_matrices,
        flux_type &characteristic_gradient) const;
      /** @} */
    };
  } // namespace MHD
} // namespace sapphirepp
#endif
