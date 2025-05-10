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
 * @file mhd-postprocessor.h
 * @author Florian Schulze (florian.schulze@mpi-hd.mpg.de)
 * @brief Define @ref sapphirepp::MHD::MHDPostprocessor
 */

#ifndef MHD_POSTPROCESSOR_H
#define MHD_POSTPROCESSOR_H

#include <deal.II/base/exceptions.h>

#include <deal.II/lac/vector.h>

#include <deal.II/numerics/data_component_interpretation.h>
#include <deal.II/numerics/data_postprocessor.h>

#include <string>
#include <vector>

#include "mhd-equations.h"

namespace sapphirepp
{
  namespace MHD
  {
    using namespace dealii;



    /**
     * @brief Postprocessor unit for the MHD equations.
     *
     * @tparam dim Dimension of the configuration space \f$ (\mathbf{x}) \f$,
     *         `dim`
     * @tparam divergence_cleaning Use Lagrange multiplier \f$ \psi \f$
     *         for hyperbolic divergence cleaning.
     * @see @dealref{DataPostprocessor}
     */
    template <unsigned int dim, bool divergence_cleaning>
    class MHDPostprocessor : public DataPostprocessor<dim>
    {
    public:
      /** @ref MHDEquations::n_components */
      static constexpr unsigned int n_components =
        MHDEquations<dim, divergence_cleaning>::n_components;
      /** @ref MHDEquations::n_vec_components */
      static constexpr unsigned int n_vec_components =
        MHDEquations<dim, divergence_cleaning>::n_vec_components;
      /** @ref MHDEquations::density_component */
      static constexpr unsigned int density_component =
        MHDEquations<dim, divergence_cleaning>::density_component;
      /** @ref MHDEquations::first_momentum_component */
      static constexpr unsigned int first_momentum_component =
        MHDEquations<dim, divergence_cleaning>::first_momentum_component;

      /** @ref MHDEquations::state_type */
      using state_type =
        typename MHDEquations<dim, divergence_cleaning>::state_type;

      /** Number of components of the postprocessor output */
      static constexpr unsigned int n_components_out = 1 + n_vec_components;
      /** Outputindex of the pressure component \f$ P \f$. */
      static constexpr unsigned int pressure_component_out = 0;
      /**
       * Starting outputindex of the velocity components \f$ \mathbf{u} \f$
       * (@ref n_vec_components components).
       */
      static constexpr unsigned int first_velocity_component_out = 1;



      /**
       * @brief Constructor
       *
       * @param mhd_equations Instance of the underlying @ref MHDEquations.
       * @param prefix Prefix for the variable names.
       */
      MHDPostprocessor(
        const MHDEquations<dim, divergence_cleaning> &mhd_equations,
        const std::string                            &prefix = "");



      /**
       * @brief Create a list of the variable names
       *
       * @return std::vector<std::string> variable_names
       */
      std::vector<std::string>
      get_names() const override;



      /**
       * @brief Create a list of component interpretation
       *        for the output variables.
       *
       * @return std::vector<DataComponentInterpretation>
       *         data_component_interpretation.
       */
      std::vector<
        dealii::DataComponentInterpretation::DataComponentInterpretation>
      get_data_component_interpretation() const override;



      /**
       * @brief Return, which data has to be provided
       *        to compute the derived quantities.
       *
       * @return dealii::UpdateFlags
       * @see @dealref{DataPostprocessor.get_needed_update_flags(),classDataPostprocessor,aadecdd040447b395164397ea1196f721}
       */
      dealii::UpdateFlags
      get_needed_update_flags() const override;



      /**
       * @brief Evaluates the derived quantities.
       *
       * @param inputs @dealref{DataPostprocessorInputs,structDataPostprocessorInputs_1_1Vector}
       * @param computed_quantities Derived quantities
       * @see @dealref{DataPostprocessor.evaluate_vector_field(),classDataPostprocessor,ac907e98f8f03ea7e6ac25237271dc7b7}
       */
      void
      evaluate_vector_field(
        const DataPostprocessorInputs::Vector<dim> &inputs,
        std::vector<Vector<double>> &computed_quantities) const override;



    private:
      /** @ref MHDEquations */
      const MHDEquations<dim, divergence_cleaning> &mhd_equations;

      /** Prefix for the variable names */
      const std::string prefix;
    };
  } // namespace MHD
} // namespace sapphirepp
#endif
