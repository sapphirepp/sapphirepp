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
 * @file debug-input-functions.h
 * @author Nils Schween (nils.schween@mpi-hd.mpg.de)
 * @author Florian Schulze (florian.schulze@mpi-hd.mpg.de)
 * @brief Define and implement @ref sapphirepp::VFP::DebugInputFunctions
 */

#ifndef VFP_DEBUGINPUTFUNCTIONS_H
#define VFP_DEBUGINPUTFUNCTIONS_H

#include <deal.II/base/exceptions.h>
#include <deal.II/base/function.h>
#include <deal.II/base/point.h>

#include <string>
#include <vector>

#include "config.h"
#include "pde-system.h"

namespace sapphirepp
{
  namespace VFP
  {

    /**
     * @brief This class collects the user-defined input functions
     *         to output them as debug information.
     *
     * @tparam dim Dimension of the reduced phase space \f$ (\mathbf{x}, p) \f$,
     *         `dim_ps`
     */
    template <unsigned int dim>
    class DebugInputFunctions : public dealii::Function<dim>
    {
    private:
      /**
       * Number of components in a vector,
       * e.g. for magnetic field.
       */
      static constexpr unsigned int n_vec_components = 3;



    public:
      /**
       * @brief Constructor
       *
       * @param scattering_frequency Scattering frequency \f$ \nu \f$
       * @param source_function Source term \f$ S \f$
       * @param magnetic_field Magnetic field \f$ \mathbf{B} \f$
       * @param background_velocity_field Background velocity field
       *                                  \f$ \mathbf{u} \f$
       */
      DebugInputFunctions(
        const ScatteringFrequency<dim>     &scattering_frequency,
        const Source<dim>                  &source_function,
        const MagneticField<dim>           &magnetic_field,
        const BackgroundVelocityField<dim> &background_velocity_field)
        : dealii::Function<dim>(
            scattering_frequency.n_components + source_function.n_components +
            magnetic_field.n_components +
            background_velocity_field.n_components + 1 + n_vec_components +
            n_vec_components * n_vec_components)
        , scattering_frequency{scattering_frequency}
        , source_function{source_function}
        , magnetic_field{magnetic_field}
        , background_velocity_field{background_velocity_field}
      {}



      /**
       * @brief Create a list of component names
       *
       * @param prefix Prefix for the component names.
       * @return std::vector<std::string> component_names
       */
      std::vector<std::string>
      create_component_name_list(const std::string &prefix = "func_")
      {
        std::vector<std::string> component_names(this->n_components);
        const char   vec_component_name[n_vec_components] = {'x', 'y', 'z'};
        unsigned int i_component                          = 0;

        component_names[i_component] = prefix + "nu";
        ++i_component;

        std::vector<std::string> source_components_names =
          PDESystem::create_component_name_list(source_function.n_components,
                                                prefix + "S_");
        for (unsigned int i = 0; i < source_function.n_components; ++i)
          {
            component_names[i_component] = source_components_names[i];
            ++i_component;
          }

        for (unsigned int d = 0; d < n_vec_components; ++d)
          {
            component_names[i_component] =
              prefix + "B_" + vec_component_name[d];
            ++i_component;
          }

        for (unsigned int d = 0; d < n_vec_components; ++d)
          {
            component_names[i_component] =
              prefix + "u_" + vec_component_name[d];
            ++i_component;
          }
        component_names[i_component] = prefix + "div_u";
        ++i_component;
        for (unsigned int d = 0; d < n_vec_components; ++d)
          {
            component_names[i_component] =
              prefix + "DT_u_" + vec_component_name[d];
            ++i_component;
          }
        for (unsigned int i = 0; i < n_vec_components; ++i)
          for (unsigned int j = 0; j < n_vec_components; ++j)
            {
              component_names[i_component] = prefix + "du" +
                                             vec_component_name[i] + "_d" +
                                             vec_component_name[j];
              ++i_component;
            }

        AssertDimension(i_component, this->n_components);

        return component_names;
      }



      /**
       * @brief Evaluate the input functions
       *
       * @param point Point in reduced phase space
       * @param values Return vector of the input functions
       *               in the same order as in @ref create_component_name_list()
       */
      void
      vector_value(const dealii::Point<dim> &point,
                   dealii::Vector<double>   &values) const override
      {
        AssertDimension(values.size(), this->n_components);
        std::vector<dealii::Point<dim>> points(1);
        points[0]                = point;
        unsigned int i_component = 0;


        std::vector<double> scattering_frequency_value(1);
        scattering_frequency.value_list(points, scattering_frequency_value);
        values[i_component] = scattering_frequency_value[0];
        ++i_component;


        dealii::Vector<double> source_function_value(
          source_function.n_components);
        source_function.vector_value(point, source_function_value);
        for (unsigned int d = 0; d < source_function.n_components; ++d)
          {
            values[i_component] = source_function_value[d];
            ++i_component;
          }


        dealii::Vector<double> magnetic_field_value(n_vec_components);
        magnetic_field.vector_value(point, magnetic_field_value);
        for (unsigned int d = 0; d < n_vec_components; ++d)
          {
            values[i_component] = magnetic_field_value[d];
            ++i_component;
          }


        dealii::Vector<double> velocity_value(n_vec_components);
        background_velocity_field.vector_value(point, velocity_value);
        for (unsigned int d = 0; d < n_vec_components; ++d)
          {
            values[i_component] = velocity_value[d];
            ++i_component;
          }

        std::vector<double> velocity_divergence(1);
        background_velocity_field.divergence_list(points, velocity_divergence);
        values[i_component] = velocity_divergence[0];
        ++i_component;

        std::vector<dealii::Vector<double>> velocity_material_derivatives(
          1, dealii::Vector<double>(n_vec_components));
        background_velocity_field.material_derivative_list(
          points, velocity_material_derivatives);
        for (unsigned int d = 0; d < n_vec_components; ++d)
          {
            values[i_component] = velocity_material_derivatives[0][d];
            ++i_component;
          }

        std::vector<dealii::FullMatrix<double>> velocity_jacobians(
          1, dealii::FullMatrix<double>(n_vec_components));
        background_velocity_field.jacobian_list(points, velocity_jacobians);
        for (unsigned int i = 0; i < n_vec_components; ++i)
          for (unsigned int j = 0; j < n_vec_components; ++j)
            {
              values[i_component] = velocity_jacobians[0][i][j];
              ++i_component;
            }


        AssertDimension(i_component, this->n_components);
      }



    private:
      /** Scattering frequency \f$ \nu \f$  */
      const ScatteringFrequency<dim> &scattering_frequency;
      /** Source term \f$ S \f$ */
      const Source<dim> &source_function;
      /** Magnetic field \f$ \mathbf{B} \f$ */
      const MagneticField<dim> &magnetic_field;
      /** Background velocity field \f$ \mathbf{u} \f$ */
      const BackgroundVelocityField<dim> &background_velocity_field;
    };
  } // namespace VFP
} // namespace sapphirepp
#endif
