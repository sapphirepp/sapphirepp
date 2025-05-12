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
 * @file slope-limiter.cpp
 * @author Florian Schulze (florian.schulze@mpi-hd.mpg.de)
 * @brief Implement utility functions related to slope limiting
 */

#include "slope-limiter.h"

#include <deal.II/base/exceptions.h>

#include <cmath>



template <unsigned int dim, bool divergence_cleaning>
sapphirepp::MHD::SlopeLimiter<dim, divergence_cleaning>::SlopeLimiter(
  const MHDParameters<dim> &mhd_parameters)
  : minmod_threshold{mhd_parameters.minmod_threshold}
  , minmod_beta{mhd_parameters.minmod_beta}
{}



template <unsigned int dim, bool divergence_cleaning>
double
sapphirepp::MHD::SlopeLimiter<dim, divergence_cleaning>::minmod(
  const std::vector<double> &values)
{
  auto [min_it, max_it] = std::minmax_element(values.begin(), values.end());
  if ((*min_it) * (*max_it) < 0.0)
    return 0.0;
  else if (std::abs(*min_it) < std::abs(*max_it))
    return *min_it;
  else
    return *max_it;
}



template <unsigned int dim, bool divergence_cleaning>
double
sapphirepp::MHD::SlopeLimiter<dim, divergence_cleaning>::minmod_gradients(
  const flux_type              &cell_gradient,
  const std::vector<flux_type> &neighbor_gradients,
  flux_type                    &limited_gradient,
  const double                  dx)
{
  double              difference = 0.;
  std::vector<double> values;
  values.reserve(neighbor_gradients.size() + 1);
  for (unsigned int c = 0; c < n_components; ++c)
    {
      for (unsigned int d = 0; d < dim; ++d)
        {
          AssertIsFinite(cell_gradient[c][d]);
          // if (std::abs(cell_gradient[c][d]) < minmod_threshold * dx * dx)
          if (std::abs(cell_gradient[c][d]) < minmod_threshold * dx)
            {
              limited_gradient[c][d] = cell_gradient[c][d];
              continue;
            }

          values.push_back(cell_gradient[c][d]);

          for (const auto &tmp : neighbor_gradients)
            if (std::isfinite(tmp[c][d]))
              values.push_back(minmod_beta * tmp[c][d]);

          limited_gradient[c][d] = minmod(values);
          AssertIsFinite(limited_gradient[c][d]);

          difference += std::fabs(cell_gradient[c][d] - limited_gradient[c][d]);

          values.clear();
        }
    }

  difference /= static_cast<double>(n_components);
  return difference;
}



template <unsigned int dim, bool divergence_cleaning>
void
sapphirepp::MHD::SlopeLimiter<dim, divergence_cleaning>::
  enforce_divergence_free_limited_gradient(flux_type &limited_gradient)
{
  double delta_p = 0;
  double delta_m = 0;

  for (unsigned int d = 0; d < dim; ++d)
    {
      delta_p +=
        std::max(limited_gradient[first_magnetic_component + d][d], 0.);
      delta_m +=
        std::max(-limited_gradient[first_magnetic_component + d][d], 0.);
    }

  const double delta = delta_p - delta_m;

  if (delta > 0)
    {
      for (unsigned int d = 0; d < dim; ++d)
        if (limited_gradient[first_magnetic_component + d][d] > 0)
          limited_gradient[first_magnetic_component + d][d] *=
            delta_m / delta_p;
    }
  else
    {
      for (unsigned int d = 0; d < dim; ++d)
        if (limited_gradient[first_magnetic_component + d][d] < 0)
          limited_gradient[first_magnetic_component + d][d] *=
            delta_p / delta_m;
    }
}



template <unsigned int dim, bool divergence_cleaning>
void
sapphirepp::MHD::SlopeLimiter<dim, divergence_cleaning>::
  limited_solution_to_dof_values(
    const state_type                  &cell_avg,
    const flux_type                   &limited_gradient,
    const std::vector<Vector<double>> &support_point_values,
    const FESystem<dim>               &fe,
    std::vector<double>               &cell_dof_values) const
{
  static_cast<void>(cell_avg);
  static_cast<void>(limited_gradient);

  fe.convert_generalized_support_point_values_to_dof_values(
    support_point_values, cell_dof_values);
}



template <unsigned int dim, bool divergence_cleaning>
void
sapphirepp::MHD::SlopeLimiter<dim, divergence_cleaning>::
  to_dof_values_using_primitive_support_points(
    const state_type              &cell_avg,
    const flux_type               &limited_gradient,
    const Point<dim>              &cell_center,
    const std::vector<Point<dim>> &support_points,
    const FESystem<dim>           &fe,
    std::vector<double>           &cell_dof_values) const
{
  AssertDimension(cell_dof_values.size(), fe.n_dofs_per_cell());

  std::vector<double>         base_dof_values;
  std::vector<Vector<double>> base_point_values;

  // loop over all base elements (respecting multiplicity) and let
  // them do the work on their share of the input argument

  unsigned int current_vector_component = 0;
  for (unsigned int base = 0; base < fe.n_base_elements(); ++base)
    {
      // We need access to the base_element, its multiplicity, the
      // number of generalized support points (n_base_points) and
      // the number of components we're dealing with.
      const FiniteElement<dim> &base_element = fe.base_element(base);
      const unsigned int        multiplicity = fe.element_multiplicity(base);
      const unsigned int        n_base_dofs  = base_element.n_dofs_per_cell();
      const unsigned int        n_base_components = base_element.n_components();

      // If the number of base degrees of freedom is zero, there is
      // nothing to do, skip the rest of the body in this case and
      // continue with the next element
      if (n_base_dofs == 0)
        {
          current_vector_component += multiplicity * n_base_components;
          continue;
        }

      Assert(base_element.is_primitive() &&
               base_element.has_generalized_support_points(),
             ExcMessage("Expect primitive FESystem with "
                        "generalized support points."));

      const size_t n_base_points =
        base_element.get_generalized_support_points().size();
      AssertDimension(support_points.size(), n_base_points);

      base_dof_values.resize(n_base_dofs);
      base_point_values.resize(n_base_points);

      for (unsigned int m = 0; m < multiplicity;
           ++m, current_vector_component += n_base_components)
        {
          for (unsigned int q_index = 0; q_index < n_base_points; ++q_index)
            {
              base_point_values[q_index].reinit(n_base_components, false);

              // Compute limited solution values at support points
              base_point_values[q_index] =
                cell_avg[current_vector_component] +
                limited_gradient[current_vector_component] *
                  (support_points[q_index] - cell_center);
            }

          base_element.convert_generalized_support_point_values_to_dof_values(
            base_point_values, base_dof_values);

          // Finally put these dof values back into cell dof
          // values vector.

          // To do this, we could really use a
          // base_to_system_index() function, but that doesn't
          // exist -- just do it by using the reverse table -- the
          // amount of work done here is not worth trying to
          // optimizing this.
          for (unsigned int i = 0; i < fe.n_dofs_per_cell(); ++i)
            if (fe.system_to_base_index(i).first == std::make_pair(base, m))
              cell_dof_values[i] =
                base_dof_values[fe.system_to_base_index(i).second];
        }
    }
}



// Explicit instantiations
template class sapphirepp::MHD::SlopeLimiter<1, false>;
template class sapphirepp::MHD::SlopeLimiter<1, true>;
template class sapphirepp::MHD::SlopeLimiter<2, false>;
template class sapphirepp::MHD::SlopeLimiter<2, true>;
template class sapphirepp::MHD::SlopeLimiter<3, false>;
template class sapphirepp::MHD::SlopeLimiter<3, true>;
