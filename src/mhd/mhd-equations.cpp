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
 * @file mhd-equations.cpp
 * @author Florian Schulze (florian.schulze@mpi-hd.mpg.de)
 * @brief Implement @ref sapphirepp::MHD::MHDEquations
 */

#include "mhd-equations.h"

#include <deal.II/base/exceptions.h>
#include <deal.II/base/utilities.h>

#include <algorithm>
#include <cmath>
#include <limits>

#include "mhd-parameters.h"
#include "sapphirepp-logstream.h"



/** @ref sapphirepp::MHD::MHDParameters<1>::epsilon_d */
static constexpr double epsilon_d =
  sapphirepp::MHD::MHDParameters<1>::epsilon_d;



template <unsigned int dim, bool divergence_cleaning>
sapphirepp::MHD::MHDEquations<dim, divergence_cleaning>::MHDEquations(
  const MHDParameters<dim> &mhd_parameters)
  : adiabatic_index{mhd_parameters.adiabatic_index}
  , divergence_cleaning_Ch{mhd_parameters.divergence_cleaning_Ch}
  , divergence_cleaning_Cr{mhd_parameters.divergence_cleaning_Cr}
  , divergence_cleaning_speed{1.}
  , divergence_cleaning_damping{0.}
{
  AssertThrow(adiabatic_index >= 1.0,
              dealii::ExcMessage("Adiabatic index must be larger than 1."));
  AssertThrow(adiabatic_index <= 2.0,
              dealii::ExcMessage("Adiabatic index must be smaller than 2."));

  AssertThrow(divergence_cleaning_Ch >= 0.,
              dealii::ExcMessage("C_h must be larger than 0."));
  AssertThrow(divergence_cleaning_Ch <= 1.,
              dealii::ExcMessage("C_h must be smaller than 1."));
  AssertThrow(divergence_cleaning_Cr >= 0.,
              dealii::ExcMessage("C_r must be positive."));
}



template <unsigned int dim, bool divergence_cleaning>
std::vector<std::string>
sapphirepp::MHD::MHDEquations<dim, divergence_cleaning>::
  create_component_name_list(const std::string &prefix)
{
  std::vector<std::string> component_names(n_components);
  const char vec_component_name[n_vec_components] = {'x', 'y', 'z'};

  component_names[density_component] = prefix + "rho";
  component_names[energy_component]  = prefix + "E";
  if constexpr (divergence_cleaning)
    component_names[divergence_cleaning_component] = prefix + "psi";
  for (unsigned int d = 0; d < dim; ++d)
    {
      component_names[first_momentum_component + d] = prefix + "p";
      component_names[first_magnetic_component + d] = prefix + "b";
    }
  for (unsigned int d = dim; d < n_vec_components; ++d)
    {
      component_names[first_momentum_component + d] =
        prefix + "p_" + vec_component_name[d];
      component_names[first_magnetic_component + d] =
        prefix + "b_" + vec_component_name[d];
    }

  return component_names;
}



template <unsigned int dim, bool divergence_cleaning>
std::vector<dealii::DataComponentInterpretation::DataComponentInterpretation>
sapphirepp::MHD::MHDEquations<dim, divergence_cleaning>::
  create_component_interpretation_list()
{
  std::vector<dealii::DataComponentInterpretation::DataComponentInterpretation>
    data_component_interpretation(n_components);

  data_component_interpretation[density_component] =
    dealii::DataComponentInterpretation::component_is_scalar;
  data_component_interpretation[energy_component] =
    dealii::DataComponentInterpretation::component_is_scalar;
  if constexpr (divergence_cleaning)
    data_component_interpretation[divergence_cleaning_component] =
      dealii::DataComponentInterpretation::component_is_scalar;
  // Interpret components in the plain as vector components
  for (unsigned int d = 0; d < dim; ++d)
    {
      data_component_interpretation[first_momentum_component + d] =
        dealii::DataComponentInterpretation::component_is_part_of_vector;
      data_component_interpretation[first_magnetic_component + d] =
        dealii::DataComponentInterpretation::component_is_part_of_vector;
    }
  // Interpret components out of the plain as scalars
  for (unsigned int d = dim; d < n_vec_components; ++d)
    {
      data_component_interpretation[first_momentum_component + d] =
        dealii::DataComponentInterpretation::component_is_scalar;
      data_component_interpretation[first_magnetic_component + d] =
        dealii::DataComponentInterpretation::component_is_scalar;
    }

  return data_component_interpretation;
}



template <unsigned int dim, bool divergence_cleaning>
void
sapphirepp::MHD::MHDEquations<dim, divergence_cleaning>::compute_flux_matrix(
  const state_type &state,
  flux_type        &flux_matrix) const
{
  AssertDimension(state.size(), n_components);

  const double pressure = compute_pressure(state);
  double       b2       = 0.;
  double       pb       = 0.;
  for (unsigned int d = 0; d < n_vec_components; ++d)
    {
      b2 += state[first_magnetic_component + d] *
            state[first_magnetic_component + d];
      pb += state[first_momentum_component + d] *
            state[first_magnetic_component + d];
    }

  for (unsigned int j = 0; j < dim; ++j)
    {
      flux_matrix[density_component][j] = state[first_momentum_component + j];

      flux_matrix[energy_component][j] =
        state[first_momentum_component + j] / state[density_component] *
          (state[energy_component] + pressure + 0.5 * b2) -
        1. / (state[density_component]) * pb *
          state[first_magnetic_component + j];

      for (unsigned int i = 0; i < n_vec_components; ++i)
        {
          flux_matrix[first_momentum_component + i][j] =
            state[first_momentum_component + j] *
              state[first_momentum_component + i] / state[density_component] -
            state[first_magnetic_component + j] *
              state[first_magnetic_component + i];

          flux_matrix[first_magnetic_component + i][j] =
            (state[first_momentum_component + j] *
               state[first_magnetic_component + i] -
             state[first_momentum_component + i] *
               state[first_magnetic_component + j]) /
            state[density_component];
        }

      flux_matrix[first_momentum_component + j][j] += pressure + 0.5 * b2;

      if constexpr (divergence_cleaning)
        {
          Assert(
            divergence_cleaning_speed > 0.,
            dealii::ExcMessage(
              "Speed for hyperbolic divergence cleaning must be positive."));

          flux_matrix[first_magnetic_component + j][j] =
            state[divergence_cleaning_component];

          flux_matrix[divergence_cleaning_component][j] =
            divergence_cleaning_speed * divergence_cleaning_speed *
            state[first_magnetic_component + j];
        }
    }
}



template <unsigned int dim, bool divergence_cleaning>
void
sapphirepp::MHD::MHDEquations<dim, divergence_cleaning>::
  add_source_divergence_cleaning(const state_type &state,
                                 state_type       &source) const
{
  if constexpr (!divergence_cleaning)
    {
      Assert(false,
             dealii::ExcMessage("You try to add the divergence cleaning source "
                                "but divergence cleaning is deactivated."));
      return;
    }
  AssertDimension(state.size(), n_components);
  AssertDimension(source.size(), n_components);

  source[divergence_cleaning_component] +=
    -divergence_cleaning_damping * state[divergence_cleaning_component];
}



template <unsigned int dim, bool divergence_cleaning>
double
sapphirepp::MHD::MHDEquations<dim, divergence_cleaning>::
  compute_maximum_normal_eigenvalue(const state_type             &state,
                                    const dealii::Tensor<1, dim> &normal) const
{
  AssertDimension(state.size(), n_components);
  Assert(std::abs(normal.norm() - 1) < epsilon_d,
         dealii::ExcMessage("Normal vector must be normalized."));

  const double pressure = compute_pressure(state);
  double       b2       = 0.;
  double       nb       = 0.;
  double       nu       = 0.;
  for (unsigned int d = 0; d < dim; ++d)
    {
      b2 += state[first_magnetic_component + d] *
            state[first_magnetic_component + d];
      nu += normal[d] * state[first_momentum_component + d] /
            state[density_component];
      nb += normal[d] * state[first_magnetic_component + d];
    }
  for (unsigned int d = dim; d < n_vec_components; ++d)
    {
      b2 += state[first_magnetic_component + d] *
            state[first_magnetic_component + d];
    }

  const double a_s2 = adiabatic_index * pressure / state[density_component];
  const double c_a2 = nb * nb / state[density_component];
  const double d_n  = (a_s2 + b2 / state[density_component]) *
                       (a_s2 + b2 / state[density_component]) -
                     4. * a_s2 * c_a2;
  Assert(a_s2 >= epsilon_d,
         ExcNonAdmissibleState(
           state,
           "Expect positive squared adiabatic sound speed, a_s^2 = " +
             dealii::Utilities::to_string(a_s2)));
  Assert(c_a2 >= 0.,
         ExcNonAdmissibleState(state,
                               "Negative squared alfven speed, c_a^2 = " +
                                 dealii::Utilities::to_string(c_a2)));
  Assert(d_n >= 0.,
         ExcNonAdmissibleState(state,
                               "Negative d_n = " +
                                 dealii::Utilities::to_string(d_n)));
  const double c_f2 =
    0.5 * (a_s2 + b2 / state[density_component] + std::sqrt(d_n));
  Assert(c_f2 >= epsilon_d,
         ExcNonAdmissibleState(state,
                               "Expect positive squared fast speed, c_f^2 = " +
                                 dealii::Utilities::to_string(c_f2)));

  return std::fabs(nu) + std::sqrt(c_f2);
  /** @todo Do we need to take divergence_cleaning_speed into account? */
}



template <unsigned int dim, bool divergence_cleaning>
double
sapphirepp::MHD::MHDEquations<dim, divergence_cleaning>::
  compute_maximum_eigenvalue(const state_type &state) const
{
  AssertDimension(state.size(), n_components);

  const double pressure = compute_pressure(state);
  double       p2       = 0.;
  double       b2       = 0.;
  for (unsigned int d = 0; d < n_vec_components; ++d)
    {
      p2 += state[first_momentum_component + d] *
            state[first_momentum_component + d];
      b2 += state[first_magnetic_component + d] *
            state[first_magnetic_component + d];
    }
  const double max_u = std::sqrt(p2) / state[density_component];

  const double a_s2 = adiabatic_index * pressure / state[density_component];
  Assert(a_s2 >= epsilon_d,
         ExcNonAdmissibleState(
           state,
           "Expect positive squared adiabatic sound speed, a_s^2 = " +
             dealii::Utilities::to_string(a_s2)));
  const double max_c_f2 = a_s2 + b2 / state[density_component];
  Assert(max_c_f2 >= epsilon_d,
         ExcNonAdmissibleState(state,
                               "Expect positive squared fast speed, c_f^2 = " +
                                 dealii::Utilities::to_string(max_c_f2)));

  return max_u + std::sqrt(max_c_f2);
}



template <unsigned int dim, bool divergence_cleaning>
double
sapphirepp::MHD::MHDEquations<dim, divergence_cleaning>::
  compute_pressure_unsafe(const state_type &state) const
{
  AssertDimension(state.size(), n_components);

  double p2 = 0.;
  double b2 = 0.;
  for (unsigned int d = 0; d < n_vec_components; ++d)
    {
      p2 += state[first_momentum_component + d] *
            state[first_momentum_component + d];
      b2 += state[first_magnetic_component + d] *
            state[first_magnetic_component + d];
    }

  const double pressure =
    (adiabatic_index - 1.) *
    (state[energy_component] - 1. / (2. * state[density_component]) * p2 -
     0.5 * b2);

  return pressure;
}



template <unsigned int dim, bool divergence_cleaning>
double
sapphirepp::MHD::MHDEquations<dim, divergence_cleaning>::compute_pressure(
  const state_type &state) const
{
  AssertDimension(state.size(), n_components);
  Assert(state[density_component] >= epsilon_d, ExcNonAdmissibleState(state));
  Assert(state[energy_component] >= epsilon_d, ExcNonAdmissibleState(state));

  const double pressure = compute_pressure_unsafe(state);

  Assert(pressure >= epsilon_d, ExcNonAdmissibleState(state));
  return pressure;
}



template <unsigned int dim, bool divergence_cleaning>
void
sapphirepp::MHD::MHDEquations<dim, divergence_cleaning>::
  compute_normale_eigenvalues(const state_type             &state,
                              const dealii::Tensor<1, dim> &normal,
                              dealii::Vector<double>       &eigenvalues) const
{
  AssertDimension(state.size(), n_components);
  AssertDimension(eigenvalues.size(), n_components);
  Assert(std::abs(normal.norm() - 1) < epsilon_d,
         dealii::ExcMessage("Normal vector must be normalized."));

  const double pressure = compute_pressure(state);
  double       b2       = 0.;
  double       nu       = 0.;
  double       nb       = 0.;
  for (unsigned int d = 0; d < dim; ++d)
    {
      b2 += state[first_magnetic_component + d] *
            state[first_magnetic_component + d];
      nu += normal[d] * state[first_momentum_component + d] /
            state[density_component];
      nb += normal[d] * state[first_magnetic_component + d];
    }
  for (unsigned int d = dim; d < n_vec_components; ++d)
    {
      b2 += state[first_magnetic_component + d] *
            state[first_magnetic_component + d];
    }

  const double a_s2 = adiabatic_index * pressure / state[density_component];
  const double c_a2 = nb * nb / state[density_component];
  const double d_n  = (a_s2 + b2 / state[density_component]) *
                       (a_s2 + b2 / state[density_component]) -
                     4. * a_s2 * c_a2;
  Assert(a_s2 >= epsilon_d,
         ExcNonAdmissibleState(
           state,
           "Expect positive squared adiabatic sound speed, a_s^2 = " +
             dealii::Utilities::to_string(a_s2)));
  Assert(c_a2 >= 0.,
         ExcNonAdmissibleState(state,
                               "Negative squared alfven speed, c_a^2 = " +
                                 dealii::Utilities::to_string(c_a2)));
  Assert(d_n >= 0.,
         ExcNonAdmissibleState(state,
                               "Negative d_n = " +
                                 dealii::Utilities::to_string(d_n)));
  const double c_s2 =
    0.5 * (a_s2 + b2 / state[density_component] - std::sqrt(d_n));
  const double c_f2 =
    0.5 * (a_s2 + b2 / state[density_component] + std::sqrt(d_n));
  Assert(c_s2 >= 0.,
         ExcNonAdmissibleState(state,
                               "Negative squared slow speed, c_s^2 = " +
                                 dealii::Utilities::to_string(c_s2)));
  Assert(c_f2 >= epsilon_d,
         ExcNonAdmissibleState(state,
                               "Expect positive squared fast speed, c_f^2 = " +
                                 dealii::Utilities::to_string(c_f2)));
  if constexpr (divergence_cleaning)
    Assert(divergence_cleaning_speed > 0.,
           dealii::ExcMessage(
             "Speed for hyperbolic divergence cleaning must be positive."));

  eigenvalues[0] = nu - std::sqrt(c_f2);
  eigenvalues[1] = nu - std::sqrt(c_a2);
  eigenvalues[2] = nu - std::sqrt(c_s2);
  eigenvalues[3] = nu;
  if constexpr (divergence_cleaning)
    eigenvalues[4] = -divergence_cleaning_speed;
  else
    eigenvalues[4] = 0.;
  eigenvalues[5] = nu + std::sqrt(c_s2);
  eigenvalues[6] = nu + std::sqrt(c_a2);
  eigenvalues[7] = nu + std::sqrt(c_f2);
  if constexpr (divergence_cleaning)
    eigenvalues[8] = divergence_cleaning_speed;
}



template <unsigned int dim, bool divergence_cleaning>
void
sapphirepp::MHD::MHDEquations<dim, divergence_cleaning>::
  compute_right_eigenvector_matrix(
    const state_type             &state,
    const dealii::Tensor<1, dim> &normal,
    dealii::FullMatrix<double>   &eigenvectors) const

{
  AssertDimension(state.size(), n_components);
  AssertDimension(eigenvectors.n(), n_components);
  AssertDimension(eigenvectors.m(), n_components);
  Assert(std::abs(normal.norm() - 1) < epsilon_d,
         dealii::ExcMessage("Normal vector must be normalized."));

  const double                        pressure = compute_pressure(state);
  dealii::Tensor<1, n_vec_components> normal_3d;
  dealii::Tensor<1, n_vec_components> u;
  double                              b2 = 0.;
  double                              u2 = 0.;
  double                              nu = 0.;
  double                              nb = 0.;
  for (unsigned int d = 0; d < dim; ++d)
    {
      normal_3d[d] = normal[d];
      u[d] = state[first_momentum_component + d] / state[density_component];
      b2 += state[first_magnetic_component + d] *
            state[first_magnetic_component + d];
      u2 += state[first_momentum_component + d] / state[density_component] *
            state[first_momentum_component + d] / state[density_component];
      nu += normal[d] * state[first_momentum_component + d] /
            state[density_component];
      nb += normal[d] * state[first_magnetic_component + d];
    }
  for (unsigned int d = dim; d < n_vec_components; ++d)
    {
      normal_3d[d] = 0.;
      u[d] = state[first_momentum_component + d] / state[density_component];
      b2 += state[first_magnetic_component + d] *
            state[first_magnetic_component + d];
      u2 += state[first_momentum_component + d] / state[density_component] *
            state[first_momentum_component + d] / state[density_component];
    }
  const double b_perp = std::sqrt(b2 - nb * nb);

  const double a_s2 = adiabatic_index * pressure / state[density_component];
  const double c_a2 = nb * nb / state[density_component];
  const double d_n  = (a_s2 + b2 / state[density_component]) *
                       (a_s2 + b2 / state[density_component]) -
                     4. * a_s2 * c_a2;
  Assert(a_s2 >= epsilon_d,
         ExcNonAdmissibleState(
           state,
           "Expect positive squared adiabatic sound speed, a_s^2 = " +
             dealii::Utilities::to_string(a_s2)));
  Assert(c_a2 >= 0.,
         ExcNonAdmissibleState(state,
                               "Negative squared alfven speed, c_a^2 = " +
                                 dealii::Utilities::to_string(c_a2)));
  Assert(d_n >= 0.,
         ExcNonAdmissibleState(state,
                               "Negative d_n = " +
                                 dealii::Utilities::to_string(d_n)));
  const double c_s2 =
    0.5 * (a_s2 + b2 / state[density_component] - std::sqrt(d_n));
  const double c_f2 =
    0.5 * (a_s2 + b2 / state[density_component] + std::sqrt(d_n));
  Assert(c_s2 >= 0.,
         ExcNonAdmissibleState(state,
                               "Negative squared slow speed, c_s^2 = " +
                                 dealii::Utilities::to_string(c_s2)));
  Assert(c_f2 >= epsilon_d,
         ExcNonAdmissibleState(state,
                               "Expect positive squared fast speed, c_f^2 = " +
                                 dealii::Utilities::to_string(c_f2)));
  const double a_s = std::sqrt(a_s2);
  const double c_a = std::sqrt(c_a2);
  const double c_s = std::sqrt(c_s2);
  const double c_f = std::sqrt(c_f2);

  double alp_s, alp_f;
  if ((c_f2 - c_s2) < epsilon_d)
    {
      alp_s = 1.;
      alp_f = 1.;
    }
  else
    {
      alp_s = std::sqrt((c_f2 - a_s2) / (c_f2 - c_s2));
      alp_f = std::sqrt((c_f2 - c_a2) / (c_f2 - c_s2));
    }
  Assert(alp_s >= 0.,
         ExcNonAdmissibleState(state,
                               "Expect non-negative value for alpha_s = " +
                                 dealii::Utilities::to_string(alp_s)));
  Assert(alp_f >= 0.,
         ExcNonAdmissibleState(state,
                               "Expect non-negative value for alpha_f = " +
                                 dealii::Utilities::to_string(alp_f)));
  if constexpr (divergence_cleaning)
    Assert(divergence_cleaning_speed > 0.,
           dealii::ExcMessage(
             "Speed for hyperbolic divergence cleaning must be positive."));

  // Construct perpendicular normal vector n_perp
  dealii::Tensor<1, n_vec_components> n_perp;
  if (b_perp >= epsilon_d)
    {
      for (unsigned int d = 0; d < n_vec_components; ++d)
        n_perp[d] = (state[first_magnetic_component + d] - nb * normal_3d[d]);
    }
  else
    {
      n_perp         = 0; // tmp for e_j
      unsigned int j = 0; // j = arg_min(normal_3d)
      for (unsigned int d = 1; d < n_vec_components; ++d)
        if (std::abs(normal_3d[d]) < std::abs(normal_3d[j]))
          j = d;
      n_perp[j] = 1.;
      n_perp    = dealii::cross_product_3d<n_vec_components>(n_perp, normal_3d);
    }
  n_perp /= n_perp.norm();

  Assert(std::abs(n_perp.norm() - 1) < epsilon_d,
         dealii::ExcMessage("n_perp is not normalized."));
  Assert(std::abs(normal_3d * n_perp) < epsilon_d,
         dealii::ExcMessage("n_perp is not perpendicular to normal"));

  const double                              n_perp_u = n_perp * u;
  const dealii::Tensor<1, n_vec_components> n_perp_cross_u =
    dealii::cross_product_3d<n_vec_components>(n_perp, u);
  const double                              sgn_nb = (nb >= 0.) ? 1. : -1.;
  const dealii::Tensor<1, n_vec_components> n_perp_cross_normal =
    dealii::cross_product_3d<n_vec_components>(n_perp, normal_3d);


  // Left fast magnetosonic mode r_1
  eigenvectors[density_component][0] = alp_f;
  eigenvectors[energy_component][0] =
    alp_f * u2 / 2. + alp_f * c_f2 / (adiabatic_index - 1.) +
    sgn_nb * alp_s * c_a * n_perp_u - alp_f * c_f * nu -
    (2. - adiabatic_index) / (adiabatic_index - 1.) * alp_f * (c_f2 - a_s2);
  for (unsigned int d = 0; d < n_vec_components; ++d)
    {
      eigenvectors[first_momentum_component + d][0] =
        alp_f * (u[d] - c_f * normal_3d[d]) + sgn_nb * alp_s * c_a * n_perp[d];
      eigenvectors[first_magnetic_component + d][0] =
        alp_s * c_f / std::sqrt(state[density_component]) * n_perp[d];
    }
  if constexpr (divergence_cleaning)
    eigenvectors[divergence_cleaning_component][0] = 0.;



  // Left alfven mode r_2
  eigenvectors[density_component][1] = 0;
  eigenvectors[energy_component][1]  = -sgn_nb * normal_3d * n_perp_cross_u;
  for (unsigned int d = 0; d < n_vec_components; ++d)
    {
      eigenvectors[first_momentum_component + d][1] =
        sgn_nb * n_perp_cross_normal[d];
      eigenvectors[first_magnetic_component + d][1] =
        n_perp_cross_normal[d] / std::sqrt(state[density_component]);
    }
  if constexpr (divergence_cleaning)
    eigenvectors[divergence_cleaning_component][1] = 0.;


  // Left slow magnetosonic mode r_3
  eigenvectors[density_component][2] = alp_s;
  eigenvectors[energy_component][2] =
    alp_s * u2 / 2. + alp_s * c_s2 / (adiabatic_index - 1.) -
    sgn_nb * alp_f * a_s * n_perp_u - alp_s * c_s * nu -
    (2. - adiabatic_index) / (adiabatic_index - 1.) * alp_s * (c_s2 - a_s2);
  for (unsigned int d = 0; d < n_vec_components; ++d)
    {
      eigenvectors[first_momentum_component + d][2] =
        alp_s * (u[d] - c_s * normal_3d[d]) - sgn_nb * alp_f * a_s * n_perp[d];
      eigenvectors[first_magnetic_component + d][2] =
        -alp_f * a_s2 / (std::sqrt(state[density_component]) * c_f) * n_perp[d];
    }
  if constexpr (divergence_cleaning)
    eigenvectors[divergence_cleaning_component][2] = 0.;


  // Density entropy mode r_4
  eigenvectors[density_component][3] = 1.;
  eigenvectors[energy_component][3]  = u2 / 2.;
  for (unsigned int d = 0; d < n_vec_components; ++d)
    {
      eigenvectors[first_momentum_component + d][3] = u[d];
      eigenvectors[first_magnetic_component + d][3] = 0;
    }
  if constexpr (divergence_cleaning)
    eigenvectors[divergence_cleaning_component][3] = 0.;


  // Magnetic divergence mode r_5
  /**
   * @note The eigenvector \f$ \mathbf{r}_{5} \f$
   * is not an eigenvector of the conserved system.
   * Instead it is choosen such that left and right eigenvectors
   * give an orthonormal system,
   * \f$ \mathbf{L} \mathbf{R} = \mathbf{1} \f$.
   */
  eigenvectors[density_component][4] = 0.;
  eigenvectors[energy_component][4]  = 0.;
  for (unsigned int d = 0; d < n_vec_components; ++d)
    {
      eigenvectors[first_momentum_component + d][4] = 0;
      eigenvectors[first_magnetic_component + d][4] = normal_3d[d];
    }
  if constexpr (divergence_cleaning)
    eigenvectors[divergence_cleaning_component][4] = -divergence_cleaning_speed;


  // Right slow magnetosonic mode r_6
  eigenvectors[density_component][5] = alp_s;
  eigenvectors[energy_component][5] =
    alp_s * u2 / 2. + alp_s * c_s2 / (adiabatic_index - 1.) +
    sgn_nb * alp_f * a_s * n_perp_u + alp_s * c_s * nu -
    (2. - adiabatic_index) / (adiabatic_index - 1.) * alp_s * (c_s2 - a_s2);
  for (unsigned int d = 0; d < n_vec_components; ++d)
    {
      eigenvectors[first_momentum_component + d][5] =
        alp_s * (u[d] + c_s * normal_3d[d]) + sgn_nb * alp_f * a_s * n_perp[d];
      eigenvectors[first_magnetic_component + d][5] =
        -alp_f * a_s2 / (std::sqrt(state[density_component]) * c_f) * n_perp[d];
    }
  if constexpr (divergence_cleaning)
    eigenvectors[divergence_cleaning_component][5] = 0.;


  // Right alfven mode r_7
  eigenvectors[density_component][6] = 0;
  eigenvectors[energy_component][6]  = sgn_nb * normal_3d * n_perp_cross_u;
  for (unsigned int d = 0; d < n_vec_components; ++d)
    {
      eigenvectors[first_momentum_component + d][6] =
        -sgn_nb * n_perp_cross_normal[d];
      eigenvectors[first_magnetic_component + d][6] =
        n_perp_cross_normal[d] / std::sqrt(state[density_component]);
    }
  if constexpr (divergence_cleaning)
    eigenvectors[divergence_cleaning_component][6] = 0.;


  // Right fast magnetosonic mode r_8
  eigenvectors[density_component][7] = alp_f;
  eigenvectors[energy_component][7] =
    alp_f * u2 / 2. + alp_f * c_f2 / (adiabatic_index - 1.) -
    sgn_nb * alp_s * c_a * n_perp_u + alp_f * c_f * nu -
    (2. - adiabatic_index) / (adiabatic_index - 1.) * alp_f * (c_f2 - a_s2);
  for (unsigned int d = 0; d < n_vec_components; ++d)
    {
      eigenvectors[first_momentum_component + d][7] =
        alp_f * (u[d] + c_f * normal_3d[d]) - sgn_nb * alp_s * c_a * n_perp[d];
      eigenvectors[first_magnetic_component + d][7] =
        alp_s * c_f / std::sqrt(state[density_component]) * n_perp[d];
    }
  if constexpr (divergence_cleaning)
    eigenvectors[divergence_cleaning_component][7] = 0.;


  if constexpr (divergence_cleaning)
    {
      // Magnetic divergence mode r_9
      /**
       * @note The eigenvector \f$ \mathbf{r}_{9} \f$
       * only exists with hyperbolic divergence cleaning.
       * Furthermore, it is not an eigenvector of the conserved system.
       * Instead it is choosen such that left and right eigenvectors
       * give an orthonormal system,
       * \f$ \mathbf{L} \mathbf{R} = \mathbf{1} \f$.
       */
      eigenvectors[density_component][8] = 0.;
      eigenvectors[energy_component][8]  = 0.;
      for (unsigned int d = 0; d < n_vec_components; ++d)
        {
          eigenvectors[first_momentum_component + d][8] = 0.;
          eigenvectors[first_magnetic_component + d][8] = normal_3d[d];
        }
      eigenvectors[divergence_cleaning_component][8] =
        divergence_cleaning_speed;
    }
}



template <unsigned int dim, bool divergence_cleaning>
void
sapphirepp::MHD::MHDEquations<dim, divergence_cleaning>::
  compute_left_eigenvector_matrix(
    const state_type             &state,
    const dealii::Tensor<1, dim> &normal,
    dealii::FullMatrix<double>   &eigenvectors) const

{
  AssertDimension(state.size(), n_components);
  AssertDimension(eigenvectors.n(), n_components);
  AssertDimension(eigenvectors.m(), n_components);
  Assert(std::abs(normal.norm() - 1) < epsilon_d,
         dealii::ExcMessage("Normal vector must be normalized."));

  const double                        pressure = compute_pressure(state);
  dealii::Tensor<1, n_vec_components> normal_3d;
  dealii::Tensor<1, n_vec_components> u;
  double                              b2 = 0.;
  double                              u2 = 0.;
  double                              nu = 0.;
  double                              nb = 0.;
  for (unsigned int d = 0; d < dim; ++d)
    {
      normal_3d[d] = normal[d];
      u[d] = state[first_momentum_component + d] / state[density_component];
      b2 += state[first_magnetic_component + d] *
            state[first_magnetic_component + d];
      u2 += state[first_momentum_component + d] / state[density_component] *
            state[first_momentum_component + d] / state[density_component];
      nu += normal[d] * state[first_momentum_component + d] /
            state[density_component];
      nb += normal[d] * state[first_magnetic_component + d];
    }
  for (unsigned int d = dim; d < n_vec_components; ++d)
    {
      normal_3d[d] = 0.;
      u[d] = state[first_momentum_component + d] / state[density_component];
      b2 += state[first_magnetic_component + d] *
            state[first_magnetic_component + d];
      u2 += state[first_momentum_component + d] / state[density_component] *
            state[first_momentum_component + d] / state[density_component];
    }
  const double b_perp = std::sqrt(b2 - nb * nb);

  const double a_s2 = adiabatic_index * pressure / state[density_component];
  const double c_a2 = nb * nb / state[density_component];
  const double d_n  = (a_s2 + b2 / state[density_component]) *
                       (a_s2 + b2 / state[density_component]) -
                     4. * a_s2 * c_a2;
  Assert(a_s2 >= epsilon_d,
         ExcNonAdmissibleState(
           state,
           "Expect positive squared adiabatic sound speed, a_s^2 = " +
             dealii::Utilities::to_string(a_s2)));
  Assert(c_a2 >= 0.,
         ExcNonAdmissibleState(state,
                               "Negative squared alfven speed, c_a^2 = " +
                                 dealii::Utilities::to_string(c_a2)));
  Assert(d_n >= 0.,
         ExcNonAdmissibleState(state,
                               "Negative d_n = " +
                                 dealii::Utilities::to_string(d_n)));
  const double c_s2 =
    0.5 * (a_s2 + b2 / state[density_component] - std::sqrt(d_n));
  const double c_f2 =
    0.5 * (a_s2 + b2 / state[density_component] + std::sqrt(d_n));
  Assert(c_s2 >= 0.,
         ExcNonAdmissibleState(state,
                               "Negative squared slow speed, c_s^2 = " +
                                 dealii::Utilities::to_string(c_s2)));
  Assert(c_f2 >= epsilon_d,
         ExcNonAdmissibleState(state,
                               "Expect positive squared fast speed, c_f^2 = " +
                                 dealii::Utilities::to_string(c_f2)));
  const double a_s = std::sqrt(a_s2);
  const double c_a = std::sqrt(c_a2);
  const double c_s = std::sqrt(c_s2);
  const double c_f = std::sqrt(c_f2);

  double alp_s, alp_f;
  if ((c_f2 - c_s2) < epsilon_d)
    {
      alp_s = 1.;
      alp_f = 1.;
    }
  else
    {
      alp_s = std::sqrt((c_f2 - a_s2) / (c_f2 - c_s2));
      alp_f = std::sqrt((c_f2 - c_a2) / (c_f2 - c_s2));
    }
  Assert(alp_s >= 0.,
         ExcNonAdmissibleState(state,
                               "Expect non-negative value for alpha_s = " +
                                 dealii::Utilities::to_string(alp_s)));
  Assert(alp_f >= 0.,
         ExcNonAdmissibleState(state,
                               "Expect non-negative value for alpha_f = " +
                                 dealii::Utilities::to_string(alp_f)));

  // Construct perpendicular normal vector n_perp
  dealii::Tensor<1, n_vec_components> n_perp;
  if (b_perp >= epsilon_d)
    {
      for (unsigned int d = 0; d < n_vec_components; ++d)
        n_perp[d] = (state[first_magnetic_component + d] - nb * normal_3d[d]);
    }
  else
    {
      n_perp         = 0; // tmp for e_j
      unsigned int j = 0; // j = arg_min(normal_3d)
      for (unsigned int d = 1; d < n_vec_components; ++d)
        if (std::abs(normal_3d[d]) < std::abs(normal_3d[j]))
          j = d;
      n_perp[j] = 1.;
      n_perp    = dealii::cross_product_3d<n_vec_components>(n_perp, normal_3d);
    }
  n_perp /= n_perp.norm();

  Assert(std::abs(n_perp.norm() - 1) < epsilon_d,
         dealii::ExcMessage("n_perp is not normalized."));
  Assert(std::abs(normal_3d * n_perp) < epsilon_d,
         dealii::ExcMessage("n_perp is not perpendicular to normal"));

  const double                              n_perp_u = n_perp * u;
  const dealii::Tensor<1, n_vec_components> n_perp_cross_u =
    dealii::cross_product_3d<n_vec_components>(n_perp, u);
  const double                              sgn_nb = (nb >= 0.) ? 1. : -1.;
  const dealii::Tensor<1, n_vec_components> n_perp_cross_normal =
    dealii::cross_product_3d<n_vec_components>(n_perp, normal_3d);

  const double theta_1 =
    alp_f * alp_f * a_s2 *
      (c_f2 + (2. - adiabatic_index) / (adiabatic_index - 1.) * a_s2) +
    alp_s * alp_s * c_f2 *
      (c_s2 + (2. - adiabatic_index) / (adiabatic_index - 1.) * a_s2);
  const double theta_2 =
    (sgn_nb * alp_f * alp_f * c_f * a_s + sgn_nb * alp_s * alp_s * c_s * c_a);
  /**
   * @note The left eigenvectors have only a numerical accuracy of
   * \f$ \sqrt{\epsilon} \f$,
   * where \f$ \epsilon \f$ is @ref MHDParameters::epsilon_d.
   * This is because of divisions with the parameter
   * \f$ \theta_1 \propto c_a^4 \f$
   * which scales like \f$ \theta_1 \sim \epsilon^2 \f$
   * for small \f$ c_a^2 \sim \epsilon \f$.
   * This is relevant for states with \f$ P \ll \rho \f$.
   */
  Assert(theta_1 >= epsilon_d * epsilon_d,
         ExcNonAdmissibleState(state,
                               "Expect positive value for theta_1 = " +
                                 dealii::Utilities::to_string(theta_1)));
  Assert(std::fabs(theta_2) >= epsilon_d,
         ExcNonAdmissibleState(state,
                               "Expect non-zero value for theta_2 = " +
                                 dealii::Utilities::to_string(theta_2)));
  if constexpr (divergence_cleaning)
    Assert(divergence_cleaning_speed > 0.,
           dealii::ExcMessage(
             "Speed for hyperbolic divergence cleaning must be positive."));


  // Left fast magnetosonic mode l_1
  eigenvectors[0][density_component] =
    alp_f * a_s2 * u2 / (4. * theta_1) +
    (sgn_nb * alp_f * a_s * nu - alp_s * c_s * n_perp_u) / (2. * theta_2);
  eigenvectors[0][energy_component] = alp_f * a_s2 / (2. * theta_1);
  for (unsigned int d = 0; d < n_vec_components; ++d)
    {
      eigenvectors[0][first_momentum_component + d] =
        -alp_f * a_s2 / (2. * theta_1) * u[d] -
        sgn_nb * alp_f * a_s / (2. * theta_2) * normal_3d[d] +
        alp_s * c_s / (2. * theta_2) * n_perp[d];
      eigenvectors[0][first_magnetic_component + d] =
        alp_s * c_f * std::sqrt(state[density_component]) *
        (c_s2 + (2. - adiabatic_index) / (adiabatic_index - 1.) * a_s2) /
        (2. * theta_1) * n_perp[d];
    }
  if constexpr (divergence_cleaning)
    eigenvectors[0][divergence_cleaning_component] = 0.;


  // Left alfven mode l_2
  eigenvectors[1][density_component] = sgn_nb / 2. * normal_3d * n_perp_cross_u;
  eigenvectors[1][energy_component]  = 0.;
  for (unsigned int d = 0; d < n_vec_components; ++d)
    {
      eigenvectors[1][first_momentum_component + d] =
        sgn_nb / 2. * n_perp_cross_normal[d];
      eigenvectors[1][first_magnetic_component + d] =
        std::sqrt(state[density_component]) / 2. * n_perp_cross_normal[d];
    }
  if constexpr (divergence_cleaning)
    eigenvectors[1][divergence_cleaning_component] = 0.;


  // Left slow magnetosonic mode l_3
  eigenvectors[2][density_component] =
    alp_s * c_f2 * u2 / (4. * theta_1) +
    (sgn_nb * alp_s * c_a * nu + alp_f * c_f * n_perp_u) / (2. * theta_2);
  eigenvectors[2][energy_component] = alp_s * c_f2 / (2. * theta_1);
  for (unsigned int d = 0; d < n_vec_components; ++d)
    {
      eigenvectors[2][first_momentum_component + d] =
        -alp_s * c_f2 / (2. * theta_1) * u[d] -
        sgn_nb * alp_s * c_a / (2. * theta_2) * normal_3d[d] -
        alp_f * c_f / (2. * theta_2) * n_perp[d];
      eigenvectors[2][first_magnetic_component + d] =
        -alp_f * c_f * std::sqrt(state[density_component]) *
        (c_f2 + (2. - adiabatic_index) / (adiabatic_index - 1.) * a_s2) /
        (2. * theta_1) * n_perp[d];
    }
  if constexpr (divergence_cleaning)
    eigenvectors[2][divergence_cleaning_component] = 0.;


  // Density entropy mode l_4
  eigenvectors[3][density_component] =
    1. - (alp_f * alp_f * a_s2 + alp_s * alp_s * c_f2) * u2 / (2. * theta_1);
  eigenvectors[3][energy_component] =
    -(alp_f * alp_f * a_s2 + alp_s * alp_s * c_f2) / theta_1;
  for (unsigned int d = 0; d < n_vec_components; ++d)
    {
      eigenvectors[3][first_momentum_component + d] =
        (alp_f * alp_f * a_s2 + alp_s * alp_s * c_f2) / theta_1 * u[d];
      eigenvectors[3][first_magnetic_component + d] =
        alp_f * alp_s * c_f * std::sqrt(state[density_component]) *
        (c_f2 - c_s2) / theta_1 * n_perp[d];
    }
  if constexpr (divergence_cleaning)
    eigenvectors[3][divergence_cleaning_component] = 0.;


  // Magnetic divergence mode l_5
  eigenvectors[4][density_component] = 0;
  eigenvectors[4][energy_component]  = 0;
  for (unsigned int d = 0; d < n_vec_components; ++d)
    {
      eigenvectors[4][first_momentum_component + d] = 0;
      if constexpr (divergence_cleaning)
        eigenvectors[4][first_magnetic_component + d] = 0.5 * normal_3d[d];
      else
        eigenvectors[4][first_magnetic_component + d] = normal_3d[d];
    }
  if constexpr (divergence_cleaning)
    eigenvectors[4][divergence_cleaning_component] =
      -0.5 / divergence_cleaning_speed;


  // Right slow magnetosonic mode l_6
  eigenvectors[5][density_component] =
    alp_s * c_f2 * u2 / (4. * theta_1) -
    (sgn_nb * alp_s * c_a * nu + alp_f * c_f * n_perp_u) / (2. * theta_2);
  eigenvectors[5][energy_component] = alp_s * c_f2 / (2. * theta_1);
  for (unsigned int d = 0; d < n_vec_components; ++d)
    {
      eigenvectors[5][first_momentum_component + d] =
        -alp_s * c_f2 / (2. * theta_1) * u[d] +
        sgn_nb * alp_s * c_a / (2. * theta_2) * normal_3d[d] +
        alp_f * c_f / (2. * theta_2) * n_perp[d];
      eigenvectors[5][first_magnetic_component + d] =
        -alp_f * c_f * std::sqrt(state[density_component]) *
        (c_f2 + (2. - adiabatic_index) / (adiabatic_index - 1.) * a_s2) /
        (2. * theta_1) * n_perp[d];
    }
  if constexpr (divergence_cleaning)
    eigenvectors[5][divergence_cleaning_component] = 0.;


  // Right alfven mode l_7
  eigenvectors[6][density_component] =
    -sgn_nb / 2. * normal_3d * n_perp_cross_u;
  eigenvectors[6][energy_component] = 0.;
  for (unsigned int d = 0; d < n_vec_components; ++d)
    {
      eigenvectors[6][first_momentum_component + d] =
        -sgn_nb / 2. * n_perp_cross_normal[d];
      eigenvectors[6][first_magnetic_component + d] =
        std::sqrt(state[density_component]) / 2. * n_perp_cross_normal[d];
    }
  if constexpr (divergence_cleaning)
    eigenvectors[6][divergence_cleaning_component] = 0.;


  // Right fast magnetosonic mode l_8
  eigenvectors[7][density_component] =
    alp_f * a_s2 * u2 / (4. * theta_1) -
    (sgn_nb * alp_f * a_s * nu - alp_s * c_s * n_perp_u) / (2. * theta_2);
  eigenvectors[7][energy_component] = alp_f * a_s2 / (2. * theta_1);
  for (unsigned int d = 0; d < n_vec_components; ++d)
    {
      eigenvectors[7][first_momentum_component + d] =
        -alp_f * a_s2 / (2. * theta_1) * u[d] +
        sgn_nb * alp_f * a_s / (2. * theta_2) * normal_3d[d] -
        alp_s * c_s / (2. * theta_2) * n_perp[d];
      eigenvectors[7][first_magnetic_component + d] =
        alp_s * c_f * std::sqrt(state[density_component]) *
        (c_s2 + (2. - adiabatic_index) / (adiabatic_index - 1.) * a_s2) /
        (2. * theta_1) * n_perp[d];
    }
  if constexpr (divergence_cleaning)
    eigenvectors[7][divergence_cleaning_component] = 0.;


  if constexpr (divergence_cleaning)
    {
      // Magnetic divergence mode l_9
      /**
       * @note The left eigenvector \f$ \mathbf{l}_{9} \f$
       * only exists with hyperbolic divergence cleaning.
       */
      eigenvectors[8][density_component] = 0.;
      eigenvectors[8][energy_component]  = 0.;
      for (unsigned int d = 0; d < n_vec_components; ++d)
        {
          eigenvectors[8][first_momentum_component + d] = 0.;
          eigenvectors[8][first_magnetic_component + d] = 0.5 * normal_3d[d];
        }
      eigenvectors[8][divergence_cleaning_component] =
        0.5 / divergence_cleaning_speed;
    }
}



template <unsigned int dim, bool divergence_cleaning>
void
sapphirepp::MHD::MHDEquations<dim, divergence_cleaning>::
  compute_transformation_matrices(
    const state_type                            &state,
    std::array<dealii::FullMatrix<double>, dim> &left_matrices,
    std::array<dealii::FullMatrix<double>, dim> &right_matrices) const
{
  // TODO
  AssertDimension(state.size(), n_components);

  dealii::Tensor<1, dim> direction;

  for (unsigned int d = 0; d < dim; ++d)
    {
      AssertDimension(left_matrices[d].n(), n_components);
      AssertDimension(left_matrices[d].m(), n_components);
      AssertDimension(right_matrices[d].n(), n_components);
      AssertDimension(right_matrices[d].m(), n_components);

      direction    = 0.;
      direction[d] = 1.;

      compute_left_eigenvector_matrix(state, direction, left_matrices[d]);
      compute_right_eigenvector_matrix(state, direction, right_matrices[d]);
    }
}



template <unsigned int dim, bool divergence_cleaning>
void
sapphirepp::MHD::MHDEquations<dim, divergence_cleaning>::
  convert_primitive_to_conserved(const state_type &primitive_state,
                                 state_type       &conserved_state) const
{
  AssertDimension(primitive_state.size(), n_components);
  AssertDimension(conserved_state.size(), n_components);
  Assert(primitive_state[density_component] >= epsilon_d,
         ExcNonAdmissibleState(primitive_state,
                               "Non-admissible primitive state."));
  Assert(primitive_state[pressure_component] >= epsilon_d,
         ExcNonAdmissibleState(primitive_state,
                               "Non-admissible primitive state."));

  double u2 = 0.;
  double b2 = 0.;
  for (unsigned int d = 0; d < n_vec_components; ++d)
    {
      u2 += primitive_state[first_velocity_component + d] *
            primitive_state[first_velocity_component + d];
      b2 += primitive_state[first_magnetic_component + d] *
            primitive_state[first_magnetic_component + d];
    }

  const double energy =
    0.5 * primitive_state[density_component] * u2 +
    primitive_state[pressure_component] / (adiabatic_index - 1.) + 0.5 * b2;
  Assert(energy >= epsilon_d,
         ExcNonAdmissibleState(primitive_state,
                               "Non-admissible primitive state."));

  conserved_state[density_component] = primitive_state[density_component];
  conserved_state[energy_component]  = energy;
  for (unsigned int d = 0; d < n_vec_components; ++d)
    {
      conserved_state[first_momentum_component + d] =
        primitive_state[density_component] *
        primitive_state[first_velocity_component + d];
      conserved_state[first_magnetic_component + d] =
        primitive_state[first_magnetic_component + d];
    }
  if constexpr (divergence_cleaning)
    conserved_state[divergence_cleaning_component] =
      primitive_state[divergence_cleaning_component];
}



template <unsigned int dim, bool divergence_cleaning>
void
sapphirepp::MHD::MHDEquations<dim, divergence_cleaning>::
  convert_conserved_to_primitive(const state_type &conserved_state,
                                 state_type       &primitive_state) const
{
  AssertDimension(conserved_state.size(), n_components);
  AssertDimension(primitive_state.size(), n_components);

  primitive_state[density_component]  = conserved_state[density_component];
  primitive_state[pressure_component] = compute_pressure(conserved_state);
  for (unsigned int d = 0; d < n_vec_components; ++d)
    {
      primitive_state[first_velocity_component + d] =
        conserved_state[first_momentum_component + d] /
        conserved_state[density_component];
      primitive_state[first_magnetic_component + d] =
        conserved_state[first_magnetic_component + d];
    }
  if constexpr (divergence_cleaning)
    primitive_state[divergence_cleaning_component] =
      conserved_state[divergence_cleaning_component];
}



template <unsigned int dim, bool divergence_cleaning>
void
sapphirepp::MHD::MHDEquations<dim, divergence_cleaning>::
  convert_characteristic_to_conserved(const state_type &characteristic_state,
                                      const dealii::Tensor<1, dim> &normal,
                                      state_type &conserved_state) const
{
  dealii::FullMatrix<double> right_matrix(n_components);
  compute_right_eigenvector_matrix(conserved_state, normal, right_matrix);

  right_matrix.vmult(conserved_state, characteristic_state, false);
}



template <unsigned int dim, bool divergence_cleaning>
void
sapphirepp::MHD::MHDEquations<dim, divergence_cleaning>::
  convert_conserved_to_characteristic(const state_type &conserved_state,
                                      const dealii::Tensor<1, dim> &normal,
                                      state_type &characteristic_state) const
{
  dealii::FullMatrix<double> left_matrix(n_components);
  compute_left_eigenvector_matrix(conserved_state, normal, left_matrix);

  left_matrix.vmult(characteristic_state, conserved_state, false);
}



template <unsigned int dim, bool divergence_cleaning>
void
sapphirepp::MHD::MHDEquations<dim, divergence_cleaning>::
  convert_gradient_characteristic_to_conserved(
    const flux_type                             &characteristic_gradient,
    std::array<dealii::FullMatrix<double>, dim> &right_matrices,
    flux_type                                   &conserved_gradient) const
{
  for (unsigned int d = 0; d < dim; ++d)
    {
      AssertDimension(right_matrices[d].n(), n_components);
      AssertDimension(right_matrices[d].m(), n_components);

      for (unsigned int c1 = 0; c1 < n_components; ++c1)
        {
          conserved_gradient[c1][d] = 0.;
          for (unsigned int c2 = 0; c2 < n_components; ++c2)
            {
              conserved_gradient[c1][d] +=
                right_matrices[d][c1][c2] * characteristic_gradient[c2][d];
            }
        }
    }
}



template <unsigned int dim, bool divergence_cleaning>
void
sapphirepp::MHD::MHDEquations<dim, divergence_cleaning>::
  convert_gradient_conserved_to_characteristic(
    const flux_type                             &conserved_gradient,
    std::array<dealii::FullMatrix<double>, dim> &left_matrices,
    flux_type                                   &characteristic_gradient) const
{
  for (unsigned int d = 0; d < dim; ++d)
    {
      AssertDimension(left_matrices[d].n(), n_components);
      AssertDimension(left_matrices[d].m(), n_components);

      for (unsigned int c1 = 0; c1 < n_components; ++c1)
        {
          characteristic_gradient[c1][d] = 0;
          for (unsigned int c2 = 0; c2 < n_components; ++c2)
            {
              characteristic_gradient[c1][d] +=
                left_matrices[d][c1][c2] * conserved_gradient[c2][d];
            }
        }
    }
}



template <unsigned int dim, bool divergence_cleaning>
void
sapphirepp::MHD::MHDEquations<dim, divergence_cleaning>::
  compute_hyperbolic_divergence_cleaning_speed(const double       dt_cfl,
                                               const double       dx_min,
                                               const unsigned int fe_degree)
{
  if constexpr (!divergence_cleaning)
    return;

  dealii::LogStream::Prefix p("MHDEquations", saplog);
  Assert(dt_cfl > 0, dealii::ExcMessage("CFL time step must be positive."));
  Assert(dx_min > 0, dealii::ExcMessage("Cell size must be positive."));


  // c_h = C_h dx / ((2k+1) dt)
  divergence_cleaning_speed =
    divergence_cleaning_Ch * dx_min / ((2 * fe_degree + 1) * dt_cfl);

  // c_p = sqrt(C_r * c_h)
  // tau^{-1} = c_h^2 / c_p^2 = c_h / C_r
  divergence_cleaning_damping =
    divergence_cleaning_speed / divergence_cleaning_Cr;

  saplog << "divergence_cleaning_speed=" << divergence_cleaning_speed
         << ", divergence_cleaning_damping=" << divergence_cleaning_damping
         << std::endl;

  Assert(divergence_cleaning_speed > 0.,
         dealii::ExcMessage(
           "Speed for hyperbolic divergence cleaning must be positive, c_h=" +
           dealii::Utilities::to_string(divergence_cleaning_speed)));
}



// Explicit instantiations
template class sapphirepp::MHD::MHDEquations<1, false>;
template class sapphirepp::MHD::MHDEquations<1, true>;
template class sapphirepp::MHD::MHDEquations<2, false>;
template class sapphirepp::MHD::MHDEquations<2, true>;
template class sapphirepp::MHD::MHDEquations<3, false>;
template class sapphirepp::MHD::MHDEquations<3, true>;
