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

#include <algorithm>
#include <cmath>

#include "sapphirepp-logstream.h"



/** Precision for safe double comparision */
const double epsilon_d = 1e-6;



sapphirepp::MHD::MHDEquations::MHDEquations(const double adiabatic_index)
  : adiabatic_index{adiabatic_index}
{
  AssertThrow(adiabatic_index > 1.0,
              dealii::ExcMessage("Adiabatic index must be larger than 1.0."));
}



std::vector<std::string>
sapphirepp::MHD::MHDEquations::create_component_name_list(
  const std::string &prefix)
{
  std::vector<std::string> component_names(n_components);

  component_names[density_component] = prefix + "rho";
  component_names[energy_component]  = prefix + "E";
  for (unsigned int d = 0; d < spacedim; ++d)
    {
      component_names[first_momentum_component + d] = prefix + "p";
      component_names[first_magnetic_component + d] = prefix + "b";
    }

  return component_names;
}



std::vector<dealii::DataComponentInterpretation::DataComponentInterpretation>
sapphirepp::MHD::MHDEquations::create_component_interpretation_list()
{
  std::vector<dealii::DataComponentInterpretation::DataComponentInterpretation>
    data_component_interpretation(n_components);

  data_component_interpretation[density_component] =
    dealii::DataComponentInterpretation::component_is_scalar;
  data_component_interpretation[energy_component] =
    dealii::DataComponentInterpretation::component_is_scalar;
  for (unsigned int d = 0; d < spacedim; ++d)
    {
      data_component_interpretation[first_momentum_component + d] =
        dealii::DataComponentInterpretation::component_is_part_of_vector;
      data_component_interpretation[first_magnetic_component + d] =
        dealii::DataComponentInterpretation::component_is_part_of_vector;
    }

  return data_component_interpretation;
}



void
sapphirepp::MHD::MHDEquations::compute_flux_matrix(const state_type &state,
                                                   flux_type &flux_matrix) const
{
  AssertDimension(state.size(), n_components);

  const double pressure = compute_pressure(state);
  double       b2       = 0.;
  double       pb       = 0.;
  for (unsigned int d = 0; d < spacedim; ++d)
    {
      b2 += state[first_magnetic_component + d] *
            state[first_magnetic_component + d];
      pb += state[first_momentum_component + d] *
            state[first_magnetic_component + d];
    }

  for (unsigned int j = 0; j < spacedim; ++j)
    {
      flux_matrix[density_component][j] = state[first_momentum_component + j];

      flux_matrix[energy_component][j] =
        state[first_momentum_component + j] / state[density_component] *
          (state[energy_component] + pressure + 0.5 * b2) -
        1. / (state[density_component]) * pb *
          state[first_magnetic_component + j];

      for (unsigned int i = 0; i < spacedim; ++i)
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
    }
}



double
sapphirepp::MHD::MHDEquations::compute_maximum_normal_eigenvalue(
  const state_type                  &state,
  const dealii::Tensor<1, spacedim> &normal) const
{
  AssertDimension(state.size(), n_components);
  Assert(std::abs(normal.norm() - 1) < epsilon_d,
         dealii::ExcMessage("Normal vector must be normalized."));

  const double pressure = compute_pressure(state);
  double       b2       = 0.;
  double       nb       = 0.;
  double       nu       = 0.;
  for (unsigned int d = 0; d < spacedim; ++d)
    {
      b2 += state[first_magnetic_component + d] *
            state[first_magnetic_component + d];
      nu += normal[d] * state[first_momentum_component + d] /
            state[density_component];
      nb += normal[d] * state[first_magnetic_component + d];
    }

  const double a_s2 = adiabatic_index * pressure / state[density_component];
  const double c_a2 = nb * nb / state[density_component];
  const double d_n  = (a_s2 + b2 / state[density_component]) *
                       (a_s2 + b2 / state[density_component]) -
                     4. * a_s2 * c_a2;
  Assert(a_s2 >= 0.,
         dealii::ExcMessage("Negative squared adiabatic sound speed warning."));
  Assert(c_a2 >= 0.,
         dealii::ExcMessage("Negative squared alfven speed warning."));
  Assert(d_n >= 0., dealii::ExcMessage("Negative squared value warning."));
  const double c_f2 =
    0.5 * (a_s2 + b2 / state[density_component] + std::sqrt(d_n));
  Assert(c_f2 >= 0.,
         dealii::ExcMessage("Negative squared fast speed warning."));

  return std::fabs(nu) + std::sqrt(c_f2);
}



double
sapphirepp::MHD::MHDEquations::compute_pressure(const state_type &state) const
{
  AssertDimension(state.size(), n_components);
  Assert(state[density_component] > 0.,
         dealii::ExcMessage("Negative density warning."));
  Assert(state[energy_component] > 0.,
         dealii::ExcMessage("Negative energy warning."));

  double p2 = 0.;
  double b2 = 0.;
  for (unsigned int d = 0; d < spacedim; ++d)
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

  Assert(pressure > 0., dealii::ExcMessage("Negative pressure warning."));
  return pressure;
}



void
sapphirepp::MHD::MHDEquations::compute_normale_eigenvalues(
  const state_type                  &state,
  const dealii::Tensor<1, spacedim> &normal,
  dealii::Vector<double>            &eigenvalues) const
{
  AssertDimension(state.size(), n_components);
  AssertDimension(eigenvalues.size(), n_components);
  Assert(std::abs(normal.norm() - 1) < epsilon_d,
         dealii::ExcMessage("Normal vector must be normalized."));

  const double pressure = compute_pressure(state);
  double       b2       = 0.;
  double       nu       = 0.;
  double       nb       = 0.;
  for (unsigned int d = 0; d < spacedim; ++d)
    {
      b2 += state[first_magnetic_component + d] *
            state[first_magnetic_component + d];
      nu += normal[d] * state[first_momentum_component + d] /
            state[density_component];
      nb += normal[d] * state[first_magnetic_component + d];
    }

  const double a_s2 = adiabatic_index * pressure / state[density_component];
  const double c_a2 = nb * nb / state[density_component];
  const double d_n  = (a_s2 + b2 / state[density_component]) *
                       (a_s2 + b2 / state[density_component]) -
                     4. * a_s2 * c_a2;
  Assert(a_s2 >= 0.,
         dealii::ExcMessage("Negative squared adiabatic sound speed warning."));
  Assert(c_a2 >= 0.,
         dealii::ExcMessage("Negative squared alfven speed warning."));
  Assert(d_n >= 0., dealii::ExcMessage("Negative squared value warning."));
  const double c_s2 =
    0.5 * (a_s2 + b2 / state[density_component] - std::sqrt(d_n));
  const double c_f2 =
    0.5 * (a_s2 + b2 / state[density_component] + std::sqrt(d_n));
  Assert(c_s2 >= 0.,
         dealii::ExcMessage("Negative squared slow speed warning."));
  Assert(c_f2 >= 0.,
         dealii::ExcMessage("Negative squared fast speed warning."));

  eigenvalues[0] = nu - std::sqrt(c_f2);
  eigenvalues[1] = nu - std::sqrt(c_a2);
  eigenvalues[2] = nu - std::sqrt(c_s2);
  eigenvalues[3] = nu;
  eigenvalues[4] = nu;
  eigenvalues[5] = nu + std::sqrt(c_s2);
  eigenvalues[6] = nu + std::sqrt(c_a2);
  eigenvalues[7] = nu + std::sqrt(c_f2);
}



void
sapphirepp::MHD::MHDEquations::compute_right_eigenvector_matrix(
  const state_type                  &state,
  const dealii::Tensor<1, spacedim> &normal,
  dealii::FullMatrix<double>        &eigenvectors) const

{
  dealii::LogStream::Prefix p("MHDEquations", saplog);
  saplog << "Calculate right eigenvector matrix:"
         << "\n state = " << state << "\n normal = " << normal << std::endl;
  AssertDimension(state.size(), n_components);
  AssertDimension(eigenvectors.n(), n_components);
  AssertDimension(eigenvectors.m(), n_components);
  Assert(std::abs(normal.norm() - 1) < epsilon_d,
         dealii::ExcMessage("Normal vector must be normalized."));

  const double                pressure = compute_pressure(state);
  dealii::Tensor<1, spacedim> u;
  double                      b2 = 0.;
  double                      u2 = 0.;
  double                      nu = 0.;
  double                      nb = 0.;
  for (unsigned int d = 0; d < spacedim; ++d)
    {
      u[d] = state[first_momentum_component + d] / state[density_component];
      b2 += state[first_magnetic_component + d] *
            state[first_magnetic_component + d];
      u2 += state[first_momentum_component + d] / state[density_component] *
            state[first_momentum_component + d] / state[density_component];
      nu += normal[d] * state[first_momentum_component + d] /
            state[density_component];
      nb += normal[d] * state[first_magnetic_component + d];
    }
  const double b_perp = std::sqrt(b2 - nb * nb);
  saplog << "pressure = " << pressure << ", u = " << u << ", b2 = " << b2
         << ", u2 = " << u2 << ", nu = " << nu << ", nb = " << nb
         << ", b_perp = " << b_perp << std::endl;

  const double a_s2 = adiabatic_index * pressure / state[density_component];
  const double c_a2 = nb * nb / state[density_component];
  const double d_n  = (a_s2 + b2 / state[density_component]) *
                       (a_s2 + b2 / state[density_component]) -
                     4. * a_s2 * c_a2;
  saplog << "a_s2 = " << a_s2 << ", c_a2 = " << c_a2 << ", d_n = " << d_n
         << std::endl;
  Assert(a_s2 >= 0.,
         dealii::ExcMessage("Negative squared adiabatic sound speed warning."));
  Assert(c_a2 >= 0.,
         dealii::ExcMessage("Negative squared alfven speed warning."));
  Assert(d_n >= 0., dealii::ExcMessage("Negative squared value warning."));
  const double c_s2 =
    0.5 * (a_s2 + b2 / state[density_component] - std::sqrt(d_n));
  const double c_f2 =
    0.5 * (a_s2 + b2 / state[density_component] + std::sqrt(d_n));
  saplog << "c_s2 = " << c_s2 << ", c_f2 = " << c_f2 << std::endl;
  Assert(c_s2 >= 0.,
         dealii::ExcMessage("Negative squared slow speed warning."));
  Assert(c_f2 >= 0.,
         dealii::ExcMessage("Negative squared fast speed warning."));
  const double a_s = std::sqrt(a_s2);
  const double c_a = std::sqrt(c_a2);
  const double c_s = std::sqrt(c_s2);
  const double c_f = std::sqrt(c_f2);
  saplog << "a_s = " << a_s << ", c_a = " << c_a << ", c_s = " << c_s
         << ", c_f = " << c_f << std::endl;

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
  saplog << "alp_s = " << alp_s << ", alp_f = " << alp_f << std::endl;
  Assert(alp_s >= 0., dealii::ExcMessage("Expect non-negative value."));
  Assert(alp_f >= 0., dealii::ExcMessage("Expect non-negative value."));

  // Construct perpendicular normal vector n_perp
  dealii::Tensor<1, spacedim> n_perp;
  if (b_perp > 0)
    {
      for (unsigned int d = 0; d < spacedim; ++d)
        n_perp[d] = (state[first_magnetic_component + d] - nb * normal[d]);
    }
  else
    {
      n_perp         = 0; // tmp for e_j
      unsigned int j = 0; // j = arg_min(normal)
      for (unsigned int d = 1; d < spacedim; ++d)
        if (std::abs(normal[d]) < std::abs(normal[j]))
          j = d;
      n_perp[j] = 1.;
      saplog << "e_j = " << n_perp << ", ";
      n_perp = dealii::cross_product_3d<spacedim>(n_perp, normal);
    }
  n_perp /= n_perp.norm();
  saplog << "normal = " << normal << ", n_perp = " << n_perp
         << ", |n_perp| = " << n_perp.norm()
         << ", n*n_perp = " << normal * n_perp << std::endl;

  Assert(std::abs(n_perp.norm() - 1) < epsilon_d,
         dealii::ExcMessage("n_perp is not normalized."));
  Assert(std::abs(normal * n_perp) < epsilon_d,
         dealii::ExcMessage("n_perp is not perpendicular to normal"));

  const double                      n_perp_u = n_perp * u;
  const dealii::Tensor<1, spacedim> n_perp_cross_u =
    dealii::cross_product_3d<spacedim>(n_perp, u);
  const double                      sgn_nb = (nb >= 0.) ? 1. : -1.;
  const dealii::Tensor<1, spacedim> n_perp_cross_normal =
    dealii::cross_product_3d<spacedim>(n_perp, normal);
  saplog << "n_perp_u = " << n_perp_u << ", n_perp_cross_u = " << n_perp_cross_u
         << ", sgn_nb = " << sgn_nb
         << ", n_perp_cross_normal = " << n_perp_cross_normal << std::endl;


  // Left fast magnetosonic mode r_1
  eigenvectors[density_component][0] = alp_f;
  eigenvectors[energy_component][0] =
    alp_f * u2 / 2. + alp_f * c_f2 / (adiabatic_index - 1.) +
    sgn_nb * alp_s * c_a * n_perp_u - alp_f * c_f * nu -
    (2. - adiabatic_index) / (adiabatic_index - 1.) * alp_f * (c_f2 - a_s2);
  for (unsigned int d = 0; d < spacedim; ++d)
    {
      eigenvectors[first_momentum_component + d][0] =
        alp_f * (u[d] - c_f * normal[d]) + sgn_nb * alp_s * c_a * n_perp[d];
      eigenvectors[first_magnetic_component + d][0] =
        alp_s * c_f / std::sqrt(state[density_component]) * n_perp[d];
    }


  // Left alfven mode r_2
  eigenvectors[density_component][1] = 0;
  eigenvectors[energy_component][1]  = -sgn_nb * normal * n_perp_cross_u;
  for (unsigned int d = 0; d < spacedim; ++d)
    {
      eigenvectors[first_momentum_component + d][1] =
        sgn_nb * n_perp_cross_normal[d];
      eigenvectors[first_magnetic_component + d][1] =
        n_perp_cross_normal[d] / state[density_component];
    }


  // Left slow magnetosonic mode r_3
  eigenvectors[density_component][2] = alp_s;
  eigenvectors[energy_component][2] =
    alp_s * u2 / 2. + alp_s * c_s2 / (adiabatic_index - 1.) -
    sgn_nb * alp_f * a_s * n_perp_u - alp_s * c_s * nu -
    (2. - adiabatic_index) / (adiabatic_index - 1.) * alp_s * (c_s2 - a_s2);
  for (unsigned int d = 0; d < spacedim; ++d)
    {
      eigenvectors[first_momentum_component + d][2] =
        alp_s * (u[d] - c_s * normal[d]) - sgn_nb * alp_f * a_s * n_perp[d];
      eigenvectors[first_magnetic_component + d][2] =
        -alp_f * a_s2 / (std::sqrt(state[density_component]) * c_f) * n_perp[d];
    }


  // Density entropy mode r_4
  eigenvectors[density_component][3] = 1.;
  eigenvectors[energy_component][3]  = u2 / 2.;
  for (unsigned int d = 0; d < spacedim; ++d)
    {
      eigenvectors[first_momentum_component + d][3] = u[d];
      eigenvectors[first_magnetic_component + d][3] = 0;
    }


  // Magnetic entropy mode r_5
  eigenvectors[density_component][4] = 0.;
  eigenvectors[energy_component][4]  = 0.;
  for (unsigned int d = 0; d < spacedim; ++d)
    {
      eigenvectors[first_momentum_component + d][4] = 0;
      eigenvectors[first_magnetic_component + d][4] = normal[d];
    }


  // Right slow magnetosonic mode r_6
  eigenvectors[density_component][5] = alp_s;
  eigenvectors[energy_component][5] =
    alp_s * u2 / 2. + alp_s * c_s2 / (adiabatic_index - 1.) +
    sgn_nb * alp_f * a_s * n_perp_u + alp_s * c_s * nu -
    (2. - adiabatic_index) / (adiabatic_index - 1.) * alp_s * (c_s2 - a_s2);
  for (unsigned int d = 0; d < spacedim; ++d)
    {
      eigenvectors[first_momentum_component + d][5] =
        alp_s * (u[d] + c_s * normal[d]) + sgn_nb * alp_f * a_s * n_perp[d];
      eigenvectors[first_magnetic_component + d][5] =
        -alp_f * a_s2 / (std::sqrt(state[density_component]) * c_f) * n_perp[d];
    }


  // Right alfven mode r_7
  eigenvectors[density_component][6] = 0;
  eigenvectors[energy_component][6]  = sgn_nb * normal * n_perp_cross_u;
  for (unsigned int d = 0; d < spacedim; ++d)
    {
      eigenvectors[first_momentum_component + d][6] =
        -sgn_nb * n_perp_cross_normal[d];
      eigenvectors[first_magnetic_component + d][6] =
        n_perp_cross_normal[d] / state[density_component];
    }


  // Right fast magnetosonic mode r_8
  eigenvectors[density_component][7] = alp_f;
  eigenvectors[energy_component][7] =
    alp_f * u2 / 2. + alp_f * c_f2 / (adiabatic_index - 1.) -
    sgn_nb * alp_s * c_a * n_perp_u + alp_f * c_f * nu -
    (2. - adiabatic_index) / (adiabatic_index - 1.) * alp_f * (c_f2 - a_s2);
  for (unsigned int d = 0; d < spacedim; ++d)
    {
      eigenvectors[first_momentum_component + d][7] =
        alp_f * (u[d] + c_f * normal[d]) - sgn_nb * alp_s * c_a * n_perp[d];
      eigenvectors[first_magnetic_component + d][7] =
        alp_s * c_f / std::sqrt(state[density_component]) * n_perp[d];
    }
}



void
sapphirepp::MHD::MHDEquations::compute_left_eigenvector_matrix(
  const state_type                  &state,
  const dealii::Tensor<1, spacedim> &normal,
  dealii::FullMatrix<double>        &eigenvectors) const

{
  dealii::LogStream::Prefix p("MHDEquations", saplog);
  saplog << "Calculate left eigenvector matrix:"
         << "\n state = " << state << "\n normal = " << normal << std::endl;
  AssertDimension(state.size(), n_components);
  AssertDimension(eigenvectors.n(), n_components);
  AssertDimension(eigenvectors.m(), n_components);
  Assert(std::abs(normal.norm() - 1) < epsilon_d,
         dealii::ExcMessage("Normal vector must be normalized."));

  const double                pressure = compute_pressure(state);
  dealii::Tensor<1, spacedim> u;
  double                      b2 = 0.;
  double                      u2 = 0.;
  double                      nu = 0.;
  double                      nb = 0.;
  for (unsigned int d = 0; d < spacedim; ++d)
    {
      u[d] = state[first_momentum_component + d] / state[density_component];
      b2 += state[first_magnetic_component + d] *
            state[first_magnetic_component + d];
      u2 += state[first_momentum_component + d] / state[density_component] *
            state[first_momentum_component + d] / state[density_component];
      nu += normal[d] * state[first_momentum_component + d] /
            state[density_component];
      nb += normal[d] * state[first_magnetic_component + d];
    }
  const double b_perp = std::sqrt(b2 - nb * nb);
  saplog << "pressure = " << pressure << ", u = " << u << ", b2 = " << b2
         << ", u2 = " << u2 << ", nu = " << nu << ", nb = " << nb
         << ", b_perp = " << b_perp << std::endl;

  const double a_s2 = adiabatic_index * pressure / state[density_component];
  const double c_a2 = nb * nb / state[density_component];
  const double d_n  = (a_s2 + b2 / state[density_component]) *
                       (a_s2 + b2 / state[density_component]) -
                     4. * a_s2 * c_a2;
  saplog << "a_s2 = " << a_s2 << ", c_a2 = " << c_a2 << ", d_n = " << d_n
         << std::endl;
  Assert(a_s2 >= 0.,
         dealii::ExcMessage("Negative squared adiabatic sound speed warning."));
  Assert(c_a2 >= 0.,
         dealii::ExcMessage("Negative squared alfven speed warning."));
  Assert(d_n >= 0., dealii::ExcMessage("Negative squared value warning."));
  const double c_s2 =
    0.5 * (a_s2 + b2 / state[density_component] - std::sqrt(d_n));
  const double c_f2 =
    0.5 * (a_s2 + b2 / state[density_component] + std::sqrt(d_n));
  saplog << "c_s2 = " << c_s2 << ", c_f2 = " << c_f2 << std::endl;
  Assert(c_s2 >= 0.,
         dealii::ExcMessage("Negative squared slow speed warning."));
  Assert(c_f2 >= 0.,
         dealii::ExcMessage("Negative squared fast speed warning."));
  const double a_s = std::sqrt(a_s2);
  const double c_a = std::sqrt(c_a2);
  const double c_s = std::sqrt(c_s2);
  const double c_f = std::sqrt(c_f2);
  saplog << "a_s = " << a_s << ", c_a = " << c_a << ", c_s = " << c_s
         << ", c_f = " << c_f << std::endl;

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
  saplog << "alp_s = " << alp_s << ", alp_f = " << alp_f << std::endl;
  Assert(alp_s >= 0., dealii::ExcMessage("Expect non-negative value."));
  Assert(alp_f >= 0., dealii::ExcMessage("Expect non-negative value."));

  // Construct perpendicular normal vector n_perp
  dealii::Tensor<1, spacedim> n_perp;
  if (b_perp > 0)
    {
      for (unsigned int d = 0; d < spacedim; ++d)
        n_perp[d] = (state[first_magnetic_component + d] - nb * normal[d]);
    }
  else
    {
      n_perp         = 0; // tmp for e_j
      unsigned int j = 0; // j = arg_min(normal)
      for (unsigned int d = 1; d < spacedim; ++d)
        if (std::abs(normal[d]) < std::abs(normal[j]))
          j = d;
      n_perp[j] = 1.;
      saplog << "e_j = " << n_perp << ", ";
      n_perp = dealii::cross_product_3d<spacedim>(n_perp, normal);
    }
  n_perp /= n_perp.norm();
  saplog << "normal = " << normal << ", n_perp = " << n_perp
         << ", |n_perp| = " << n_perp.norm()
         << ", n*n_perp = " << normal * n_perp << std::endl;

  Assert(std::abs(n_perp.norm() - 1) < epsilon_d,
         dealii::ExcMessage("n_perp is not normalized."));
  Assert(std::abs(normal * n_perp) < epsilon_d,
         dealii::ExcMessage("n_perp is not perpendicular to normal"));

  const double                      n_perp_u = n_perp * u;
  const dealii::Tensor<1, spacedim> n_perp_cross_u =
    dealii::cross_product_3d<spacedim>(n_perp, u);
  const double                      sgn_nb = (nb >= 0.) ? 1. : -1.;
  const dealii::Tensor<1, spacedim> n_perp_cross_normal =
    dealii::cross_product_3d<spacedim>(n_perp, normal);
  saplog << "n_perp_u = " << n_perp_u << ", n_perp_cross_u = " << n_perp_cross_u
         << ", sgn_nb = " << sgn_nb
         << ", n_perp_cross_normal = " << n_perp_cross_normal << std::endl;

  const double theta_1 =
    alp_f * alp_f * a_s2 *
      (c_f2 + (2. - adiabatic_index) / (adiabatic_index - 1.) * a_s2) +
    alp_s * alp_s * c_f2 *
      (c_s2 + (2. - adiabatic_index) / (adiabatic_index - 1.) * a_s2);
  const double theta_2 =
    (sgn_nb * alp_f * alp_f * c_f * a_s + sgn_nb * alp_s * alp_s * c_s * c_a);
  saplog << "theta_1 = " << theta_1 << ", theta_2 = " << theta_2 << std::endl;
  Assert(theta_1 > 0., dealii::ExcMessage("Expect positive value."));
  Assert(theta_2 > 0., dealii::ExcMessage("Expect positive value."));


  // Left fast magnetosonic mode l_1
  eigenvectors[0][density_component] =
    alp_f * a_s2 * u2 / (4. * theta_1) +
    (sgn_nb * alp_f * a_s * nu - alp_s * c_s * n_perp_u) / (2. * theta_2);
  eigenvectors[0][energy_component] = alp_s * a_s2 / (2. * theta_1);
  for (unsigned int d = 0; d < spacedim; ++d)
    {
      eigenvectors[0][first_momentum_component + d] =
        -alp_f * a_s2 / (2. * theta_1) * u[d] -
        sgn_nb * alp_f * a_s / (2. * theta_2) * normal[d] +
        alp_s * c_s / (2. * theta_2) * n_perp[d];
      eigenvectors[0][first_magnetic_component + d] =
        alp_s * c_f * std::sqrt(state[density_component]) *
        (c_s2 + (2. - adiabatic_index) / (adiabatic_index - 1.) * a_s2) /
        (2. * theta_1) * n_perp[d];
    }


  // Left alfven mode l_2
  eigenvectors[1][density_component] = sgn_nb / 2. * normal * n_perp_cross_u;
  eigenvectors[1][energy_component]  = 0.;
  for (unsigned int d = 0; d < spacedim; ++d)
    {
      eigenvectors[1][first_momentum_component + d] =
        sgn_nb / 2. * n_perp_cross_normal[d];
      eigenvectors[1][first_magnetic_component + d] =
        std::sqrt(state[density_component]) / 2. * n_perp_cross_normal[d];
    }


  // Left slow magnetosonic mode l_3
  eigenvectors[2][density_component] =
    alp_s * c_f2 * u2 / (4. * theta_1) +
    (sgn_nb * alp_s * c_a * nu + alp_f * c_f * n_perp_u) / (2. * theta_2);
  eigenvectors[2][energy_component] = alp_s * c_f2 / (2. * theta_1);
  for (unsigned int d = 0; d < spacedim; ++d)
    {
      eigenvectors[2][first_momentum_component + d] =
        -alp_s * c_f2 / (2. * theta_1) * u[d] -
        sgn_nb * alp_s * c_a / (2. * theta_2) * normal[d] -
        alp_f * c_f / (2. * theta_2) * n_perp[d];
      eigenvectors[2][first_magnetic_component + d] =
        -alp_f * c_f * std::sqrt(state[density_component]) *
        (c_f2 + (2. - adiabatic_index) / (adiabatic_index - 1.) * a_s2) /
        (2. * theta_1) * n_perp[d];
    }


  // Density entropy mode l_4
  eigenvectors[3][density_component] =
    1. - (alp_f * alp_f * a_s2 + alp_s * alp_s * c_f2) * u2 / (2. * theta_1);
  eigenvectors[3][energy_component] =
    -(alp_f * alp_f * a_s2 + alp_s * alp_s * c_f2) / theta_1;
  for (unsigned int d = 0; d < spacedim; ++d)
    {
      eigenvectors[3][first_momentum_component + d] =
        (alp_f * alp_f * a_s2 + alp_s * alp_s * c_f2) / theta_1 * u[d];
      eigenvectors[3][first_magnetic_component + d] =
        alp_f * alp_s * c_f * std::sqrt(state[density_component]) *
        (c_f2 - c_s2) / theta_1 * n_perp[d];
    }


  // Magnetic entropy mode l_5
  eigenvectors[4][density_component] = 0;
  eigenvectors[4][energy_component]  = 0;
  for (unsigned int d = 0; d < spacedim; ++d)
    {
      eigenvectors[4][first_momentum_component + d] = 0;
      eigenvectors[4][first_magnetic_component + d] = normal[d];
    }


  // Right slow magnetosonic mode l_6
  eigenvectors[5][density_component] =
    alp_s * c_f2 * u2 / (4. * theta_1) -
    (sgn_nb * alp_s * c_a * nu + alp_f * c_f * n_perp_u) / (2. * theta_2);
  eigenvectors[5][energy_component] = alp_s * c_f2 / (2. * theta_1);
  for (unsigned int d = 0; d < spacedim; ++d)
    {
      eigenvectors[5][first_momentum_component + d] =
        -alp_s * c_f2 / (2. * theta_1) * u[d] +
        sgn_nb * alp_s * c_a / (2. * theta_2) * normal[d] +
        alp_f * c_f / (2. * theta_2) * n_perp[d];
      eigenvectors[5][first_magnetic_component + d] =
        -alp_f * c_f * std::sqrt(state[density_component]) *
        (c_f2 + (2. - adiabatic_index) / (adiabatic_index - 1.) * a_s2) /
        (2. * theta_1) * n_perp[d];
    }


  // Right alfven mode l_7
  eigenvectors[6][density_component] = -sgn_nb / 2. * normal * n_perp_cross_u;
  eigenvectors[6][energy_component]  = 0.;
  for (unsigned int d = 0; d < spacedim; ++d)
    {
      eigenvectors[6][first_momentum_component + d] =
        -sgn_nb / 2. * n_perp_cross_normal[d];
      eigenvectors[6][first_magnetic_component + d] =
        std::sqrt(state[density_component]) / 2. * n_perp_cross_normal[d];
    }


  // Right fast magnetosonic mode l_8
  eigenvectors[7][density_component] =
    alp_f * a_s2 * u2 / (4. * theta_1) -
    (sgn_nb * alp_f * a_s * nu - alp_s * c_s * n_perp_u) / (2. * theta_2);
  eigenvectors[7][energy_component] = alp_s * a_s2 / (2. * theta_1);
  for (unsigned int d = 0; d < spacedim; ++d)
    {
      eigenvectors[7][first_momentum_component + d] =
        -alp_f * a_s2 / (2. * theta_1) * u[d] +
        sgn_nb * alp_f * a_s / (2. * theta_2) * normal[d] -
        alp_s * c_s / (2. * theta_2) * n_perp[d];
      eigenvectors[7][first_magnetic_component + d] =
        alp_s * c_f * std::sqrt(state[density_component]) *
        (c_s2 + (2. - adiabatic_index) / (adiabatic_index - 1.) * a_s2) /
        (2. * theta_1) * n_perp[d];
    }
}



void
sapphirepp::MHD::MHDEquations::convert_primitive_to_conserved(
  const state_type &primitive_state,
  state_type       &conserved_state) const
{
  AssertDimension(primitive_state.size(), n_components);
  AssertDimension(conserved_state.size(), n_components);
  Assert(primitive_state[density_component] > 0.,
         dealii::ExcMessage("Negative density warning."));
  Assert(primitive_state[pressure_component] > 0.,
         dealii::ExcMessage("Negative pressure warning."));

  double u2 = 0.;
  double b2 = 0.;
  for (unsigned int d = 0; d < spacedim; ++d)
    {
      u2 += primitive_state[first_velocity_component + d] *
            primitive_state[first_velocity_component + d];
      b2 += primitive_state[first_magnetic_component + d] *
            primitive_state[first_magnetic_component + d];
    }

  const double energy =
    0.5 * primitive_state[density_component] * u2 +
    primitive_state[pressure_component] / (adiabatic_index - 1.) + 0.5 * b2;
  Assert(energy > 0., dealii::ExcMessage("Negative energy warning."));

  conserved_state[density_component] = primitive_state[density_component];
  conserved_state[energy_component]  = energy;
  for (unsigned int d = 0; d < spacedim; ++d)
    {
      conserved_state[first_momentum_component + d] =
        primitive_state[density_component] *
        primitive_state[first_velocity_component + d];
      conserved_state[first_magnetic_component + d] =
        primitive_state[first_magnetic_component + d];
    }
}



void
sapphirepp::MHD::MHDEquations::convert_conserved_to_primitive(
  const state_type &conserved_state,
  state_type       &primitive_state) const
{
  AssertDimension(conserved_state.size(), n_components);
  AssertDimension(primitive_state.size(), n_components);

  primitive_state[density_component]  = conserved_state[density_component];
  primitive_state[pressure_component] = compute_pressure(conserved_state);
  for (unsigned int d = 0; d < spacedim; ++d)
    {
      primitive_state[first_velocity_component + d] =
        conserved_state[first_momentum_component + d] /
        conserved_state[density_component];
      primitive_state[first_magnetic_component + d] =
        conserved_state[first_magnetic_component + d];
    }
}
