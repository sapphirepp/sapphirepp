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
 * @file test-mhd-equations.cpp
 * @author Florian Schulze (florian.schulze@mpi-hd.mpg.de)
 * @brief Tests for @ref sapphirepp::MHD::MHDEquations
 */

#include <deal.II/base/exceptions.h>
#include <deal.II/base/utilities.h>

#include <cmath>
#include <limits>
#include <random>

#include "mhd-equations.h"
#include "mhd-parameters.h"
#include "sapphirepp-logstream.h"

using namespace sapphirepp;
using namespace MHD;



/** @ref sapphirepp::MHD::MHDParameters<1>::epsilon_d */
static constexpr double epsilon_d =
  sapphirepp::MHD::MHDParameters<1>::epsilon_d;



class RandomNumber
{
public:
  RandomNumber(double mean = 0., double stddev = 1.)
    : dist_normal{mean, stddev}
    , dist_uniform{0., 1.}
  {
    seed = std::random_device{}();
    saplog << "RandomNumber(" << mean << ", " << stddev << ")"
           << " with seed: " << seed << std::endl;
    rng.seed(seed);
  }



  double
  operator()()
  {
    return dist_normal(rng);
  }



  double
  rand_uniform(double low = 0., double high = 1.)
  {
    return low + (high - low) * dist_uniform(rng);
  }



  void
  set_seed(unsigned int s)
  {
    seed = s;
    rng.seed(seed);
    saplog << "Random seed set to: " << seed << std::endl;
  }



private:
  std::default_random_engine             rng;
  std::normal_distribution<double>       dist_normal;
  std::uniform_real_distribution<double> dist_uniform;
  unsigned int                           seed;
};



template <unsigned int dim>
dealii::Tensor<1, dim>
random_normal(RandomNumber &rnd)
{
  dealii::LogStream::Prefix p("random_normal", saplog);
  dealii::Tensor<1, dim>    normal;

  for (unsigned int d = 0; d < dim; ++d)
    normal[d] = rnd();
  normal /= normal.norm();

  saplog << "Generated normal: " << normal << ", norm = " << normal.norm()
         << std::endl;
  Assert(std::abs(normal.norm() - 1) < epsilon_d,
         dealii::ExcMessage("Normal vector is not normalized."));

  return normal;
}



template <unsigned int dim, bool divergence_cleaning>
void
random_state(const MHDEquations<dim, divergence_cleaning> &mhd_equations,
             RandomNumber                                 &rnd,
             const bool                                    magnetic_field,
             typename MHDEquations<dim, divergence_cleaning>::state_type &state)
{
  using MHDEqs           = MHDEquations<dim, divergence_cleaning>;
  const double mhd_floor = 10 * epsilon_d;

  dealii::LogStream::Prefix p("random_state", saplog);

  AssertDimension(state.size(), MHDEqs::n_components);

  for (unsigned int c = 0; c < MHDEqs::n_components; ++c)
    state[c] = rnd();
  state[MHDEqs::density_component] = std::abs(state[MHDEqs::density_component]);
  state[MHDEqs::density_component] =
    state[MHDEqs::density_component] < mhd_floor ?
      mhd_floor :
      state[MHDEqs::density_component];
  state[MHDEqs::energy_component] = std::abs(state[MHDEqs::energy_component]);
  state[MHDEqs::energy_component] =
    state[MHDEqs::energy_component] < mhd_floor ?
      mhd_floor :
      state[MHDEqs::energy_component];
  if (!magnetic_field)
    {
      for (unsigned int d = 0; d < MHDEqs::n_vec_components; ++d)
        state[MHDEqs::first_magnetic_component + d] = 0.;
      if (divergence_cleaning)
        state[MHDEqs::divergence_cleaning_component] = 0.;
    }
  // Use velocity as random variable for less vacuum states
  for (unsigned int d = 0; d < MHDEqs::n_vec_components; ++d)
    state[MHDEqs::first_momentum_component + d] =
      state[MHDEqs::first_momentum_component + d] *
      state[MHDEqs::density_component];


  double pressure = mhd_equations.compute_pressure_unsafe(state);

  // Test specific state
  // const std::string tmp =
  //   "3.17207e-05 -2.29659 1.42733 -1.66069 158721. 0.00000 0.00000 0.00000";
  // state                               = 0;
  // std::vector<double> hardcoded_state = dealii::Utilities::string_to_double(
  //   dealii::Utilities::split_string_list(tmp, " "));
  // for (unsigned int c = 0;
  //      c < MHDEqs::n_components && c < hardcoded_state.size();
  //      ++c)
  //   state[c] = hardcoded_state[c];
  // pressure = mhd_equations.compute_pressure_unsafe(state);

  if (pressure < mhd_floor)
    {
      saplog << "P_old=" << pressure << std::endl;
      // Ensure positive pressure by setting E' = E + dP/(gamma-1)
      state[MHDEqs::energy_component] +=
        (mhd_floor - pressure) / (mhd_equations.adiabatic_index - 1.);
      // TODO: Due to double precision the new pressure can be below mhd_floor
      pressure = mhd_equations.compute_pressure_unsafe(state);
    }

  saplog << "Generated state: " << state << ", P=" << pressure << std::endl;
  AssertThrow(state[MHDEqs::density_component] >= epsilon_d,
              ExcNonAdmissibleState<dim>(state));
  AssertThrow(state[MHDEqs::energy_component] >= epsilon_d,
              ExcNonAdmissibleState<dim>(state));
  AssertThrow(pressure >= epsilon_d, ExcNonAdmissibleState<dim>(state));
}



template <unsigned int dim, bool divergence_cleaning>
void
test_mhd_equation(RandomNumber      &rnd,
                  const unsigned int num            = 1,
                  const bool         magnetic_field = divergence_cleaning)
{
  using MHDEqs     = MHDEquations<dim, divergence_cleaning>;
  using state_type = typename MHDEqs::state_type;
  using flux_type  = typename MHDEqs::flux_type;

  const double adiabatic_index = rnd.rand_uniform(1., 2.);
  saplog << "Test MHDEquations with dim=" << dim
         << ", divergence_cleaning=" << divergence_cleaning
         << ", adiabatic_index=" << adiabatic_index << std::endl;
  saplog.push("MHDEquations<" + std::to_string(dim) + "," +
              std::to_string(divergence_cleaning) + ">");

  // Declare variables
  dealii::Tensor<1, dim> normal;
  state_type state(MHDEqs::n_components), tmp_1(MHDEqs::n_components),
    tmp_2(MHDEqs::n_components);
  flux_type                  flux_matrix;
  dealii::Vector<double>     eigenvalues(MHDEqs::n_components);
  dealii::FullMatrix<double> right_eigenvectors(MHDEqs::n_components),
    left_eigenvectors(MHDEqs::n_components), tmp_matrix(MHDEqs::n_components);
  dealii::IdentityMatrix                      id_matrix(MHDEqs::n_components);
  dealii::FullMatrix<double>                  identity(id_matrix);
  std::array<dealii::FullMatrix<double>, dim> left_matrices;
  std::array<dealii::FullMatrix<double>, dim> right_matrices;
  for (unsigned int d = 0; d < dim; ++d)
    {
      left_matrices[d]  = dealii::FullMatrix<double>(MHDEqs::n_components);
      right_matrices[d] = dealii::FullMatrix<double>(MHDEqs::n_components);
    }


  // Setup MHD equations
  MHDEqs mhd_equations(adiabatic_index);


  unsigned int skipped_tests = 0;
  for (unsigned int n = 0; n < num; ++n)
    {
      // Generate a state and normal
      random_state(mhd_equations, rnd, magnetic_field, state);
      normal                = random_normal<dim>(rnd);
      const double pressure = mhd_equations.compute_pressure(state);
      saplog << "state[" << n << "]:\t" << state << ", P=" << pressure
             << std::endl;
      saplog << "normal: \t" << normal << std::endl;
      dealii::LogStream::Prefix p("state[" + std::to_string(n) + "]", saplog);

      // Test conversion to primitive state
      mhd_equations.convert_conserved_to_primitive(state, tmp_1);
      saplog << "primitive: " << tmp_1 << std::endl;
      mhd_equations.convert_primitive_to_conserved(tmp_1, tmp_2);
      saplog << "conserved: " << tmp_2 << std::endl;
      for (unsigned int c = 0; c < MHDEqs::n_components; ++c)
        AssertThrow(std::abs(state[c] - tmp_2[c]) < epsilon_d,
                    dealii::ExcMessage("Problem with conversion to "
                                       "primitive state"));

      const double a_s2 =
        adiabatic_index * pressure / state[MHDEqs::density_component];
      saplog << "adiabatic sound speed a_s2=" << a_s2 << std::endl;
      if (a_s2 < std::sqrt(epsilon_d))
        {
          saplog << "The eigensystem has problems, "
                    "if the adiabatic sound speed is too small. "
                    "Skip remaining tests!"
                 << std::endl;
          skipped_tests++;
          continue;
        }

      // Setup hyperbolic divergence cleaning
      const double       dx        = 1.;
      const unsigned int fe_degree = 1;
      const double       max_eigenvalue =
        mhd_equations.compute_maximum_eigenvalue(state);
      const double dt = dx / ((2. * fe_degree + 1.) * max_eigenvalue);
      saplog << "dx=" << dx << ", fe_degree=" << fe_degree
             << ", max_eigenvalue=" << max_eigenvalue << ", dt=" << dt
             << std::endl;
      mhd_equations.compute_hyperbolic_divergence_cleaning_speed(dt,
                                                                 dx,
                                                                 fe_degree);

      if (max_eigenvalue > 1. / std::sqrt(epsilon_d))
        {
          saplog << "The eigensystem has problems, "
                    "if the max eigenvalue is too big. "
                    "Skip remaining tests!"
                 << std::endl;
          skipped_tests++;
          continue;
        }

      // Test flux and sources
      mhd_equations.compute_flux_matrix(state, flux_matrix);
      saplog << "Flux matrix: ";
      for (unsigned int c = 0; c < MHDEqs::n_components; ++c)
        saplog << flux_matrix[c] << " ";
      saplog << std::endl;
      if constexpr (divergence_cleaning)
        {
          mhd_equations.add_source_divergence_cleaning(state, tmp_1);
          saplog << "HDC source: " << tmp_1 << std::endl;
        }

      // Test eigenvalue
      const double maximum_normal_eigenvalue =
        mhd_equations.compute_maximum_normal_eigenvalue(state, normal);
      saplog << "maximum_normal_eigenvalue=" << maximum_normal_eigenvalue
             << std::endl;
      mhd_equations.compute_normale_eigenvalues(state, normal, eigenvalues);
      saplog << "Eigenvalues: " << eigenvalues << std::endl;


      // Test eigensystem
      saplog << "Test left and right eigenvector matrices:" << std::endl;
      saplog << "Right eigenvector matrix R:" << std::endl;
      mhd_equations.compute_right_eigenvector_matrix(state,
                                                     normal,
                                                     right_eigenvectors);
      right_eigenvectors.print(saplog, 12, 5);
      saplog << "Left eigenvector matrix L:" << std::endl;
      mhd_equations.compute_left_eigenvector_matrix(state,
                                                    normal,
                                                    left_eigenvectors);
      left_eigenvectors.print(saplog, 12, 5);

      // The rest of the tests have only a precision sqrt(epsilon_d)
      left_eigenvectors.mmult(tmp_matrix, right_eigenvectors);
      saplog << "L*R:" << std::endl;
      tmp_matrix.print(saplog, 12, 5);
      for (unsigned int i = 0; i < MHDEqs::n_components; ++i)
        for (unsigned int j = 0; j < MHDEqs::n_components; ++j)
          AssertThrow(std::abs(tmp_matrix[i][j] - identity[i][j]) <
                        std::sqrt(epsilon_d),
                      dealii::ExcMessage("Problem with eigenvectors: l_" +
                                         std::to_string(i + 1) + " * r_" +
                                         std::to_string(j + 1) +
                                         " != delta_ij"));

      right_eigenvectors.mmult(tmp_matrix, left_eigenvectors);
      saplog << "R*L:" << std::endl;
      tmp_matrix.print(saplog, 12, 5);
      if (state[MHDEqs::energy_component] / state[MHDEqs::density_component] >
          1. / std::sqrt(epsilon_d))
        {
          saplog << "The eigensystem has problems, "
                    "if the energy is much larger the density. "
                    "Skip remaining tests!"
                 << std::endl;
          skipped_tests++;
          continue;
        }
      for (unsigned int i = 0; i < MHDEqs::n_components; ++i)
        for (unsigned int j = 0; j < MHDEqs::n_components; ++j)
          AssertThrow(std::abs(tmp_matrix[i][j] - identity[i][j]) <
                        std::sqrt(epsilon_d),
                      dealii::ExcMessage("Problem with eigenvectors: (R*L)_" +
                                         std::to_string(i + 1) +
                                         std::to_string(j + 1) +
                                         " != delta_ij"));


      // Test transformation matrices
      mhd_equations.compute_transformation_matrices(state,
                                                    left_matrices,
                                                    right_matrices);

      // Test conversion to characteristic state
      mhd_equations.convert_conserved_to_characteristic(state, normal, tmp_1);
      saplog << "characteristic: " << tmp_1 << std::endl;
      mhd_equations.convert_characteristic_to_conserved(tmp_1, normal, tmp_2);
      saplog << "conserved: " << tmp_2 << std::endl;
      for (unsigned int c = 0; c < MHDEqs::n_components; ++c)
        AssertThrow(std::abs(state[c] - tmp_2[c]) < std::sqrt(epsilon_d),
                    dealii::ExcMessage("Problem with conversion to "
                                       "characteristic state"));
    }
  saplog.pop();
  saplog << "Needed to skip " << skipped_tests << "/" << num
         << " due to vacuum states." << std::endl;
}



int
main(int argc, char *argv[])
{
  try
    {
      dealii::Utilities::MPI::MPI_InitFinalize mpi_initialization(argc,
                                                                  argv,
                                                                  1);
      saplog.init(std::numeric_limits<unsigned int>::max(), true);
      saplog.depth_console(1);

      RandomNumber rnd{};
      // rnd.set_seed(1083480086);

      const unsigned int num = 10000;
      test_mhd_equation<1, false>(rnd, num);
      test_mhd_equation<2, false>(rnd, num);
      test_mhd_equation<3, false>(rnd, num);
      test_mhd_equation<1, true>(rnd, num);
      test_mhd_equation<2, true>(rnd, num);
      test_mhd_equation<3, true>(rnd, num);
    }
  catch (std::exception &exc)
    {
      sapphirepp::saplog.print_error(exc);
      return 1;
    }
  catch (...)
    {
      std::cerr << std::endl;
      std::cerr << "\n"
                << "----------------------------------------------------"
                << "\n"
                << "Unknown exception!" << "\n"
                << "Aborting!" << "\n"
                << "----------------------------------------------------"
                << std::endl;
      return 1;
    }
  return 0;
}
