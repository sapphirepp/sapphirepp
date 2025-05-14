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

#include <cmath>
#include <limits>
#include <random>

#include "mhd-equations.h"
#include "sapphirepp-logstream.h"

using namespace sapphirepp;
using namespace MHD;



const double epsilon_d = 1e-6;



class RandomNumber
{
public:
  RandomNumber(double mean = 0., double stddev = 1.)
    : dist_normal{mean, stddev}
    , dist_uniform{0., 1.}
  {
    seed = std::random_device{}();
    saplog << "RandomNumber(" << mean << ", " << stddev << ")"
           << " with seed: \n"
           << seed << std::endl;
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



template <unsigned int dim, bool hdc>
void
random_state(const MHDEquations<dim, hdc>                &mhd_equations,
             RandomNumber                                &rnd,
             const bool                                   magnetic_field,
             typename MHDEquations<dim, hdc>::state_type &state)
{
  constexpr unsigned int n_components = MHDEquations<dim, hdc>::n_components;
  constexpr unsigned int n_vec_components =
    MHDEquations<dim, hdc>::n_vec_components;
  constexpr unsigned int density_component =
    MHDEquations<dim, hdc>::density_component;
  constexpr unsigned int energy_component =
    MHDEquations<dim, hdc>::energy_component;
  constexpr unsigned int first_magnetic_component =
    MHDEquations<dim, hdc>::first_magnetic_component;
  constexpr unsigned int divergence_cleaning_component =
    MHDEquations<dim, hdc>::divergence_cleaning_component;
  using state_type = typename MHDEquations<dim, hdc>::state_type;

  dealii::LogStream::Prefix p("random_state", saplog);

  AssertDimension(state.size(), n_components);

  for (unsigned int c = 0; c < n_components; ++c)
    state[c] = rnd();
  state[density_component] = std::abs(state[density_component]);
  state[energy_component]  = std::abs(state[energy_component]);
  if (!magnetic_field)
    {
      for (unsigned int d = 0; d < n_vec_components; ++d)
        state[first_magnetic_component + d] = 0.;
      if (hdc)
        state[divergence_cleaning_component] = 0.;
    }

  // Test specific state
  // state = 0;
  // state[0] = 0.25472;
  // state[1] = 0.220437;
  // state[2] = 0.0386911;
  // state[3] = 0.;
  // state[4] = 0.815913;
  // state[5] = 3.87308e-09;
  // state[6] = -1.89152e-11;
  // state[7] = 0.;

  double pressure = mhd_equations.compute_pressure_unsafe(state);
  if (pressure < 0)
    {
      // Ensure positive pressure by setting: E' = 2E - P/(gamma-1)
      state[energy_component] = 2 * state[energy_component] -
                                pressure / (mhd_equations.adiabatic_index - 1.);
      pressure = mhd_equations.compute_pressure_unsafe(state);
    }

  saplog << "Generated state: " << state << ", P = " << pressure << std::endl;
  AssertThrow(state[density_component] > 0., ExcNonAdmissibleState<dim>(state));
  AssertThrow(state[energy_component] > 0., ExcNonAdmissibleState<dim>(state));
  AssertThrow(pressure > 0., ExcNonAdmissibleState<dim>(state));
}



template <unsigned int dim, bool hdc>
void
test_mhd_equation(RandomNumber      &rnd,
                  const unsigned int num            = 1,
                  const bool         magnetic_field = hdc)
{
  constexpr unsigned int n_components = MHDEquations<dim, hdc>::n_components;
  using state_type = typename MHDEquations<dim, hdc>::state_type;
  using flux_type  = typename MHDEquations<dim, hdc>::flux_type;

  dealii::LogStream::Prefix p("MHDEquations<" + std::to_string(dim) + "," +
                                std::to_string(hdc) + ">",
                              saplog);

  const double adiabatic_index = rnd.rand_uniform(1., 2.);
  saplog << "Test MHDEquations with dim=" << dim
         << ", divergence_cleaning=" << hdc
         << ", adiabatic_index=" << adiabatic_index << std::endl;

  // Declare variables
  dealii::Tensor<1, dim> normal;
  state_type state(n_components), tmp_1(n_components), tmp_2(n_components);
  flux_type  flux_matrix;
  dealii::Vector<double>     eigenvalues(n_components);
  dealii::FullMatrix<double> right_eigenvectors(n_components),
    left_eigenvectors(n_components), tmp_matrix(n_components);
  dealii::IdentityMatrix                      id_matrix(n_components);
  dealii::FullMatrix<double>                  identity(id_matrix);
  std::array<dealii::FullMatrix<double>, dim> left_matrices;
  std::array<dealii::FullMatrix<double>, dim> right_matrices;
  for (unsigned int d = 0; d < dim; ++d)
    {
      left_matrices[d]  = dealii::FullMatrix<double>(n_components);
      right_matrices[d] = dealii::FullMatrix<double>(n_components);
    }


  // Setup MHD equations
  MHDEquations<dim, hdc> mhd_equations(adiabatic_index);


  for (unsigned int n = 0; n < num; ++n)
    {
      // Generate a state and normal
      random_state(mhd_equations, rnd, magnetic_field, state);
      normal = random_normal<dim>(rnd);
      saplog << "state[" << n << "]:\t" << state << std::endl;
      saplog << "normal: \t" << normal << std::endl;
      dealii::LogStream::Prefix p("state[" + std::to_string(n) + "]", saplog);

      // Test conversion to primitive state
      mhd_equations.convert_conserved_to_primitive(state, tmp_1);
      saplog << "primitive: " << tmp_1 << std::endl;
      mhd_equations.convert_primitive_to_conserved(tmp_1, tmp_2);
      saplog << "conserved: " << tmp_2 << std::endl;
      for (unsigned int c = 0; c < n_components; ++c)
        AssertThrow(std::abs(state[c] - tmp_2[c]) < epsilon_d,
                    dealii::ExcMessage("Problem with conversion to "
                                       "primitive state"));

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

      // Test flux and sources
      mhd_equations.compute_flux_matrix(state, flux_matrix);
      saplog << "Flux matrix: ";
      for (unsigned int c = 0; c < n_components; ++c)
        saplog << "\n " << flux_matrix[c];
      saplog << std::endl;
      if constexpr (hdc)
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

      left_eigenvectors.mmult(tmp_matrix, right_eigenvectors);
      saplog << "L*R:" << std::endl;
      tmp_matrix.print(saplog, 12, 5);
      for (unsigned int i = 0; i < n_components; ++i)
        for (unsigned int j = 0; j < n_components; ++j)
          AssertThrow(std::abs(tmp_matrix[i][j] - identity[i][j]) < epsilon_d,
                      dealii::ExcMessage("Problem with eigenvectors: l_" +
                                         std::to_string(i + 1) + " * r_" +
                                         std::to_string(j + 1) +
                                         " != delta_ij"));

      right_eigenvectors.mmult(tmp_matrix, left_eigenvectors);
      saplog << "R*L:" << std::endl;
      tmp_matrix.print(saplog, 12, 5);
      for (unsigned int i = 0; i < n_components; ++i)
        for (unsigned int j = 0; j < n_components; ++j)
          AssertThrow(std::abs(tmp_matrix[i][j] - identity[i][j]) < epsilon_d,
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
      for (unsigned int c = 0; c < n_components; ++c)
        AssertThrow(std::abs(state[c] - tmp_2[c]) < epsilon_d,
                    dealii::ExcMessage("Problem with conversion to "
                                       "characteristic state"));
    }
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
      saplog.depth_console(2);

      RandomNumber rnd{};
      // rnd.set_seed(1083480086);

      const unsigned int num = 1000;
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
