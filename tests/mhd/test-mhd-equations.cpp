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

#include "mhd-equations.h"
#include "sapphirepp-logstream.h"

using namespace sapphirepp;
using namespace MHD;



const double epsilon_d = 1e-6;



template <unsigned int dim, bool hdc>
typename MHDEquations<dim, hdc>::state_type
generate_state(const MHDEquations<dim, hdc> &mhd_equations)
{
  constexpr unsigned int n_components = MHDEquations<dim, hdc>::n_components;
  constexpr unsigned int density_component =
    MHDEquations<dim, hdc>::density_component;
  constexpr unsigned int energy_component =
    MHDEquations<dim, hdc>::energy_component;
  using state_type = typename MHDEquations<dim, hdc>::state_type;

  state_type state(n_components);

  //  0.25472 0.220437 0.0386911 0 0.815913 3.87308e-09 -1.89152e-11 0
  // state[0] = 0.25472;
  // state[1] = 0.220437;
  // state[2] = 0.0386911;
  // state[3] = 0.;
  // state[4] = 0.815913;
  // state[5] = 3.87308e-09;
  // state[6] = -1.89152e-11;
  // state[7] = 0.;

  state[0] = 1.;
  state[1] = 0.;
  state[2] = 0.;
  state[3] = 0.;
  state[4] = 0.6;
  state[5] = 0.;
  state[6] = 0.;
  state[7] = 0.;

  const double pressure = mhd_equations.compute_pressure_unsafe(state);
  saplog << "state = " << state << ", P = " << pressure << std::endl;
  AssertThrow(state[density_component] > 0., ExcNonAdmissibleState<dim>(state));
  AssertThrow(state[energy_component] > 0., ExcNonAdmissibleState<dim>(state));
  AssertThrow(pressure > 0., ExcNonAdmissibleState<dim>(state));

  return state;
}



template <unsigned int dim, bool hdc>
void
test_mhd_equation(const double adiabatic_index)
{
  constexpr unsigned int n_components = MHDEquations<dim, hdc>::n_components;
  // constexpr unsigned int density_component =
  //   MHDEquations<dim, hdc>::density_component;
  // constexpr unsigned int energy_component =
  //   MHDEquations<dim, hdc>::energy_component;
  using state_type = typename MHDEquations<dim, hdc>::state_type;

  saplog << "Test MHDEquations<" << dim << "," << hdc << ">(" << adiabatic_index
         << ")" << std::endl;
  MHDEquations<dim, hdc> mhd_equations(adiabatic_index);

  state_type conserved_state = generate_state(mhd_equations);
  saplog << "conserved_state: " << conserved_state << std::endl;

  // Test conversion to primitive state
  state_type primitive_state(n_components);
  state_type conserved_state_2(n_components);
  mhd_equations.convert_conserved_to_primitive(conserved_state,
                                               primitive_state);
  saplog << "primitive_state: " << primitive_state << std::endl;
  mhd_equations.convert_primitive_to_conserved(primitive_state,
                                               conserved_state_2);
  saplog << "conserved_state: " << conserved_state_2 << std::endl;
  // Compare conserved states to numerical precision
  for (unsigned int c = 0; c < n_components; ++c)
    AssertThrow(std::abs(conserved_state[c] - conserved_state_2[c]) < epsilon_d,
                dealii::ExcMessage("Problem with conversion"));
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

      constexpr unsigned int dim             = 1;
      constexpr bool         hdc             = false;
      const double           adiabatic_index = 1.4;

      test_mhd_equation<dim, hdc>(adiabatic_index);
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
