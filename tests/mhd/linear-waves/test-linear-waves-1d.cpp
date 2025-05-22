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
 * @file test/mhd/linear-waves/test-linear-waves-1d.cpp
 * @author Florian Schulze (florian.schulze@mpi-hd.mpg.de)
 * @brief Implement 1d tests for linear-waves example
 */

#include <deal.II/base/mpi.h>
#include <deal.II/base/parameter_handler.h>

#include <mpi.h>

#include "config.h"
#include "mhd-parameters.h"
#include "output-parameters.h"
#include "sapphirepp-logstream.h"
#include "test-run-mhd.h"



const unsigned int dim = 1;



int
main(int argc, char *argv[])
{
  try
    {
      /** [Main function setup] */
      using namespace sapphirepp;
      using namespace MHD;
      dealii::Utilities::MPI::MPI_InitFinalize mpi_initialization(argc,
                                                                  argv,
                                                                  1);
      saplog.init(argc, argv);

      std::string parameter_filename = "parameter.prm";
      if (argc > 1)
        parameter_filename = argv[1];

      double max_L2_error = 1e-10;
      if (argc > 2)
        max_L2_error = std::stod(argv[2]);

      dealii::ParameterHandler prm;
      PhysicalParameters       physical_parameters(dim);
      Utils::OutputParameters  output_parameters;
      MHDParameters<dim>       mhd_parameters;

      physical_parameters.declare_parameters(prm);
      output_parameters.declare_parameters(prm);
      mhd_parameters.declare_parameters(prm);

      prm.parse_input(parameter_filename);

      physical_parameters.parse_parameters(prm);
      output_parameters.parse_parameters(prm);
      mhd_parameters.parse_parameters(prm);
      /** [Main function setup] */


      /** [Copy MHD parameter] */
      physical_parameters.box_length = std::vector<double>(dim);
      for (unsigned int d = 0; d < dim; ++d)
        {
          physical_parameters.box_length[d] =
            std::abs(mhd_parameters.p1[d] - mhd_parameters.p2[d]);

          AssertThrow((mhd_parameters.boundary_conditions[2 * d + 0] ==
                       MHD::BoundaryConditionsMHD::periodic) &&
                        (mhd_parameters.boundary_conditions[2 * d + 1] ==
                         MHD::BoundaryConditionsMHD::periodic),
                      dealii::ExcMessage("This example assumes periodic BC."));
        }
      /** [Copy MHD parameter] */

      /** [Setup analytic solution] */
      MHDEquations<dim, divergence_cleaning> mhd_equations(
        mhd_parameters.adiabatic_index);
      // Set divergence cleaning speed to arbitrary values
      mhd_equations.compute_hyperbolic_divergence_cleaning_speed(1., 1., 1);

      InitialConditionMHD<dim, divergence_cleaning> analytic_solution(
        physical_parameters, mhd_equations);
      /** [Setup analytic solution] */

      return test_run_mhd<dim>(mhd_parameters,
                               physical_parameters,
                               output_parameters,
                               analytic_solution,
                               max_L2_error);
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
