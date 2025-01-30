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
 * @file test/mhd/shock-tube/test-shock-tube.cpp
 * @author Florian Schulze (florian.schulze@mpi-hd.mpg.de)
 * @brief Implement tests for shock-tube example
 */

#include <deal.II/base/mpi.h>
#include <deal.II/base/parameter_handler.h>

#include <mpi.h>

#include "config.h"
#include "grid-data-function.h"
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

      std::string comparison_path = "athena-results";
      if (argc > 2)
        comparison_path = argv[2];

      std::string comparison_base_filename = "Sod";
      if (argc > 3)
        comparison_base_filename = argv[3];

      double max_L2_error = 1e-10;
      if (argc > 4)
        max_L2_error = std::stod(argv[4]);

      dealii::ParameterHandler prm;
      PhysicalParameters       physical_parameters;
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

      /** [Setup exact solution] */
      Utils::GridDataFunction<dim_mhd, MHDEquations<dim>::spacedim>
        exact_solution(comparison_path,
                       comparison_base_filename,
                       MHDEquations<dim>::n_components,
                       0.,
                       5);
      /** [Setup exact solution] */

      return test_run_mhd(mhd_parameters,
                          physical_parameters,
                          output_parameters,
                          exact_solution,
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
