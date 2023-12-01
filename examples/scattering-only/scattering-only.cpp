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
 * @file examples/scattering-only/scattering-only.cpp
 * @author Florian Schulze (florian.schulze@mpi-hd.mpg.de)
 * @brief Implement main function for scattering-only example
 */

/** [Includes] */
#include <deal.II/base/mpi.h>
#include <deal.II/base/parameter_handler.h>

#include <mpi.h>

#include "config.h"
#include "output-parameters.h"
#include "sapphirepp-logstream.h"
#include "vfp-parameters.h"
#include "vfp-solver.h"
/** [Includes] */

/** [Main function] */
int
main(int argc, char *argv[])
{
  /** [Main function] */
  /** [Try-Catch begin] */
  try
    {
      using namespace sapphirepp;
      using namespace VFP;
      /** [Try-Catch begin] */
      /** [MPI initialization] */
      dealii::Utilities::MPI::MPI_InitFinalize mpi_initialization(argc,
                                                                  argv,
                                                                  1);
      /** [MPI initialization] */
      /** [Saplog] */
      saplog.init(2);
      /** [Saplog] */
      /** [Command line argument] */
      std::string parameter_filename = "parameter.prm";
      if (argc > 1)
        parameter_filename = argv[1];
      saplog << "Start example scattering-only with parameter file \""
             << parameter_filename << "\"" << std::endl;
      /** [Command line argument] */

      /** [Run time parameters] */
      dealii::ParameterHandler prm;
      PhysicalParameters       physical_parameters;
      Utils::OutputParameters  output_parameters;
      VFPParameters<dimension> vfp_parameters;
      /** [Run time parameters] */

      /** [Declare parameters] */
      physical_parameters.declare_parameters(prm);
      output_parameters.declare_parameters(prm);
      vfp_parameters.declare_parameters(prm);
      /** [Declare parameters] */
      /** [Parse parameters] */
      prm.parse_input(parameter_filename);

      physical_parameters.parse_parameters(prm);
      output_parameters.parse_parameters(prm);
      vfp_parameters.parse_parameters(prm);
      /** [Parse parameters] */

      /** [VFP Solver] */
      VFPSolver<dimension> vfp_solver(vfp_parameters,
                                      physical_parameters,
                                      output_parameters);
      vfp_solver.run();
      /** [VFP Solver] */

      /** [L2 error] */
      InitialValueFunction<dimension> analytic_solution(
        physical_parameters, vfp_parameters.expansion_order);
      analytic_solution.set_time(vfp_parameters.final_time);

      const double L2_error =
        vfp_solver.compute_global_error(analytic_solution,
                                        VectorTools::L2_norm,
                                        VectorTools::L2_norm);

      saplog << "L2 error: " << L2_error << std::endl;
    }
  /** [L2 error] */
  /** [Try-Catch end] */
  catch (std::exception &exc)
    {
      std::cerr << std::endl
                << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Exception on processing: " << std::endl
                << exc.what() << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      return 1;
    }
  catch (...)
    {
      std::cerr << std::endl
                << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Unknown exception!" << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      return 1;
    }
  return 0;
}
/** [Try-Catch end] */
