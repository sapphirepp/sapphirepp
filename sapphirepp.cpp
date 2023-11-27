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
 * @file sapphirepp.cpp
 * @author Nils Schween (nils.schween@mpi-hd.mpg.de)
 * @author Florian Schulze (florian.schulze@mpi-hd.mpg.de)
 * @brief Implement main function for @sapphire
 */

#include <deal.II/base/mpi.h>

#include <mpi.h>

#include "output-parameters.h"
#include "vfp-solver.h"

// Deactivating the spatial advection term is equivalent to assuming a
// homogeneous distribution function (i.e. a distribution function which
// does not depend on x,y z). In this program this is equivalent to set
// dimension of the configuration to zero.
const int dim = 2;

int
main(int argc, char *argv[])
{
  try
    {
      using namespace sapphirepp;
      using namespace VFP;
      Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

      saplog.depth_console(2);
      const unsigned int mpi_rank =
        Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);
      if (mpi_rank > 0)
        {
          saplog.depth_console(0);
          char buffer[4];
          std::snprintf(buffer, 4, "%03d", mpi_rank);
          saplog.push("mpi" + std::string(buffer));
        }

      std::string parameter_filename = "parameter-template.prm";
      if (argc > 1)
        parameter_filename = argv[1];

      unsigned int mpi_size = Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD);
      saplog << "Start sapphire with parameter file \"" << parameter_filename
             << "\" on " << mpi_size << " processor(s) ["
             << Utilities::System::get_date() << " "
             << Utilities::System::get_time() << "]" << std::endl;

      saplog.push("main");
      ParameterHandler        prm;
      VFPParameters<dim>      vfp_parameters;
      PhysicalParameters      physical_parameters;
      Utils::OutputParameters output_parameters;

      vfp_parameters.declare_parameters(prm);
      physical_parameters.declare_parameters(prm);
      output_parameters.declare_parameters(prm);

      if (mpi_rank == 0)
        {
          saplog << "Writing template parameter file \t["
                 << Utilities::System::get_time() << "]" << std::endl;
          prm.print_parameters("parameter-template.prm", ParameterHandler::PRM);
        }

      prm.parse_input(parameter_filename);

      vfp_parameters.parse_parameters(prm);
      physical_parameters.parse_parameters(prm);
      output_parameters.parse_parameters(prm);

      saplog.pop();
      VFPSolver<dim> vfp_solver(vfp_parameters,
                                physical_parameters,
                                output_parameters);
      vfp_solver.run();
    }
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
