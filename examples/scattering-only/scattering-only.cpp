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

/// [Includes]
#include <deal.II/base/mpi.h>

#include <mpi.h>

#include "config.h"
#include "output-module.h"
#include "vfp-equation-solver.h"
/// [Includes]

/// [main func]
int
main(int argc, char *argv[])
{
  /// [main func]
  /// [try catch block begin]
  try
    {
      /// [try catch block begin]
      /// [MPI init]
      using namespace Sapphire;
      using namespace VFP;
      dealii::Utilities::MPI::MPI_InitFinalize mpi_initialization(argc,
                                                                  argv,
                                                                  1);
      /// [MPI init]
      /// [Saplog]
      saplog.depth_console(2);
      const unsigned int mpi_rank =
        dealii::Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);
      if (mpi_rank > 0)
        {
          saplog.depth_console(0);
          saplog.push("mpi" + dealii::Utilities::to_string(mpi_rank, 3));
        }
      /// [Saplog]
      /// [Command line argument]
      std::string parameter_filename = "parameter.prm";
      if (argc > 1)
        parameter_filename = argv[1];

      int mpi_size = dealii::Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD);
      saplog << "Start example scattering-only with parameter file \""
             << parameter_filename << "\" on " << mpi_size << " processor(s) ["
             << dealii::Utilities::System::get_date() << " "
             << dealii::Utilities::System::get_time() << "]" << std::endl;
      /// [Command line argument]

      /// [Run time params]
      saplog.push("main");
      ParameterHandler               prm;
      VFPSolverControl<dimension>    vfp_solver_control;
      PhysicalProperties             physical_properties;
      Utils::OutputModule<dimension> output_module;
      /// [Run time params]

      /// [Declare params]
      vfp_solver_control.declare_parameters(prm);
      physical_properties.declare_parameters(prm);
      output_module.declare_parameters(prm);
      /// [Declare params]
      /// [Parse params]
      prm.parse_input(parameter_filename);

      vfp_solver_control.parse_parameters(prm);
      physical_properties.parse_parameters(prm);
      output_module.parse_parameters(prm);

      saplog.pop(); // pop "main" prefix
      /// [Parse params]

      /// [VFP Solver]
      VFPEquationSolver<dimension> vfp_equation_solver(vfp_solver_control,
                                                       physical_properties,
                                                       output_module);
      vfp_equation_solver.run();
      /// [VFP Solver]

      /// [L2 error]
      InitialValueFunction<dimension> analytic_solution(
        physical_properties, vfp_solver_control.expansion_order);
      analytic_solution.set_time(vfp_solver_control.final_time);

      const double L2_error =
        vfp_equation_solver.compute_global_error(analytic_solution,
                                                 VectorTools::L2_norm,
                                                 VectorTools::L2_norm);

      saplog << "L2 error: " << L2_error << std::endl;
    }
  /// [L2 error]
  /// [try catch block end]
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
/// [try catch block end]
