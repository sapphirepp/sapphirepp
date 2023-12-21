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
 * @file tests/gyro-motion-f0/test-gyro-motion-f0.cpp
 * @author Florian Schulze (florian.schulze@mpi-hd.mpg.de)
 * @brief Implement tests for gyro-motion-f0 example
 */

#include <deal.II/base/mpi.h>
#include <deal.II/base/parameter_handler.h>

#include <mpi.h>

#include "config.h"
#include "output-parameters.h"
#include "sapphirepp-logstream.h"
#include "vfp-parameters.h"
#include "vfp-solver.h"



int
main(int argc, char *argv[])
{
  try
    {
      using namespace sapphirepp;
      using namespace VFP;

      dealii::Utilities::MPI::MPI_InitFinalize mpi_initialization(argc,
                                                                  argv,
                                                                  1);

      saplog.init(2);

      std::string parameter_filename = "parameter.prm";
      if (argc > 1)
        parameter_filename = argv[1];

      double max_L2_error = 1e-10;
      if (argc > 2)
        max_L2_error = std::stod(argv[2]);

      dealii::Timer            timer;
      ParameterHandler         prm;
      VFPParameters<dimension> vfp_parameters;
      PhysicalParameters       physical_parameters;
      Utils::OutputParameters  output_parameters;

      vfp_parameters.declare_parameters(prm);
      physical_parameters.declare_parameters(prm);
      output_parameters.declare_parameters(prm);

      prm.parse_input(parameter_filename);

      vfp_parameters.parse_parameters(prm);
      physical_parameters.parse_parameters(prm);
      output_parameters.parse_parameters(prm);

      timer.start();
      VFPSolver<dimension> vfp_solver(vfp_parameters,
                                      physical_parameters,
                                      output_parameters);
      vfp_solver.run();
      timer.stop();


      saplog << "Calculate analytic solution" << std::endl;
      const double gyroperiod =
        2. * M_PI * vfp_parameters.gamma * vfp_parameters.mass /
        (physical_parameters.B0 * vfp_parameters.charge);
      AssertThrow(std::fmod(vfp_parameters.final_time, gyroperiod) == 0.,
                  dealii::ExcMessage(
                    "Final time is not a multiple of the gyroperiod.\n\t" +
                    dealii::Utilities::to_string(vfp_parameters.final_time) +
                    " != n * " + dealii::Utilities::to_string(gyroperiod)));


      InitialValueFunction<dimension> analytic_solution(
        physical_parameters, vfp_parameters.expansion_order);

      const ComponentSelectFunction<dimension> weight(
        0,
        (vfp_parameters.expansion_order + 1) *
          (vfp_parameters.expansion_order + 1));


      const double L2_error = vfp_solver.compute_global_error(
        analytic_solution, VectorTools::L2_norm, VectorTools::L2_norm, &weight);
      const double L2_norm =
        vfp_solver.compute_weighted_norm(VectorTools::L2_norm,
                                         VectorTools::L2_norm,
                                         &weight);


      saplog << "L2_error = " << L2_error << ", L2_norm = " << L2_norm
             << ", rel error = " << L2_error / L2_norm
             << ", CPU/wall time = " << timer.cpu_time() << "/"
             << timer.wall_time() << " s" << std::endl;

      AssertThrow((L2_error / L2_norm) < max_L2_error,
                  dealii::ExcMessage(
                    "L2 error is too large! (" +
                    dealii::Utilities::to_string(L2_error / L2_norm) + " > " +
                    dealii::Utilities::to_string(max_L2_error) + ")"));
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
