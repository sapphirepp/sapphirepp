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
 * @file tests/vfpgyro-advection/test-advection.cpp
 * @author Florian Schulze (florian.schulze@mpi-hd.mpg.de)
 * @brief Implement tests for gyro-advection example
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
      saplog.init(argc, argv);

      std::string parameter_filename = "parameter.prm";
      if (argc > 1)
        parameter_filename = argv[1];

      double max_L2_error = VFPParameters<dimension>::epsilon_d;
      if (argc > 2)
        max_L2_error = std::stod(argv[2]);

      dealii::Timer            timer;
      ParameterHandler         prm;
      VFPParameters<dimension> vfp_parameters(vfp_flags);
      PhysicalParameters       physical_parameters;
      Utils::OutputParameters  output_parameters;

      vfp_parameters.declare_parameters(prm);
      physical_parameters.declare_parameters(prm);
      output_parameters.declare_parameters(prm);

      prm.parse_input(parameter_filename);

      vfp_parameters.parse_parameters(prm);
      physical_parameters.parse_parameters(prm);
      output_parameters.parse_parameters(prm);

      const double gyroperiod =
        vfp_parameters.gamma * vfp_parameters.mass /
        (physical_parameters.B0 / (2 * M_PI) * vfp_parameters.charge);
      const double gyroradius =
        vfp_parameters.gamma * vfp_parameters.mass * vfp_parameters.velocity /
        (vfp_parameters.charge * physical_parameters.B0);
      const double box_length =
        std::abs(vfp_parameters.p1[0] - vfp_parameters.p2[0]);
      const double crossing_time = box_length / physical_parameters.u0;
      saplog << "particle_velocity = " << vfp_parameters.velocity
             << ", gyroperiod = " << gyroperiod
             << ", gyroradius = " << gyroradius
             << ", box_length = " << box_length
             << ", crossing_time = " << crossing_time
             << ", final_time = " << vfp_parameters.final_time << std::endl;

      AssertThrow(std::fmod(vfp_parameters.final_time, gyroperiod) == 0.,
                  dealii::ExcMessage(
                    "Final time is not a multiple of the gyroperiod.\n\t" +
                    dealii::Utilities::to_string(vfp_parameters.final_time) +
                    " != n * " + dealii::Utilities::to_string(gyroperiod)));

      AssertThrow((std::fmod(vfp_parameters.final_time, crossing_time) == 0.) ||
                    (physical_parameters.u0 == 0.),
                  dealii::ExcMessage(
                    "Final time is not a multiple of the crossing time.\n\t" +
                    dealii::Utilities::to_string(vfp_parameters.final_time) +
                    " != n * " + dealii::Utilities::to_string(crossing_time)));

      timer.start();
      VFPSolver<dimension> vfp_solver(vfp_parameters,
                                      physical_parameters,
                                      output_parameters);
      vfp_solver.run();
      timer.stop();


      saplog << "Calculate analytic solution" << std::endl;
      InitialValueFunction<dimension> analytic_solution(
        physical_parameters, vfp_solver.get_pde_system().system_size);

      const dealii::ComponentSelectFunction<dimension> weight(
        0, vfp_solver.get_pde_system().system_size);

      const double L2_error =
        vfp_solver.compute_global_error(analytic_solution,
                                        dealii::VectorTools::L2_norm,
                                        dealii::VectorTools::L2_norm,
                                        &weight);
      const double L2_norm =
        vfp_solver.compute_weighted_norm(dealii::VectorTools::L2_norm,
                                         dealii::VectorTools::L2_norm,
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
