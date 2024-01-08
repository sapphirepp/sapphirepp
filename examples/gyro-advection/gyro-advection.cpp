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
 * @file examples/gyro-advection/gyro-advection.cpp
 * @author Florian Schulze (florian.schulze@mpi-hd.mpg.de)
 * @brief Implement main function for gyro-advection example
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

      saplog.init(5);

      std::string parameter_filename = "parameter.prm";
      if (argc > 1)
        parameter_filename = argv[1];

      dealii::ParameterHandler prm;
      PhysicalParameters       physical_parameters;
      Utils::OutputParameters  output_parameters;
      VFPParameters<dimension> vfp_parameters;

      physical_parameters.declare_parameters(prm);
      output_parameters.declare_parameters(prm);
      vfp_parameters.declare_parameters(prm);

      prm.parse_input(parameter_filename);

      physical_parameters.parse_parameters(prm);
      output_parameters.parse_parameters(prm);
      vfp_parameters.parse_parameters(prm);

      const double gyroperiod =
        vfp_parameters.gamma * vfp_parameters.mass /
        (physical_parameters.B0_2pi * vfp_parameters.charge);
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

      VFPSolver<dimension> vfp_solver(vfp_parameters,
                                      physical_parameters,
                                      output_parameters);
      vfp_solver.run();

      saplog << "Calculate analytic solution" << std::endl;

      if (std::fmod(vfp_parameters.final_time, gyroperiod) != 0.)
        saplog << "Warning: Final time is not a multiple of the gyroperiod."
               << std::endl;

      if ((std::fmod(vfp_parameters.final_time, crossing_time) != 0.) &&
          (physical_parameters.u0 != 0.))
        saplog << "Warning: Final time is not a multiple of the crossing time."
               << std::endl;

      InitialValueFunction<dimension> analytic_solution(
        physical_parameters, vfp_parameters.expansion_order);

      const dealii::ComponentSelectFunction<dimension> weight(
        0,
        (vfp_parameters.expansion_order + 1) *
          (vfp_parameters.expansion_order + 1));

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
             << ", rel error = " << L2_error / L2_norm << std::endl;
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
