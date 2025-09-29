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
 * @file tests/vfp/convergence-study/test-convergence-study.cpp
 * @author Florian Schulze (florian.schulze@mpi-hd.mpg.de)
 * @brief Implement tests for convergence-study example
 */

#include <deal.II/base/discrete_time.h>
#include <deal.II/base/mpi.h>
#include <deal.II/base/parameter_handler.h>

#include <mpi.h>

#include <fstream>
#include <iostream>

#include "config.h"
#include "output-parameters.h"
#include "sapphirepp-logstream.h"
#include "test-run-vfp.h"
#include "vfp-parameters.h"



const unsigned int dim = sapphirepp::VFP::dimension;



int
main(int argc, char *argv[])
{
  try
    {
      /** [Main function setup] */
      using namespace sapphirepp;
      using namespace VFP;
      dealii::Utilities::MPI::MPI_InitFinalize mpi_initialization(argc,
                                                                  argv,
                                                                  1);
      saplog.init(argc, argv);

      std::string parameter_filename = "parameter.prm";
      if (argc > 1)
        parameter_filename = argv[1];

      double max_L2_error = VFPParameters<dim>::epsilon_d;
      if (argc > 2)
        max_L2_error = std::stod(argv[2]);

      dealii::ParameterHandler prm;
      PhysicalParameters       physical_parameters;
      Utils::OutputParameters  output_parameters;
      VFPParameters<dim>       vfp_parameters(vfp_flags);

      physical_parameters.declare_parameters(prm);
      output_parameters.declare_parameters(prm);
      vfp_parameters.declare_parameters(prm);

      prm.parse_input(parameter_filename);

      physical_parameters.parse_parameters(prm);
      output_parameters.parse_parameters(prm);
      vfp_parameters.parse_parameters(prm);
      /** [Main function setup] */


      /** [Copy VFP parameter] */
      physical_parameters.box_length =
        std::abs(vfp_parameters.p1[0] - vfp_parameters.p2[0]);
      physical_parameters.velocity = vfp_parameters.velocity;
      physical_parameters.gamma    = vfp_parameters.gamma;
      physical_parameters.charge   = vfp_parameters.charge;
      physical_parameters.mass     = vfp_parameters.mass;

      AssertThrow(dim == 1,
                  dealii::ExcMessage("This example assumes 'dim = 1'."));
      AssertThrow(vfp_parameters.expansion_order == 1,
                  dealii::ExcMessage(
                    "This example assumes 'Expansion order = 1'."));
      AssertThrow((vfp_parameters.boundary_conditions[0] ==
                   VFP::BoundaryConditions::periodic) &&
                    (vfp_parameters.boundary_conditions[1] ==
                     VFP::BoundaryConditions::periodic),
                  dealii::ExcMessage("This example assumes periodic BC."));
      /** [Copy VFP parameter] */


      /** [Setup exact solution] */
      const unsigned int system_size = (vfp_parameters.expansion_order + 1) *
                                       (vfp_parameters.expansion_order + 1);
      InitialValueFunction<dim> exact_solution(physical_parameters,
                                               system_size);
      const dealii::ComponentSelectFunction<dim> weight(0, system_size);
      /** [Setup exact solution] */

      /** [Start test run] */
      return test_run_vfp<dim>(vfp_parameters,
                               physical_parameters,
                               output_parameters,
                               exact_solution,
                               max_L2_error,
                               &weight);
      /** [Start test run] */
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
