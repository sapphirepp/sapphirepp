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
 * @file examples/convergence-study/convergence-study.cpp
 * @author Florian Schulze (florian.schulze@mpi-hd.mpg.de)
 * @brief Implement main function for convergence-study example
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
#include "vfp-parameters.h"
#include "vfp-solver.h"



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
      /** [Main function setup] */


      /** [Copy VFP parameter] */
      physical_parameters.box_length =
        std::abs(vfp_parameters.p1[0] - vfp_parameters.p2[0]);
      physical_parameters.velocity = vfp_parameters.velocity;
      physical_parameters.gamma    = vfp_parameters.gamma;
      physical_parameters.charge   = vfp_parameters.charge;
      physical_parameters.mass     = vfp_parameters.mass;

      AssertThrow(dimension == 1,
                  dealii::ExcMessage("This example assumes 'dimension = 1'."));
      AssertThrow(vfp_parameters.expansion_order == 1,
                  dealii::ExcMessage(
                    "This example assumes 'Expansion order = 1'."));
      AssertThrow((vfp_parameters.boundary_conditions[0] ==
                   VFP::BoundaryConditions::periodic) &&
                    (vfp_parameters.boundary_conditions[1] ==
                     VFP::BoundaryConditions::periodic),
                  dealii::ExcMessage("This example assumes periodic BC."));
      /** [Copy VFP parameter] */


      /** [Create error file] */
      std::ofstream error_file(output_parameters.output_path / "error.csv");
      AssertThrow(!error_file.fail(),
                  dealii::ExcFileNotOpen(output_parameters.output_path /
                                         "error.csv"));
      error_file << "timestep"
                 << "; "
                 << "time"
                 << "; "
                 << "L2_norm"
                 << "; "
                 << "L2_error"
                 << "; "
                 << "relative_error" << std::endl;
      /** [Create error file] */


      /** [Setup vfp_solver] */
      VFPSolver<dimension> vfp_solver(vfp_parameters,
                                      physical_parameters,
                                      output_parameters);
      vfp_solver.setup();
      /** [Setup vfp_solver] */


      /** [Setup analytic solution] */
      const unsigned int system_size = vfp_solver.get_pde_system().system_size;

      InitialValueFunction<dimension> analytic_solution(physical_parameters,
                                                        system_size);
      const dealii::ComponentSelectFunction<dimension> weight(0, system_size);

      PETScWrappers::MPI::Vector analytic_solution_vector;
      analytic_solution_vector.reinit(
        vfp_solver.get_dof_handler().locally_owned_dofs(), MPI_COMM_WORLD);
      /** [Setup analytic solution] */


      /** [Time loop] */
      DiscreteTime discrete_time(0,
                                 vfp_parameters.final_time,
                                 vfp_parameters.time_step);
      for (; discrete_time.is_at_end() == false; discrete_time.advance_time())
        {
          saplog << "Time step " << std::setw(6) << std::right
                 << discrete_time.get_step_number()
                 << " at t = " << discrete_time.get_current_time() << " \t["
                 << Utilities::System::get_time() << "]" << std::endl;

          analytic_solution.set_time(discrete_time.get_current_time());
          /** [Time loop] */


          /** [Output solution] */
          if ((discrete_time.get_step_number() %
               output_parameters.output_frequency) == 0)
            {
              LogStream::Prefix prefix("Output", saplog);
              saplog << "Output solution" << std::endl;

              dealii::DataOut<dimension> data_out;
              data_out.attach_dof_handler(vfp_solver.get_dof_handler());

              // Output numeric solution
              data_out.add_data_vector(vfp_solver.get_current_solution(),
                                       PDESystem::create_component_name_list(
                                         system_size));

              // Output projected analytic solution
              vfp_solver.project(analytic_solution, analytic_solution_vector);
              data_out.add_data_vector(analytic_solution_vector,
                                       PDESystem::create_component_name_list(
                                         system_size, "project_f_"));

              // Output interpolated analytic solution
              dealii::VectorTools::interpolate(vfp_solver.get_dof_handler(),
                                               analytic_solution,
                                               analytic_solution_vector);
              data_out.add_data_vector(analytic_solution_vector,
                                       PDESystem::create_component_name_list(
                                         system_size, "interpol_f_"));

              data_out.build_patches(vfp_parameters.polynomial_degree);
              output_parameters.write_results<dimension>(
                data_out,
                discrete_time.get_step_number(),
                discrete_time.get_current_time());
            }
          /** [Output solution] */


          /** [Calculate error] */
          {
            LogStream::Prefix prefix2("Error", saplog);
            saplog << "Calculate error" << std::endl;

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

            error_file << discrete_time.get_step_number() << "; "
                       << discrete_time.get_current_time() << "; " << L2_norm
                       << "; " << L2_error << "; " << L2_error / L2_norm
                       << std::endl;
          }
          /** [Calculate error] */


          /** [Time step] */
          switch (vfp_parameters.time_stepping_method)
            {
              case TimeSteppingMethod::forward_euler:
              case TimeSteppingMethod::backward_euler:
              case TimeSteppingMethod::crank_nicolson:
                vfp_solver.theta_method(discrete_time.get_current_time(),
                                        discrete_time.get_next_step_size());
                break;
              case TimeSteppingMethod::erk4:
                vfp_solver.explicit_runge_kutta(
                  discrete_time.get_current_time(),
                  discrete_time.get_next_step_size());
                break;
              case TimeSteppingMethod::lserk:
                vfp_solver.low_storage_explicit_runge_kutta(
                  discrete_time.get_current_time(),
                  discrete_time.get_next_step_size());
                break;
              default:
                AssertThrow(false, ExcNotImplemented());
            }
        }
      /** [Time step] */


      /** [Last time step] */
      analytic_solution.set_time(discrete_time.get_current_time());

      {
        LogStream::Prefix prefix("Output", saplog);
        saplog << "Output results" << std::endl;

        dealii::DataOut<dimension> data_out;
        data_out.attach_dof_handler(vfp_solver.get_dof_handler());

        // Output numeric solution
        data_out.add_data_vector(vfp_solver.get_current_solution(),
                                 PDESystem::create_component_name_list(
                                   system_size));

        // Output projected analytic solution
        vfp_solver.project(analytic_solution, analytic_solution_vector);
        data_out.add_data_vector(
          analytic_solution_vector,
          PDESystem::create_component_name_list(system_size, "project_f_"));

        // Output interpolated analytic solution
        dealii::VectorTools::interpolate(vfp_solver.get_dof_handler(),
                                         analytic_solution,
                                         analytic_solution_vector);
        data_out.add_data_vector(
          analytic_solution_vector,
          PDESystem::create_component_name_list(system_size, "interpol_f_"));

        data_out.build_patches(vfp_parameters.polynomial_degree);
        output_parameters.write_results<dimension>(
          data_out,
          discrete_time.get_step_number(),
          discrete_time.get_current_time());
      }

      {
        LogStream::Prefix prefix2("Error", saplog);
        saplog << "Calculate L2 error" << std::endl;

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

        error_file << discrete_time.get_step_number() << "; "
                   << discrete_time.get_current_time() << "; " << L2_norm
                   << "; " << L2_error << "; " << L2_error / L2_norm
                   << std::endl;
      }
      /** [Last time step] */


      /** [End simulation] */
      saplog << "Simulation ended at t = " << discrete_time.get_current_time()
             << " \t[" << Utilities::System::get_time() << "]" << std::endl;

      error_file.close();
      /** [End simulation] */
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
