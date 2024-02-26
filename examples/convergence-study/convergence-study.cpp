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
      using namespace sapphirepp;
      using namespace VFP;
      dealii::Utilities::MPI::MPI_InitFinalize mpi_initialization(argc,
                                                                  argv,
                                                                  1);

      saplog.init(3);

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

      AssertThrow(vfp_parameters.expansion_order == 1,
                  dealii::ExcMessage("Expansion order be one."));

      // Copy vfp_parameters to physical parameters for inital function
      physical_parameters.velocity = vfp_parameters.velocity;
      physical_parameters.gamma    = vfp_parameters.gamma;
      physical_parameters.charge   = vfp_parameters.charge;
      physical_parameters.mass     = vfp_parameters.mass;
      physical_parameters.box_length =
        std::abs(vfp_parameters.p1[0] - vfp_parameters.p2[0]);

      const double gyroperiod =
        vfp_parameters.gamma * vfp_parameters.mass /
        (physical_parameters.B0_2pi * vfp_parameters.charge);
      const double gyroradius =
        vfp_parameters.gamma * vfp_parameters.mass * vfp_parameters.velocity /
        (vfp_parameters.charge * physical_parameters.B0);
      const double box_length =
        std::abs(vfp_parameters.p1[0] - vfp_parameters.p2[0]);
      saplog << "particle_velocity = " << vfp_parameters.velocity
             << ", gyroperiod = " << gyroperiod
             << ", gyroradius = " << gyroradius
             << ", box_length = " << box_length
             << ", final_time = " << vfp_parameters.final_time << std::endl;

      if (std::fmod(vfp_parameters.final_time, gyroperiod) != 0.)
        saplog << "Warning: Final time is not a multiple of the gyroperiod."
               << std::endl;

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

      VFPSolver<dimension> vfp_solver(vfp_parameters,
                                      physical_parameters,
                                      output_parameters);
      vfp_solver.setup();

      const unsigned int num_exp_coefficients =
        (vfp_parameters.expansion_order + 1) *
        (vfp_parameters.expansion_order + 1);

      InitialValueFunction<dimension> analytic_solution(
        physical_parameters, vfp_parameters.expansion_order);
      const dealii::ComponentSelectFunction<dimension> weight(
        0, num_exp_coefficients);

      PETScWrappers::MPI::Vector analytic_solution_vector;
      analytic_solution_vector.reinit(
        vfp_solver.get_dof_handler().locally_owned_dofs(), MPI_COMM_WORLD);

      std::vector<std::string> component_names(num_exp_coefficients);
      std::vector<std::string> component_names_interpol(num_exp_coefficients);
      const std::vector<std::array<unsigned int, 3>> lms_indices =
        PDESystem::create_lms_indices(num_exp_coefficients);
      for (unsigned int i = 0; i < num_exp_coefficients; ++i)
        {
          const std::array<unsigned int, 3> &lms = lms_indices[i];
          component_names[i] = "analytic_f_" + std::to_string(lms[0]) +
                               std::to_string(lms[1]) + std::to_string(lms[2]);
          component_names_interpol[i] =
            "analytic_interpol_f_" + std::to_string(lms[0]) +
            std::to_string(lms[1]) + std::to_string(lms[2]);
        }

      saplog << "Calculate analytic solution" << std::endl;

      dealii::DataOut<dimension> data_out;
      data_out.attach_dof_handler(vfp_solver.get_dof_handler());
      analytic_solution.set_time(0.);


      vfp_solver.project(analytic_solution, analytic_solution_vector);
      data_out.add_data_vector(analytic_solution_vector, component_names);

      dealii::VectorTools::interpolate(vfp_solver.get_dof_handler(),
                                       analytic_solution,
                                       analytic_solution_vector);
      data_out.add_data_vector(analytic_solution_vector,
                               component_names_interpol);

      data_out.build_patches(vfp_parameters.polynomial_degree);
      output_parameters.write_results<dimension>(data_out,
                                                 0,
                                                 0.,
                                                 "analytic_solution");

      saplog << "Compute error" << std::endl;
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

      error_file << 0 << "; " << 0.0 << "; " << L2_norm << "; " << L2_error
                 << "; " << L2_error / L2_norm << std::endl;

      DiscreteTime discrete_time(0,
                                 vfp_parameters.final_time,
                                 vfp_parameters.time_step);
      for (; discrete_time.is_at_end() == false; discrete_time.advance_time())
        {
          {
            LogStream::Prefix prefix("VFP", saplog);
            saplog << "Time step " << std::setw(6) << std::right
                   << discrete_time.get_step_number()
                   << " at t = " << discrete_time.get_current_time() << " \t["
                   << Utilities::System::get_time() << "]" << std::endl;

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

            vfp_solver.output_results(discrete_time.get_step_number() + 1,
                                      discrete_time.get_next_time());
          }

          LogStream::Prefix prefix2("analytic", saplog);
          saplog << "Calculate analytic solution" << std::endl;
          analytic_solution.set_time(discrete_time.get_next_time());

          if ((discrete_time.get_step_number() + 1) %
                output_parameters.output_frequency ==
              0)
            {
              dealii::DataOut<dimension> data_out;
              data_out.attach_dof_handler(vfp_solver.get_dof_handler());
              vfp_solver.project(analytic_solution, analytic_solution_vector);
              data_out.add_data_vector(analytic_solution_vector,
                                       component_names);

              dealii::VectorTools::interpolate(vfp_solver.get_dof_handler(),
                                               analytic_solution,
                                               analytic_solution_vector);
              data_out.add_data_vector(analytic_solution_vector,
                                       component_names_interpol);

              data_out.build_patches(vfp_parameters.polynomial_degree);
              output_parameters.write_results<dimension>(
                data_out,
                discrete_time.get_step_number() + 1,
                discrete_time.get_next_time(),
                "analytic_solution");
            }

          saplog << "Compute error" << std::endl;
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

          error_file << discrete_time.get_step_number() + 1 << "; "
                     << discrete_time.get_next_time() << "; " << L2_norm << "; "
                     << L2_error << "; " << L2_error / L2_norm << std::endl;
        }
      saplog << "Simulation ended at t = " << discrete_time.get_current_time()
             << " \t[" << Utilities::System::get_time() << "]" << std::endl;

      error_file.close();
      vfp_solver.get_timer().print_wall_time_statistics(MPI_COMM_WORLD);
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
