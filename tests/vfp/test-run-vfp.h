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
 * @file tests/vfp/test-run-vfp.h
 * @author Florian Schulze (florian.schulze@mpi-hd.mpg.de)
 * @brief Implement functions for vfp test runs.
 */

#ifndef TESTS_TESTRUNVFP_H
#define TESTS_TESTRUNVFP_H

#include <deal.II/base/discrete_time.h>
#include <deal.II/base/function.h>
#include <deal.II/base/mpi.h>
#include <deal.II/base/parameter_handler.h>

#include <deal.II/numerics/data_out.h>

#include <mpi.h>

#include <fstream>
#include <iostream>

#include "output-parameters.h"
#include "sapphirepp-logstream.h"
#include "vfp-parameters.h"
#include "vfp-solver.h"



/**
 * @brief Output solution and exact solution for the VFP test run.
 *
 * @tparam dim Dimension of the reduced phase space \f$ (\mathbf{x}, p) \f$,
 *         `dim_ps`
 * @param vfp_solver @ref sapphirepp::VFP::VFPSolver
 * @param vfp_parameters Parameters for the VFP equation
 * @param output_parameters Parameters for the output
 * @param exact_solution Exact solution to compare to
 * @param time_step_number Time step number
 * @param cur_time Simulation time
 * @param filename Overwrite for the @ref base_file_name
 */
template <unsigned int dim>
void
test_run_vfp_output(const sapphirepp::VFP::VFPSolver<dim>     &vfp_solver,
                    const sapphirepp::VFP::VFPParameters<dim> &vfp_parameters,
                    sapphirepp::Utils::OutputParameters &output_parameters,
                    const dealii::Function<dim>         &exact_solution,
                    const unsigned int                   time_step_number,
                    const double                         cur_time,
                    const std::string                   &filename = "")
{
  using namespace sapphirepp;
  using namespace VFP;

  // dealii::TimerOutput::Scope timer_section(vfp_solver.get_timer(),
  //                                          "VFP - Output");
  dealii::LogStream::Prefix prefix("Output", saplog);
  saplog << "Output results at t = " << cur_time << std::endl;
  dealii::DataOut<dim> data_out;
  data_out.attach_dof_handler(vfp_solver.get_dof_handler());

  AssertDimension(exact_solution.n_components,
                  vfp_solver.get_pde_system().system_size);
  dealii::PETScWrappers::MPI::Vector exact_solution_vector;
  exact_solution_vector.reinit(
    vfp_solver.get_dof_handler().locally_owned_dofs(), MPI_COMM_WORLD);

  const std::string function_symbol =
    ((vfp_flags & VFPFlags::scaled_distribution_function) != VFPFlags::none) ?
      "g" :
      "f";

  // Output numeric solution
  data_out.add_data_vector(vfp_solver.get_current_solution(),
                           PDESystem::create_component_name_list(
                             vfp_solver.get_pde_system().system_size,
                             "numeric_" + function_symbol + "_"));

  // Output projected exact solution (skip in steady-state)
  if ((vfp_flags & VFPFlags::time_evolution) != VFPFlags::none)
    {
      vfp_solver.project(exact_solution, exact_solution_vector);
      data_out.add_data_vector(exact_solution_vector,
                               PDESystem::create_component_name_list(
                                 vfp_solver.get_pde_system().system_size,
                                 "project_" + function_symbol + "_"));
    }

  // Output interpolated exact solution
  dealii::VectorTools::interpolate(vfp_solver.get_dof_handler(),
                                   exact_solution,
                                   exact_solution_vector);
  data_out.add_data_vector(exact_solution_vector,
                           PDESystem::create_component_name_list(
                             vfp_solver.get_pde_system().system_size,
                             "interpol_" + function_symbol + "_"));

  // Output the partition of the mesh
  const auto           &triangulation = vfp_solver.get_triangulation();
  dealii::Vector<float> subdomain(triangulation.n_active_cells());
  for (unsigned int i = 0; i < subdomain.size(); ++i)
    subdomain(i) = static_cast<float>(triangulation.locally_owned_subdomain());
  data_out.add_data_vector(subdomain, "subdomain");

  // if (output_parameters.debug_input_functions)
  //   {
  //     // Add input functions
  //     dealii::PETScWrappers::MPI::Vector temp(
  //       debug_input_functions_dof_handler.locally_owned_dofs(),
  //       mpi_communicator);
  //     dealii::VectorTools::interpolate(debug_input_functions_dof_handler,
  //                                      debug_input_functions,
  //                                      temp);
  //     data_out.add_data_vector(
  //       debug_input_functions_dof_handler,
  //       temp,
  //       debug_input_functions.create_component_name_list());
  //   }

  data_out.build_patches(vfp_parameters.polynomial_degree);
  output_parameters.write_results<dim>(data_out,
                                       time_step_number,
                                       cur_time,
                                       filename);

  // probe_location.probe_all_points(
  //   dof_handler,
  //   mapping,
  //   locally_relevant_current_solution,
  //   time_step_number,
  //   cur_time,
  //   function_symbol);
}



/**
 * @brief Calculate the L2 error for the VFP test run.
 *
 * @tparam dim Dimension of the reduced phase space \f$ (\mathbf{x}, p) \f$,
 *         `dim_ps`
 * @param vfp_solver @ref sapphirepp::VFP::VFPSolver
 * @param exact_solution Exact solution to compare to
 * @param error_file Output stream to save the error
 * @param max_L2_error Maximum expected L2 error.
 *        Do not compare for `max_L2_error = 0`.
 * @param weight Weights for the error
 * @param time_step_number Time step number
 * @param cur_time Simulation time
 */
template <unsigned int dim>
void
test_run_vfp_error(const sapphirepp::VFP::VFPSolver<dim> &vfp_solver,
                   dealii::Function<dim>                 &exact_solution,
                   std::ostream                          &error_file,
                   const double                           max_L2_error = 0.,
                   const dealii::Function<dim, double>   *weight = nullptr,
                   const unsigned int                     time_step_number = 0,
                   const double                           cur_time         = 0.)
{
  using namespace sapphirepp;

  // dealii::TimerOutput::Scope timer_section(vfp_solver.get_timer(),
  //                                          "VFP - Error");
  dealii::LogStream::Prefix prefix_error("Error", saplog);
  saplog << "Calculate L2 error" << std::endl;

  const double L2_error =
    vfp_solver.compute_global_error(exact_solution,
                                    dealii::VectorTools::L2_norm,
                                    dealii::VectorTools::L2_norm,
                                    weight);
  const double L2_norm =
    vfp_solver.compute_weighted_norm(dealii::VectorTools::L2_norm,
                                     dealii::VectorTools::L2_norm,
                                     weight);

  saplog << "L2_error = " << L2_error << ", L2_norm = " << L2_norm
         << ", rel error = " << L2_error / L2_norm << std::endl;

  error_file << time_step_number << "; " //
             << cur_time << "; "         //
             << L2_norm << "; "          //
             << L2_error << "; "         //
             << L2_error / L2_norm << std::endl;


  if (max_L2_error > 0.)
    AssertThrow((L2_error / L2_norm) < max_L2_error,
                dealii::ExcMessage(
                  "L2 error is too large! (" +
                  dealii::Utilities::to_string(L2_error / L2_norm) + " > " +
                  dealii::Utilities::to_string(max_L2_error) + ")"));
}



/**
 * @brief Test run for the VFP module.
 *
 * @tparam dim Dimension of the reduced phase space \f$ (\mathbf{x}, p) \f$,
 *         `dim_ps`
 * @param vfp_parameters Parameters for the VFP equation
 * @param physical_parameters User defined parameters of the problem
 * @param output_parameters Parameters for the output
 * @param exact_solution Exact solution to compare to
 * @param max_L2_error Maximum expected L2 error.
 *        Do not compare for `max_L2_error = 0`.
 * @param weight Weights for the error
 * @param compare_each_time_step Compare against `max_L2_error` each time step?
 * @return Returns `0` on success,
 *         or `1` if an error or exception occurred during the test run.
 */
template <unsigned int dim>
int
test_run_vfp(const sapphirepp::VFP::VFPParameters<dim> &vfp_parameters,
             const sapphirepp::PhysicalParameters      &physical_parameters,
             sapphirepp::Utils::OutputParameters       &output_parameters,
             dealii::Function<dim>                     &exact_solution,
             const double                               max_L2_error = 0.,
             const dealii::Function<dim, double>       *weight       = nullptr,
             const bool compare_each_time_step                       = true)
{
  sapphirepp::saplog << "Start test run VFP" << std::endl;
  try
    {
      using namespace sapphirepp;
      using namespace VFP;
      dealii::LogStream::Prefix prefix_test("VFPTest", saplog);

      if constexpr ((vfp_flags & VFPFlags::time_evolution) == VFPFlags::none)
        AssertThrow(false, dealii::ExcNotImplemented());

      /** [Create error file] */
      std::ofstream error_file(output_parameters.output_path / "error.csv");
      AssertThrow(!error_file.fail(),
                  dealii::ExcFileNotOpen(output_parameters.output_path /
                                         "error.csv"));
      error_file << "timestep" << "; " //
                 << "time" << "; "     //
                 << "L2_norm" << "; "  //
                 << "L2_error" << "; " //
                 << "relative_error" << std::endl;
      /** [Create error file] */


      /** [Setup vfp_solver] */
      VFPSolver<dim> vfp_solver(vfp_parameters,
                                physical_parameters,
                                output_parameters);
      vfp_solver.setup();
      /** [Setup vfp_solver] */


      /** [Time loop] */
      dealii::DiscreteTime discrete_time(0,
                                         vfp_parameters.final_time,
                                         vfp_parameters.time_step);
      for (; discrete_time.is_at_end() == false; discrete_time.advance_time())
        {
          saplog << "Time step " << std::setw(6) << std::right
                 << discrete_time.get_step_number()
                 << " at t = " << discrete_time.get_current_time() << " \t["
                 << dealii::Utilities::System::get_time() << "]" << std::endl;

          exact_solution.set_time(discrete_time.get_current_time());
          /** [Time loop] */


          /** [Output solution] */
          if ((discrete_time.get_step_number() %
               output_parameters.output_frequency) == 0)
            {
              test_run_vfp_output<dim>(vfp_solver,
                                       vfp_parameters,
                                       output_parameters,
                                       exact_solution,
                                       discrete_time.get_step_number(),
                                       discrete_time.get_current_time());
              /** [Output solution] */


              /** [Calculate error] */
              test_run_vfp_error<dim>(vfp_solver,
                                      exact_solution,
                                      error_file,
                                      compare_each_time_step ? max_L2_error :
                                                               0.,
                                      weight,
                                      discrete_time.get_step_number(),
                                      discrete_time.get_current_time());
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
              case TimeSteppingMethod::lserk4:
                vfp_solver.low_storage_explicit_runge_kutta(
                  discrete_time.get_current_time(),
                  discrete_time.get_next_step_size());
                break;
              default:
                AssertThrow(false, dealii::ExcNotImplemented());
            }
        }
      /** [Time step] */


      /** [Last time step] */
      exact_solution.set_time(discrete_time.get_current_time());

      test_run_vfp_output<dim>(vfp_solver,
                               vfp_parameters,
                               output_parameters,
                               exact_solution,
                               discrete_time.get_step_number(),
                               discrete_time.get_current_time());

      test_run_vfp_error<dim>(vfp_solver,
                              exact_solution,
                              error_file,
                              max_L2_error,
                              weight,
                              discrete_time.get_step_number(),
                              discrete_time.get_current_time());
      /** [Last time step] */


      /** [End simulation] */
      error_file.close();

      saplog << "Simulation ended at t = " << discrete_time.get_current_time()
             << " \t[" << dealii::Utilities::System::get_time() << "]"
             << std::endl;

      {
        dealii::LogStream::Prefix              prefix("Summary", saplog);
        dealii::Utilities::System::MemoryStats memory_stats;
        dealii::Utilities::System::get_memory_stats(memory_stats);
        saplog << "Peak (local) resident memory size (HWM):    \t" //
               << memory_stats.VmHWM << " KiB"                     //
               << " \t= " << (memory_stats.VmHWM >> 10) << " MiB " //
               << " \t= " << (memory_stats.VmHWM >> 20) << " GiB " //
               << std::endl;
        saplog << "Peak (local) virtual/available memory size: \t"  //
               << memory_stats.VmPeak << " KiB"                     //
               << " \t= " << (memory_stats.VmPeak >> 10) << " MiB " //
               << " \t= " << (memory_stats.VmPeak >> 20) << " GiB " //
               << std::endl;

        vfp_solver.get_timer().print_wall_time_statistics(MPI_COMM_WORLD);
      }
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

  sapphirepp::saplog << "Succeeded test run VFP." << std::endl;
  return 0;
}

#endif
