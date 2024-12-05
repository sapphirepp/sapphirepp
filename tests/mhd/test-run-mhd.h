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
 * @file test-run-mhd.h
 * @author Florian Schulze (florian.schulze@mpi-hd.mpg.de)
 * @brief Implement functions for test runs.
 */

#ifndef TESTS_TESTRUNMHD_H
#define TESTS_TESTRUNMHD_H

#include <deal.II/base/discrete_time.h>
#include <deal.II/base/function.h>
#include <deal.II/base/mpi.h>
#include <deal.II/base/parameter_handler.h>

#include <deal.II/numerics/data_out.h>

#include <mpi.h>

#include <fstream>
#include <iostream>

#include "mhd-parameters.h"
#include "mhd-solver.h"
#include "output-parameters.h"
#include "sapphirepp-logstream.h"



template <unsigned int dim>
int
test_run_mhd(const sapphirepp::MHD::MHDParameters<dim> &mhd_parameters,
             const sapphirepp::PhysicalParameters      &physical_parameters,
             sapphirepp::Utils::OutputParameters       &output_parameters,
             dealii::Function<3>                       &exact_solution,
             const double                               max_L2_error = 1e-10)
{
  try
    {
      using namespace sapphirepp;
      using namespace MHD;

      saplog << "Start test run mhd" << std::endl;
      LogStream::Prefix  p("Test", saplog);
      const unsigned int spacedim = MHDEquations::spacedim;

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


      /** [Setup mhd_solver] */
      MHDSolver<dim> mhd_solver(mhd_parameters,
                                physical_parameters,
                                output_parameters);
      mhd_solver.setup();
      /** [Setup mhd_solver] */


      /** [Setup analytic solution] */
      AssertDimension(exact_solution.n_components, MHDEquations::n_components)
        PETScWrappers::MPI::Vector analytic_solution_vector;
      analytic_solution_vector.reinit(
        mhd_solver.get_dof_handler().locally_owned_dofs(), MPI_COMM_WORLD);
      /** [Setup analytic solution] */


      /** [Time loop] */
      DiscreteTime discrete_time(0,
                                 mhd_parameters.final_time,
                                 mhd_parameters.time_step);
      for (; discrete_time.is_at_end() == false; discrete_time.advance_time())
        {
          saplog << "Time step " << std::setw(6) << std::right
                 << discrete_time.get_step_number()
                 << " at t = " << discrete_time.get_current_time() << " \t["
                 << Utilities::System::get_time() << "]" << std::endl;

          exact_solution.set_time(discrete_time.get_current_time());
          /** [Time loop] */


          /** [Output solution] */
          if ((discrete_time.get_step_number() %
               output_parameters.output_frequency) == 0)
            {
              LogStream::Prefix prefix("Output", saplog);
              saplog << "Output solution" << std::endl;

              dealii::DataOut<dim, spacedim> data_out;
              data_out.attach_dof_handler(mhd_solver.get_dof_handler());

              // Output numeric solution
              data_out.add_data_vector(
                mhd_solver.get_current_solution(),
                MHDEquations::create_component_name_list("numeric_"),
                dealii::DataOut<dim, spacedim>::type_dof_data,
                MHDEquations::create_component_interpretation_list());

              // Output projected analytic solution
              mhd_solver.project(exact_solution, analytic_solution_vector);
              data_out.add_data_vector(
                analytic_solution_vector,
                MHDEquations::create_component_name_list("project_"),
                dealii::DataOut<dim, spacedim>::type_dof_data,
                MHDEquations::create_component_interpretation_list());

              // Output interpolated analytic solution
              dealii::VectorTools::interpolate(mhd_solver.get_dof_handler(),
                                               exact_solution,
                                               analytic_solution_vector);
              data_out.add_data_vector(
                analytic_solution_vector,
                MHDEquations::create_component_name_list("interpol_"),
                dealii::DataOut<dim, spacedim>::type_dof_data,
                MHDEquations::create_component_interpretation_list());

              // Output cell average and shock indicator
              data_out.add_data_vector(mhd_solver.get_cell_average_component(
                                         MHDEquations::density_component),
                                       "average_roh",
                                       DataOut<dim, spacedim>::type_cell_data);
              data_out.add_data_vector(mhd_solver.get_shock_indicator(),
                                       "shock_indicator",
                                       DataOut<dim, spacedim>::type_cell_data);
              data_out.add_data_vector(
                mhd_solver.get_positivity_limiter_indicator(),
                "positivity_limiter",
                DataOut<dim, spacedim>::type_cell_data);

              // Output the partition of the mesh
              const auto &triangulation = mhd_solver.get_triangulation();
              dealii::Vector<float> subdomain(triangulation.n_active_cells());
              for (unsigned int i = 0; i < subdomain.size(); ++i)
                subdomain(i) =
                  static_cast<float>(triangulation.locally_owned_subdomain());
              data_out.add_data_vector(subdomain, "subdomain");

              data_out.build_patches(mhd_parameters.polynomial_degree);
              output_parameters.write_results<dim, spacedim>(
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
              mhd_solver.compute_global_error(exact_solution,
                                              dealii::VectorTools::L2_norm,
                                              dealii::VectorTools::L2_norm);
            const double L2_norm =
              mhd_solver.compute_weighted_norm(dealii::VectorTools::L2_norm,
                                               dealii::VectorTools::L2_norm);

            saplog << "L2_error = " << L2_error << ", L2_norm = " << L2_norm
                   << ", rel error = " << L2_error / L2_norm << std::endl;

            error_file << discrete_time.get_step_number() << "; "
                       << discrete_time.get_current_time() << "; " << L2_norm
                       << "; " << L2_error << "; " << L2_error / L2_norm
                       << std::endl;
          }
          /** [Calculate error] */


          /** [Time step] */
          switch (mhd_parameters.time_stepping_method)
            {
              case TimeSteppingMethodMHD::forward_euler:
                mhd_solver.forward_euler_method(
                  discrete_time.get_current_time(),
                  discrete_time.get_next_step_size());
                break;
              case TimeSteppingMethodMHD::erk2:
              case TimeSteppingMethodMHD::erk4:
                mhd_solver.explicit_runge_kutta(
                  discrete_time.get_current_time(),
                  discrete_time.get_next_step_size());
                break;
              default:
                AssertThrow(false, ExcNotImplemented());
            }
        }
      /** [Time step] */


      /** [Last time step] */
      exact_solution.set_time(discrete_time.get_current_time());

      {
        LogStream::Prefix prefix("Output", saplog);
        saplog << "Output results" << std::endl;

        dealii::DataOut<dim, spacedim> data_out;
        data_out.attach_dof_handler(mhd_solver.get_dof_handler());

        // Output numeric solution
        data_out.add_data_vector(
          mhd_solver.get_current_solution(),
          MHDEquations::create_component_name_list("numeric_"),
          dealii::DataOut<dim, spacedim>::type_dof_data,
          MHDEquations::create_component_interpretation_list());

        // Output projected analytic solution
        mhd_solver.project(exact_solution, analytic_solution_vector);
        data_out.add_data_vector(
          analytic_solution_vector,
          MHDEquations::create_component_name_list("project_"),
          dealii::DataOut<dim, spacedim>::type_dof_data,
          MHDEquations::create_component_interpretation_list());

        // Output interpolated analytic solution
        dealii::VectorTools::interpolate(mhd_solver.get_dof_handler(),
                                         exact_solution,
                                         analytic_solution_vector);
        data_out.add_data_vector(
          analytic_solution_vector,
          MHDEquations::create_component_name_list("interpol_"),
          dealii::DataOut<dim, spacedim>::type_dof_data,
          MHDEquations::create_component_interpretation_list());

        // Output cell average and shock indicator
        data_out.add_data_vector(mhd_solver.get_cell_average_component(
                                   MHDEquations::density_component),
                                 "average_roh",
                                 DataOut<dim, spacedim>::type_cell_data);
        data_out.add_data_vector(mhd_solver.get_shock_indicator(),
                                 "shock_indicator",
                                 DataOut<dim, spacedim>::type_cell_data);
        data_out.add_data_vector(mhd_solver.get_positivity_limiter_indicator(),
                                 "positivity_limiter",
                                 DataOut<dim, spacedim>::type_cell_data);

        // Output the partition of the mesh
        const auto           &triangulation = mhd_solver.get_triangulation();
        dealii::Vector<float> subdomain(triangulation.n_active_cells());
        for (unsigned int i = 0; i < subdomain.size(); ++i)
          subdomain(i) =
            static_cast<float>(triangulation.locally_owned_subdomain());
        data_out.add_data_vector(subdomain, "subdomain");

        data_out.build_patches(mhd_parameters.polynomial_degree);
        output_parameters.write_results<dim, spacedim>(
          data_out,
          discrete_time.get_step_number(),
          discrete_time.get_current_time());
      }

      {
        LogStream::Prefix prefix2("Error", saplog);
        saplog << "Calculate L2 error" << std::endl;

        const double L2_error =
          mhd_solver.compute_global_error(exact_solution,
                                          dealii::VectorTools::L2_norm,
                                          dealii::VectorTools::L2_norm);
        const double L2_norm =
          mhd_solver.compute_weighted_norm(dealii::VectorTools::L2_norm,
                                           dealii::VectorTools::L2_norm);

        saplog << "L2_error = " << L2_error << ", L2_norm = " << L2_norm
               << ", rel error = " << L2_error / L2_norm << std::endl;

        error_file << discrete_time.get_step_number() << "; "
                   << discrete_time.get_current_time() << "; " << L2_norm
                   << "; " << L2_error << "; " << L2_error / L2_norm
                   << std::endl;

        AssertThrow((L2_error / L2_norm) < max_L2_error,
                    dealii::ExcMessage(
                      "L2 error is too large! (" +
                      dealii::Utilities::to_string(L2_error / L2_norm) + " > " +
                      dealii::Utilities::to_string(max_L2_error) + ")"));
      }
      /** [Last time step] */


      /** [End simulation] */
      error_file.close();

      saplog << "Simulation ended at t = " << discrete_time.get_current_time()
             << " \t[" << Utilities::System::get_time() << "]" << std::endl;

      mhd_solver.get_timer().print_wall_time_statistics(MPI_COMM_WORLD);
      /** [End simulation] */
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

  sapphirepp::saplog << "Succeeded test run mhd." << std::endl;
  return 0;
}

#endif
