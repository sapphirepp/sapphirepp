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
 * @file examples/mhd/linear-waves/linear-waves.cpp
 * @author Florian Schulze (florian.schulze@mpi-hd.mpg.de)
 * @brief Implement main function for linear-waves example
 */

#include <deal.II/base/discrete_time.h>
#include <deal.II/base/mpi.h>
#include <deal.II/base/parameter_handler.h>

#include <deal.II/numerics/data_out.h>

#include <mpi.h>

#include <fstream>
#include <iostream>

#include "config.h"
#include "mhd-parameters.h"
#include "mhd-solver.h"
#include "output-parameters.h"
#include "sapphirepp-logstream.h"



int
main(int argc, char *argv[])
{
  try
    {
      /** [Main function setup] */
      using namespace sapphirepp;
      using namespace MHD;
      dealii::Utilities::MPI::MPI_InitFinalize mpi_initialization(argc,
                                                                  argv,
                                                                  1);

      // saplog.init(2);
      saplog.init(100);

      std::string parameter_filename = "parameter.prm";
      if (argc > 1)
        parameter_filename = argv[1];

      dealii::ParameterHandler prm;
      PhysicalParameters       physical_parameters(dim_mhd);
      Utils::OutputParameters  output_parameters;
      MHDParameters<dim_mhd>   mhd_parameters;

      physical_parameters.declare_parameters(prm);
      output_parameters.declare_parameters(prm);
      mhd_parameters.declare_parameters(prm);

      prm.parse_input(parameter_filename);

      physical_parameters.parse_parameters(prm);
      output_parameters.parse_parameters(prm);
      mhd_parameters.parse_parameters(prm);
      /** [Main function setup] */


      /** [Copy MHD parameter] */
      const unsigned int spacedim    = MHDEquations::spacedim;
      const unsigned int dimension   = dim_mhd;
      physical_parameters.box_length = std::vector<double>(dimension, 1.);
      for (unsigned int d = 0; d < dimension; ++d)
        {
          physical_parameters.box_length[d] =
            std::abs(mhd_parameters.p1[d] - mhd_parameters.p2[d]);

          AssertThrow((mhd_parameters.boundary_conditions[2 * d + 0] ==
                       MHD::BoundaryConditionsMHD::periodic) &&
                        (mhd_parameters.boundary_conditions[2 * d + 1] ==
                         MHD::BoundaryConditionsMHD::periodic),
                      dealii::ExcMessage("This example assumes periodic BC."));
        }
      /** [Copy MHD parameter] */


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
      MHDSolver<dim_mhd> mhd_solver(mhd_parameters,
                                    physical_parameters,
                                    output_parameters);
      mhd_solver.setup();
      /** [Setup mhd_solver] */


      /** [Setup analytic solution] */
      InitialConditionMHD<spacedim> analytic_solution(
        physical_parameters, mhd_parameters.adiabatic_index);

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

          analytic_solution.set_time(discrete_time.get_current_time());
          /** [Time loop] */


          /** [Output solution] */
          if ((discrete_time.get_step_number() %
               output_parameters.output_frequency) == 0)
            {
              LogStream::Prefix prefix("Output", saplog);
              saplog << "Output solution" << std::endl;

              dealii::DataOut<dim_mhd, spacedim> data_out;
              data_out.attach_dof_handler(mhd_solver.get_dof_handler());

              // Output numeric solution
              data_out.add_data_vector(
                mhd_solver.get_current_solution(),
                MHDEquations::create_component_name_list("numeric_"),
                dealii::DataOut<dim_mhd, spacedim>::type_dof_data,
                MHDEquations::create_component_interpretation_list());

              // Output projected analytic solution
              mhd_solver.project(analytic_solution, analytic_solution_vector);
              data_out.add_data_vector(
                analytic_solution_vector,
                MHDEquations::create_component_name_list("project_"),
                dealii::DataOut<dim_mhd, spacedim>::type_dof_data,
                MHDEquations::create_component_interpretation_list());

              // Output interpolated analytic solution
              dealii::VectorTools::interpolate(mhd_solver.get_dof_handler(),
                                               analytic_solution,
                                               analytic_solution_vector);
              data_out.add_data_vector(
                analytic_solution_vector,
                MHDEquations::create_component_name_list("interpol_"),
                dealii::DataOut<dim_mhd, spacedim>::type_dof_data,
                MHDEquations::create_component_interpretation_list());

              data_out.build_patches(mhd_parameters.polynomial_degree);
              output_parameters.write_results<dim_mhd, spacedim>(
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
              mhd_solver.compute_global_error(analytic_solution,
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
      analytic_solution.set_time(discrete_time.get_current_time());

      {
        LogStream::Prefix prefix("Output", saplog);
        saplog << "Output results" << std::endl;

        dealii::DataOut<dim_mhd, spacedim> data_out;
        data_out.attach_dof_handler(mhd_solver.get_dof_handler());

        // Output numeric solution
        data_out.add_data_vector(
          mhd_solver.get_current_solution(),
          MHDEquations::create_component_name_list("numeric_"),
          dealii::DataOut<dim_mhd, spacedim>::type_dof_data,
          MHDEquations::create_component_interpretation_list());

        // Output projected analytic solution
        mhd_solver.project(analytic_solution, analytic_solution_vector);
        data_out.add_data_vector(
          analytic_solution_vector,
          MHDEquations::create_component_name_list("project_"),
          dealii::DataOut<dim_mhd, spacedim>::type_dof_data,
          MHDEquations::create_component_interpretation_list());

        // Output interpolated analytic solution
        dealii::VectorTools::interpolate(mhd_solver.get_dof_handler(),
                                         analytic_solution,
                                         analytic_solution_vector);
        data_out.add_data_vector(
          analytic_solution_vector,
          MHDEquations::create_component_name_list("interpol_"),
          dealii::DataOut<dim_mhd, spacedim>::type_dof_data,
          MHDEquations::create_component_interpretation_list());

        data_out.build_patches(mhd_parameters.polynomial_degree);
        output_parameters.write_results<dim_mhd, spacedim>(
          data_out,
          discrete_time.get_step_number(),
          discrete_time.get_current_time());
      }

      {
        LogStream::Prefix prefix2("Error", saplog);
        saplog << "Calculate L2 error" << std::endl;

        const double L2_error =
          mhd_solver.compute_global_error(analytic_solution,
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
      /** [Last time step] */


      /** [End simulation] */
      saplog << "Simulation ended at t = " << discrete_time.get_current_time()
             << " \t[" << Utilities::System::get_time() << "]" << std::endl;

      error_file.close();
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
  return 0;
}
