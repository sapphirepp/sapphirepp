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
///
#include <deal.II/base/mpi.h>

#include <mpi.h>

#include <fstream>

#include "output-module.h"
#include "vfp-equation-solver.h"

enum class TestParameter
{
  final_time, //< Test if gyroperiod is correct
  expansion_order,
  time_step,
  num_cells,
  polynomial_degree,
  mpi_processes
};

const unsigned int dim = 2;

int
error_with_parameter(std::string               parameter_filename,
                     const TestParameter       test_parameter,
                     const std::vector<double> values,
                     const std::vector<double> max_expected_errors)
{
  using namespace Sapphire;
  using namespace VFP;

  if (max_expected_errors.size() > 0)
    AssertThrow(max_expected_errors.size() == values.size(),
                ExcMessage(
                  "Expected errors and values have different lengths"));

  std::string parameter_name;

  switch (test_parameter)
    {
      case TestParameter::final_time:
        {
          parameter_name = "final_time";
          break;
        }
      case TestParameter::expansion_order:
        {
          parameter_name = "expansion_order";
          break;
        }
      case TestParameter::time_step:
        {
          parameter_name = "time_step";
          break;
        }
      case TestParameter::num_cells:
        {
          parameter_name = "num_cells";
          break;
        }
      case TestParameter::polynomial_degree:
        {
          parameter_name = "polynomial_degree";
          break;
        }
      case TestParameter::mpi_processes:
        {
          parameter_name = "mpi_processes";
          break;
        }
      default:
        Assert(false, ExcNotImplemented());
    }

  saplog.push("Gyro");
  saplog << "Test for " << parameter_name << std::endl;

  Timer                    timer;
  ParameterHandler         prm;
  VFPSolverControl<dim>    vfp_solver_control;
  PhysicalProperties       physical_properties;
  Utils::OutputModule<dim> output_module;

  vfp_solver_control.declare_parameters(prm);
  physical_properties.declare_parameters(prm);
  output_module.declare_parameters(prm);

  prm.parse_input(parameter_filename);

  vfp_solver_control.parse_parameters(prm);
  physical_properties.parse_parameters(prm);
  output_module.parse_parameters(prm);

  saplog.pop();

  std::ofstream log_file(output_module.output_path /
                           ("error_with_" + parameter_name + ".csv"),
                         std::ios::app);
  saplog.attach(log_file, false);

  saplog << "# " << parameter_name
         << "\t L2\t Linfty\t CPU time [s]\t Wall time [s]\t n_dof"
         << std::endl;
  for (unsigned int i = 0; i < values.size(); i++)
    {
      switch (test_parameter)
        {
          case TestParameter::final_time:
            vfp_solver_control.final_time = values[i];
            break;
          case TestParameter::expansion_order:
            vfp_solver_control.expansion_order = uint(values[i]);
            break;
          case TestParameter::time_step:
            vfp_solver_control.time_step = values[i];
            break;
          case TestParameter::num_cells:
            if (dim == 1)
              vfp_solver_control.n_cells = {uint(values[i])};
            else if (dim == 2)
              vfp_solver_control.n_cells = {uint(values[i]), uint(values[i])};
            else if (dim == 3)
              vfp_solver_control.n_cells = {uint(values[i]),
                                            uint(values[i]),
                                            uint(values[i])};
            break;
          case TestParameter::polynomial_degree:
            vfp_solver_control.polynomial_degree = uint(values[i]);
            break;
          case TestParameter::mpi_processes:
            if (values[i] == Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD))
              break;
            else
              continue;
          default:
            Assert(false, ExcNotImplemented());
        }
      saplog.push("Gyro");
      saplog << parameter_name << "=" << values[i] << std::endl;
      saplog.pop();

      output_module.base_file_name =
        "solution_" + parameter_name + "_" + std::to_string(values[i]) + "_";

      timer.start();
      VFPEquationSolver<dim> vfp_equation_solver(vfp_solver_control,
                                                 physical_properties,
                                                 output_module);
      vfp_equation_solver.run();
      timer.stop();

      InitialValueFunction<dim> cylinder(physical_properties,
                                         vfp_solver_control.expansion_order);


      ComponentSelectFunction<dim> mask(
        0,
        (vfp_solver_control.expansion_order + 1) *
          (vfp_solver_control.expansion_order + 1));

      const double L2_error = vfp_equation_solver.compute_global_error(
        cylinder, VectorTools::L2_norm, VectorTools::L2_norm, &mask);
      const double Linfty_error = vfp_equation_solver.compute_global_error(
        cylinder, VectorTools::L2_norm, VectorTools::Linfty_norm, &mask);

      saplog << values[i] << "\t" << L2_error << "\t" << Linfty_error << "\t"
             << timer.cpu_time() << "\t" << timer.wall_time() << "\t"
             << vfp_equation_solver.get_n_dofs() << std::endl;

      if (max_expected_errors.size() > 0)
        {
          if (test_parameter == TestParameter::mpi_processes)
            {
              AssertThrow(
                timer.cpu_time() < max_expected_errors[i],
                ExcMessage("Runtime longer than expected for " +
                           std::to_string(int(values[i])) + " mpi processes: " +
                           std::to_string(timer.cpu_time()) + "s > " +
                           std::to_string(max_expected_errors[i]) + "s"));
            }
          else
            {
              AssertThrow(L2_error < max_expected_errors[i],
                          ExcMessage(
                            "L2 error does not match expected value for " +
                            parameter_name + "=" + std::to_string(values[i]) +
                            ": " + std::to_string(L2_error) + " > " +
                            std::to_string(max_expected_errors[i])));
            }
        }
    }

  saplog.detach();
  log_file.close();

  return 0;
}

int
calculate_gyro(std::string parameter_filename, double max_expected_error = 0)
{
  using namespace Sapphire;
  using namespace VFP;

  saplog.push("Gyro");
  Timer                  timer;
  ParameterHandler       prm;
  VFPSolverControl<dim>  vfp_solver_control;
  PhysicalProperties     physical_properties;
  Utils::OutputModule<2> output_module;

  vfp_solver_control.declare_parameters(prm);
  physical_properties.declare_parameters(prm);
  output_module.declare_parameters(prm);

  prm.parse_input(parameter_filename);

  vfp_solver_control.parse_parameters(prm);
  physical_properties.parse_parameters(prm);
  output_module.parse_parameters(prm);

  TransportOnly particle_properties;
  saplog << "Gyroperiod: "
         << 2 * M_PI * particle_properties.gamma / physical_properties.B0
         << std::endl;
  saplog << "Gyroradius: " << particle_properties.gamma / physical_properties.B0
         << std::endl;

  saplog.pop();

  timer.start();
  VFPEquationSolver<dim> vfp_equation_solver(vfp_solver_control,
                                             physical_properties,
                                             output_module);
  vfp_equation_solver.run();
  timer.stop();

  InitialValueFunction<dim> cylinder(physical_properties,
                                     vfp_solver_control.expansion_order);
  double                    error =
    vfp_equation_solver.compute_global_error(cylinder,
                                             VectorTools::L2_norm,
                                             //  VectorTools::Linfty_norm);
                                             VectorTools::L2_norm);

  saplog << "Error = " << error << " CPU/wall time = " << timer.cpu_time()
         << "/" << timer.wall_time() << " s" << std::endl;

  if (max_expected_error > 0)
    AssertThrow(error < max_expected_error,
                ExcMessage("Error does not match expected value: " +
                           std::to_string(error) + " > " +
                           std::to_string(max_expected_error)));

  return 0;
}


int
main(int argc, char *argv[])
{
  try
    {
      using namespace Sapphire;
      using namespace dealii;
      Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

      saplog.pop();
      saplog.depth_console(1);
      saplog.depth_file(0);
      // saplog.depth_console(-1);
      const unsigned int mpi_rank =
        Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);
      if (mpi_rank > 0)
        {
          saplog.depth_console(0);
          saplog.depth_file(0);
          char buffer[4];
          std::snprintf(buffer, 4, "%03d", mpi_rank);
          saplog.push("mpi" + std::string(buffer));
        }

      std::string parameter_filename = "parameter.prm";
      if (argc > 1)
        parameter_filename = argv[1];

      std::string parameter_name = "gyro";
      if (argc > 2)
        parameter_name = argv[2];

      int mpi_size = Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD);
      saplog << "Start gyroradius test with parameter file \""
             << parameter_filename << "\" on " << mpi_size << " processor(s) ["
             << Utilities::System::get_date() << " "
             << Utilities::System::get_time() << "]" << std::endl;

      if (parameter_name == "gyro")
        {
          double max_expected_error = 0;
          if (argc > 3)
            max_expected_error = std::stod(argv[3]);

          return calculate_gyro(parameter_filename, max_expected_error);
        }

      TestParameter test_parameter;
      if (parameter_name == "expansion_order")
        test_parameter = TestParameter::expansion_order;
      else if (parameter_name == "final_time")
        test_parameter = TestParameter::final_time;
      else if (parameter_name == "time_step")
        test_parameter = TestParameter::time_step;
      else if (parameter_name == "num_cells")
        test_parameter = TestParameter::num_cells;
      else if (parameter_name == "polynomial_degree")
        test_parameter = TestParameter::polynomial_degree;
      else if (parameter_name == "mpi_processes")
        test_parameter = TestParameter::mpi_processes;
      else
        AssertThrow(false,
                    ExcMessage("Unknown parameter name: " + parameter_name));

      saplog << "Test parameter is " << parameter_name << std::endl;

      std::vector<double> values;
      std::vector<double> max_expected_errors;
      if (argc > 3)
        {
          const std::string expected_output_file = argv[3];
          saplog << "Read expected output from \"" << expected_output_file
                 << "\"" << std::endl;
          std::ifstream expected_output(expected_output_file);
          AssertThrow(expected_output.is_open(),
                      ExcMessage("Could not open expected output file"));
          std::string line;
          double      value, max_expected_error;
          std::getline(expected_output, line);
          while (std::getline(expected_output, line))
            {
              std::istringstream iss(line);
              iss >> value;
              values.push_back(value);
              iss >> max_expected_error;
              max_expected_errors.push_back(max_expected_error);
            }
          expected_output.close();
        }
      else
        {
          switch (test_parameter)
            {
              case TestParameter::final_time:
                values = {0.5, 1, 1.5, 2, 3, 4};
                break;
              case TestParameter::expansion_order:
                values = {1, 2, 3, 4, 5};
                break;
              case TestParameter::time_step:
                values = {0.01, 0.05, 0.1, 0.2, 0.5};
                break;
              case TestParameter::num_cells:
                values = {8, 16, 32, 64};
                break;
              case TestParameter::polynomial_degree:
                values = {1, 2};
                break;
              case TestParameter::mpi_processes:
                values = {double(mpi_size)};
                break;
              default:
                Assert(false, ExcNotImplemented());
            }
        }

      return error_with_parameter(parameter_filename,
                                  test_parameter,
                                  values,
                                  max_expected_errors);
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
