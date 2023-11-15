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

#include <deal.II/base/mpi.h>

#include <mpi.h>

#include "config.h"
#include "output-parameters.h"
#include "vfp-equation-solver.h"


const unsigned int dim = 1;

int
convergence_with_expansion_order(const std::string         &parameter_filename,
                                 const std::vector<double> &values)
{
  using namespace Sapphire;
  using namespace VFP;

  saplog.push("Tests");
  saplog << "Compute convergence with expansion_order" << std::endl;

  Timer                        timer;
  ParameterHandler             prm;
  VFPParameters<dim>           vfp_parameters;
  PhysicalProperties           physical_properties;
  Utils::OutputParameters<dim> output_parameters;

  vfp_parameters.declare_parameters(prm);
  physical_properties.declare_parameters(prm);
  output_parameters.declare_parameters(prm);

  prm.parse_input(parameter_filename);

  vfp_parameters.parse_parameters(prm);
  physical_properties.parse_parameters(prm);
  output_parameters.parse_parameters(prm);

  std::ofstream log_file(output_parameters.output_path /
                           ("convergence_expansion_order.csv"),
                         std::ios::app);
  saplog.attach(log_file, false);

  saplog.pop();
  saplog << "# "
         << "expansion_order"
         << "\t L2\t Linfty\t CPU time [s]\t Wall time [s]\t n_dof"
         << std::endl;
  saplog.push("Tests");
  for (unsigned int i = 0; i < values.size(); i++)
    {
      saplog << "expansion_order"
             << "=" << values[i] << std::endl;
      vfp_parameters.expansion_order = uint(values[i]);
      output_parameters.base_file_name =
        "expansion_order_" + dealii::Utilities::to_string(values[i]);

      VFPEquationSolver<dim> vfp_equation_solver(vfp_parameters,
                                                 physical_properties,
                                                 output_parameters);
      vfp_equation_solver.run();
      timer.stop();

      InitialValueFunction<dim> analytic_solution(
        physical_properties, vfp_parameters.expansion_order);
      analytic_solution.set_time(vfp_parameters.final_time);

      const double L2_error =
        vfp_equation_solver.compute_global_error(analytic_solution,
                                                 VectorTools::L2_norm,
                                                 VectorTools::L2_norm);
      const double Linfty_error =
        vfp_equation_solver.compute_global_error(analytic_solution,
                                                 VectorTools::L2_norm,
                                                 VectorTools::Linfty_norm);
      saplog.pop();
      saplog << values[i] << "\t" << L2_error << "\t" << Linfty_error << "\t"
             << timer.cpu_time() << "\t" << timer.wall_time() << "\t"
             << vfp_equation_solver.get_n_dofs() << std::endl;
      saplog.push("Tests");
    }

  saplog.detach();
  log_file.close();

  return 0;
}


int
test_run(const std::string &parameter_filename, const double max_L2_error)
{
  using namespace Sapphire;
  using namespace VFP;

  saplog.push("Tests");
  saplog << "Test run with parameter file \"" << parameter_filename
         << "\" and maximal L2 error of " << max_L2_error << std::endl;

  dealii::Timer                timer;
  ParameterHandler             prm;
  VFPParameters<dim>           vfp_parameters;
  PhysicalProperties           physical_properties;
  Utils::OutputParameters<dim> output_parameters;

  vfp_parameters.declare_parameters(prm);
  physical_properties.declare_parameters(prm);
  output_parameters.declare_parameters(prm);
  prm.parse_input(parameter_filename);

  vfp_parameters.parse_parameters(prm);
  physical_properties.parse_parameters(prm);
  output_parameters.parse_parameters(prm);

  timer.start();
  VFPEquationSolver<dim> vfp_equation_solver(vfp_parameters,
                                             physical_properties,
                                             output_parameters);
  vfp_equation_solver.run();
  timer.stop();

  InitialValueFunction<dim> analytic_solution(physical_properties,
                                              vfp_parameters.expansion_order);
  analytic_solution.set_time(vfp_parameters.final_time);

  const double L2_error =
    vfp_equation_solver.compute_global_error(analytic_solution,
                                             VectorTools::L2_norm,
                                             VectorTools::L2_norm);

  saplog << "L2_error = " << L2_error
         << ", CPU/wall time = " << timer.cpu_time() << "/" << timer.wall_time()
         << " s" << std::endl;

  AssertThrow(L2_error < max_L2_error,
              dealii::ExcMessage(
                "L2 error is too large! (" +
                dealii::Utilities::to_string(L2_error) + " > " +
                dealii::Utilities::to_string(max_L2_error) + ")"));
  saplog.pop();
  return 0;
}


int
main(int argc, char *argv[])
{
  try
    {
      using namespace Sapphire;
      using namespace VFP;
      dealii::Utilities::MPI::MPI_InitFinalize mpi_initialization(argc,
                                                                  argv,
                                                                  1);

      saplog.init();
      saplog.pop();
      saplog.depth_console(2);
      saplog.depth_file(0);
      saplog << "Start test-scattering-only-1d" << std::endl;

      std::string parameter_filename = "parameter.prm";
      if (argc > 1)
        parameter_filename = argv[1];

      std::string run_type = "test_run";
      if (argc > 2)
        run_type = argv[2];

      if (run_type == "test_run")
        {
          double max_L2_error = 1e-10;
          if (argc > 3)
            max_L2_error = std::stod(argv[3]);
          test_run(parameter_filename, max_L2_error);
        }
      else if (run_type == "expansion_order")
        {
          std::vector<double> values = {1, 2, 3, 4, 5, 6, 7, 8};
          if (argc > 3)
            {
              values.clear();
              for (int i = 3; i < argc; i++)
                values.push_back(std::stod(argv[i]));
            }
          convergence_with_expansion_order(parameter_filename, values);
        }
      else
        AssertThrow(false,
                    dealii::ExcMessage("Unknown run type \"" + run_type +
                                       "\"!"));
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
