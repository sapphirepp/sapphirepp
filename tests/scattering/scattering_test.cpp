#include <deal.II/base/mpi.h>

#include <mpi.h>

#include "config.h"
#include "output-module.h"
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
                           ("convergence_expansion_order.csv"),
                         std::ios::app);
  saplog.attach(log_file, false);

  saplog << "# "
         << "expansion_order"
         << "\t L2\t Linfty\t CPU time [s]\t Wall time [s]\t n_dof"
         << std::endl;
  for (unsigned int i = 0; i < values.size(); i++)
    {
      vfp_solver_control.expansion_order = uint(values[i]);

      saplog.push("Tests");
      saplog << "expansion_order"
             << "=" << values[i] << std::endl;
      saplog.pop();

      output_module.base_file_name =
        "expansion_order_" + dealii::Utilities::to_string(values[i]);

      timer.start();
      VFPEquationSolver<dim> vfp_equation_solver(vfp_solver_control,
                                                 physical_properties,
                                                 output_module);
      vfp_equation_solver.run();
      timer.stop();

      InitialValueFunction<dim> analytic_solution(
        physical_properties, vfp_solver_control.expansion_order);
      analytic_solution.set_time(vfp_solver_control.final_time);

      const double L2_error =
        vfp_equation_solver.compute_global_error(analytic_solution,
                                                 VectorTools::L2_norm,
                                                 VectorTools::L2_norm);
      const double Linfty_error =
        vfp_equation_solver.compute_global_error(analytic_solution,
                                                 VectorTools::L2_norm,
                                                 VectorTools::Linfty_norm);

      saplog << values[i] << "\t" << L2_error << "\t" << Linfty_error << "\t"
             << timer.cpu_time() << "\t" << timer.wall_time() << "\t"
             << vfp_equation_solver.get_n_dofs() << std::endl;
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

  dealii::Timer            timer;
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

  timer.start();
  VFPEquationSolver<dim> vfp_equation_solver(vfp_solver_control,
                                             physical_properties,
                                             output_module);
  vfp_equation_solver.run();
  timer.stop();

  saplog.push("Tests");
  InitialValueFunction<dim> analytic_solution(
    physical_properties, vfp_solver_control.expansion_order);
  analytic_solution.set_time(vfp_solver_control.final_time);

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

      saplog.pop();
      saplog.depth_console(1);
      saplog.depth_file(0);
      const unsigned int mpi_rank =
        dealii::Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);
      if (mpi_rank > 0)
        {
          saplog.depth_console(0);
          saplog.push("mpi" + dealii::Utilities::to_string(mpi_rank, 3));
        }
      int mpi_size = dealii::Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD);
      saplog << "Start scattering_test on " << mpi_size << " processor(s) ["
             << dealii::Utilities::System::get_date() << " "
             << dealii::Utilities::System::get_time() << "]" << std::endl;

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
