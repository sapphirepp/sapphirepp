#include <deal.II/base/mpi.h>

#include <mpi.h>

#include "config.h"
#include "output-module.h"
#include "vfp-equation-solver.h"


const unsigned int dim = 1;

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
      saplog.depth_console(2);
      const unsigned int mpi_rank =
        dealii::Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);
      if (mpi_rank > 0)
        {
          saplog.depth_console(0);
          saplog.push("mpi" + dealii::Utilities::to_string(mpi_rank, 3));
        }

      std::string parameter_filename = "parameter.prm";
      if (argc > 1)
        parameter_filename = argv[1];

      double max_L2_error = 1e-10;
      if (argc > 2)
        max_L2_error = std::stod(argv[2]);

      int mpi_size = dealii::Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD);
      saplog << "Start test scattering with parameter file \""
             << parameter_filename << "\" on " << mpi_size << " processor(s) ["
             << dealii::Utilities::System::get_date() << " "
             << dealii::Utilities::System::get_time() << "]" << std::endl;
      saplog << "Max L2 error: " << max_L2_error << std::endl;

      saplog.push("main");
      dealii::Timer                  timer;
      ParameterHandler               prm;
      VFPSolverControl<dimension>    vfp_solver_control;
      PhysicalProperties             physical_properties;
      Utils::OutputModule<dimension> output_module;

      vfp_solver_control.declare_parameters(prm);
      physical_properties.declare_parameters(prm);
      output_module.declare_parameters(prm);
      prm.parse_input(parameter_filename);

      vfp_solver_control.parse_parameters(prm);
      physical_properties.parse_parameters(prm);
      output_module.parse_parameters(prm);

      saplog.pop(); // pop "main" prefix

      VFPEquationSolver<dimension> vfp_equation_solver(vfp_solver_control,
                                                       physical_properties,
                                                       output_module);
      timer.start();
      vfp_equation_solver.run();
      timer.stop();

      InitialValueFunction<dimension> analytic_solution(
        physical_properties, vfp_solver_control.expansion_order);
      analytic_solution.set_time(vfp_solver_control.final_time);

      const double L2_error =
        vfp_equation_solver.compute_global_error(analytic_solution,
                                                 VectorTools::L2_norm,
                                                 VectorTools::L2_norm);

      saplog << "L2_error=" << L2_error
             << " CPU/wall time = " << timer.cpu_time() << "/"
             << timer.wall_time() << " s" << std::endl;

      AssertThrow(L2_error < max_L2_error,
                  dealii::ExcMessage("L2 error is too large!"));
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
