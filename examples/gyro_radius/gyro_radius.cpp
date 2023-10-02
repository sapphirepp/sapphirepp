#include <deal.II/base/mpi.h>

#include <mpi.h>

#include "output-module.h"
#include "vfp-equation-solver.h"

int
main(int argc, char *argv[])
{
  try
    {
      using namespace Sapphire;
      using namespace VFP;
      Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

      saplog.pop();
      saplog.depth_console(1);
      // saplog.depth_console(-1);
      const unsigned int mpi_rank =
        Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);
      if (mpi_rank > 0)
        {
          saplog.depth_console(0);
          char buffer[4];
          std::snprintf(buffer, 4, "%03d", mpi_rank);
          saplog.push("mpi" + std::string(buffer));
        }

      std::string parameter_filename = "parameter.prm";
      if (argc > 1)
        parameter_filename = argv[1];

      int mpi_size = Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD);
      saplog << "Start gyro_radius test with parameter file \""
             << parameter_filename << "\" on " << mpi_size << " processor(s) ["
             << Utilities::System::get_date() << " "
             << Utilities::System::get_time() << "]" << std::endl;


      saplog.push("Gyro");
      Timer                  timer;
      ParameterHandler       prm;
      VFPSolverControl       vfp_solver_control;
      PhysicalProperties     physical_properties;
      Utils::OutputModule<2> output_module;

      vfp_solver_control.declare_parameters(prm);
      physical_properties.declare_parameters(prm);
      output_module.declare_parameters(prm);

      prm.parse_input(parameter_filename);

      vfp_solver_control.parse_parameters(prm);
      physical_properties.parse_parameters(prm);
      output_module.parse_parameters(prm);

      output_module.init(prm);

      TransportOnly particle_properties;
      saplog << "Gyroperiod: "
             << 2 * M_PI * particle_properties.gamma / physical_properties.B0
             << std::endl;
      saplog << "Gyroradius: "
             << particle_properties.gamma / physical_properties.B0 << std::endl;

      saplog.pop();

      timer.start();
      VFPEquationSolver vfp_equation_solver(vfp_solver_control,
                                            physical_properties,
                                            output_module);
      vfp_equation_solver.run();
      timer.stop();

      InitialValueFunction<VFPSolverControl::dim> cylinder(
        physical_properties, vfp_solver_control.expansion_order);
      double error =
        vfp_equation_solver.compute_global_error(cylinder,
                                                 VectorTools::L2_norm,
                                                 //  VectorTools::Linfty_norm);
                                                 VectorTools::L2_norm);

      saplog << "Error = " << error << " CPU/wall time = " << timer.cpu_time()
             << "/" << timer.wall_time() << " s" << std::endl;
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
