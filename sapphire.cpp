#include <deal.II/base/mpi.h>

#include <mpi.h>

#include "vfp-equation-solver.h"

int
main(int argc, char *argv[])
{
  try
    {
      using namespace Sapphire;
      using namespace VFP;
      Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

      Sapphire::saplog.depth_console(10);

      std::string parameter_filename = "../vfp-equation.json";
      // std::string parameter_filename = "../parameter-template.json";
      if (argc > 1)
        parameter_filename = argv[1];

      Utils::ParameterParser parameter_parser(parameter_filename);

      parameter_parser.write_template_parameters("../parameter-template.json");

      VFPEquationSolver      vfp_equation_solver(parameter_parser);
      vfp_equation_solver.run();
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
