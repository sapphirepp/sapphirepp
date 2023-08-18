#include <deal.II/base/mpi.h>
#include <deal.II/base/multithread_info.h>
#include <deal.II/base/parameter_handler.h>

#include <mpi.h>

#include <iostream>

#include "hd-solver.h"
#include "output-module.h"
#include "parameter-parser.h"
#include "sapphire-logstream.h"

int
main(int argc, char **argv)
{
  using namespace Sapphire;
  using namespace Hydro;
  using namespace dealii;

  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

  try
    {
      deallog.depth_console(0);
      Sapphire::saplog.depth_console(100);

      std::string parameter_filename = "../parameter.prm";
      if (argc > 1)
        parameter_filename = argv[1];
      Sapphire::saplog << "parameter_filename = " << parameter_filename
                       << std::endl;

      ParameterParser prm(parameter_filename);
      prm.write_template_parameters("../parameter-template.json");

      OutputModule<dimension> output_module(prm);

      HDSolver<dimension> hd_solver(prm, output_module);
      hd_solver.run();
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
