#include <deal.II/base/mpi.h>
#include <deal.II/base/multithread_info.h>
#include <deal.II/base/parameter_handler.h>

#include <mpi.h>

#include <iostream>

#include "burgers-eq.h"
#include "conservation-eq.h"
#include "numerics.h"
#include "output-module.h"
#include "parameter-parser.h"
#include "sapphire-logstream.h"

int
main(int argc, char *argv[])
{
  try
    {
      using namespace Sapphire::Hydro;
      using namespace Sapphire::Utils;

      Sapphire::saplog.depth_console(10);

      DEBUG_PRINT(std::cout, 3, argc);
      for (int i = 3; i < argc; ++i)
        DEBUG_PRINT(std::cout, 3, argv[i]);

      // dealii::Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv,
      // 1); //Only use MPI on one core
      dealii::Utilities::MPI::MPI_InitFinalize mpi_initialization(
        argc, argv); // Use TBB multithreading

      DEBUG_PRINT(std::cout,
                  3,
                  "n_cores = " << dealii::MultithreadInfo::n_cores());
      DEBUG_PRINT(std::cout,
                  3,
                  "n_threads = " << dealii::MultithreadInfo::n_threads());

      const unsigned int                        dim  = 1;
      const double                              beta = 1.0;

      std::string parameter_filename = "../parameter.prm";
      if (argc > 1)
        parameter_filename = argv[1];
      DEBUG_PRINT(std::cout, 0, parameter_filename);

      ParameterParser prm(parameter_filename);
      prm.write_template_parameters("../parameter-template.json");

      OutputModule<dim> output_module(prm);

      BurgersEq<dim> burgers_eq(prm, output_module, beta);
      burgers_eq.run();
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
