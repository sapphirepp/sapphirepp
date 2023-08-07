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

      // dealii::Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv,
      // 1); //Only use MPI on one core
      dealii::Utilities::MPI::MPI_InitFinalize mpi_initialization(
        argc, argv); // Use TBB multithreading

      Sapphire::saplog << "n_cores = " << dealii::MultithreadInfo::n_cores()
                       << std::endl;
      Sapphire::saplog << "n_threads = " << dealii::MultithreadInfo::n_threads()
                       << std::endl;

      const unsigned int dim  = 1;
      const double       beta = 1.0;

      std::string parameter_filename = "../../examples/parameter.json";
      if (argc > 1)
        parameter_filename = argv[1];
      Sapphire::saplog << "parameter_filename = " << parameter_filename
                       << std::endl;

      ParameterParser   prm(parameter_filename);
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
