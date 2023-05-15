#include "burgers-eq.h"
#include <deal.II/base/mpi.h>
#include <mpi.h>
#include <iostream>

int main(int argc, char *argv[])
{
    try
    {
        using namespace Sapphire::Hydro;
        dealii::Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

        BurgersEq burgers_eq;
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
