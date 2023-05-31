#include "conservation-eq.h"

#include <deal.II/base/mpi.h>

#include <iostream>
#include <mpi.h>

int main(int argc, char *argv[]) {
  try {
    using namespace Sapphire::Hydro;
    dealii::Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

    // const unsigned int dim = 1;
    const unsigned int dim = 2;
    const Tensor<1, dim> beta({-0.5});
    // const Tensor<1, dim> beta({+0.5, 0.0});
    // const Tensor<1, dim> beta({0, +0.5});

    // TODO: Activate MPI/OPM
    InitialCondition<dim> initial_condition(beta);
    BoundaryValues<dim> boundary_values(beta);
    ExactSolution<dim> exact_solution(beta);

    ConservationEq<dim> conservation_eq(beta, &initial_condition,
                                        &boundary_values, &exact_solution);
    conservation_eq.run();
  } catch (std::exception &exc) {
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
  } catch (...) {
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
