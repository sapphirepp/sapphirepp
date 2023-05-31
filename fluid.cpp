#include "conservation-eq.h"

#include <deal.II/base/mpi.h>

#include <iostream>
#include <mpi.h>

namespace Sapphire {
namespace Hydro {
using namespace dealii;

/**
 * @brief Exact analytical solution of the conservation equation.
 *
 *  \( u(x, t) = u_0(x - a \cdot t) \)
 */
template <int dim> class ExactSolution : public Function<dim> {
public:
  ExactSolution(const Tensor<1, dim> &beta, const double time = 0.0)
      : Function<dim>(1, time), beta(beta) {}
  double value(const Point<dim> &p,
               const unsigned int component = 0) const override {
    (void)component; // supress unused variable warning

    // return 1.0;
    // Tensor<1, dim> normal;
    // normal[0] = 1.0;
    // normal[1] = 1.0;
    // normal /= normal.norm();
    // return std::sin(numbers::PI * normal * (p - beta * this->get_time()));
    const double sigma = 0.1;
    return std::exp(-(p - beta * this->get_time()) *
                    (p - beta * this->get_time()) / (2.0 * sigma * sigma));
  }

private:
  const Tensor<1, dim> beta;
};

/**
 * @brief Iniitial condition for the conservation equation.
 *
 * \( u_0(x) = sin(x) \)
 */
template <int dim> class InitialCondition : public Function<dim> {
public:
  InitialCondition(const Tensor<1, dim> &beta) : Function<dim>(1), beta(beta) {}
  double value(const Point<dim> &p,
               const unsigned int component = 0) const override {
    return ExactSolution<dim>(beta, 0.0).value(p, component);
  }

private:
  const Tensor<1, dim> beta;
};

template <int dim> class BoundaryValues : public Function<dim> {
public:
  BoundaryValues(const Tensor<1, dim> &beta, const double time = 0.0)
      : Function<dim>(1, time), beta(beta) {}
  double value(const Point<dim> &p,
               const unsigned int component = 0) const override {
    return ExactSolution<dim>(beta, this->get_time()).value(p, component);
  }

private:
  const Tensor<1, dim> beta;
};

} // namespace Hydro
} // namespace Sapphire

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
