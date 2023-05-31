#include "conservation-eq.h"

#include <deal.II/base/mpi.h>

#include <iostream>
#include <mpi.h>

namespace PhysicalSetup {
using namespace dealii;

/**
 * @brief Exact analytical solution of constant linear advection equation.
 *
 * \f$ u(x, t) = u_0(\vec{x} - \vec{\beta} t) \F$
 * with
 * \f$ u_0(\vec{x}) =  ... \f$
 *
 * \tparam dim: space dimension
 */
template <int dim> class ExactSolution : public Function<dim> {
public:
  ExactSolution(const Tensor<1, dim> &beta, const double time = 0.0)
      : Function<dim>(1, time), beta(beta) {}

  double value(const Point<dim> &p,
               const unsigned int component = 0) const override {
    (void)component; // supress unused variable warning

    // return 1.0;
    Tensor<1, dim> normal;
    normal[0] = 1.0;
    normal[1] = 1.0;
    normal /= normal.norm();
    return std::sin(numbers::PI * normal * (p - beta * this->get_time()));
    const double sigma = 0.1;
    return std::exp(-(p - beta * this->get_time()) *
                    (p - beta * this->get_time()) / (2.0 * sigma * sigma));
  }

private:
  const Tensor<1, dim> beta;
};

/**
 * @brief Iniitial condition extracted from an exact solution.
 *
 * \f$ u_0(\vec{x}) = u(\vec{x}, t = 0) \f$
 *
 * \tparam dim: space dimension
 */
template <int dim> class InitialCondition : public Function<dim> {
public:
  InitialCondition(Function<dim> *exact_solution)
      : Function<dim>(1), exact_solution(exact_solution) {}

  double value(const Point<dim> &p,
               const unsigned int component = 0) const override {
    exact_solution->set_time(0.0);
    return exact_solution->value(p, component);
  }

private:
  const SmartPointer<Function<dim>> exact_solution;
};

/**
 * @brief Boundary values extracted from an exact solution.
 *
 * \f$ u_b(\vec{x}, t) = u(\vec{x}, t) \f$
 *
 * @tparam dim: space dimension
 */
template <int dim> class BoundaryValues : public Function<dim> {
public:
  BoundaryValues(Function<dim> *exact_solution, const double time = 0.0)
      : Function<dim>(1, time), exact_solution(exact_solution) {}

  double value(const Point<dim> &p,
               const unsigned int component = 0) const override {
    exact_solution->set_time(this->get_time());
    return exact_solution->value(p, component);
  }

private:
  const SmartPointer<Function<dim>> exact_solution;
};

} // namespace PhysicalSetup

int main(int argc, char *argv[]) {
  try {
    using namespace Sapphire::Hydro;
    dealii::Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

    // const unsigned int dim = 1;
    const unsigned int dim = 2;
    // const unsigned int dim = 3;
    const Tensor<1, dim> beta({-0.5});
    // const Tensor<1, dim> beta({+0.5, 0.0});
    // const Tensor<1, dim> beta({0, +0.5});

    PhysicalSetup::ExactSolution<dim> exact_solution(beta);
    SmartPointer<Function<dim>> exact_solution_ptr(&exact_solution);

    PhysicalSetup::BoundaryValues<dim> boundary_values(&exact_solution);
    // PhysicalSetup::BoundaryValues<dim> boundary_values(exact_solution_ptr);
    PhysicalSetup::InitialCondition<dim> initial_condition(&exact_solution);

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
