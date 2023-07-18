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

namespace PhysicalSetup
{
  using namespace dealii;

  /**
   * @brief Velocity field for a rigid rotator.
   *
   * \f$ \mathbf{\beta}(\v{x}, t) = \omega (y \mathbf{e}_x - x \mathbf{e}_y) \f$
   *
   * \tparam dim: space dimension
   */
  template <int dim>
  class VelocityFieldRigidRotator : public TensorFunction<1, dim, double>
  {
  public:
    VelocityFieldRigidRotator(const double &omega, const double time = 0.0)
      : TensorFunction<1, dim, double>(time)
      , omega(omega)
    {}

    Tensor<1, dim>
    value(const Point<dim> &p) const override
    {
      AssertDimension(dim, 2);

      Tensor<1, dim> values;
      values[0] = omega * p[1];
      values[1] = -omega * p[0];
      return values;
    }

    Tensor<2, dim>
    gradient(const Point<dim> &p) const override
    {
      (void)p; // suppress unused parameter warning
      AssertDimension(dim, 2);

      Tensor<2, dim> values;
      values[0][0] = 0.0;
      values[0][1] = omega;
      values[1][0] = -omega;
      values[1][1] = 0.0;
      return values;
    }

  private:
    const double omega;
  };

  /**
   * @brief Exact analytical solution of rigid rotator.
   *
   * \tparam dim: space dimension
   */
  template <int dim>
  class ExactSolutionRigidRotator : public Function<dim>
  {
  public:
    ExactSolutionRigidRotator(const double &omega, const double time = 0.0)
      : Function<dim>(1, time)
      , omega(omega)
    {}

    double
    value(const Point<dim> &p, const unsigned int component = 0) const override
    {
      AssertIndexRange(component, 1);
      AssertDimension(dim, 2);

      Point<dim> x;
      x[0] = p[0] * std::cos(omega * this->get_time()) -
             p[1] * std::sin(omega * this->get_time());
      x[1] = p[0] * std::sin(omega * this->get_time()) +
             p[1] * std::cos(omega * this->get_time());

      const double width  = 0.2;
      const double length = 0.7;

      if ((std::abs(x[0]) < width) && (std::abs(x[1]) < length))
        return 1.0;
      if ((std::abs(x[1]) < width) && (std::abs(x[0]) < length))
        return 1.0;

      return 0.0;
    }

  private:
    const double omega;
  };

  /**
   * @brief Constant velocity field.
   *
   * \tparam dim: space dimension
   */
  template <int dim>
  class ConstantVelocityField : public TensorFunction<1, dim, double>
  {
  public:
    ConstantVelocityField(const Tensor<1, dim> &beta)
      : TensorFunction<1, dim, double>()
      , beta(beta)
    {}

    Tensor<1, dim>
    value(const Point<dim> &p) const override
    {
      (void)p; // suppress unused parameter warning

      return beta;
    }

    Tensor<2, dim>
    gradient(const Point<dim> &p) const override
    {
      (void)p; // suppress unused parameter warning

      Tensor<2, dim> values;
      values[0][0] = 0.0;
      values[0][1] = 0.0;
      values[1][0] = 0.0;
      values[1][1] = 0.0;
      return values;
    }

  private:
    const Tensor<1, dim> beta;
  };

  /**
   * @brief Exact analytical solution of constant linear advection equation.
   *
   * \f$ u(x, t) = u_0(x - \beta t) \f$
   * with
   * <!-- \f$ u_0(\mathbf{x}) =  1 \f$ -->
   * <!-- \f$ u_0(\mathbf{x}) =  \sin(\pi * \hat{n} \cdot \mathbf{x}) \f$ -->
   * \f$ u_0(\mathbf{x}) =  \exp(-\mathbf{x} \cdot \mathbf{x} / (2 \sigma^2))
   * \f$
   *
   * \tparam dim: space dimension
   */
  template <int dim>
  class ExactSolutionConstantAdvection : public Function<dim>
  {
  public:
    ExactSolutionConstantAdvection(const Tensor<1, dim> &beta,
                                   const double          time = 0.0)
      : Function<dim>(1, time)
      , beta(beta)
    {}

    double
    value(const Point<dim> &p, const unsigned int component = 0) const override
    {
      AssertIndexRange(component, 1);

      const Point<dim> x = p - beta * this->get_time();

      // return 1.0;
      // Tensor<1, dim> normal;
      // normal[0] = 1.0;
      // normal[1] = 1.0;
      // normal /= normal.norm();
      // return std::sin(numbers::PI * normal * x);
      const double sigma = 0.1;
      return std::exp(-(x * x) / (2.0 * sigma * sigma));
    }

  private:
    const Tensor<1, dim> beta;
  };

  /**
   * @brief Exact solution of Burgers' equation.
   *
   * \f$ u(x, t) = \Theta(x - t/2) \f$
   * with \f$ \Theta(x) \f$ being the Heaviside step function.
   *
   * \tparam dim: space dimension (must be 1)
   */
  template <int dim>
  class ExactSolutionBurgerEq : public Function<dim>
  {
  public:
    ExactSolutionBurgerEq(const double time = 0.0)
      : Function<dim>(1, time)
    {
      AssertDimension(dim, 1);
    }

    double
    value(const Point<dim> &p, const unsigned int component = 0) const override
    {
      (void)component; // suppress unused parameter warning
      Point<dim> x;

      // return -std::sin(numbers::PI * p[0]);

      // x[0] = p[0] - this->get_time() / 2.0;
      // if (x[0] > 0.0)
      //   return 1.0;
      // else
      //   return 0.0;

      x[0] = p[0] - this->get_time() * 3.0;
      if (x[0] > -0.5)
        return 1.0;
      else
        return 2.0;
    }
  };

  /**
   * @brief Initial condition extracted from an exact solution.
   *
   * \f$ u_0(\mathbf{x}) = u(\mathbf{x}, t = 0) \f$
   *
   * \tparam dim: space dimension
   */
  template <int dim>
  class InitialCondition : public Function<dim>
  {
  public:
    InitialCondition(Function<dim> *exact_solution)
      : Function<dim>(1)
      , exact_solution(exact_solution)
    {}

    double
    value(const Point<dim> &p, const unsigned int component = 0) const override
    {
      exact_solution->set_time(0.0);
      return exact_solution->value(p, component);
    }

  private:
    const SmartPointer<Function<dim>> exact_solution;
  };

  /**
   * @brief Boundary values extracted from an exact solution.
   *
   * \f$ u_b(\mathbf{x}, t) = u(\mathbf{x}, t) \f$
   *
   * @tparam dim: space dimension
   */
  template <int dim>
  class BoundaryValues : public Function<dim>
  {
  public:
    BoundaryValues(Function<dim> *exact_solution, const double time = 0.0)
      : Function<dim>(1, time)
      , exact_solution(exact_solution)
    {}

    double
    value(const Point<dim> &p, const unsigned int component = 0) const override
    {
      exact_solution->set_time(this->get_time());
      return exact_solution->value(p, component);
    }

  private:
    const SmartPointer<Function<dim>> exact_solution;
  };

} // namespace PhysicalSetup

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

      // const unsigned int dim = 1;
      // const unsigned int dim = 2;
      // const unsigned int dim = 3;

      // const double omega = 0.2;
      // // const double omega = 2 * numbers::PI; //Full circle
      // // const double omega = numbers::PI; // Half circle
      // PhysicalSetup::VelocityFieldRigidRotator<dim> beta(omega);
      // PhysicalSetup::ExactSolutionRigidRotator<dim> exact_solution(omega);

      // // const Tensor<1, dim> v({-0.5});
      // // const Tensor<1, dim> v({+0.5, 0.0});
      // // const Tensor<1, dim> v({0, +0.5});
      // // PhysicalSetup::ConstantVelocityField<dim> beta(v);
      // // PhysicalSetup::ExactSolutionConstantAdvection<dim>
      // exact_solution(v);

      // PhysicalSetup::BoundaryValues<dim> boundary_values(&exact_solution);
      // // SmartPointer<Function<dim>> exact_solution_ptr(&exact_solution);
      // // PhysicalSetup::BoundaryValues<dim>
      // boundary_values(exact_solution_ptr);
      // PhysicalSetup::InitialCondition<dim>
      // initial_condition(&exact_solution);

      // ConservationEq<dim> conservation_eq(&beta, &initial_condition,
      //                                     &boundary_values, &exact_solution);

      // conservation_eq.run();

      const unsigned int                        dim  = 1;
      const double                              beta = 1.0;
      PhysicalSetup::ExactSolutionBurgerEq<dim> exact_solution;
      PhysicalSetup::BoundaryValues<dim>   boundary_values(&exact_solution);
      PhysicalSetup::InitialCondition<dim> initial_condition(&exact_solution);

      std::string parameter_filename = "../parameter.prm";
      if (argc > 1)
        parameter_filename = argv[1];
      DEBUG_PRINT(std::cout, 0, parameter_filename);

      ParameterParser prm(parameter_filename);
      prm.write_template_parameters("../parameter-template.prm");

      OutputModule<dim> output_module(prm);

      BurgersEq<dim> burgers_eq(&initial_condition,
                                &boundary_values,
                                &exact_solution,
                                prm,
                                output_module,
                                beta);
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
