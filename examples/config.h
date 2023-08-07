#include <deal.II/base/function.h>
#include <deal.II/base/tensor_function.h>


namespace Sapphire
{
  namespace Hydro
  {
    using namespace dealii;
    template <int dim>
    class ExactSolutionBurgersEq : public Function<dim>
    {
    public:
      ExactSolutionBurgersEq(const double time = 0.0)
        : Function<dim>(1, time)
      {
        AssertDimension(dim, 1);
      }

      double
      value(const Point<dim>  &p,
            const unsigned int component = 0) const override
      {
        (void)component; // suppress unused parameter warning
        Point<dim> x;

        return -std::sin(numbers::PI * p[0]);

        // x[0] = p[0];
        // return 0.25 + 0.5 * std::sin(numbers::PI * (2 * x[0] - 1));

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

    template <int dim>
    class InitialConditionBurgersEq : public Function<dim>
    {
    public:
      InitialConditionBurgersEq()
        : Function<dim>(1)
        , exact_solution()
      {
        exact_solution.set_time(0.0);
      }

      double
      value(const Point<dim>  &p,
            const unsigned int component = 0) const override
      {
        return exact_solution.value(p, component);
      }

    private:
      ExactSolutionBurgersEq<dim> exact_solution;
    };

    template <int dim>
    class BoundaryValuesBurgersEq : public Function<dim>
    {
    public:
      BoundaryValuesBurgersEq(const double time = 0.0)
        : Function<dim>(1, time)
        , exact_solution()
      {}

      double
      value(const Point<dim>  &p,
            const unsigned int component = 0) const override
      {
        return exact_solution.value(p, component);
      }

      void
      set_time(const double new_time) override
      {
        exact_solution.set_time(new_time);
      }



    private:
      ExactSolutionBurgersEq<dim> exact_solution;
    };

  } // namespace Hydro
} // namespace Sapphire
