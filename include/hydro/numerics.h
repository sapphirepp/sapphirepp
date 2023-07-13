/**
 * @file numerics.h
 * @author Florian Schulze (florian.schulze@mpi-hd.mpg.de)
 * @brief Implement numerical functions to solve the hydrodynamics equations
 * @version 0.1
 * @date 2023-07-12
 */

#ifndef HYDROSOLVER_NUMERICS_H
#define HYDROSOLVER_NUMERICS_H

#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_control.h>

#include <deal.II/lac/vector.h>

#include <deal.II/base/parameter_handler.h>

namespace Sapphire {
namespace Hydro {
using namespace dealii;

enum class TimeSteppingScheme { ForwardEuler, ExplicitRK };
enum class FluxType { Central, Upwind, LaxFriedrich };
enum class SlopeLimiter {
  NoLimiter,
  LinearReconstruction,
  MinMod,
  MUSCL,
  GerneralizedSlopeLimiter
};

class HDSolverControl {
public:
  HDSolverControl(const TimeSteppingScheme scheme, const FluxType flux_type,
                  const SlopeLimiter limiter, const unsigned int fe_degree,
                  const double time_step, const double end_time,
                  const unsigned int refinement_level,
                  const unsigned int max_iterations, const double tolerance)
      : scheme(scheme), flux_type(flux_type), limiter(limiter),
        fe_degree(fe_degree), time_step(time_step), end_time(end_time),
        refinement_level(refinement_level), max_iterations(max_iterations),
        tolerance(tolerance){};
  HDSolverControl(const HDSolverControl &hd_solver_control) = default;
  ~HDSolverControl() = default;

  static void declare_parameters(ParameterHandler &prm) {
    prm.enter_subsection("Hydrodynamics");

    prm.enter_subsection("Time stepping");
    prm.declare_entry("Scheme", "Forward Euler",
                      Patterns::Selection("Forward Euler|Explicit RK"),
                      "Time stepping scheme");
    prm.declare_entry("Time step", "0.1", Patterns::Double(0.0),
                      "Time step size");
    prm.declare_entry("End time", "1.0", Patterns::Double(0.0),
                      "End time of simulation");
    prm.leave_subsection(); // subsection Time stepping

    prm.enter_subsection("Finite element");
    prm.declare_entry("Polynomial degree", "1", Patterns::Integer(1),
                      "Polynomial degree of finite element");
    prm.leave_subsection(); // subsection Finite element

    prm.enter_subsection("Mesh");
    prm.declare_entry("Refinement level", "1", Patterns::Integer(0),
                      "Refinement level of mesh");
    prm.leave_subsection(); // subsection Mesh

    prm.enter_subsection("Numerical flux");
    prm.declare_entry("Numerical flux", "Lax-Friedrichs",
                      Patterns::Selection("Central|Upwind|Lax-Friedrichs"),
                      "Numerical flux");
    prm.leave_subsection(); // subsection Numerical flux

    prm.enter_subsection("Solpe limiter");
    prm.declare_entry(
        "Solpe limiter", "No limiter",
        Patterns::Selection(
            "No limiter|Linear reconstruction|MinMod|MUSCL|Generalised slope "
            "limiter"),
        "Slope limiter");
    prm.leave_subsection(); // subsection Slope limiter

    prm.enter_subsection("Linear solver");
    prm.declare_entry("Max iterations", "1000", Patterns::Integer(0),
                      "Maximum number of iterations");
    prm.declare_entry("Tolerance", "1e-12", Patterns::Double(0.0),
                      "Tolerance of solver");
    prm.leave_subsection(); // subsection Linear solver

    prm.leave_subsection(); // subsection Hydrodynamics
  };

  static HDSolverControl parse_parameters(ParameterHandler &prm) {
    std::string s;
    // TODO_BE: understand why this path does not work
    //  std::vector<std::string> path = {"Hydrodynamics", "Time stepping"};
    //  s = prm.get(path, "Solver");

    prm.enter_subsection("Hydrodynamics");

    prm.enter_subsection("Time stepping");
    TimeSteppingScheme scheme;
    s = prm.get("Scheme");
    if (s == "Forward Euler")
      scheme = TimeSteppingScheme::ForwardEuler;
    else if (s == "Explicit RK")
      scheme = TimeSteppingScheme::ExplicitRK;
    else
      AssertThrow(false, ExcNotImplemented());
    double time_step = prm.get_double("Time step");
    double end_time = prm.get_double("End time");
    prm.leave_subsection(); // subsection Time stepping

    prm.enter_subsection("Finite element");
    unsigned int fe_degree = prm.get_integer("Polynomial degree");
    prm.leave_subsection(); // subsection Finite element

    prm.enter_subsection("Mesh");
    unsigned int refinement_level = prm.get_integer("Refinement level");
    prm.leave_subsection(); // subsection Mesh

    prm.enter_subsection("Numerical flux");
    FluxType flux_type;
    s = prm.get("Numerical flux");
    if (s == "Central")
      flux_type = FluxType::Central;
    else if (s == "Upwind")
      flux_type = FluxType::Upwind;
    else if (s == "Lax-Friedrichs")
      flux_type = FluxType::LaxFriedrich;
    else
      AssertThrow(false, ExcNotImplemented());
    prm.leave_subsection(); // subsection Numerical flux

    prm.enter_subsection("Solpe limiter");
    SlopeLimiter limiter;
    s = prm.get("Solpe limiter");
    if (s == "No limiter")
      limiter = SlopeLimiter::NoLimiter;
    else if (s == "Linear reconstruction")
      limiter = SlopeLimiter::LinearReconstruction;
    else if (s == "MinMod")
      limiter = SlopeLimiter::MinMod;
    else if (s == "MUSCL")
      limiter = SlopeLimiter::MUSCL;
    else if (s == "Generalised slope limiter")
      limiter = SlopeLimiter::GerneralizedSlopeLimiter;
    else
      AssertThrow(false, ExcNotImplemented());
    prm.leave_subsection(); // subsection Solpe limiter

    prm.enter_subsection("Linear solver");
    unsigned int max_iterations = prm.get_integer("Max iterations");
    double tolerance = prm.get_double("Tolerance");
    prm.leave_subsection(); // subsection Linear solver

    prm.leave_subsection(); // subsection Hydrodynamics

    return HDSolverControl(scheme, flux_type, limiter, fe_degree, time_step,
                           end_time, refinement_level, max_iterations,
                           tolerance);
  };

  const TimeSteppingScheme scheme;
  const FluxType flux_type;
  const SlopeLimiter limiter;

  const unsigned int fe_degree;
  const double time_step;
  const double end_time;
  const unsigned int refinement_level;

  const unsigned int max_iterations;
  const double tolerance;
};

double minmod(const std::vector<double> &values);

template <int dim>
void minmod(const std::vector<Tensor<1, dim>> &values, const unsigned int n,
            Tensor<1, dim> &return_value);

/**
 * @brief Compute the numerical flux
 *
 * @param flux_1 Flux on the first side of the interface \f$ \mathbf{f}_+ \f$
 * @param flux_2 Flux on the second/neighbor side of the interface \f$
 * \mathbf{f}_- \f$
 * @param n Nommal vector of the interface \f$ \hat{n} \f$
 * @param value_1 Value of the function on the first side of the interface \f$
 * u_+ \f$
 * @param value_2 Value of the function on the second/neighbor side of the
 * interface \f$ u_- \f$
 * @return double Flux across the interface \f$ \hat{n} \cdot \mathbf{f} \f$
 */
template <int dim>
double compute_numerical_flux(const Tensor<1, dim> &flux_1,
                              const Tensor<1, dim> &flux_2,
                              const Tensor<1, dim> &n, const double &value_1,
                              const double &value_2,
                              const HDSolverControl &hd_solver_control);

template <int dim>
void compute_limited_slope(const double &cell_average,
                           const Tensor<1, dim> cell_average_grad,
                           const std::vector<double> &neighbor_cell_averages,
                           const std::vector<Tensor<1, dim>> &neighbor_distance,
                           const unsigned int n_neighbors,
                           Tensor<1, dim> &limited_slope,
                           const HDSolverControl &hd_solver_control);

} // namespace Hydro
} // namespace Sapphire
#endif