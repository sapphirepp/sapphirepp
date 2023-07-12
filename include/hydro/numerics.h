/**
 * @file numerics.h
 * @author Florian Schulze (florian.schulze@mpi-hd.mpg.de)
 * @brief Implement numerical functions to solve the hydrodynamics equations
 * @version 0.1
 * @date 2023-07-12
 *
 */

#ifndef HYDROSOLVER_NUMERICS_H
#define HYDROSOLVER_NUMERICS_H

#include <deal.II/lac/vector.h>

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
  HDSolverControl(
      const TimeSteppingScheme scheme = TimeSteppingScheme::ForwardEuler,
      const FluxType flux_type = FluxType::Central,
      const SlopeLimiter limiter = SlopeLimiter::NoLimiter,
      const unsigned int fe_degree = 1, const double time_step = 0.1,
      const double end_time = 1.0, const unsigned int refinement_level = 1)
      : scheme(scheme), flux_type(flux_type), limiter(limiter),
        fe_degree(fe_degree), time_step(time_step), end_time(end_time),
        refinement_level(refinement_level){};

  // Copy constructor
  HDSolverControl(const HDSolverControl &solver_control)
      : scheme(solver_control.scheme), flux_type(solver_control.flux_type),
        limiter(solver_control.limiter), fe_degree(solver_control.fe_degree),
        time_step(solver_control.time_step), end_time(solver_control.end_time),
        refinement_level(solver_control.refinement_level){};

  const TimeSteppingScheme scheme;
  const FluxType flux_type;
  const SlopeLimiter limiter;

  const unsigned int fe_degree;
  const double time_step;
  const double end_time;
  const unsigned int refinement_level;
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
                              const HDSolverControl &solver_control);

template <int dim>
void compute_limited_slope(const double &cell_average,
                           const Tensor<1, dim> cell_average_grad,
                           const std::vector<double> &neighbor_cell_averages,
                           const std::vector<Tensor<1, dim>> &neighbor_distance,
                           const unsigned int n_neighbors,
                           Tensor<1, dim> &limited_slope,
                           const HDSolverControl &solver_control);

} // namespace Hydro
} // namespace Sapphire
#endif