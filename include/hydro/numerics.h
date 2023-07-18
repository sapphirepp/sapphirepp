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

#include "parameter_parser.h"

namespace Sapphire
{
  namespace Hydro
  {
    using namespace dealii;

    class HDSolverControl
    {
    public:
      HDSolverControl(const Sapphire::Utils::ParameterParser &prm);

      const TimeSteppingScheme    scheme;
      const FluxType              flux_type;
      const SlopeLimiter          limiter;
      const SlopeLimiterCriterion limiter_criterion;

      const unsigned int fe_degree;
      const double       time_step;
      const double       end_time;
      const unsigned int refinement_level;

      const unsigned int max_iterations;
      const double       tolerance;

      double Upwind_eta = 1.0; //< Parameter controling upwind or central flux
      double LaxFriedrichs_C;
    };

    double
    minmod(const std::vector<double> &values);

    template <int dim>
    void
    minmod(const std::vector<Tensor<1, dim>> &values,
           const unsigned int                 n,
           Tensor<1, dim>                    &return_value);

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
    double
    compute_numerical_flux(const Tensor<1, dim>  &flux_1,
                           const Tensor<1, dim>  &flux_2,
                           const Tensor<1, dim>  &n,
                           const double          &value_1,
                           const double          &value_2,
                           const HDSolverControl &hd_solver_control);

    template <int dim>
    void
    compute_limited_slope(const double              &cell_average,
                          const Tensor<1, dim>       cell_average_grad,
                          const std::vector<double> &neighbor_cell_averages,
                          const std::vector<Tensor<1, dim>> &neighbor_distance,
                          const unsigned int                 n_neighbors,
                          Tensor<1, dim>                    &limited_slope,
                          const HDSolverControl             &hd_solver_control);

  } // namespace Hydro
} // namespace Sapphire
#endif