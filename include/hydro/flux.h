/**
 * @file flux.h
 * @author Florian Schulze (florian.schulze@mpi-hd.mpg.de)
 * @brief
 * @version 0.1
 * @date 2023-07-28
 */


#ifndef HYDROSOLVER_FLUX_H
#define HYDROSOLVER_FLUX_H

#include <deal.II/lac/vector.h>

#include "config.h"
#include "parameter-parser.h"

namespace Sapphire
{
  namespace Hydro
  {
    using namespace dealii;

    template <int dim, typename Number>
    inline DEAL_II_ALWAYS_INLINE Tensor<1, dim, Number>
    euler_velocity(const Tensor<1, dim + 2, Number> &conserved_variables)
    {
      const Number inverse_density = Number(1.) / conserved_variables[0];

      Tensor<1, dim, Number> velocity;
      for (unsigned int d = 0; d < dim; ++d)
        velocity[d] = conserved_variables[1 + d] * inverse_density;

      return velocity;
    }

    template <int dim, typename Number>
    inline DEAL_II_ALWAYS_INLINE Number
    euler_pressure(const Tensor<1, dim + 2, Number> &conserved_variables)
    {
      const Tensor<1, dim, Number> velocity =
        euler_velocity<dim>(conserved_variables);

      Number rho_u_dot_u = conserved_variables[1] * velocity[0];
      for (unsigned int d = 1; d < dim; ++d)
        rho_u_dot_u += conserved_variables[1 + d] * velocity[d];

      return (gamma - 1.) * (conserved_variables[dim + 1] - 0.5 * rho_u_dot_u);
    }

    template <int dim, typename Number>
    inline DEAL_II_ALWAYS_INLINE Tensor<1, dim + 2, Tensor<1, dim, Number>>
    euler_flux(const Tensor<1, dim + 2, Number> &conserved_variables)
    {
      const Tensor<1, dim, Number> velocity =
        euler_velocity<dim>(conserved_variables);
      const Number pressure = euler_pressure<dim>(conserved_variables);

      Tensor<1, dim + 2, Tensor<1, dim, Number>> flux;
      for (unsigned int d = 0; d < dim; ++d)
        {
          flux[0][d] = conserved_variables[1 + d];
          for (unsigned int e = 0; e < dim; ++e)
            flux[e + 1][d] = conserved_variables[e + 1] * velocity[d];
          flux[d + 1][d] += pressure;
          flux[dim + 1][d] =
            velocity[d] * (conserved_variables[dim + 1] + pressure);
        }

      return flux;
    }

    template <int n_components, int dim, typename Number>
    inline DEAL_II_ALWAYS_INLINE Tensor<1, n_components, Number>
    operator*(const Tensor<1, n_components, Tensor<1, dim, Number>> &matrix,
              const Tensor<1, dim, Number>                          &vector)
    {
      Tensor<1, n_components, Number> result;
      for (unsigned int d = 0; d < n_components; ++d)
        result[d] = matrix[d] * vector;
      return result;
    }

    template <int dim, typename Number>
    inline DEAL_II_ALWAYS_INLINE Tensor<1, dim + 2, Number>
    euler_numerical_flux(const Tensor<1, dim + 2, Number> &u_m,
                         const Tensor<1, dim + 2, Number> &u_p,
                         const Tensor<1, dim, Number>     &normal)
    {
      const auto velocity_m = euler_velocity<dim>(u_m);
      const auto velocity_p = euler_velocity<dim>(u_p);

      const auto pressure_m = euler_pressure<dim>(u_m);
      const auto pressure_p = euler_pressure<dim>(u_p);

      const auto flux_m = euler_flux<dim>(u_m);
      const auto flux_p = euler_flux<dim>(u_p);

      switch (numerical_flux_type)
        {
          case lax_friedrichs_modified:
            {
              const auto lambda =
                0.5 * std::sqrt(std::max(velocity_p.norm_square() +
                                           gamma * pressure_p * (1. / u_p[0]),
                                         velocity_m.norm_square() +
                                           gamma * pressure_m * (1. / u_m[0])));

              return 0.5 * (flux_m * normal + flux_p * normal) +
                     0.5 * lambda * (u_m - u_p);
            }

          case harten_lax_vanleer:
            {
              const auto avg_velocity_normal =
                0.5 * ((velocity_m + velocity_p) * normal);
              const auto   avg_c = std::sqrt(std::abs(
                0.5 * gamma *
                (pressure_p * (1. / u_p[0]) + pressure_m * (1. / u_m[0]))));
              const Number s_pos =
                std::max(Number(), avg_velocity_normal + avg_c);
              const Number s_neg =
                std::min(Number(), avg_velocity_normal - avg_c);
              const Number inverse_s = Number(1.) / (s_pos - s_neg);

              return inverse_s *
                     ((s_pos * (flux_m * normal) - s_neg * (flux_p * normal)) -
                      s_pos * s_neg * (u_m - u_p));
            }

          default:
            {
              Assert(false, ExcNotImplemented());
              return {};
            }
        }
    }

    template <int dim, typename Number>
    VectorizedArray<Number>
    evaluate_function(const Function<dim>                       &function,
                      const Point<dim, VectorizedArray<Number>> &p_vectorized,
                      const unsigned int                         component)
    {
      VectorizedArray<Number> result;
      for (unsigned int v = 0; v < VectorizedArray<Number>::size(); ++v)
        {
          Point<dim> p;
          for (unsigned int d = 0; d < dim; ++d)
            p[d] = p_vectorized[d][v];
          result[v] = function.value(p, component);
        }
      return result;
    }


    template <int dim, typename Number, int n_components = dim + 2>
    Tensor<1, n_components, VectorizedArray<Number>>
    evaluate_function(const Function<dim>                       &function,
                      const Point<dim, VectorizedArray<Number>> &p_vectorized)
    {
      AssertDimension(function.n_components, n_components);
      Tensor<1, n_components, VectorizedArray<Number>> result;
      for (unsigned int v = 0; v < VectorizedArray<Number>::size(); ++v)
        {
          Point<dim> p;
          for (unsigned int d = 0; d < dim; ++d)
            p[d] = p_vectorized[d][v];
          for (unsigned int d = 0; d < n_components; ++d)
            result[d][v] = function.value(p, d);
        }
      return result;
    }

  } // namespace Hydro
} // namespace Sapphire

#endif