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

#include "parameter-parser.h"

namespace Sapphire
{
  namespace Hydro
  {
    using namespace dealii;

    template <int dim>
    class Flux
    {
    public:
      static constexpr unsigned int n_components             = dim + 2;
      static constexpr unsigned int first_momentum_component = 0;
      static constexpr unsigned int density_component        = dim;
      static constexpr unsigned int energy_component         = dim + 1;

      Flux(const Sapphire::Utils::ParameterParser &prm, const double &beta);

      void
      flux(const double      &u,
           Tensor<1, dim>    &flux_value,
           const unsigned int component) const;
      void
      flux(const Vector<double>        &u,
           std::vector<Tensor<1, dim>> &flux_value) const;

      double
      numerical_flux(const double         &u_1,
                     const double         &u_2,
                     const Tensor<1, dim> &n,
                     const unsigned int    component); // TODO_HD: const


    private:
      const FluxType flux_type;
      const double   beta;

      double Upwind_eta      = 1.0;
      double LaxFriedrichs_C = 3.0; // TODO_HD: Calulate C


      Tensor<1, dim> flux_1, flux_2;
    };

  } // namespace Hydro
} // namespace Sapphire

#endif