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
      Flux(const Sapphire::Utils::ParameterParser &prm, const double &beta);

      void
      flux(const double &u, Tensor<1, dim> &flux_value) const;

      double
      numerical_flux(const double         &u_1,
                     const double         &u_2,
                     const Tensor<1, dim> &n);

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