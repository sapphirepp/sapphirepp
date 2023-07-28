#include "flux.h"

#include <iostream>

template <int dim>
Sapphire::Hydro::Flux<dim>::Flux(const Sapphire::Utils::ParameterParser &prm,
                                 const double                           &beta)
  : flux_type(prm.hdsolver_flux_type)
  , beta(beta)
{}


template <int dim>
void
Sapphire::Hydro::Flux<dim>::flux(const double   &u,
                                 Tensor<1, dim> &flux_value) const
{
  flux_value[0] = beta * u * u;
}

template <int dim>
double
Sapphire::Hydro::Flux<dim>::numerical_flux(const double         &u_1,
                                           const double         &u_2,
                                           const Tensor<1, dim> &n)
{
  double numerical_flux = 0;
  flux(u_1, flux_1);
  flux(u_2, flux_2);
  switch (flux_type)
    {
      case FluxType::Central:
        {
          numerical_flux += 0.5 * (flux_1 + flux_2) * n;
          break;
        }

      case FluxType::Upwind:
        {
          numerical_flux += 0.5 * (flux_1 + flux_2) * n;
          numerical_flux +=
            0.5 * Upwind_eta * (std::abs(flux_1 * n) - std::abs(flux_2 * n));
          break;
        }

      case FluxType::LaxFriedrichs:
        {
          numerical_flux += 0.5 * (flux_1 + flux_2) * n;
          numerical_flux += 0.5 * LaxFriedrichs_C * (u_1 - u_2);
          break;
        }

      default:
        Assert(false, ExcNotImplemented());
        break;
    }

  return numerical_flux;
}

// explicit instantiation
template class Sapphire::Hydro::Flux<1>;
template class Sapphire::Hydro::Flux<2>;
template class Sapphire::Hydro::Flux<3>;
