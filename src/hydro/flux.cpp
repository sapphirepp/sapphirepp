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
Sapphire::Hydro::Flux<dim>::flux(const double      &u,
                                 Tensor<1, dim>    &flux_value,
                                 const unsigned int component) const
{
  AssertIndexRange(component, n_components);
  AssertDimension(dim, 1);

  // TODO_HD: Change to HD fluxes
  flux_value[0] = beta * u * u;
}

template <int dim>
void
Sapphire::Hydro::Flux<dim>::flux(const Vector<double>        &u,
                                 std::vector<Tensor<1, dim>> &flux_value) const
{
  AssertDimension(u.size(), n_components);
  AssertDimension(flux_value.size(), n_components);
  AssertDimension(dim, 1);

  for (unsigned int i = 0; i < n_components; ++i)
    flux(u[i], flux_value[i], i);
}

template <int dim>
double
Sapphire::Hydro::Flux<dim>::numerical_flux(const double         &u_1,
                                           const double         &u_2,
                                           const Tensor<1, dim> &n,
                                           const unsigned int    component)
{
  double numerical_flux = 0;
  flux(u_1, flux_1, component);
  flux(u_2, flux_2, component);
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
