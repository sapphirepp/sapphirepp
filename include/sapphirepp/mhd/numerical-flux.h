// -----------------------------------------------------------------------------
//
// Copyright (C) 2023 by the Sapphire++ authors
//
// This file is part of Sapphire++.
//
// Sapphire++ is free software: you can redistribute it and/or modify it under
// the terms of the GNU Lesser General Public License as published by the Free
// Software Foundation, either version 3 of the License, or (at your option) any
// later version.
//
// Sapphire++ is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
// A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more
// details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with Sapphire++. If not, see <https://www.gnu.org/licenses/>.
//
// -----------------------------------------------------------------------------

/**
 * @file numerical-flux.h
 * @author Florian Schulze (florian.schulze@mpi-hd.mpg.de)
 * @brief Define @ref sapphirepp::MHD::NumericalFlux
 */

#ifndef MHD_NUMERICALFLUX_H
#define MHD_NUMERICALFLUX_H

#include <deal.II/base/tensor.h>

#include <deal.II/lac/vector.h>

#include <vector>

#include "mhd-equations.h"
#include "mhd-flags.h"

namespace sapphirepp
{
  namespace MHD
  {
    template <unsigned int dim>
    class NumericalFlux
    {
    public:
      NumericalFlux(const MHDEquations<dim> &mhd_equations);



      void
      compute_numerical_normal_flux(
        const dealii::Tensor<1, dim>                 &normal,
        const typename MHDEquations<dim>::state_type &state_1,
        const typename MHDEquations<dim>::state_type &state_2,
        typename MHDEquations<dim>::state_type &numerical_normal_flux) const;



    private:
      const MHDEquations<dim> mhd_equations;
    };
  } // namespace MHD
} // namespace sapphirepp
#endif
