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
 * @file mhd-equations.h
 * @author Florian Schulze (florian.schulze@mpi-hd.mpg.de)
 * @brief Define @ref sapphirepp::MHD::MHDEquations
 */

#ifndef MHD_MHDEQUATIONS_H
#define MHD_MHDEQUATIONS_H

#include <deal.II/base/tensor.h>

#include <string>
#include <vector>

namespace sapphirepp
{
  namespace MHD
  {
    using namespace dealii;



    /**
     * @brief Define functions related to MHD equations.
     *
     * This class calculates everything related to the PDE system resulting from
     * the MHD equations, including the mapping between the components and the
     * corresponding physical quantities.
     */
    template <unsigned int dim>
    class MHDEquations
    {
    public:
      /** @{ */
      /** Dimension of the magnetic field */
      static constexpr unsigned int dim_B = 3;
      /** Number of components */
      static constexpr unsigned int n_components             = 2 + dim + dim_B;
      static constexpr unsigned int density_component        = 0;
      static constexpr unsigned int first_momentum_component = 1;
      static constexpr unsigned int energy_component         = dim + 1;
      static constexpr unsigned int first_magnetic_component = dim + 2;
      /** @} */

      using state_type = Tensor<1, n_components>;
      using flux_type  = Tensor<1, n_components, Tensor<1, dim>>;



      /** Constructor */
      MHDEquations();


      /** @{ */
      /**
       * @brief Create a list of component names, `rho`, `u_i`, `E`, `B_i`.
       *
       * @param prefix Prefix for the component names.
       * @return std::vector<std::string> component_names.
       */
      static std::vector<std::string>
      create_component_name_list(const std::string &prefix = "");
      /** @} */
    };
  } // namespace MHD
} // namespace sapphirepp
#endif
