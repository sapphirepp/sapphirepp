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
 * @file mhd-parameters.h
 * @author Florian Schulze (florian.schulze@mpi-hd.mpg.de)
 * @brief Define @ref sapphirepp::MHD::MHDParameters
 */

#ifndef MHD_MHDPARAMETERS_H
#define MHD_MHDPARAMETERS_H

#include <deal.II/base/parameter_handler.h>

#include <string>
#include <vector>

#include "mhd-flags.h"

namespace sapphirepp
{
  namespace MHD
  {
    using namespace dealii;

    /**
     * @brief Parameters class for the MHD module
     *
     * Defines and reads all parameters that are needed for the VFP module.
     *
     * @tparam dim Total dimension of the problem in configuration space
     *         \f$ (\mathbf{x}) \f$ (`dim_cs`)
     */
    template <unsigned int dim>
    class MHDParameters
    {
    public:
      /** @{ */
      /** The @ref GridTypeMHD for the grid generation */
      GridTypeMHD grid_type;

      /**
       * First corner point for
       * @dealref{subdivided_hyper_rectangle(),namespaceGridGenerator,
       * ac76417d7404b75cf53c732f456e6e971}
       * (only for @ref GridTypeMHD::hypercube)
       */
      dealii::Point<dim> p1;
      /**
       * Second corner point opposite to p1 for
       * @dealref{subdivided_hyper_rectangle(),namespaceGridGenerator,
       * ac76417d7404b75cf53c732f456e6e971}
       * (only for @ref GridTypeMHD::hypercube)
       */
      dealii::Point<dim> p2;
      /**
       * A vector of dim positive values denoting the number of cells to
       * generate in that direction for
       * @dealref{subdivided_hyper_rectangle(),namespaceGridGenerator,
       * ac76417d7404b75cf53c732f456e6e971}
       * (only for @ref GridTypeMHD::hypercube)
       */
      std::vector<unsigned int> n_cells;

      /** File for grid creation (only for @ref GridTypeMHD::file) */
      std::string grid_file;

      /**
       * A vector of dimension `dim*2` containing the boundary indicators
       * @ref BoundaryConditionsMHD for each boundary of the hypercube
       */
      std::vector<BoundaryConditionsMHD> boundary_conditions;
      /** @} */



      /** @{ */
      /** Time stepping method */
      TimeSteppingMethodMHD time_stepping_method;
      /** Courant number/CFL number to infer time step from CFL condition */
      double courant_number;
      /** Maximum time step size for the simulation in dimensionless units */
      double max_time_step;
      /** End time for the simulation in dimensionless units */
      double final_time;
      /** @} */



      /** @{ */
      /** Polynomial degree of the DG shape functions */
      unsigned int polynomial_degree;
      /** @} */



      /** @{ */
      /** Adiabatic index \f$ \gamma \f$ */
      double adiabatic_index;
      /** @} */



      /**
       * @addtogroup numerical-parameters Numerical parameters
       * Parameters for the numerical algorithms.
       * @note !DO NOT CHANGE unless you are aware what the parameters do!
       * @{
       */
      /** Precision for double / zero comparision. */
      static constexpr double epsilon_d = 1e-8;
      /** Floor value for pressure, density and energy. */
      double mhd_floor = 10. * epsilon_d;

      /** Threshold for KXRCF shock indicator. */
      double indicator_threshold = 0.1;

      /** minmod threshold parameter \f$ M \f$ for slope limiter. */
      double minmod_threshold = 0.;
      /** minmod limiter parameter \f$ \beta \f$ for slope limiter. */
      double minmod_beta = 2.;

      /** Constant \f$ C_h \f$ for hyperbolic divergence cleaning. */
      double divergence_cleaning_Ch = 0.8;
      /** Constant \f$ C_r \f$ for hyperbolic divergence cleaning. */
      double divergence_cleaning_Cr = 0.18;
      /** @} */



      /** @brief Constructor */
      MHDParameters();

      /**
       * @brief Delcare parameters in parameter file
       *
       * @param prm @dealref{ParameterHandler}
       */
      void
      declare_parameters(ParameterHandler &prm);

      /**
       * @brief Parse parameters from parameter file
       *
       * @param prm @dealref{ParameterHandler}
       */
      void
      parse_parameters(ParameterHandler &prm);
    };
  } // namespace MHD
} // namespace sapphirepp
#endif
