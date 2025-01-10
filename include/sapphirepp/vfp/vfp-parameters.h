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
 * @file vfp-parameters.h
 * @author Nils Schween (nils.schween@mpi-hd.mpg.de)
 * @author Florian Schulze (florian.schulze@mpi-hd.mpg.de)
 * @brief Define @ref sapphirepp::VFP::VFPParameters
 */

#ifndef VFP_PARAMETERS_H
#define VFP_PARAMETERS_H

#include <deal.II/base/parameter_handler.h>

#include <string>
#include <vector>

#include "config.h"
#include "vfp-flags.h"

namespace sapphirepp
{
  namespace VFP
  {
    using namespace dealii;

    /**
     * @brief Parameters class for the VFP module
     *
     * Defines and reads all parameters that are needed for the VFP module.
     *
     * @tparam dim Total dimension of the problem in reduced phase space \f$
     *         (\mathbf{x}, p) \f$ (`dim_ps`)
     */
    template <unsigned int dim>
    class VFPParameters
    {
    private:
      /** Is the momentum term activated? */
      static constexpr bool momentum =
        (vfp_flags & VFPFlags::momentum) != VFPFlags::none ? true : false;
      /** Dimension in reduced phase space */
      static constexpr int dim_ps = dim;
      /** Dimension of the configuration space */
      static constexpr int dim_cs = dim - momentum;



    public:
      /** @{ */
      /** The @ref GridType for the grid generation */
      GridType grid_type;

      /**
       * First corner point for
       * @dealref{subdivided_hyper_rectangle(),namespaceGridGenerator,
       * ac76417d7404b75cf53c732f456e6e971}
       * (only for @ref GridType::hypercube)
       */
      dealii::Point<dim> p1;
      /**
       * Second corner point opposite to p1 for
       * @dealref{subdivided_hyper_rectangle(),namespaceGridGenerator,
       * ac76417d7404b75cf53c732f456e6e971}
       * (only for @ref GridType::hypercube)
       */
      dealii::Point<dim> p2;
      /**
       * A vector of dim positive values denoting the number of cells to
       * generate in that direction for
       * @dealref{subdivided_hyper_rectangle(),namespaceGridGenerator,
       * ac76417d7404b75cf53c732f456e6e971}
       * (only for @ref GridType::hypercube)
       */
      std::vector<unsigned int> n_cells;

      /**
       * Number of cells in the shock region (only for @ref GridType::shock)
       */
      unsigned int n_cells_shock;
      /**
       * Width of the shock region (only for @ref GridType::shock)
       */
      double shock_width;
      /**
       * Scaling factor for the cell width in x-direction outside the shock
       * region (only for @ref GridType::shock)
       */
      double scaling_factor_shock;

      /** File for grid creation (only for @ref GridType::file) */
      std::string grid_file;

      /**
       * A vector of dimension `dim*2` containing the boundary indicators
       * @ref BoundaryConditions for each boundary of the hypercube
       */
      std::vector<BoundaryConditions> boundary_conditions;
      /** @} */


      /**  @{ */
      /** Time stepping method */
      TimeSteppingMethod time_stepping_method;
      /** Time step size for the simulation in dimensionless units */
      double time_step;
      /** End time for the simulation in dimensionless units */
      double final_time;
      /** Crank-Nicolson parameter \f$ \theta \f$ */
      double theta;
      /** @} */


      /** @{ */
      /** Expansion order \f$ l_{\rm max} \f$ */
      unsigned int expansion_order;
      /** @} */


      /** @{ */
      /** Polynomial degree of the DG shape functions */
      unsigned int polynomial_degree;
      /** @} */


      /** @{ */
      /**
       * Mass of the particles in dimensionless units. The reference values are
       * defined in @ref VFP::ReferenceValues.
       */
      double mass;
      /**
       * Charge of the particles in dimensionless units. The reference values
       * are defined in @ref VFP::ReferenceValues.
       */
      double charge;
      /** @} */


      /** @{ */
      /**
       * Lorentz factor of the particles.
       * Only specified in the transport-only (i.e. p-independent) case.
       */
      double gamma;
      /**
       * Velocity of the particles in dimensionless units.
       * Only specified in the transport-only (i.e. p-independent) case.
       */
      double velocity;
      /** @} */


      /** @{ */
      /**
       * Postprocess to preform phase space reconstruction?
       */
      bool perform_phase_space_reconstruction;
      /**
       * Points to perform phase space reconstruction
       */
      std::vector<dealii::Point<dim>> reconstruction_points;
      /**
       * Number of theta points for phase space reconstruction
       */
      unsigned int n_theta;
      /**
       * Number of phi points for phase space reconstruction
       */
      unsigned int n_phi;
      /** @} */



      /** @brief Constructor */
      VFPParameters();

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
  } // namespace VFP
} // namespace sapphirepp
#endif
