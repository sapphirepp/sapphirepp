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

#include <ostream>
#include <string>
#include <vector>

#include "config.h"
#include "vfp-flags.h"

namespace sapphirepp
{
  namespace VFP
  {
    using namespace dealii;

    template <unsigned int dim>
    class VFPParameters
    {
    private:
      static constexpr bool momentum =
        (vfp_flags & VFPFlags::momentum) != VFPFlags::none ? true : false;
      static constexpr int dim_ps = dim;
      static constexpr int dim_cs = dim - momentum;

    public:
      VFPParameters();

      void
      declare_parameters(ParameterHandler &prm);
      void
      parse_parameters(ParameterHandler &prm);

      // Runtime settings
      // These settings are read from a parameter file
      unsigned int expansion_order;
      // Mesh
      GridType                        grid_type;
      std::string                     grid_file;
      dealii::Point<dim>              p1;
      dealii::Point<dim>              p2;
      std::vector<unsigned int>       n_cells;
      std::vector<BoundaryConditions> boundary_conditions;
      // Finite element
      unsigned int polynomial_degree;
      // Time stepping
      TimeSteppingMethod time_stepping_method;
      double             theta;
      double             time_step;
      double             final_time;
      // Particle properties
      //  NOTE: All physical quantities are dimensionless. The reference values
      //  are defined in the reference-values.h header.
      double mass;
      double charge;
      // TransportOnly
      //  In the transport-only case (i.e. no dependence on p) the energy of the
      //  particles has to be given.
      double gamma;
      double velocity;
    };
  } // namespace VFP
} // namespace sapphirepp
#endif
