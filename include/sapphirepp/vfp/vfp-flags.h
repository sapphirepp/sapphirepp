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
 * @file vfp-flags.h
 * @author Nils Schween (nils.schween@mpi-hd.mpg.de),
 * @author Florian Schulze (florian.schulze@mpi-hd.mpg.de)
 * @brief Define @ref sapphirepp::VFP::VFPFlags and other `enums` used by the
 *        @ref sapphirepp::VFP module
 */

#ifndef VFP_VFPFLAGS_H
#define VFP_VFPFLAGS_H

namespace sapphirepp
{

  /**
   * @namespace sapphirepp::VFP
   * @brief Namespace for the Vlasov-Fokker-Planck module
   */
  namespace VFP
  {
    enum class VFPFlags
    {
      none = 0,

      /// Activate the spatial advection term
      /// \f$ (\mathbf{u} + \mathbf{v}) \cdot \nabla_x f \f$
      spatial_advection = 1 << 0,

      /// Activate the collision term
      /// \f$ \frac{\nu}{2} \Delta_{\theta, \varphi} f \f$
      collision = 1 << 1,

      /// Activate the magnetic field
      /// \f$ q \mathbf{v} \cdot \left( \mathbf{B} \times \nabla_{p} f \right)
      /// \f$
      magnetic = 1 << 2,

      /// If the fields \f$ \mathbf{u} \f$ and \f$\mathbf{B} \f$ are time
      /// independent, this flag should be used. \n
      /// It results in a big performance gain, because the `dg_matrix` only has
      /// to be assembled once. \n
      /// By default the fields are assumed to be time dependent.
      time_independent_fields = 1 << 3,

      /// Activate the momentum term
      /// \f$ \left( \gamma m \frac{\mathrm{D} \mathbf{u}}{\mathrm{D} t} +
      /// \mathbf{p} \cdot\nabla_{x} \mathbf{u} \right) \cdot \nabla_{p} f \f$
      momentum = 1 << 4,

      /// Use a linear momentum variable, \f$ p \f$. \n
      /// By default a logarithmic momentum variable, \f$ \ln(p) \f$, is used.
      linear_p = 1 << 5,

      /// Activate the source term
      /// \f$ S(\mathbf{x}, \mathbf{p}, t) \f$
      source = 1 << 6,

      /// If the source \f$ S \f$ is time independent, this flag should be used.
      /// It results in a big performance gain, because the `system_rhs` only
      /// has to be assembled once. \n
      /// By default the source is assumed to be time dependent.
      time_independent_source = 1 << 7
    };

    constexpr VFPFlags
    operator|(VFPFlags f1, VFPFlags f2)
    {
      return static_cast<VFPFlags>(static_cast<int>(f1) | static_cast<int>(f2));
    }

    constexpr VFPFlags
    operator&(VFPFlags f1, VFPFlags f2)
    {
      return static_cast<VFPFlags>(static_cast<int>(f1) & static_cast<int>(f2));
    }

    template <typename StreamType>
    inline StreamType &
    operator<<(StreamType &os, VFPFlags f)
    {
      os << "VFP flags: \n";
      if ((f & VFPFlags::spatial_advection) != VFPFlags::none)
        os << "	 - Spatial Advection\n";
      if ((f & VFPFlags::collision) != VFPFlags::none)
        os << "	 - Collision\n";
      if ((f & VFPFlags::magnetic) != VFPFlags::none)
        os << "	 - Magnetic\n";
      if ((f & VFPFlags::time_independent_fields) != VFPFlags::none)
        os << "	 - Time Independent Fields\n";
      else
        os << "	 - Time Dependent Fields\n";
      if ((f & VFPFlags::momentum) != VFPFlags::none)
        {
          os << "	 - Momentum";
          if ((f & VFPFlags::linear_p) != VFPFlags::none)
            os << "	(linear)\n";
          else
            os << "	(logarithmic)\n";
        }
      if ((f & VFPFlags::source) != VFPFlags::none)
        {
          os << "	 - Source";
          if ((f & VFPFlags::time_independent_source) != VFPFlags::none)
            os << "	(time independent)\n";
          else
            os << " (time dependent)\n";
        }
      return os;
    }

    enum class BoundaryConditions
    {
      continuous_gradients,
      zero_inflow,
      periodic
    };

    enum class TimeSteppingMethod
    {
      crank_nicolson,
      forward_euler,
      backward_euler,
      erk4,
      lserk
    };

    enum class GridType
    {
      hypercube,
      file
    };
  } // namespace VFP
} // namespace sapphirepp
#endif
