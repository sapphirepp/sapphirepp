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

    /**
     * @brief Flags to activate the different terms of the VFP equation.
     *
     * We split the VFP equation into different terms,
     * \f{align}{
     *   \frac{\partial f}{\partial t} & \quad & \text{(time-evolution term)} \\
     *   & + (\mathbf{u} + \mathbf{v}) \cdot \nabla_{x} f
     *   & \text{(spatial advection term)} \\
     *   & - \gamma m \frac{\mathrm{D} \mathbf{u}}{\mathrm{D} t}
     *     \cdot \nabla_{p}f
     *     - \mathbf{p} \cdot\nabla_{x} \mathbf{u}\cdot \nabla_{p} f
     *   & \text{(momentum term)} \\
     *   & + q \mathbf{v} \cdot \left( \mathbf{B} \times \nabla_{p} f \right)
     *   & \text{(rotation term)} \\
     *   = & \frac{\nu}{2} \Delta_{\theta, \varphi} f
     *   & \text{(collision term)} \\
     *   & + S \,. & \text{(source term)} \\
     * \f}
     * These terms can be individually activated or deactivated.
     */
    enum class VFPFlags
    {
      none = 0,
      /**
       * Activate the time-evolution term
       * \f$ \frac{\partial f}{\partial t} \f$
       */
      time_evolution = 1 << 0,

      /**
       * Activate the spatial advection term
       * \f$ (\mathbf{u} + \mathbf{v}) \cdot \nabla_x f \f$
       */
      spatial_advection = 1 << 1,

      /**
       * Activate the collision term
       * \f$ \frac{\nu}{2} \Delta_{\theta, \varphi} f \f$
       */
      collision = 1 << 2,

      /**
       * Activate the rotation term, i.e. the magnetic field
       * \f$ q \mathbf{v} \cdot \left( \mathbf{B} \times \nabla_{p} f \right)
       * \f$
       */
      rotation = 1 << 3,

      /**
       * If the fields \f$ \mathbf{u} \f$ and \f$\mathbf{B} \f$ are time
       * independent, this flag should be used. \n
       * This significantly improves performance, as it requires the `dg_matrix`
       * to be assembled only once. \n
       * By default the fields are assumed to be time dependent.
       */
      time_independent_fields = 1 << 4,

      /**
       * Activate the momentum term
       * \f$ \left( \gamma m \frac{\mathrm{D} \mathbf{u}}{\mathrm{D} t} +
       * \mathbf{p} \cdot\nabla_{x} \mathbf{u} \right) \cdot \nabla_{p} f \f$
       */
      momentum = 1 << 5,

      /**
       * Use a linear momentum variable, \f$ p \f$. \n
       * By default a logarithmic momentum variable, \f$ \ln(p) \f$, is used.
       */
      linear_p = 1 << 6,

      /**
       * Activate the source term
       * \f$ S(\mathbf{x}, \mathbf{p}, t) \f$
       */
      source = 1 << 7,

      /**
       * If the source \f$ S \f$ is time independent, this flag should be used.
       * \n
       * This significantly improves performance, as it requires the
       * `system_rhs` to be assembled only once.
       * \n By default the source is assumed to be time dependent.
       */
      time_independent_source = 1 << 8,
    };



    /**
     * @brief Implement the bitwise OR operator for the VFPFlags
     *
     * @param f1 First VFPFlags
     * @param f2 Second VFPFlags
     * @return constexpr VFPFlags f1 | f2
     */
    constexpr VFPFlags
    operator|(VFPFlags f1, VFPFlags f2)
    {
      return static_cast<VFPFlags>(static_cast<int>(f1) | static_cast<int>(f2));
    }



    /**
     * @brief Implement the bitwise AND operator for the VFPFlags
     *
     * @param f1 First VFPFlags
     * @param f2 Second VFPFlags
     * @return constexpr VFPFlags f1 & f2
     */
    constexpr VFPFlags
    operator&(VFPFlags f1, VFPFlags f2)
    {
      return static_cast<VFPFlags>(static_cast<int>(f1) & static_cast<int>(f2));
    }



    /**
     * @brief Print the VFPFlags to a stream
     *
     * @tparam StreamType Type of the output stream
     * @param os Output stream
     * @param f VFPFlags
     * @return StreamType& os
     */
    template <typename StreamType>
    inline StreamType &
    operator<<(StreamType &os, VFPFlags f)
    {
      os << "VFP flags: \n";
      if ((f & VFPFlags::time_evolution) != VFPFlags::none)
        os << "	 - Time_evolution term\n";
      if ((f & VFPFlags::spatial_advection) != VFPFlags::none)
        os << "	 - Spatial Advection\n";
      if ((f & VFPFlags::collision) != VFPFlags::none)
        os << "	 - Collision\n";
      if ((f & VFPFlags::rotation) != VFPFlags::none)
        os << "	 - Rotation\n";
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



    /**
     * @brief Boundary conditions for the VFP equation
     */
    enum class BoundaryConditions
    {
      /**
       * Continuous. The distribution function at the boundary is continuously
       * extrapolated from the interior.
       */
      continuous,

      /**
       * Zero inflow. Only outflow, no inflow at the boundary.
       */
      zero_inflow,

      /**
       * Periodic boundary conditions. Has to be set on both sides of the
       * boundary.
       */
      periodic
    };



    /**
     * @brief Time stepping methods for the VFP equation
     */
    enum class TimeSteppingMethod
    {
      /**
       * Implicit Crank-Nicolson method.
       */
      crank_nicolson,

      /**
       * Explicit Euler method.
       */
      forward_euler,

      /**
       * Implicit Euler method.
       */
      backward_euler,

      /**
       * Fourth order Runge-Kutta method.
       */
      erk4,

      /**
       * Low storage fourth order Runge-Kutta method.
       */
      lserk
    };



    /**
     * @brief Grid generation types for the VFP equation
     */
    enum class GridType
    {
      /**
       * Create a hypercube grid.
       */
      hypercube,

      /**
       * Create a grid for a shock problem in x-direction.
       */
      shock,

      /**
       * Use a grid that is read from a file.
       */
      file
    };
  } // namespace VFP
} // namespace sapphirepp
#endif
