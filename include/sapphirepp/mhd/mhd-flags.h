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
 * @file mhd-flags.h
 * @author Florian Schulze (florian.schulze@mpi-hd.mpg.de)
 * @brief Define @ref sapphirepp::MHD::MHDFlags and other `enums` used by the
 *        @ref sapphirepp::MHD module
 */

#ifndef MHD_MHDFLAGS_H
#define MHD_MHDFLAGS_H

namespace sapphirepp
{

  /**
   * @namespace sapphirepp::MHD
   * @brief Namespace for the magnetohydrodynamics module
   */
  namespace MHD
  {

    /**
     * @brief Flags for the MHD equation.
     */
    enum class MHDFlags
    {
      none = 0,
    };



    /**
     * @brief Implement the bitwise OR operator for the MHDFlags
     *
     * @param f1 First MHDFlags
     * @param f2 Second MHDFlags
     * @return constexpr MHDFlags f1 | f2
     */
    constexpr MHDFlags
    operator|(MHDFlags f1, MHDFlags f2)
    {
      return static_cast<MHDFlags>(static_cast<int>(f1) | static_cast<int>(f2));
    }



    /**
     * @brief Implement the bitwise AND operator for the MHDFlags
     *
     * @param f1 First MHDFlags
     * @param f2 Second MHDFlags
     * @return constexpr MHDFlags f1 & f2
     */
    constexpr MHDFlags
    operator&(MHDFlags f1, MHDFlags f2)
    {
      return static_cast<MHDFlags>(static_cast<int>(f1) & static_cast<int>(f2));
    }



    /**
     * @brief Print the MHDFlags to a stream
     *
     * @tparam StreamType Type of the output stream
     * @param os Output stream
     * @param f MHDFlags
     * @return StreamType& os
     */
    template <typename StreamType>
    inline StreamType &
    operator<<(StreamType &os, MHDFlags f)
    {
      os << "MHD flags: \n";
      if ((f & MHDFlags::none) != MHDFlags::none)
        os << "	 - None\n";
      return os;
    }



    /**
     * @brief Boundary conditions for the MHD equations
     */
    enum class BoundaryConditionsMHD
    {
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
     * @brief Time stepping methods for the MHD equation
     */
    enum class TimeSteppingMethodMHD
    {
      /**
       * Explicit Euler method.
       */
      forward_euler,

      /**
       * Second order Runge-Kutta method.
       */
      erk2,

      /**
       * Fourth order Runge-Kutta method.
       */
      erk4
    };



    /**
     * @brief Grid generation types for the MHD equations
     */
    enum class GridTypeMHD
    {
      /**
       * Create a hypercube grid.
       */
      hypercube,

      /**
       * Use a grid that is read from a file.
       */
      file
    };
  } // namespace MHD
} // namespace sapphirepp
#endif
