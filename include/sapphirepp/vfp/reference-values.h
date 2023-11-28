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
 * @file reference-values.h
 * @author Nils Schween (nils.schween@mpi-hd.mpg.de)
 * @author Florian Schulze (florian.schulze@mpi-hd.mpg.de)
 * @brief Define @ref sapphirepp::VFP::ReferenceValues
 */

#ifndef VFP_REFERENCEVALUES_H
#define VFP_REFERENCEVALUES_H

namespace sapphirepp
{
  namespace VFP
  {

    /**
     * @brief Reference values for the VFP equation
     *
     * @sapphire uses dimensionless units \f$ x^{*}, t^{*}, p^{*}, v^{*} \f$ to
     * solve the VFP equation. These dimensionless are defined using reference
     * values, \f$ m_{0}, q_{0}, c, B_{0} \f$, in the following way:
     *
     * | Type     | Dimensionless unit             | Reference value                            |
     * |:---------|:-------------------------------|:-------------------------------------------|
     * | Length   | \f$ x^{*} = x/r_{g,0} \f$      | \f$ r_{g,0} =  m_{0} c/q_{0} B_{0} \f$     |
     * | Time     | \f$ t^{*} = t \omega_{g,0} \f$ | \f$ \omega_{g,0} = q_{0} B_{0}/  m_{0} \f$ |
     * | Momentum | \f$ p^{*} = p/p_{0} \f$        | \f$ p_{0} =  m_{0} c \f$                   |
     * | Velocity | \f$ v^{*} = v/v_{0} \f$        | \f$ v_0 = c\f$                             |
     *
     * Note the reference values are **not** used when solving the VFP equation.
     * They are useful in post-processing the solution to produce physical
     * results.
     */
    struct ReferenceValues
    {
      /** \f$ m_{0} \f$ in kg (proton mass) */
      const double mass = 1.672621923e-27;

      /** \f$ c \f$ in m/s (speed of light) */
      const double velocity = 299792458;

      /** \f$ B_{0} \f$ in Tesla (1 microGauss) */
      const double magnetic_field_strength = 1.e-10;

      /** \f$ q_{0} \f$ in Columb (elementary charge) */
      const double charge = 1.602176634e-19;
    };



    template <typename StreamType>
    StreamType &
    operator<<(StreamType &os, const ReferenceValues &reference_values)
    {
      os << "Reference Values: \n"
         << "	Mass: " << reference_values.mass << " kg \n"
         << "	Velocity: " << reference_values.velocity << " m/s \n"
         << "	Magnetic field strength: "
         << reference_values.magnetic_field_strength << "T \n"
         << "	Charge: " << reference_values.charge << "C \n";
      return os;
    }
  } // namespace VFP
} // namespace sapphirepp

#endif
