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
     * @brief Reference values for the VFP equation.
     *
     * @sapphire uses dimensionless units
     * \f$ x^{*}, t^{*}, p^{*}, v^{*} \f$
     * to solve the VFP equation.
     * These dimensionless are defined using reference values,
     * \f$ \underline{m}, \underline{q}, c, \underline{B} \f$,
     * in the following way:
     *
     * | Type     | Dimensionless unit                       | Reference value                                                              |
     * |:---------|:-----------------------------------------|:-----------------------------------------------------------------------------|
     * | Length   | \f$ x^{*} = x/\underline{r}_{g} \f$      | \f$ \underline{r}_{g} =  \underline{m} c/\underline{q} \underline{B} \f$     |
     * | Time     | \f$ t^{*} = t \underline{\omega}_{g} \f$ | \f$ \underline{\omega}_{g} = \underline{q} \underline{B}/  \underline{m} \f$ |
     * | Momentum | \f$ p^{*} = p/\underline{p} \f$          | \f$ \underline{p} =  \underline{m} c \f$                                     |
     * | Velocity | \f$ v^{*} = v/c \f$                      | Speed of light \f$ c \f$                                                     |
     *
     * Note the reference values are **only** used for the
     * @ref VFPFlags::radiation_reaction "radiation reaction" term.
     * If the radiation reaction term is deactivated,
     * they are **not** used when solving the VFP equation.
     * They are useful in post-processing the solution
     * to produce physical results.
     */
    struct ReferenceValues
    {
      /**
       * Reference mass \f$ \underline{m} \f$ in kilogram.
       * The default is the proton mass,
       * i.e. \f$ \underline{m} = 1.672621923 \times 10^{-27} \, \mathrm{kg}\f$.
       */
      double mass;

      /**
       * Speed of light \f$ c \f$ in meters per second.
       * The default is \f$ c = 299792458 \, \frac{\mathrm{m}}{\mathrm{s}} \f$ .
       */
      double speed_of_light;

      /**
       * Reference magnetic field strength \f$ \underline{B} \f$ in Tesla.
       * The default is 1 microGauss,
       * i.e. \f$ \underline{B} = 10^{-10} \, \mathrm{T}\f$.
       */
      double magnetic_field_strength;

      /**
       * Reference charge \f$ \underline{q} \f$ in Columb.
       * The default is the elementary charge,
       * i.e. \f$ \underline{q} = 1.602176634 \times 10^{-19} \, \mathrm{C}\f$.
       */
      double charge;

      /**
       * The reference length \f$ \underline{r}_{g} \f$ in meter.
       * The default is
       * \f$
       *    \underline{r}_{g} = \underline{m} c/\underline{q} \underline{B}
       * \f$ .
       */
      double length;

      /**
       * The reference frequency \f$ \underline{\omega}_{g} \f$ in 1/s.
       * The default is the gyro-frequency, i.e.
       * \f$
       *    \underline{\omega}_{g} = \underline{q} \underline{B}/ \underline{m}
       * \f$ .
       */
      double frequency;

      /**
       * Reference time in seconds.
       * The default is \f$ \underline{t} = 1/\underline{\omega}_{g} \f$.
       */
      double time;

      /**
       * Reference momentum \f$ \underline{p} \f$ in kg m/s.
       * The default is \f$ \underline{p} = \underline{m} c \f$.
       */
      double momentum;

      /**
       * Vacuum permeability \f$ \mu_{0} \f$ in N/A^2.
       * The default value is
       * \f$
       *    \mu_{0}
       *    = 1.25663706127 \times 10^{-6} \, \frac{\mathrm{N}}{\mathrm{A}^2}
       *    \approx 4 \pi \times 10^{-7} \, \frac{\mathrm{N}}{\mathrm{A}^2}
       * \f$.
       */
      double vacuum_permeability;

      /**
       * Radiation reaction characteristic time \f$ \underline{\tau}_R \f$
       * in seconds,
       * \f$
       *    \underline{\tau}_R
       *    = \frac{9 \pi \underline{m}^3 c \underline{\omega}_g}{
       *      \mu_{0} \underline{q}^4 \underline{B}^2}
       * \f$ .
       */
      double radiation_reaction_characteristic_time;

      /** Conversion of parsec to meter. */
      const double parsec_to_m = 3.0857e16;
    };



    /**
     * @brief Output reference values to a stream
     *
     * @tparam StreamType Type of the stream
     * @param os Output stream
     * @param reference_values Reference values
     * @return StreamType& os
     */
    template <typename StreamType>
    StreamType &
    operator<<(StreamType &os, const ReferenceValues &reference_values)
    {
      os << "Reference Values: \n"
         << "	Mass: " << reference_values.mass << " kg \n"
         << "	Speed of light: " << reference_values.speed_of_light << " m/s \n"
         << "	Magnetic field strength: "
         << reference_values.magnetic_field_strength << " T \n"
         << "	Charge: " << reference_values.charge << " C \n"
         << "	Length: " << reference_values.length
         << " m = " << reference_values.length / reference_values.parsec_to_m
         << " pc \n"
         << "	Time: " << reference_values.time << " s \n"
         << "	Frequency: " << reference_values.frequency << " 1/s \n"
         << "	Momentum: " << reference_values.momentum
         << " kg m/s = " << reference_values.momentum / reference_values.charge
         << " eV/c \n"
         << " Vacuum permeability: " << reference_values.vacuum_permeability
         << " N/A^2 \n"
         << " Radiation reaction characteristic time: "
         << reference_values.radiation_reaction_characteristic_time << " s \n";
      return os;
    }
  } // namespace VFP
} // namespace sapphirepp

#endif
