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

#ifndef VFP_REFERENCEVALUES_H
#define VFP_REFERENCEVALUES_H

namespace sapphirepp
{
  namespace VFP
  {
    struct ReferenceValues
    {
      // reference values
      const double mass     = 1.672621923e-27;       // kg (proton mass)
      const double velocity = 299792458;             // m/s (speed of light)
      const double magnetic_field_strength = 1.e-10; // Tesla (1 microGauss)
      const double charge = 1.602176634e-19; // Coulmb ( elementary charge)
    };

    template <typename StreamType>
    StreamType &
    operator<<(StreamType &os, const ReferenceValues &reference_values);
  } // namespace VFP
} // namespace sapphirepp

#endif
