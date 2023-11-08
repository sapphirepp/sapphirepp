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

#include "reference-values.h"

#include <deal.II/base/conditional_ostream.h>

#include <ostream>

template <typename StreamType>
StreamType &
Sapphire::VFP::operator<<(
  StreamType                           &os,
  const Sapphire::VFP::ReferenceValues &reference_values)
{
  os << "Reference Values: \n"
     << "	Mass: " << reference_values.mass << " kg \n"
     << "	Velocity: " << reference_values.velocity << " m/s \n"
     << "	Magnetic field strength: " << reference_values.magnetic_field_strength
     << "T \n"
     << "	Charge: " << reference_values.charge << "C \n";
  return os;
}

// explicit instantiation
template std::ostream &
Sapphire::VFP::operator<<(std::ostream &,
                          const Sapphire::VFP::ReferenceValues &);
template dealii::ConditionalOStream &
Sapphire::VFP::operator<<(dealii::ConditionalOStream &,
                          const Sapphire::VFP::ReferenceValues &);
