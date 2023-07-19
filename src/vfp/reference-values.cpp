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
