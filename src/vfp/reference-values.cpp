#include "reference-values.h"

#include <deal.II/base/conditional_ostream.h>

#include <ostream>

template <typename StreamType>
StreamType &Sapphire::operator<<(
    StreamType &os, const Sapphire::ReferenceValues &reference_values) {
  os << "Reference Values: \n"
     << "	Mass: " << reference_values.mass << " kg \n"
     << "	Velocity: " << reference_values.velocity << " m/s \n"
     << "	Magnetic field strength: "
     << reference_values.magnetic_field_strength << "T \n"
     << "	Charge: " << reference_values.charge << "C \n";
  return os;
}

// explicit instantiation
template std::ostream &Sapphire::operator<<(std::ostream &,
                                            const Sapphire::ReferenceValues &);
template dealii::ConditionalOStream &Sapphire::operator<<(
    dealii::ConditionalOStream &, const Sapphire::ReferenceValues &);
