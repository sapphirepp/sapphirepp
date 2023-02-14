#include "reference-values.h"

#include <deal.II/base/conditional_ostream.h>

#include <ostream>

template <typename StreamType>
StreamType &VFPEquation::operator<<(
    StreamType &os, const VFPEquation::ReferenceValues &reference_values) {
  os << "Reference Values: \n"
     << "	Energy: " << reference_values.energy << " GeV \n"
     << "	Mass: " << reference_values.mass << " GeV/c^2 \n"
     << "	Gamma: " << reference_values.gamma << "\n"
     << "	Velocity: " << reference_values.velocity << " m/s \n"
     << "	Magnetic field strength: "
     << reference_values.magnetic_field_strength << "T \n"
     << "	Charge: " << reference_values.charge << "C \n"
     << "	Time: " << reference_values.time << "s \n\n";
  return os;
}

// explicit instantiation
template std::ostream &VFPEquation::operator<<(
    std::ostream &, const VFPEquation::ReferenceValues &);
template dealii::ConditionalOStream &VFPEquation::operator<<(
    dealii::ConditionalOStream &, const VFPEquation::ReferenceValues &);
