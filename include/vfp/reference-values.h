#ifndef VFPEQUATION_REFERENCEVALUES_H
#define VFPEQUATION_REFERENCEVALUES_H

namespace Sapphire {
struct ReferenceValues {
  // reference values
  const double mass = 1.672621923e-27;            // kg (proton mass)
  const double velocity = 299792458;              // m/s (speed of light)
  const double magnetic_field_strength = 1.e-10;  // Tesla (1 microGauss)
  const double charge = 1.602176634e-19;          // Coulmb ( elementary charge)
};

template <typename StreamType>
StreamType &operator<<(StreamType &os, const ReferenceValues &reference_values);
}  // namespace Sapphire

#endif
