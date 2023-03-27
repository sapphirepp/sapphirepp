#ifndef VFPEQUATION_REFERENCEVALUES_H
#define VFPEQUATION_REFERENCEVALUES_H

namespace Sapphire {
struct ReferenceValues {
  // reference values
  const double energy = 1;       // GeV
  const double mass = 0.9382721;  // GeV/c^2 (proton mass)
  const double gamma = energy / mass;
  const double velocity = 299792458;              // m/s (speed of light)
  const double magnetic_field_strength = 1.e-10;  // Tesla (1 microGauss)
  const double charge = 1.602176634e-19;          // Coulmb ( elementary charge)
  // NOTE: For the computation of the reciprocal of the gyro frequency, mass
  // needs to be converted GeV/c^2 to kg. Hence, the factor 1.e9 *
  // e/c^2
  const double time = (gamma * mass * 1.602176634e-19) /
                      (299792458. * 299792458. * charge *
                       magnetic_field_strength);  // s (1/gyro_frequency)
};

template <typename StreamType>
StreamType &operator<<(StreamType &os, const ReferenceValues &reference_values);
}  // namespace Sapphire

#endif
