#include "reference-values.h"

#include <ostream>

struct ReferenceValues {
  // reference values
  const double energy = 10;       // GeV
  const double mass = 0.9382721;  // GeV/c^2 (proton mass)
  const double gamma = energy / mass;
  const double velocity = 299792458;              // m/s (speed of light)
  const double magnetic_field_strength = 1.e-10;  // Tesla (1 microGauss)
  const double charge = 1.602176634e-19;          // Coulmb ( elementary charge)
  const double time =
      (gamma * mass) /
      (charge * magnetic_field_strength);  // s (1/gyro_frequency)
};

std::ostream &operator<<(std::ostream &os,
                         const ReferenceValues &reference_values) {
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
