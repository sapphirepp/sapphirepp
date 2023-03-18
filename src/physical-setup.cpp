#include "physical-setup.h"

// Initial values
template <int dim>
void Sapphire::InitialValueFunction<dim>::vector_value(
    const dealii::Point<dim> &p, dealii::Vector<double> &f) const {
  Assert(dim <= 3, dealii::ExcNotImplemented());
  Assert(f.size() == (expansion_order + 1) * (expansion_order + 1),
         dealii::ExcDimensionMismatch(
             f.size(), (expansion_order + 1) * (expansion_order + 1)));
  // NOTE: The zeroth component of values corresponds to f_000, the first
  // component to f_110 etc.
  //
  // NOTE: If the momentum term is activated the last component
  // of point, i.e. p[dim] is the p coordinate.
  //
  // EXAMPLES
  // DIM = 1
  // constant function
  static_cast<void>(p);
  f[0] = 0.;
  // Gaussian
  // f[0] = 1. * std::exp(-(std::pow(p[0], 2)));

  // Constant disc
  // if (std::abs(p[0]) < 1.) f[0] = 1.;

  // space/momentum dependent distribution
  // f[0] = 1. + 0.2 * p[0];

  // sinusoidal distribution
  // f[0] = std::sin((1. * 3.14159265359) / 2 * p[0]) + 1.;

  // DIM = 2
  // constant
  // static_cast<void>(p);
  // f[0] = 0.;
  // Gaussian
  // f[0] = 1. * std::exp(-((std::pow(p[0], 2) + std::pow(p[1], 2))));

  // Constant disc
  // if (p.norm() <= 1.) f[0] = 1.;

  // Constant cross (e.g. for a rigid rotator)
  // if (std::abs(p[0]) <= 3 && std::abs(p[1]) <= 0.5) f[0] = 1.;
  // if (std::abs(p[0]) <= 0.5 && std::abs(p[1]) <= 3) f[0] = 1.;

  // Initialise f_000, f_110, f_100, f111, ... with a Gaussian
  // std::fill(
  //     f.begin(), f.end(),
  //     1. * std::exp(-((std::pow(p[0] - 0.5, 2) + std::pow(p[1] - 0.5, 2)) /
  //                     0.01)));

  // DIM 3
  // f[0] = 1. * std::exp(-((std::pow(p[0], 2) + std::pow(p[1], 2) +
  //                                std::pow(p[2], 2))));
}
// explicit instantiation
template class Sapphire::InitialValueFunction<1>;
template class Sapphire::InitialValueFunction<2>;
template class Sapphire::InitialValueFunction<3>;

template <int dim>
void Sapphire::ScatteringFrequency<dim>::value_list(
    const std::vector<dealii::Point<dim>> &points,
    std::vector<double> &scattering_frequencies,
    const unsigned int component) const {
  Assert(scattering_frequencies.size() == points.size(),
         dealii::ExcDimensionMismatch(scattering_frequencies.size(),
                                      points.size()));
  static_cast<void>(component);
  // EXAMPLES:
  // Constant scattering frequency
  std::fill(scattering_frequencies.begin(), scattering_frequencies.end(), 0.5);
}
// explicit instantiation
template class Sapphire::ScatteringFrequency<1>;
template class Sapphire::ScatteringFrequency<2>;
template class Sapphire::ScatteringFrequency<3>;

// Source term implementation
template <int dim>
void Sapphire::Source<dim>::vector_value(
    const dealii::Point<dim> &p, dealii::Vector<double> &values) const {
  Assert(values.size() == (expansion_order + 1) * (expansion_order + 1),
         dealii::ExcDimensionMismatch(
             values.size(), (expansion_order + 1) * (expansion_order + 1)));
  // The zeroth component of values corresponds to isotropic part of the
  // distribution

  // EXAMPLES:
  // Transport only
  // 1D Gaussian (isotropic)
  // values[0] = 0.01 * std::exp(-(std::pow(p[0], 2)));
  // 2D Gaussian (isotropic)
  // values[0] = .01 * std::exp(-(std::pow(p[0], 2) + std::pow(p[1], 2)));
  // 2D pulsating Gaussian (isotropic, time dependent)
  values[0] = 0.01 * (std::sin(this->get_time()) + 1.) *
              std::exp(-(std::pow(p[0], 2) + std::pow(p[1], 2)));
}
// explicit instantiation
template class Sapphire::Source<1>;
template class Sapphire::Source<2>;
template class Sapphire::Source<3>;

// Magnetic field implementation
template <int dim>
void Sapphire::MagneticField<dim>::vector_value(
    const dealii::Point<dim> &point,
    dealii::Vector<double> &magnetic_field) const {
  Assert(magnetic_field.size() == 3,
         dealii::ExcDimensionMismatch(magnetic_field.size(), 3));

  // EXAMPLES:
  // constant magnetic field
  static_cast<void>(point);  // suppress compiler warning
  magnetic_field[0] = 0.;
  magnetic_field[1] = 0.;
  magnetic_field[2] = 3.;
}

// explicit instantiation
template class Sapphire::MagneticField<1>;
template class Sapphire::MagneticField<2>;
template class Sapphire::MagneticField<3>;

// Velocity field implementation
template <int dim>
void Sapphire::BackgroundVelocityField<dim>::vector_value(
    const dealii::Point<dim> &point, dealii::Vector<double> &velocity) const {
  Assert(velocity.size() == 3,
         dealii::ExcDimensionMismatch(velocity.size(), 3));
  // EXAMPLES:
  // constant velocity field
  static_cast<void>(point);  // suppress compiler warning
  velocity[0] = .0;          // u_x
  velocity[1] = .0;          // u_y
  velocity[2] = .0;          // u_z

  // space dependent velocity
  // u_x = 1/5 * x
  // velocity[0] = 1. / 5 * point[0];
  // velocity[1] = .0;
  // velocity[2] = .0;

  // u_x = 0.1 * cos(pi/2 * x)
  // velocity[0] = 0.1 * std::cos(pi / 2 * point[0]);
  // velocity[1] = .0;
  // velocity[2] = .0;

  // time-dependent velocity-field
  // u_x = 1/10 * t
  // velocity[0] = 1. / 10 * this->get_time();
  // velocity[1] = .0;
  // velocity[2] = .0;

  // u_x = 1/2 * sin(t)
  // velocity[0] = 1. / 2 * std::sin(this->get_time());
  // velocity[1] = .0;
  // velocity[2] = .0;

  // time- and space dependent velocity
  // u_x = -1/10 * sin(t) * cos(pi/2 * x)
  // velocity[0] = -0.1 * std::sin(this->get_time()) * std::cos(pi / 2 * point[0]);
  // velocity[1] = .0;
  // velocity[2] = .0;

  // u_x = 1/10 * x * (1 - e^(-t))
  // velocity[0] = 0.1 * point[0] * (1 - std::exp(-this->get_time()));
  // velocity[1] = .0;
  // velocity[2] = .0;

  // rigid rotator
  // velocity[0] = -point[1];
  // velocity[1] = point[0];
  // velocity[2] = 0.;
}

// Divergence
template <int dim>
void Sapphire::BackgroundVelocityField<dim>::divergence_list(
    const std::vector<dealii::Point<dim>> &points,
    std::vector<double> &divergence) {
  Assert(divergence.size() == points.size(),
         dealii::ExcDimensionMismatch(divergence.size(), points.size()));
  // EXAMPLES:
  // constant velocity
  std::fill(divergence.begin(), divergence.end(), 0.);

  // space-dependent velocity field
  // u_x = 1/5 * x => div u = 1/5
  // std::fill(divergence.begin(), divergence.end(), 1. / 5);

  // u_x = 0.1 * cos(pi/2 * x) => div u = -0.1 * pi/2 * sin(pi/2 * x)
  // for (unsigned int i = 0; i < points.size(); ++i)
  //   divergence[i] = -0.1 * pi / 2 * std::sin(pi / 2 * points[i][0]);

  // time-dependent velocity field
  // u_x = 1/10 * t => div u = 0
  // std::fill(divergence.begin(), divergence.end(), 0.);

  // u_x = 1/2 * sin(t) => div u = 0
  // std::fill(divergence.begin(), divergence.end(), 0.);

  // time- and space-dependent velocity field
  // u_x = -1/10 * sin(t) * cos(pi/2 * x) =>
  // div u = 0.1 * pi/2 * sin(t) * sin(pi/2 * x)
  // for (unsigned int i = 0; i < points.size(); ++i)
  //   divergence[i] = 0.1 * pi / 2 * std::sin(this->get_time()) *
  //                   std::sin(pi / 2 * points[i][0]);

  // u_x = 1/10 * x * (1 - e^(-t))  => div u = 0.1 (1  - e^(-t))
  // std::fill(divergence.begin(), divergence.end(),
  //           0.1 * (1 - std::exp(-this->get_time())));
}

// Material derivative
template <int dim>
void Sapphire::BackgroundVelocityField<dim>::material_derivative_list(
    const std::vector<dealii::Point<dim>> &points,
    std::vector<dealii::Vector<double>> &material_derivatives) {
  Assert(material_derivatives[0].size() == 3,
         dealii::ExcDimensionMismatch(material_derivatives[0].size(), 3));
  Assert(
      material_derivatives.size() == points.size(),
      dealii::ExcDimensionMismatch(material_derivatives.size(), points.size()));
  for (unsigned int i = 0; i < points.size(); ++i) {
    // EXAMPLES:
    // constant velocity
    material_derivatives[i][0] = 0.;  // D/Dt u_x
    material_derivatives[i][1] = 0.;  // D/Dt u_y
    material_derivatives[i][2] = 0.;  // D/Dt u_z

    // space dependent velocity field
    // u_x = 1/5 * x => D/Dt u_x = 1/25 * x
    // material_derivatives[i][0] = 1. / 25 * points[i][0];
    // material_derivatives[i][1] = 0.;
    // material_derivatives[i][2] = 0.;

    // u_x = 0.1 * cos(pi/2 * x)
    // => D/Dt u_x = 1/100 * pi/2 * cos(pi/2*x) * sin(pi/2 * x)
    // material_derivatives[i][0] = 1. / 100 * pi / 2 *
    //                              std::sin(pi / 2 * points[i][0]) *
    //                              std::cos(pi / 2 * points[i][0]);
    // material_derivatives[i][1] = 0.;
    // material_derivatives[i][2] = 0.;

    // time-dependent velocity field
    // u_x = 1/10 * t => D/Dt u_x = 1/10
    // material_derivatives[i][0] = 1. / 10;
    // material_derivatives[i][1] = 0.;
    // material_derivatives[i][2] = 0.;

    // u_x = 1/2 * sin(t) => D/Dt u_x = 1/2 * cos(t)
    // material_derivatives[i][0] = 1. / 2 * std::cos(this->get_time());
    // material_derivatives[i][1] = 0.;
    // material_derivatives[i][2] = 0.;

    // time- and space-dependent velocity field
    // u_x = -1/10 * sin(t) * cos(pi/2 * x)
    // => D/Dt u_x = -0.1 * cos(pi/2 * x) * (cos(t) + 0.1 * pi/2 * sin(t) *
    // sin(pi/2 * x))
    // material_derivatives[i][0] =
    //     -0.1 * std::cos((pi / 2 * points[i][0])) *
    //     (std::cos(this->get_time()) +
    //      0.1 * pi / 2 * std::sin(this->get_time()) *
    //          std::sin(this->get_time()) * std::sin(pi / 2 * points[i][0]));
    // material_derivatives[i][1] = 0.;
    // material_derivatives[i][2] = 0.;

    // u_x = 1/10 * x * (1 - e^(-t))  => D/Dt u_x = 0.1 * x * (e^(-t)  + (1 -
    // e^(-t))^2)
    // material_derivatives[i][0] = 0.1 * points[i][0] *
    //                              (std::exp(-this->get_time()) +
    //                               0.1 * (1 - std::exp(-this->get_time())) *
    //                                   (1 - std::exp(-this->get_time())));
    // material_derivatives[i][1] = 0.;
    // material_derivatives[i][2] = 0.;
  }
}

// Jacobian matrix
template <int dim>
void Sapphire::BackgroundVelocityField<dim>::jacobian_list(
    const std::vector<dealii::Point<dim>> &points,
    std::vector<std::vector<dealii::Vector<double>>> &jacobians) const {
  Assert(jacobians.size() == points.size(),
         dealii::ExcDimensionMismatch(jacobians.size(), points.size()));
  Assert(jacobians[0].size() == 3,
         dealii::ExcDimensionMismatch(jacobians[0].size(), 3));
  for (unsigned int i = 0; i < points.size(); ++i) {
    // EXAMPLES:
    // constant velocity field
    jacobians[i][0][0] = 0.;  // \partial u_x / \partial x
    jacobians[i][0][1] = 0.;  // \partial u_x / \partial y
    jacobians[i][0][2] = 0.;  // \partial u_x / \partial z

    jacobians[i][1][0] = 0.;  // \partial u_y / \partial x
    jacobians[i][1][1] = 0.;  // \partial u_y / \partial y
    jacobians[i][1][2] = 0.;  // \partial u_y / \partial z

    jacobians[i][2][0] = 0.;  // \partial u_z / \partial x
    jacobians[i][2][1] = 0.;  // \partial u_z / \partial y
    jacobians[i][2][2] = 0.;  // \partial u_z / \partial z

    // space dependent velocity field
    // u_x = 1/5 * x
    // jacobians[i][0][0] = 1. / 5;
    // jacobians[i][0][1] = 0.;
    // jacobians[i][0][2] = 0.;

    // jacobians[i][1][0] = 0.;
    // jacobians[i][1][1] = 0.;
    // jacobians[i][1][2] = 0.;

    // jacobians[i][2][0] = 0.;
    // jacobians[i][2][1] = 0.;
    // jacobians[i][2][2] = 0.;

    // u_x = 0.1 * cos(pi/2 * x)
    // jacobians[i][0][0] = 0.1 * pi/2 * std::sin(pi/2 * points[i][0]);
    // jacobians[i][0][1] = 0.;
    // jacobians[i][0][2] = 0.;

    // jacobians[i][1][0] = 0.;
    // jacobians[i][1][1] = 0.;
    // jacobians[i][1][2] = 0.;

    // jacobians[i][2][0] = 0.;
    // jacobians[i][2][1] = 0.;
    // jacobians[i][2][2] = 0.;

    // time-dependent velocity field
    // u_x = 1/10 * t and u_x = 1/2 * sin(t)
    // jacobians[i][0][0] = 0.;
    // jacobians[i][0][1] = 0.;
    // jacobians[i][0][2] = 0.;

    // jacobians[i][1][0] = 0.;
    // jacobians[i][1][1] = 0.;
    // jacobians[i][1][2] = 0.;

    // jacobians[i][2][0] = 0.;
    // jacobians[i][2][1] = 0.;
    // jacobians[i][2][2] = 0.;

    // time- and space dependent velocity field
    // u_x = -1/10 * sin(t) * cos(pi/2 * x)
    // jacobians[i][0][0] = 0.1 * pi / 2 * std::sin(this->get_time()) *
    //                      std::sin(pi / 2 * points[i][0]);
    // jacobians[i][0][1] = 0.;
    // jacobians[i][0][2] = 0.;

    // jacobians[i][1][0] = 0.;
    // jacobians[i][1][1] = 0.;
    // jacobians[i][1][2] = 0.;

    // jacobians[i][2][0] = 0.;
    // jacobians[i][2][1] = 0.;
    // jacobians[i][2][2] = 0.;

    // u_x = 1/10 * x * (1 - e^(-t))
    // jacobians[i][0][0] = 0.1 * (1 - std::exp(-this->get_time()));x
    // jacobians[i][0][1] = 0.;
    // jacobians[i][0][2] = 0.;

    // jacobians[i][1][0] = 0.;
    // jacobians[i][1][1] = 0.;
    // jacobians[i][1][2] = 0.;

    // jacobians[i][2][0] = 0.;
    // jacobians[i][2][1] = 0.;
    // jacobians[i][2][2] = 0.;
  }
}

// explicit instantiation
template class Sapphire::BackgroundVelocityField<1>;
template class Sapphire::BackgroundVelocityField<2>;
template class Sapphire::BackgroundVelocityField<3>;
