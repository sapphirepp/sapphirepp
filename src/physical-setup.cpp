#include "physical-setup.h"

// Magnetic field implementation
template <int dim>
void VFPEquation::MagneticField<dim>::vector_value(
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
template class VFPEquation::MagneticField<1>;
template class VFPEquation::MagneticField<2>;
template class VFPEquation::MagneticField<3>;

// Velocity field implementation
template <int dim>
void VFPEquation::BackgroundVelocityField<dim>::vector_value(
    const dealii::Point<dim> &point, dealii::Vector<double> &velocity) const {
  Assert(velocity.size() == 3,
         dealii::ExcDimensionMismatch(velocity.size(), 3));
  // EXAMPLES:
  // constant velocity field
  static_cast<void>(point);  // suppress compiler warning
  // velocity[0] = .4;          // u_x
  // velocity[1] = .2;          // u_y
  // velocity[2] = .0;          // u_z

  // space dependent velocity
  // u_x = 1/5 * x
  velocity[0] = 1. / 5 * point[0];
  velocity[1] = .0;
  velocity[2] = .0;

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
  // velocity[0] =
  //     -0.1 * std::sin(this->get_time()) * std::cos(pi / 2 * point[0]);
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
void VFPEquation::BackgroundVelocityField<dim>::divergence_list(
    const std::vector<dealii::Point<dim>> &points,
    std::vector<double> &divergence) {
  Assert(divergence.size() == points.size(),
         dealii::ExcDimensionMismatch(divergence.size(), points.size()));
  // EXAMPLES:
  // constant velocity
  // std::fill(divergence.begin(), divergence.end(), 0.);

  // space-dependent velocity field
  // u_x = 1/5 * x => div u = 1/5
  std::fill(divergence.begin(), divergence.end(), 1. / 5);

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
  //               std::sin(pi / 2 * points[i][0]);

  // u_x = 1/10 * x * (1 - e^(-t))  => div u = 0.1 (1  - e^(-t))
  // std::fill(divergence.begin(), divergence.end(),
  //           0.1 * (1 - std::exp(-this->get_time())));
}

// Material derivative
template <int dim>
void VFPEquation::BackgroundVelocityField<dim>::material_derivative_list(
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
    // material_derivatives[i][0] = 0.;  // D/Dt u_x
    // material_derivatives[i][1] = 0.;  // D/Dt u_y
    // material_derivatives[i][2] = 0.;  // D/Dt u_z

    // space dependent velocity field
    // u_x = 1/5 * x => D/Dt u_x = 1/25 * x
    material_derivatives[i][0] = 1. / 25 * points[i][0];
    material_derivatives[i][1] = 0.;
    material_derivatives[i][2] = 0.;

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
void VFPEquation::BackgroundVelocityField<dim>::jacobian_list(
    const std::vector<dealii::Point<dim>> &points,
    std::vector<std::vector<dealii::Vector<double>>> &jacobians) const {
  Assert(jacobians.size() == points.size(),
         dealii::ExcDimensionMismatch(jacobians.size(), points.size()));
  Assert(jacobians[0].size() == 3,
         dealii::ExcDimensionMismatch(jacobians[0].size(), 3));
  for (unsigned int i = 0; i < points.size(); ++i) {
    // EXAMPLES:
    // constant velocity field
    // jacobians[i][0][0] = 0.;  // \partial u_x / \partial x
    // jacobians[i][0][1] = 0.;  // \partial u_x / \partial y
    // jacobians[i][0][2] = 0.;  // \partial u_x / \partial z

    // jacobians[i][1][0] = 0.;  // \partial u_y / \partial x
    // jacobians[i][1][1] = 0.;  // \partial u_y / \partial y
    // jacobians[i][1][2] = 0.;  // \partial u_y / \partial z

    // jacobians[i][2][0] = 0.;  // \partial u_z / \partial x
    // jacobians[i][2][1] = 0.;  // \partial u_z / \partial y
    // jacobians[i][2][2] = 0.;  // \partial u_z / \partial z

    // space dependent velocity field
    // u_x = 1/5 * x
    jacobians[i][0][0] = 1. / 5;
    jacobians[i][0][1] = 0.;
    jacobians[i][0][2] = 0.;

    jacobians[i][1][0] = 0.;
    jacobians[i][1][1] = 0.;
    jacobians[i][1][2] = 0.;

    jacobians[i][2][0] = 0.;
    jacobians[i][2][1] = 0.;
    jacobians[i][2][2] = 0.;

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
template class VFPEquation::BackgroundVelocityField<1>;
template class VFPEquation::BackgroundVelocityField<2>;
template class VFPEquation::BackgroundVelocityField<3>;
