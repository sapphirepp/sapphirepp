<p align="center">
<img
src="https://git.mpi-hd.mpg.de/schween/sapphire/media/branch/main/logo/sapphire-logo-text.png"
width=563 height=300 > 
</p>

# About
Sapphire is an acronym and stands for "**S**imulating **a**strophysical
**p**lasmas and **p**articles with **hi**gly **r**elativistic **e**nergies". It
is made to simulate the interaction of particles with a background plasma. To
this end it solves a Vlasov-Fokker-Planck equation in mixed coordinates, namely

\[
  \frac{\partial f}{\partial t} + (\mathbf{u} + \mathbf{v}) \cdot \nabla_{x}
  f - \gamma m \frac{\mathrm{D} \mathbf{u}}{\mathrm{D} t} \cdot \nabla_{p}f
  - \mathbf{p} \cdot\nabla_{x} \mathbf{u}\cdot \nabla_{p} f
  + q \mathbf{v} \cdot \left( \mathbf{B} \times \nabla_{p} f \right) =
  \frac{\nu}{2} \Delta_{\theta, \varphi} f + S(\mathbf{x}, \mathbf{p}, t)\, .
\]

$f$ is a number density of particles in phase space. The above equation models
how this distribution evolves in time assuming that the particles interact with
the a background plasma, which is determined by means of its velocity
$\mathbf{u}$ and its magnetic field $\mathbf{B}$. The interaction is determined
by the right-hand side of the equation: It models particles whose direction
diffuses.

To solve the above equation, we expand $f$ using real spherical harmonics, namely
\[
 f(\mathbf{x}, \mathbf{p}) = \sum^{\infty}_{l = 0} \sum^{l}_{m = 0} \sum^{1}_{s
 = 0} f_{lms}(\mathbf{x}, p ,t) Y_{lms}(\theta,\varphi) \,.
\]

This turns the original equation into a system of equations determining the
coefficients $f_{lms}$ of the expansion. And they are the output of the
simulation.

The system of equations is solved applying a finite element method, namely the
discontinuous Galerkin method. Sapphire is implemented with the help of the
finite element library [deal.ii](https://dealii.org/).

# Compiling 
Clone the repository

```shell
git clone https://git.mpi-hd.mpg.de/schween/sapphire.git
```

## Personal Computer
The compilation of Sapphire requires an installed version of deal.ii with
[PETSc](https://petsc.org) and [p4est](https://p4est.org) support. Please read
the deal.ii's [readme](https://www.dealii.org/developer/readme.html) for further
details. After having installed deal.ii, the path in `CMakeLists.txt` to the
deal.ii installation needs to be adapted, i.e. change

```cmake
set(DEAL_II_DIR "/path/to/dealii")
```

After having set the `DEAL_II_DIR`, sapphire can be compiled with the following
commands

```shell
cmake -S . -B build && make -C build
``` 

compiles it.

## MPIK cluster

On the MPIK cluster deal.ii with its dependencies is installed for all members of
APT under

```
/lfs/l8/tap/tapshared/lib/dealii
```

To compile sapphire the following set of modules needs to be loaded

```
module load python
module load cmake
module load gcc/12.1.0
module load openmpi/4.1.4
module load boost/1.81.0
module load hdf5-parallel/1.12.2
```

and the following set of variables needs to be set

```
export PETSC_DIR=/lfs/l8/tap/tapshared/lib/petsc
export PETSC_ARCH=
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/lfs/l8/tap/tapshared/lib/zstd/lib
```

After having set the `DEAL_II_DIR` to the above path, loaded the modules and set
the above variables, call the above cmake command to compile sapphire.

Note that the modules need to be loaded and the variables need to be set in your
qsub script as well.

# Usage
A user of sapphire has many options. They can be categorised into compile time
options, run-time options and physical input. The compile time options change
the VFP equation. The run time options determine the computational domain, the
time stepping method, the order of the spherical harmonic expansion, the
polynomial degree of the shape functions of the finite elements and where to
store the output of simulation runs. The physical input allows a user to set the
reference values for the dimensionless units, to prescribe the initial
distribution of particles, the background velocity field \(\mathbf{u}\), the
magnetic field \(\mathbf{B}\), the scattering frequency \(\nu\) and the source
term \(S\).

## Compile time options
Sapphire is able to solve various "versions" of the above Vlasov-Fokker-Planck
equation. "Version" means that the user can decide, which terms of the equation
will be included in his/her simulation. This decision is made before sapphire is
compiled. The corresponding options are found in `include/vfp-solver-control.h`. 
The options are:
- `terms`: This option determines which terms are included in the simulation.
Possible options are:
    - `spatial_advection` and the related term is \((\mathbf{u} +\mathbf{v})\cdot\nabla_{x}f\). Not including the spatial advection term
      is tantamount to setting the dimension of the configuration space to zero,
      i.e. the expansion coefficients do not depend on \(\mathbf{x}\).
    - `momentum` and the related terms are \(-\gamma m \mathrm{D}\mathbf{u}/\mathrm{D} t \cdot \nabla_{p} - \mathbf{p} \cdot  \nabla_{x} \mathbf{u} \cdot \nabla_{p}f\). Not  including the momentum
      terms is tantamount to expansion coefficients which are independent of the
      magnitude of momentum, i.e \(p\).
    - `magnetic` and the related term is \(q\mathbf{v} \cdot (\mathbf{B} \times \nabla_{p} f)\).
    - `collision` and the related term is \(\nu/2 \Delta_{\theta, \varphi} f\).
	- `source` and the related term is \(S(\mathbf{x}, \mathbf{p}, t)\)

  An example setup is
  ```cpp
  static constexpr TermFlags terms = 
   TermFlags::spatial_advection | TermFlags::magnetic | TermFlags::collision
  ```
  Note the operator `|` means `and`.
 - `dim_configuration_space`: This option decides if \(\mathbf{x}\) has one, two
   or three components, i.e. the expansion coefficients either dependent on
   \(x\), \(x\) and \(y\) or \(x, y \) and \(z\). Note that the total dimension
   of the problem must not exceed three. In particular, if the momentum terms
   are included its maximum value is two.
 - `logarithmic_p`: If this option is `true`, \(\ln p\) instead of \(p\) is used.
   Note this option only changes the program, if the momentum terms are included.
 - `time_dependent_fields`: If the background velocity \(\mathbf{u}\) or the
   magnetic field \(\mathbf{B}\) are time dependent, this option must be set to
   `true`. Note that setting it to `true` incurs very high computational costs.
   Only set it to `true`, if necessary.
 - `time_dependent_source`: If the source term \(S\) is time
   dependent it must be set to `true`.

## Run-time options
Once the version of the VFP equation which is going to be solved is determined
by the compile time options, it is possible to change some parameters for each
start of sapphire. These parameters concern the triangulation of the simulation
domain, the time stepping method, the polynomial degree of the shape
functions of the finite elements, the expansion order of the spherical
harmonic expansion and the output of the simulations. They are stored in a file
whose name is `vfp-equation.prm`. Note there is a template file called
`parameter-template.prm`. Before the first start of sapphire, this file must be
copied (or renamed) to `vfp-equation.prm`.

- The `Mesh` (or triangulation) section offers the following options:
    - `Point 1` and `Point 2` are two diagonally opposite points which determine the
      rectangular domain of the simulation. If the total dimension of the problem is
      less then three, only the first values are considered.
    - `Number of cells` determines in how many cells the domain will be divided in each
      of the coordinate direction. For a two dimensional problem, setting it to
      `16,16,16`, will yield a grid consisting of \(16^2 = 256\) cells.
    - `Periodicity` decides if periodic boundary conditions are applied and in which
      direction. For example, `true, false, false` yields a problem which is
      periodic in the \(x\)-direction.
  
- The `Time Stepping` section offers the following options:
    - `Method` allows a user to decide which method to use for time stepping. The
      possible values are
        - [Forward Euler](https://en.wikipedia.org/wiki/Euler_method)
        - [Backward Euler](https://en.wikipedia.org/wiki/Backward_Euler_method)
        - [Crank-Nicolson](https://en.wikipedia.org/wiki/Crank%E2%80%93Nicolson_method)
        - [ERK4](https://en.wikipedia.org/wiki/Runge%E2%80%93Kutta_methods):
          Explicit Runge-Kutta method with four stages
        - [LSERK4](https://en.wikipedia.org/wiki/Runge%E2%80%93Kutta_methods) :
          Explicit low storage Runge-Kutta method with five stages
    - `Time step size` determines the time step size. This option is crucial for
      the stability of a simulation run. If it is too large, the simulation will
      become unstable.
    - `Final time` sets the end of the simulation.

- The `Expansion` section only offers one option and that is
    - `Expansion order`. It determines how many coefficients of the spherical
      harmonic expansion are computed.
  
- The `Finite element` section has only one option as well. It is
    - `Polynomial degree`. It sets the degree of the shape functions of the finite
      elements. For example, if it is \(1\), in each cell linear functions are used to
      approximate the expansion coefficients.

- The `Output` section offers
    - `Results folder` an option to set the path to the folder where the simulation
      results should be stored.
    - `Simulation identifier` is the name of the current simulation run. A
      folder with this name will be created in the results folder and the output
      of the corresponding run will be stored in this folder.
    - `Format` sets the format of the output. Currently only one format is
      implemented. It is `vtu`.


## Physical input

### Reference values

All computations are done in dimensionless units. The units are
	| Type           | Definition                           | Reference                                        |
	|:---------------|:-------------------------------------|:-------------------------------------------------|
	| Length         | \(x^{*} = x/r_{g,0}\)                | \(r_{g,0} = \gamma_{0} m_{0} v_{0}/q_{0} B_{0}\) |
	| Time           | \(t^{*} = t \omega_{g,0}\)           | \(\omega_{g,0} = q_{0} B_{0}/ \gamma_{0} m_{0}\) |
	| Velocity       | \(v^{*} = v/v_{0} \)                 | \(v_{0} = r_{g,0} \omega_{g,0} \)                |
	| Energy         | \(E^{*} = E/E_{0} \)                 | \(E_{0} = \gamma_{0} m_{0} {v_{0}}^{2} \)        |
	| Momentum       | \(p^{*} = p/p_{0} \)                 | \(p_{0} = \gamma_{0} m_{0} v_{0} \)              |
	| Lorentz factor | \(\gamma^{*} = \gamma / \gamma_{0}\) | \(\gamma_{0} =(1 - {v^{*}}^2)^{-1/2}\)           |
	| Mass           | \(m^{*} = m / m_{0} \)               |                                                  |
	| Charge         | \(q^{*} = q / q_{0} \)               |                                                  |
	| Magnetic field | \(B^{*} = B/ B_{0} \)                |                                                  |

All quantities with a zero subscript are reference values, which are set by the
user. The values \(E_{0}, p_{0}\) and \(\gamma_{0}\) are not independent. The
user only sets the reference energy, the other two quantities are computed. 

The reference values are set in `include/reference-values.h`. There default
values are

```cpp
struct ReferenceValues {
  // reference values
  const double energy = 10;       // GeV
  const double mass = 0.9382721;  // GeV/c^2 (proton mass)
  const double gamma = energy / mass;
  const double velocity = 299792458;              // m/s (speed of light)
  const double magnetic_field_strength = 1.e-10;  // Tesla (1 microGauss)
  const double charge = 1.602176634e-19;          // Coulmb ( elementary charge)
};
```
Except for the energy and the mass (which are in units of GeV and GeV/c^2)
everything is in SI units.

### Initial Condition, Fields, Source and Collisions

The initial condition, the fields, the source and the scattering frequency must
be provided as analytical expressions (and not, for example, as data in a file).
All of them are defined in the file `src/physical-setup.cpp`. The following
examples demonstrate how to implement a physical setup in sapphire.

#### Fields

A constant magnetic field \(B_{z} = 3\) is implemented as follows:

```cpp
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
  magnetic_field[0] = 0.; // B_x
  magnetic_field[1] = 0.; // B_y
  magnetic_field[2] = 3.; // B_z
}
```

This is the definition of the method `vector_value` of the object
`MagneticField`. The method has two inputs, a point and a vector, which is
called `magnetic_field`. The point
contains the coordinates at which to evaluate the magnetic field. And the vector
(is actually an output argument) and will contain the magnetic field at point. 

A \(x\)-dependent background velocity field \(u_{x} = 1/5 x\) is implemented as
follows

```cpp
// Velocity field implementation
template <int dim>
void Sapphire::BackgroundVelocityField<dim>::vector_value(
    const dealii::Point<dim> &point, dealii::Vector<double> &velocity) const {
  Assert(velocity.size() == 3,
         dealii::ExcDimensionMismatch(velocity.size(), 3));
  // EXAMPLES:
  // space dependent velocity
  // u_x = 1/5 * x
  // velocity[0] = 1. / 5 * point[0]; \\ u_x
  // velocity[1] = .0; \\ u_y
  // velocity[2] = .0; \\ u_z
}
```

Again, this is the definition of the method `vector_value` of the object
`BackgroundVelocityfield`. This time it is not a constant velocity field, but a
velocity field, which depends on \(x\). This is expressed with using the first
coordinate in `point` , i.e. `point[0]`. When sapphire needs to know the
velocity \(\mathbf{u}\), it calls `vector_value` and it hands over the
coordinates where it needs to know it, namely it gives it a `point` and en empty
vector, which then is filled as prescribed by `vector_value`.


A time dependent velocity field, say \(u_{x} = 1/10 t\), is implemented as
follows

```cpp
// Velocity field implementation
template <int dim>
void Sapphire::BackgroundVelocityField<dim>::vector_value(
    const dealii::Point<dim> &point, dealii::Vector<double> &velocity) const {
  Assert(velocity.size() == 3,
         dealii::ExcDimensionMismatch(velocity.size(), 3));
  // EXAMPLES:
  // time-dependent velocity-field
  // u_x = 1/10 * t
  velocity[0] = 1. / 10 * this->get_time(); \\ u_x
  velocity[1] = .0; \\ u_y
  velocity[2] = .0; \\ u_z
}
```

For a time dependent magnetic field it works exactly the same. Note that
independent of the total dimension of the problem, all the three components of
\(\mathbf{u}\) and \(\mathbf{B}\) are set. This allows the user to work with
"out-of-the-plane" components of the fields.

If the simulation include momentum terms, it is necessary to also define
analytical expressions for the divergence, the material derivative and the
Jacobian of the background velocity field. Follow the examples in
`src/physical-setup.cpp`.

#### Initial condition

The initial conditions determine the distribution of the particles at \(t = 0\).
The output of the simulation are the expansion coefficients of the spherical
harmonic expansion, namely \(f_{lms}(\mathbf{x},p,t)\). Hence, the initial
condition have to define what these are at the beginning of the simulation.
\(f_{000}\) represents the isotropic part of the particle distribution whereas
the coefficients with \(l > 0\) encode the its anisotropies. The initial value
function allows to set all the coefficients, but most of the time it is enough
to set the first one. The following example demonstrates how to do this:

```cpp
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
  // f[0] = 1. * std::exp(-(std::pow(p[0], 2)));
  // DIM = 2
  // f[0] = 1. * std::exp(-((std::pow(p[0], 2) + std::pow(p[1], 2))));
  // DIM 3
  f[0] = 1. * std::exp(-((std::pow(p[0], 2) + std::pow(p[1], 2) + std::pow(p[2], 2))));
}
```
This shows an isotropic Gaussian distribution of particles at \(t = 0\).
Depending on the total dimension of the problem a different line needs to be commented
out. A very important point to know is, that if momentum terms are included, the
last component of the point `p` corresponds to the \(p\)-direction For example, if
`terms` include the `momentum` flag and if `dim_configuration_space` is set to
one, then `p[2]` is \(p\). If `dim_configuration_space` was two, then
`p[3]` would be \(p\). 

#### Source term

The source term works very similar to the initial condition function. In
principal, it is possible to inject particles with anisotropies. But in many
cases, an isotropic injection is sufficient.

```cpp
template <int dim>
void Sapphire::Source<dim>::vector_value(
    const dealii::Point<dim> &p, dealii::Vector<double> &values) const {
  Assert(values.size() == (expansion_order + 1) * (expansion_order + 1),
         dealii::ExcDimensionMismatch(
             values.size(), (expansion_order + 1) * (expansion_order + 1)));
  // The zeroth component of values corresponds to isotropic part of the
  // distribution

  // EXAMPLES:
  // 2D pulsating Gaussian (isotropic, time dependent)
  values[0] = 0.01 * (std::sin(this->get_time()) + 1.) *
              std::exp(-(std::pow(p[0], 2) + std::pow(p[1], 2)));
}
```

The example shows a time dependent Gaussian injection of particles. As in the
case of the initial conditions, the last component of the point `p` represents
the momentum, if and only if momentum terms are included.

#### Scattering frequency

The scattering frequency is important if collisions are included in the
simulation. The following example shows a constant scattering frequency 

```cpp
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
``` 

The signature of this function different compared to the previous examples. The
scattering frequency is a scalar. Hence, it is not a vector which is returned,
but a value. Moreover, it is implemented such that it takes multiple points at
once (i.e. it takes a vector of points) and it returns the corresponding values
as a vector as well.

Another example is a momentum dependent scattering frequency. Assume `terms`
includes the `momentum` flag and `dim_configuration_space=1`, then 

```cpp
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
  // Momentum dependent scattering frequency
  for(unsigned int i = 0; i < points.size(); ++i) {
	scattering_frequencies[i] = 0.3/points[i][dim]
  }
}
```
gives the relation \(\nu = 0.3/p\).

### Transport only case

If the momentum terms are **not** included in the simulation, we speak of the
"transport only case". The energy of the particles cannot change and thus the
particles are only transported and **not** accelerated. An implication of this
is that the expansion coefficients do not depend on \(p\). In this case it is
necessary to set the energy of the particles. This can be done in the file
`include/physical-setup.h`. There is

```cpp
struct TransportOnly {
  // In the transport-only case (i.e. no dependence on p) the energy of the
  // particles has to be given.
  const double energy = 1.;
};

```
The Lorentz factor, the velocity and the momentum of the particles are derived
from their energy.

## Particle Properties

The particles are described by two properties: their mass and their charge. Both
can be set in the file `include/physical-setup.h`. There is a struct

```cpp
struct ParticleProperties {
  const double mass = 1.;
  const double charge = 1.;
};
```

# Help - I do not know what to do

Just write an email to nils.schween@mpi-hd.mpg.de and I will see if I can
support. 
