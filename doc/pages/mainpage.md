@mainpage notitle

<img src="sapphire-logo-text.png" alt="Sapphire++ logo" style="display:block;
margin: 0 auto; width:40%">

@tableofcontents

# About {#about}

@sapphire is an acronym and stands for \"<strong>S</strong>imulating
<strong>a</strong>strophysical <strong>p</strong>lasmas and
<strong>p</strong>articles with <strong>hi</strong>ghly
<strong>r</strong>elativistic <strong>e</strong>nergies in
C<strong>++</strong>\".

It is a code to simulate the interaction of charged particles with a background
plasma, a typical example is the propagation and acceleration of cosmic rays. To
this end it solves a Vlasov-Fokker-Planck (VFP) equation in mixed coordinates,
namely

$$
  \frac{\partial f}{\partial t} + (\mathbf{u} + \mathbf{v}) \cdot \nabla_{x} f -
  \gamma m \frac{\mathrm{D} \mathbf{u}}{\mathrm{D} t} \cdot \nabla_{p}f -
  \mathbf{p} \cdot\nabla_{x} \mathbf{u}\cdot \nabla_{p} f +
  q \mathbf{v} \cdot \left( \mathbf{B} \times \nabla_{p} f \right) =
  \frac{\nu}{2} \Delta_{\theta, \varphi} f + S\, .
$$

$f$ is the distribution function of the particles, i.e. their number density in
phase space. $\mathbf{u}$ is the velocity of the background plasma, e.g. the
velocities of a plasma around a collisionless shock as observed after supernova
explosions. Note that the momentum $\mathbf{p}$ is given in mixed coordinates,
meaning it is defined relative to the background field $\mathbf{u}$.
$\mathbf{B}$ is mean magnetic field of this plasma. $S$ is a source term
representing the injection of the charged particles. They interact with the
background plasma via isotropic and elastic scattering off fluctuations of the
plasma's mean electromagnetic fields. $\nu$ is the corresponding scattering
frequency. $\mathbf{u}$, $\mathbf{B}$, $S$ and $\nu$ are the input of the
simulation and can be analytical functions of $t, \mathbf{x}$ and $p$. The
distribution function $f$ is its output.

@sapphire tries to exploit that the distributions of energetic and charged
particles are in many physical environments almost isotropic. It, thus,
expresses $f$ as a series of spherical harmonics, i.e.

$$
 f(t, \mathbf{x}, \mathbf{p}) = \sum^{\infty}_{l = 0} \sum^{l}_{m = 0}
 \sum^{1}_{s = 0} f_{lms}(t, \mathbf{x}, p) Y_{lms}(\theta,\varphi) \, ,
$$

where $\mathbf{p} = p \, (\cos\theta, \sin\theta \cos\varphi, \sin\theta
\sin\varphi)^T$ is used. The zeroth order term is the isotropic part and higher
order terms represent the anisotropies of the distribution function. @sapphire
computes the expansion coefficients $f_{lms}$, hence it solves a system of
partial differential equations (PDEs), which is derived from the above VFP
equation.

A distinguishing feature of @sapphire is that it solves this system of PDEs by
applying the discontinuous Galerkin (dG) method. The dG method is a finite
element method which is particularly suited for the shown VFP equation.

The full implementation details can be found in @cite Schween2024a and
@cite Schween2025.

Even though primarily developed in the context of cosmic-ray physics, @sapphire
can be applied in the context of inertial confinement fusion, or for example to
compute the transport of radiation (interpreting $f$ as the intensity).

@sapphire is developed by
[Nils Schween](https://github.com/nils-schween),
[Florian Schulze](https://github.com/floschulze) and
[Brian Reville](https://github.com/brevrev)
members of the [Astrophysical Plasma
Theory](https://www.mpi-hd.mpg.de/mpi/en/research/scientific-divisions-and-groups/independent-research-groups/apt)
group located at the [Max-Planck-Institut für
Kernphysik](https://www.mpi-hd.mpg.de/mpi/en/) in Heidelberg, Germany.

# Installation {#installation}

## System requirements {#requirements}

@sapphire uses the C++ library @dealii to solve the differential equations with
finite elements. We use MPI to parallelize @sapphire. For optimized parallel
grid generation and linear algebra operations, we use the
[p4est](https://www.p4est.org) and [PETSc](https://petsc.org) through their
interface in @dealii. In order to compile and use @sapphire you need the
following programs installed:

- [CMake](http://www.cmake.org/) version @cmake_min_version or later
- [GNU make](http://www.gnu.org/software/make/)
- [Open-MPI](https://www.open-mpi.org)
- [Threading Building Blocks (TBB)](https://github.com/oneapi-src/oneTBB)
- [boost](https://www.boost.org)
- [zlib](https://zlib.net)
- [HDF5](https://www.hdfgroup.org/solutions/hdf5/) with MPI support
- [LAPACK](https://www.netlib.org/lapack/)
- [p4est](https://www.p4est.org)
- [PETSc](https://petsc.org)
- [assimp](https://www.assimp.org/)
- @dealii version @dealii_min_version or later
- For generating the documentation: [Doxygen](https://www.doxygen.nl) version
  @doxygen_min_version or later with `dot`
- For visualization we recommand to use [VisIt](http://www.llnl.gov/visit/) or
  @paraview

To install @dealii and all necessary prerequisites for @sapphire, we provide an
[installation script](https://github.com/sapphirepp/sapphirepp/blob/main/scripts/install-dealii.sh).
Download the script and it with run the script by executing:

```shell
curl -O https://raw.githubusercontent.com/sapphirepp/sapphirepp/refs/heads/main/scripts/install-dealii.sh
chmod u+x install-dealii.sh
./install-dealii.sh
```

Follow the on-screen instructions during the installation process. There are a
few common errors occurring during installation:

- If the script aborts due to failed tests from a missing Fortran compiler. This
  sometimes occurs when testing the [PETSc](https://petsc.org) installation. You
  can try rerunning the script but skip the installation
  [PETSc](https://petsc.org) (it should be installed successfully already).

  ```shell
  ./install-dealii.sh --skip-prerequisites
  ```

- If the installation process is interrupted during the compilation of @dealii,
  this is often due to compiler running out of memory. At this stage, you can
  resume the installation with the following commands:

  ```shell
  cd dealii-X.X.X/build
  make -j N install
  ```

  where `X.X.X` represents the version number of @dealii, and `N` specifies the
  number of threads to use for the compilation. To avoid memory issues, we
  recommend restarting the compilation with fewer threads.

For macOS user, @dealii offers prepackaged `.dmg` files with all dependencies
included. To install, follow the
[deal.II Mac OS X Instructions](https://github.com/dealii/dealii/wiki/MacOSX).

If you have trouble installing @dealii,
refer to the
[deal.II installation instructions](https://www.dealii.org/current/readme.html).

## Compiling @sapphire {#compilation}

You can obtain @sapphire either as a tarball from the
[release page](https://github.com/sapphirepp/sapphirepp/releases), or by
cloning the `git` repository:

```shell
git clone https://github.com/sapphirepp/sapphirepp.git
```

Afterwards you can compile @sapphire, with `cmake` and `make`:

```shell
cd sapphirepp
export DEAL_II_DIR="path/to/deal.II"
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DEXAMPLES=ON
make --directory=build
```

Developers are advised to use the following `cmake` options:

```shell
cmake -S . -B build -DCMAKE_BUILD_TYPE=Debug -DEXAMPLES=ON -DTESTS=ON -DDOC=ON -DCMAKE_CXX_FLAGS="-Werror -Wall -Wextra" -DCMAKE_CXX_CLANG_TIDY=clang-tidy
```

@note Do not forget to switch to the optimized `Release` build before starting
      extensive simulation runs.

## Getting started {#getting-started}

To run @sapphire, you need to provide a parameter file:

```shell
./build/sapphirepp parameter-template.prm
```

You enable parallel execution by using `mpirun`:

```shell
mpirun -np N ./build/sapphirepp parameter-template.prm
```

where `N` is the number of processors to use.

For a brief introduction, refer to our [quick-start guide](#quick-start).
Additional information and detailed examples can be found in the
[examples](#examples) section.

If you have any questions,
feel free to reach out to us
via the [GitHub Discussions page](https://github.com/sapphirepp/sapphirepp/discussions)
for support and community interaction.

<div class="section_buttons">

|                        Next |
|----------------------------:|
| [Quick start](#quick-start) |

</div>
