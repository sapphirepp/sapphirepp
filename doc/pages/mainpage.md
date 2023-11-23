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
phase space. $\mathbf{U}$ is the velocity of the background plasma, e.g. the
velocities of a plasma around a collisionless shock as observed after supernova
explosions. $\mathbf{B}$ is mean magnetic field of this plasma. $S$ is a source
term representing the injection of the charged particles. They interact with the
background plasma via isotropic and elastic scattering off fluctuations of the
plasma's mean electromagnetic fields. $\nu$ is the corresponding scattering
frequency. $\mathbf{U}$, $\mathbf{B}$, $S$ and $\nu$ are the input of the
simulation and can be analytical functions of $t, \mathbf{x}$ and $p$. The
distribution function $f$ is its output.

@sapphire tries to exploit that the distributions of energetic and charged
particles are in many physical environments almost isotropic. It, thus,
expresses $f$ as a series of spherical harmonics, i.e.

$$
	f(t, \mathbf{x}, \mathbf{p}) = \sum^{\infty}_{l = 0} \sum^{l}_{m = 0}
	\sum^{1}_{s = 0} f_{lms}(t, \mathbf{x}, p) Y_{lms}(\theta,\varphi) \, ,
$$

whose zeroth order term is the isotropic part and higher order terms represent
the anisotropies of the distribution function. @sapphire computes the expansion
coefficients $f_{lms}$, hence it solves a system of partial differential equations
(PDEs), which is derived from the above VFP equation.

A distinguishing feature of @sapphire is that it solves this system of PDEs by
applying the discontinuous Galerkin (dG) method. The dG method is a finite
element method which is particularly suited for the shown VFP equation.

Even though primarily developed in the context of cosmic-ray physics, @sapphire
can be applied in the context of inertial confinement fusion, or for example to
compute the transport of radiation (interpreting $f$ as the intensity).

@sapphire is developed by the [Astrophysical Plasma
Theory](https://www.mpi-hd.mpg.de/mpi/en/research/scientific-divisions-and-groups/independent-research-groups/apt)
group located at the [Max Planck Institute for Nuclear
Physics](https://www.mpi-hd.mpg.de/mpi/en/) in Heidelberg, Germany.



# Installation {#installation}


## System requirements {#requirements}


@sapphire uses the C++ library @dealii to solve the differential equations with
finite elements. We use MPI to parallelize @sapphire. For optimized parallel
grid generation and linear algebra operations, we use the
[p4est](https://www.p4est.org), [PETSc](https://petsc.org) and though their
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
- @dealii version @dealii_min_version or later
- For generating the documentation: [Doxygen](https://www.doxygen.nl) version
  @doxygen_min_version or later with `dot`
- For visualization we recommand to use [VisIt](http://www.llnl.gov/visit/) or
  [ParaView](http://www.paraview.org/)

To install @dealii and all necessary prerequisites for @sapphire, we provide an
[installation script](https://github.com/sapphirepp/sapphirepp/blob/main/scripts/install-dealii.sh).
Run the script by executing:

```shell
chmod u+x install-dealii.sh
./install_dealii.sh
```

Follow the on-screen instructions during the installation process. If the script
aborts due to failed tests from a missing Fortran compiler, try rerunning the
script and skip the installation of already installed packages. If the
installation process is interrupted during the compilation of @dealii, you can
resume the installation with the following commands:


```shell
cd dealii-X.X.X/build
make -j N install
```

where `X.X.X` is the version number of @dealii and `N` is the number of threads
to use for the compilation (should match number of processors on your system).

For macOS user, @dealii offers prepackaged `.dmg` files with all dependencies
included. To install, follow the
[deal.II Mac OSX Instructions](https://github.com/dealii/dealii/wiki/MacOSX).


## Compiling @sapphire {#compilation}

You can obtain @sapphire either as a tarball from the
[release page](https://github.com/sapphirepp/sapphirepp/releases), or by
cloning the `git` repository:

```shell
git clone https://github.com/sapphirepp/sapphirepp
```

@note Since the repository is still private, you need to use the following link:
  ```shell
  git clone https://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@github.com/sapphirepp/sapphirepp
  ```

Afterwards you can compile @sapphire, with `cmake` and `make`:

```shell
cd sapphirepp
export DEAL_II_DIR="path/to/deal.II"
cmake -S . -B build -DCMAKE_BUILD_TYPE=Debug
cd build && make
```

@note On the MPIK cluster you need to load the following modules:
  ```shell
  module load python/3.11.6 cmake gcc/13.1.0 ucx/1.15.x-zen4 openmpi/4.1.5-zen4 boost/1.83.0 hdf5/1.14.1-2-zen4 intel-oneapi-mkl/2023.1.0-zen4 
  export DEAL_II_DIR="/lfs/l8/tap/tapshared/lib/dealii"
  ```

For extensive simulation runs, it's recommended to use to the release version:

```shell
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
```

Developers are advised to use the following options:

```shell
cmake -S . -B build -DCMAKE_BUILD_TYPE=Debug -DEXAMPLES=ON -DTESTS=ON -DDOC=ON -DCMAKE_CXX_FLAGS="-Werror -Wall -Wextra" -DCMAKE_CXX_CLANG_TIDY=clang-tidy
```


## Getting started {#getting-started}

To run @sapphire, you need to provide a configuration file:

```shell
./sapphirepp parameter.prm
```

You enable parallel execution by using `mpirun`:

```shell
mpirun ./sapphirepp parameter.prm
```

For more details on using @sapphire, refer to the [examples](#examples).


# Licence

Sapphire++ is distributed under the [LGPL 3.0 license](https://github.com/sapphirepp/sapphirepp/blob/main/LICENSE).

If you use this software in your research, please cite the following paper:

> [Your Name], [Paper Title], [Journal], [Year], DOI: [DOI]

Here's the BibTeX entry for the paper:

```bibtex
@article{YourPaper,
  title={Paper Title},
  author={Your Name},
  journal={Journal},
  year={Year},
  doi={DOI}
}
```

@todo Add paper reference

