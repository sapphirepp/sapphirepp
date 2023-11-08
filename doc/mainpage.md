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

It is made to simulate the interaction of particles with a background plasma. To
this end it solves a Vlasov-Fokker-Planck equation in mixed coordinates, namely

$$
  \frac{\partial f}{\partial t} + (\mathbf{u} + \mathbf{v}) \cdot \nabla_{x} f -
  \gamma m \frac{\mathrm{D} \mathbf{u}}{\mathrm{D} t} \cdot \nabla_{p}f -
  \mathbf{p} \cdot\nabla_{x} \mathbf{u}\cdot \nabla_{p} f +
  q \mathbf{v} \cdot \left( \mathbf{B} \times \nabla_{p} f \right) =
  \frac{\nu}{2} \Delta_{\theta, \varphi} f + S(\mathbf{x}, \mathbf{p}, t)\, .
$$


# Installation {#installation}


## System requirements {#requirements}


@sapphire uses the C++ library @dealii to solve the differential equations with
finite elements. We use MPI to parallelize @sapphire. For optimized parallel
grid generation, linear algebra operations and linear system solvers, we use the
[p4est](https://www.p4est.org), [PETSc](https://petsc.org),
[SLEPc](https://slepc.upv.es) and interface of @dealii.
In order to compile and use @sapphire you need the following programs installed:

- [CMake](http://www.cmake.org/) version 3.5.0 or later
- [GNU make](http://www.gnu.org/software/make/)
- [HDF5](https://www.hdfgroup.org/solutions/hdf5/) with MPI support
- [Open-MPI](https://www.open-mpi.org)
- [boost](https://www.boost.org)
- [LAPACK](https://www.netlib.org/lapack/)
- [Threading Building Blocks (TBB)](https://github.com/oneapi-src/oneTBB)
- [p4est](https://www.p4est.org)
- [PETSc](https://petsc.org)
- [SLEPc](https://slepc.upv.es)
- @dealii
- For generating the documentation: [Doxygen](https://www.doxygen.nl) with `dot`
- For visualization we recommand to use [VisIt](http://www.llnl.gov/visit/) or
  [ParaView](http://www.paraview.org/)

To install @dealii and all necessary prerequisites for @sapphire, we provide an
[installation script](https://git.mpi-hd.mpg.de/schween/sapphire/install_dealii.sh).
Run the script by executing:

```shell
./install_dealii.sh
```

Follow the on-screen instructions during the installation process. If the script
aborts due to failed tests from a missing Fortran compiler, try rerunning the
script and skip the installation of already installed packages. If the
installation process is interrupted during the compilation of @dealii, you can
resume the installation with the following commands:


```shell
cd dealii-X.X.X/build
make -j install
```
If this process fails, it may be that your system has run out of memory. In 
such case, specifiy the number of threads after `-j`.

For macOS user, @dealii offers prepackaged `.dmg` files with all dependencies
included. To install, follow the
[deal.II Mac OSX Instructions](https://github.com/dealii/dealii/wiki/MacOSX).


## Compiling @sapphire {#compilation}

You can obtain @sapphire either as a `tar.gz` file
[here](https://git.mpi-hd.mpg.de/schween/sapphire/releases/latest), or by
cloning the `git` project:

```shell
git clone https://git.mpi-hd.mpg.de/schween/sapphire.git
```

To compile the project, use `cmake`:

```shell
cd sapphire
export DEAL_II_DIR="path/to/deal.II"
cmake -S . -B build -DCMAKE_BUILD_TYPE=Debug
cd build && make -j
```

For extensive simulation runs, it's recommended to use to the release version:

```shell
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
```

Developers are advised to use the following options:

```shell
cmake -S . -B build -DCMAKE_BUILD_TYPE=Debug -DEXAMPLES=ON -DTESTS=ON -DDOC=ON
```

To run sapphire, use:

```shell
./sapphire parameter.prm
```

or

```shell
mpirun ./sapphire parameter.prm
```

For more details on using @sapphire, refer to the [examples](#examples).

@todo Replace `Gitea` links <https://git.mpi-hd.mpg.de/schween/sapphire> with
`GitHub` links (including the correct version)

@todo Rename `sapphire` folder to `sapphirepp`
