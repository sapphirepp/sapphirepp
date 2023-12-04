# Advanced example: Scattering only {#scattering-only}

@tableofcontents


## Introduction {#introduction-scattering-only}

In this example, we want to give an introduction to @sapphire and show its basic
usage. To this end, we structure the example in three parts. First, we give an
introduction and physical motivation. Next, we show how to implement the example
in @sapphire, by going through the respective files line by line giving some
explanatory comments. In the results section, we show the expected
output of the program and explain how to interpret it.

As an introductory example, we want to study the simple case with only a
scattering term in the VFP equation,

$$
  \frac{\partial f}{\partial t} =  \frac{\nu}{2} \Delta_{\theta, \varphi} f \,.
$$

Here $f(\mathbf{x}, \mathbf{p}, t)$ is the distribution function and $\nu$ is
the scattering frequency. The Laplacian $\Delta_{\theta, \varphi}$ only acts on
the directional part of the momentum $\mathbf{p} = (p \cos\theta, p \cos\varphi
\sin\theta, p \sin\varphi \cos\theta)^{\top}$.

@sapphire uses an expansion with real spherical harmonics,
$Y_{lms}(\theta,\varphi)$, to solve the VFP equation. The distribution function
can therefore be written as

$$
  f(\mathbf{x}, \mathbf{p}, t) = \sum^{\infty}_{l = 0} \sum^{l}_{m = 0}
  \sum^{1}_{s = 0} f_{lms}(\mathbf{x}, p ,t) Y_{lms}(\theta,\varphi) \,.
$$

This expansion uses, that the distribution function is nearly isotropic, with
anisotropies described by the higher moments $f_{l>0}$. The VFP equation now
turns into a system of equations for the expansion coefficients
$f_{lms}(\mathbf{x}, p, t)$. Notice, that $f_{lms}$ depends only on the absolute
value of the momentum $p$. This reduces the dimensionality of the problem from
$(\mathrm{dim\_cs})^2$ to $\mathrm{dim\_ps} = \mathrm{dim\_cs}+1$. Here, we
introduce call the nomenclature configuration space $\mathrm{cs}$ for the
dimension of $\mathbf{x}$ and (reduced) phase space $\mathrm{ps}$ for
$(\mathbf{x}, p)$.

Collecting all expansion coefficients into a vector $\boldsymbol{f} = (f_{000},
f_{110}, f_{100}, f_{111}, f_{220}, f_{210}, f_{200}, f_{211}, f_{221}, \dots)^
{\top}$ the VFP equation now becomes

$$
  \partial_{t} \boldsymbol{f} + \nu \, \mathbf{C} \boldsymbol{f} = 0 \,.
$$

The matrix $\mathbf{C}$ represents the scattering operator and is given by (see
@cite Schween2024)

$$
  (\mathbf{C})_{i(l'm's')j(l,m,s)} = \frac{l(l + 1)}{2}\delta_{l'l}\delta_{m'm}
  \delta_{s's} \,.
$$

Focusing on the equation for a single expansion coefficient, it becomes obvious
that in this example the different expansion coefficients are independent,

$$
  \partial_{t} f_{lms} + \nu \, \frac{l(l + 1)}{2} f_{lms} = 0 \,.
$$

This differential equation describes and exponential decay,

$$
  f_{lms}(t) = f_{lms, 0} \exp\left(-\nu \frac{l(l + 1)}{2} t \right) \,,
$$

where $f_{lms, 0} = f_{lms}(t=0)$. We will use this solution to check the
correctness of the numerical solution.


## Implementation {#implementation-scattering-only}

To implement an example like this in @sapphire, we need two different files. The
first one is the header file `config.h` which implements the physical setup,
i.e. the initial condition, scattering frequency, magnetic field, etc. The
second file is a `.cpp` file that implements the `main` function. In our
example, we will call this file `scattering-only.cpp`.

@note The `config.h` file will be included in the @sapphire library. Therefore,
 it must not be renamed and has to define certain variables and functions used
 by @sapphire.

Furthermore, we need to adapt the `CMakeLists.txt` file to compile the example,
and provide a parameter file that specifies the runtime parameters.


### config.h {#config-scattering-only}

We start by investigating the `config.h` file line by line. Notice that the
complete file is given at the bottom of the page and in the
[examples/scattering-only](https://github.com/sapphirepp/sapphirepp/tree/main/examples/scattering-only)
folder.

First, we have to
make sure that the file is only included once, and then import some
dependencies:

@snippet{lineno} examples/scattering-only/config.h Includes

To avoid confusion of variable and class names, we collect everything related to
@sapphire in the @ref sapphirepp namespace:

@snippet{lineno} examples/scattering-only/config.h Namespace sapphirepp

Often we parametrize the physical setup with some runtime parameter. Since it is
set up dependent what these parameters are, they have to be specified by the
user. The @ref sapphirepp::PhysicalParameters "PhysicalParameters" class allows
for this. It uses the @dealii @dealref{ParameterHandler} to read parameters from
a parameter file.

The @ref sapphirepp::PhysicalParameters "PhysicalParameters" class consists of
**public** variables for the user defined runtime parameter, a default
constructor and two functions to **declare** and **parse** the parameter from
the parameter file. In this example, we have two parameters that we want to
specify, the scattering frequency $\nu$ (`nu`) and the initial value of the
expansion coefficients, $f_{lms,0}$ (`f0`).

@snippet{lineno} examples/scattering-only/config.h PhysicalParameters

The @ref sapphirepp::PhysicalParameters::declare_parameters
"declare_parameters()" function is using the @dealref{ParameterHandler} class to
declare parameter in the parameter file. We collect all parameter in a
subsection "Physical parameters" (`prm.enter_subsection()`). The declaration of
the parameter is straight forward, using the
@dealref{ParameterHandler.declare_entry(),classParameterHandler,a6d65f458be69e23a348221cb67fc411d}
function. It needs the name of the parameter, a default value and its
type/pattern as input. Additionally, one can give a description of the
parameter. After all parameters are declared, we have to close the subsection
again (`prm.leave_subsection()`).

@snippet{lineno} examples/scattering-only/config.h Declare parameters

The parsing of the parameter is equally simple. We use the
@dealref{ParameterHandler.get_double(),classParameterHandler,aeaf3c7846747695b1f327677e3716ec5}
function to get the value for a previously declared parameter. Again, we first
have to enter the subsection at the start and leave it at the end.

@snippet{lineno} examples/scattering-only/config.h Parse parameters

Next, we define static variables and functions related to the VFP equation. We
therefore use the namespace @ref sapphirepp::VFP "VFP" to collect them in one
place.

@snippet{lineno} examples/scattering-only/config.h Namespace VFP

First, we define the dimensionality of the problem. Since the solution does not
depend on either $\mathbf{x}$ nor $p$, we just use one space dimension. We have
to define this variable as a `static constexpr`, so it is known at compile time.

@snippet{lineno} examples/scattering-only/config.h Dimension

Next, we define which terms of the VFP equation we use. Again, we define a
`static constexpr` which the **compiler** can use to determine the active terms
of the VFP equation. Selecting only the relevant terms results in big
performance gains, even though setting the respective terms to 0 at runtime
would produce the same results. In this example, we only have a scattering term,
and no $p$ dependence. We therefore only activate the `collision` term, while
all other terms (and the $p$ dependence) are turned off by default.

@snippet{lineno} examples/scattering-only/config.h VFP flags

We define initial values as a @dealref{Function} class provided by @dealii. We
define our own class @ref sapphirepp::VFP::InitialValueFunction
"InitialValueFunction" that inherits all properties of the parent class
@dealref{Function}. Keeping the style of @dealii, we keep the dimension as a
template parameter (even though it is at the moment defined by the compile time
variable `dimension`).

As we have seen before, the initial condition depends on only one parameter,
`f0`. But in this example, we will extend the scope of this function, by using
it as the analytic solution we can compare to. Therefore, the function will also
depend on `nu` as a function of time. The function has to provide a value for
each component $f_{lms}$. Therefore, it is a **vector-valued** function with
$(l_{\rm max} +1)^2$ components. To convert between the system index $i$ and the
spherical harmonic indices $l, m, s$ we need a mapping given by `lms_indices`.

@snippet{lineno} examples/scattering-only/config.h InitialValueFunction constructor

After defining the constructor of the function, we have to define its value. For
this we override the
@dealref{vector_value(),classFunction,ae316ebc05d21989d573024f8a23c49cb}
function of the parent class. At a *point* `p` it has to provide a *value* `f`
for each component. In this example, value is given by

$$
  f_{lms}(t) = f_{lms, 0} \exp\left(-\nu \frac{l(l + 1)}{2} t\right)\,.
$$

@snippet{lineno} examples/scattering-only/config.h InitialValueFunction value

The definition of the scattering frequency works similar. The only difference
is, that the scattering frequency is a scalar function, therefore we use the
function @dealref{value_list(),classFunction,abe86ee7f7f12cf4041d1e714c0fb42f3}
to get the scattering frequency at multiple points in one function call.

@snippet{lineno} examples/scattering-only/config.h ScatteringFrequency

Since @sapphire expects a definition of the function @ref
sapphirepp::VFP::Source "Source", @ref sapphirepp::VFP::MagneticField
"MagneticField" and @ref sapphirepp::VFP::BackgroundVelocityField
"BackgroundVelocityField" in `config.h`, we have to implement them here. We can
however leave the implementation empty, since the functions won't be used in
this example.

@snippet{lineno} examples/scattering-only/config.h Source

Empty implementation of the magnetic field:

@snippet{lineno} examples/scattering-only/config.h MagneticField

Empty implementation of the background velocity field:

@snippet{lineno} examples/scattering-only/config.h BackgroundVelocityField

Last, we have to close the namespaces and the include guard.

@snippet{lineno} examples/scattering-only/config.h Close namespaces


### scattering-only.cpp {#main-scattering-only}

In the second file, we will define the `main` function for the scattering
example. Since the setup is already implemented in the `config.h` file, this
file will not change much for different examples. We will explain the file in
detail once. More advanced examples will only point out differences to this
example.

First, we include a list of header files:

- `deal.II/base/mpi.h` and `mpi.h` allow MPI parallelization
- `deal.II/base/parameter_handler.h` to read parameters from a parameter
  file
- `config.h` is the configuration file we just implemented
- `vfp-solver.h` implements the VFP equation solver. This class solves the VFP
  equation and is the main interface to @sapphire
- `vfp-parameters.h` declares the parameters related to VFP equation
  solver, like the time step size and final time
- `output-parameters.h` declares the parameters related to output, like the
  output directory and the output frequency
- `sapphirepp-logstream.h` allows output to the console using different levels
  of verbosity

@snippet{lineno} examples/scattering-only/scattering-only.cpp Includes

Then we define the `main` function, allowing for command line arguments:

@snippet{lineno} examples/scattering-only/scattering-only.cpp Main function

We will run the program inside a try-catch block to catch any exceptions and
output them. Furthermore, we will use the namespace @ref sapphirepp and @ref
sapphirepp::VFP, so we don't have to prefix the respective classes.

@snippet{lineno} examples/scattering-only/scattering-only.cpp Try-Catch begin

Before we start anything else, we need to initialize MPI for parallel runs. We
use the @dealref{MPI_InitFinalize,classUtilities_1_1MPI_1_1MPI__InitFinalize}
utility that initializes MPI, PETSc and p4est at the start and shuts them down
at the end of the program. The last argument of the constructor specifies how
many threads should be used per MPI process. We will avoid TBB parallelization
and only use one thread per process.

@snippet{lineno} examples/scattering-only/scattering-only.cpp MPI initialization

In order to output status and debug information to the console, we will use the
custom log stream @ref sapphirepp::saplog "saplog". We can specify the verbosity
via the argument of the @ref sapphirepp::Utils::SapphireppLogStream::init
"init()" function. A verbosity of `2` correspond to progress information.

@snippet{lineno} examples/scattering-only/scattering-only.cpp Saplog

We want to specify the parameter file as a command line argument. Therefore, we
need to process the command line arguments but default to the file
`parameter.prm` in case no argument is given. Furthermore, we want to confirm
the parameter file by outputting it to the console (using the `saplog` stream).

@snippet{lineno} examples/scattering-only/scattering-only.cpp Command line argument

Next, we create all objects that declare runtime parameter. Namely, the @ref
sapphirepp::VFP::VFPParameters "VFPParameters", @ref
sapphirepp::Utils::OutputParameters "OutputParameters" and @ref
sapphirepp::PhysicalParameters "PhysicalParameters" class. We will use the
@dealref{ParameterHandler} handle the parameter file.

@snippet{lineno} examples/scattering-only/scattering-only.cpp Run time parameters

Before reading the parameter file, we have to declare all parameters that the
ParameterHandler should expect.

@snippet{lineno} examples/scattering-only/scattering-only.cpp Declare parameters

After declaring the parameters, we can parse the parameter file and the
parameters of the objects.

@snippet{lineno} examples/scattering-only/scattering-only.cpp Parse parameters

Finally, we can create the VFP equation solver @ref sapphirepp::VFP::VFPSolver
"VFPSolver" and run it.

@snippet{lineno} examples/scattering-only/scattering-only.cpp VFP Solver

After the simulation ends, we can compute the error by comparing to the analytic
solution (implemented in @ref sapphirepp::VFP::InitialValueFunction
"InitialValueFunction").

@snippet{lineno} examples/scattering-only/scattering-only.cpp L2 error

In case an exception is thrown, we will catch it, output it to `std::err` and
return with an error code.

@snippet{lineno} examples/scattering-only/scattering-only.cpp Try-Catch end


### Compiling {#compiling-scattering-only}

In order to compile the example, we need to create a `CMakeLists.txt` file. This
file is used by the [CMake](https://cmake.org) build system to create makefiles
that will compile the program. Since the `config.h` file is included in the main
library, we have to link against it. As a consequence, @sapphire can not be
understood as an independent library. This is reflected in the fact, that all
examples or applications of sapphire have to be included alongside the main
library.

To simplify things a little, we define a bunch of variables in the root
`CMakeLists.txt`. As long as applications are located in the examples folder,
there are only two things to do:

1. Add the application as a subfolder to the `examples/CMakeLists.txt` file:

    ```cmake
    add_subdirectory(scattering-only)
    ```

2. Add the following `CMakeLists.txt` in the application folder:
  
@include{lineno} examples/scattering-only/CMakeLists.txt

Now calling `cmake` inside the root folder (see [Compiling
Sapphire++](#compilation)) will result in creating a
`build/examples/scattering-only` folder. To build the example executing the
following commands (staring in the `build` folder):

```shell
cmake .
cd examples/scattering-only
make -j
```


### Parameter file {#parameter-scattering-only}


The last file we need is the parameter file. This can be either a `.prm`,
`.json` or `.xml` file. More information on the different input formats can be
found in the @dealref{ParameterHandler} documentation.

The parameter file is split into three sections. First, the `Physical
parameters` section with all user defined parameters from the @ref
sapphirepp::PhysicalParameters "PhysicalParameters" class. Next, the `Output`
section, containing the parameters from the @ref
sapphirepp::Utils::OutputParameters "OutputParameters" class. Third, the `VFP`
section with all parameters related to @ref sapphirepp::VFP::VFPParameters
"VFPParameters".  
These sections can contain a number of parameters and subsections using the
following syntax:

```param
subsection Section Name
  set Parameter Name = value
  ...

  subsection Subsection Name
    set Parameter Name = value
    ...
  end
end
```

Notice that section and parameter names are case-sensitive and  can contain
spaces.

@include{lineno} examples/scattering-only/parameter.prm


## Results {#results-scattering-only}

After compiling the example program, we can run it with:

```shell
./scattering-only parameter.prm
```

This results in console output similar to:

```output
Sapphire::Start Sapphire++ v1.1.0-dev with 1 MPI process(es) [2023/11/23 18:43:53]
Sapphire::Start example scattering-only with parameter file "parameter.prm"
WARNING: spatial advection is deactivated, but dim_cs > 0
  If the spatial advection term is deactivated,
  the distribution function is assumed to be homogeneous
  i.e. the dimension of the configuration space should be zero.
Sapphire:VFP::Run VFP equation solver.          [18:43:53]
Sapphire:VFP::Time step    0 at t = 0.00000     [18:43:53]
Sapphire:VFP::Time step    1 at t = 0.100000    [18:43:54]
Sapphire:VFP::Time step    2 at t = 0.200000    [18:43:54]
Sapphire:VFP::Time step    3 at t = 0.300000    [18:43:54]
Sapphire:VFP::Time step    4 at t = 0.400000    [18:43:55]
Sapphire:VFP::Time step    5 at t = 0.500000    [18:43:55]
Sapphire:VFP::Time step    6 at t = 0.600000    [18:43:55]
Sapphire:VFP::Time step    7 at t = 0.700000    [18:43:56]
Sapphire:VFP::Time step    8 at t = 0.800000    [18:43:56]
Sapphire:VFP::Time step    9 at t = 0.900000    [18:43:56]
Sapphire:VFP::Simulation ended at t = 1.00000   [18:43:57]

+------------------------------------------+------------------+------------+------------------+
| Total wallclock time elapsed             |     3.584s     0 |     3.584s |     3.584s     0 |
|                                          |                  |                               |
| Section                      | no. calls |   min time  rank |   avg time |   max time  rank |
+------------------------------------------+------------------+------------+------------------+
| DG matrix                    |        21 |     3.233s     0 |     3.233s |     3.233s     0 |
| ERK4                         |        10 |     3.113s     0 |     3.113s |     3.113s     0 |
| FE system                    |         1 |    0.1036s     0 |    0.1036s |    0.1036s     0 |
| Grid setup                   |         1 |  0.002364s     0 |  0.002364s |  0.002364s     0 |
| Mass matrix                  |         1 |   0.07479s     0 |   0.07479s |   0.07479s     0 |
| Output                       |        11 |    0.1217s     0 |    0.1217s |    0.1217s     0 |
| Project f onto the FEM space |         1 |    0.0128s     0 |    0.0128s |    0.0128s     0 |
| Setup                        |         1 |    0.3662s     0 |    0.3662s |    0.3662s     0 |
+------------------------------------------+------------------+------------+------------------+
Sapphire:VFP::Compute the global error
Sapphire::L2 error: 0.000199350
```

We first see two star-up messages, followed by a warning. This is related to the
problem that we have no spatial dependence in this simple scattering-only
example, but still use `dimension = 1` in the code. Since this is only a toy
problem, we can ignore this warning, in real applications it would indicate that
we might have made a mistake in the setup.  
Next, we see that the @ref sapphirepp::VFP::VFPSolver "VFPSolver" starts, and
evolves the solution in 10 time steps. Last we get detailed timing information.
We see that most time is spent in the DG matrix assembly, followed by the ERK4
time stepping. This of course will change depending on the machine and build
options (Debug vs. Release).  
Last, we see that the L2 error is $\sim 10^{-4}$. Using a 4th order Runge-Kutta
scheme and a time step size of $0.1$, this is expected.

Let's now have a look at the output files saved in `results/scattering-only/`.
We see, this folder contains a `log.prm` file and a number of `solution_*.vtu`
files. The `log.prm` file is a copy of the parameter file used for the run,
allowing to reproduce/identify the results at a later stage. The
`solution_*.vtu` files contain the solution $f_{lms}(t)$ at the respective time
steps. These files can be visualized using visualization software like
[ParaView](https://www.paraview.org) or
[VisIt](https://visit-dav.github.io/visit-website/).

After opening the files as a time series, we can see that it contains all
coefficients `f_lms` and the `subdomain` for each domain (in this case we only
have one domain with id `0`). Showing only the `f_l00` solution in a 2d plot, we
get something like:

![Results GIF](https://sapphirepp.org/img/examples/scattering-only.gif)

Let's compare this with our expectation. The analytic solution is given by
$f_{lms}(t) = f_{lms, 0} \, e^{-\nu \frac{l(l + 1)}{2} t}$. As we can see, the
higher $l$ modes in deed decay faster, while $f_{000}$ is constant. We can also
check, that the solution is independent of $m$ and $s$. With
$f_{lms, 0} = 1$, $\nu = 0.1$ and $t = 1$ we expect the following values:

$$
  f_{lms}(t) = \exp\left(-0.1 \frac{l(l + 1)}{2} \right)
$$


We can use the solution at $t=1$ for a more quantitative comparison. The
expected values are:

```text
f_0ms = 1.0000
f_1ms = 0.9048
f_2ms = 0.7408
f_3ms = 0.5488
f_4ms = 0.3679
f_5ms = 0.2231
f_6ms = 0.1225
f_7ms = 0.0608
f_8ms = 0.0273
f_9ms = 0.0111
```

![Results at t=1](https://sapphirepp.org/img/examples/scattering-only-t1.png)


## The plain program {#plain-scattering-only}


### config.h {#plain-config-scattering-only}

@include{lineno} examples/scattering-only/config.h


### scattering-only.cpp {#plain-main-scattering-only}

@include{lineno} examples/scattering-only/scattering-only.cpp


<div class="section_buttons">

| Previous              |
|:----------------------|
| [Examples](#examples) |

</div>


---

@author Florian Schulze (<florian.schulze@mpi-hd.mpg.de>)
@date 2023-11-23
