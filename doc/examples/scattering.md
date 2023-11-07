# Introduction and Basic Usage: Scattering {#scattering}

@tableofcontents


## Introduction {#introduction}

In this example, we want to give an introduction to @sapphire and show the basic
usage. All examples are structured in the following way: First, we give an
introduction, showing which parts of the VFP equation we want to solve, giving a
physical motivation and show how @sapphire translates this problem to finite
elements. Next, we show how to implement the example in @sapphire. We start with
a detailed description of the `config.h`, going through the file line by line
giving some explanatory comments. Similarly, we present the `scattering.cpp`
which implements the `main` function. In the results section, we show the
expected output of the program and explain how to interpret it.

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
f_{110}, f_{100}, f_{111}, f_{220}, f_{210}, f_{200}, f_{211}, f_{221} \dots)^
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


## Implementation {#code}

The implementation of the example in @sapphire is split in two files. The first
one is the header file `config.h` which implements compile time flags and the
physical setup. The later one consist of the definition of the initial condition
$f_{lms, 0}$ (@vfpref{InitialValueFunction}), the scattering frequency $\nu$
(@vfpref{ScatteringFrequency}), the source $S$ (@vfpref{Source}), the magnetic
field $\mathbf{B}$ (@vfpref{MagneticField}) and the background velocity field
$\mathbf{u}$ (@vfpref{BackgroundVelocityField}). The second file is the
`scattering.cpp` file that implements the `main` function. The since the
physical setup is defined in the `config.h` file, the `main` will be nearly
identical for most use-cases of @sapphire. Runtime parameter will be defined in
a parameter file `parameter.prm`.


### config.h {#config}

We start by going line by line through the `config.h` file. First, we have to
make sure that the file is only included once, and then import some
dependencies:

@snippet{lineno} examples/scattering/config.h Includes

Everything implemented in @sapphire is part of the namespace @ref Sapphire.

@snippet{lineno} examples/scattering/config.h Namespace Sapphire

Often we parametrize the physical setup with some runtime parameter. Since it is
set up dependent what these parameters are, they have to be specified by the
user. The @vfpref{PhysicalProperties} class allows for this. It uses the @dealii
concept of a @dealref{ParameterHandler}.

The @vfpref{PhysicalProperties} class consists of **public** variables for the
user defined runtime parameter, a default constructor and two functions to
**declare** and **parse** the parameter from the parameter file. In this
example, we have two parameters that we want to specify, the scattering
frequency $\nu$ and the initial value of the expansion coefficients,
$f_{lms,0}$. We will call these parameters `nu` and `f0` respectively, assuming
all expansion coefficients have the same initial value.

@snippet{lineno} examples/scattering/config.h Physical prop

The @ref Sapphire::PhysicalProperties::declare_parameters "declare_parameters()"
function is using the @dealref{ParameterHandler} class to declare parameter in
the parameter file. We sort all parameter in a subsection "Physical properties".
In addition, we write a message to the custom @ref
Sapphire::Utils::SapphireLogStream "SapphireLogStream" `saplog`. The
@dealref{LogStream::Prefix,classLogStream_1_1Prefix} ensures, that the message
is prefixed and only shown, if detailed output is requested. The declaration of
the parameter is straight forward, using the
@dealref{ParameterHandler.declare_entry(),classParameterHandler,a6d65f458be69e23a348221cb67fc411d}
function. It takes the name of the parameter, a default value and its
type/pattern. Additionally, one can give a description of the parameter.

@snippet{lineno} examples/scattering/config.h Declare params

The parsing of the parameter is equally simple. We use the
@dealref{ParameterHandler.get_double(),classParameterHandler,aeaf3c7846747695b1f327677e3716ec5}
function to get the value for a previously declared parameter.

@snippet{lineno} examples/scattering/config.h Parse params

At the end, we define the runtime parameter as **public** variables of the
class, so that subsequent functions have easy access to them.

@snippet{lineno} examples/scattering/config.h Runtime params

Next, we define static variables and functions related to the VFP equation. We
therefore use the namespace @ref Sapphire::VFP "VFP" to collect them in one
place.

@snippet{lineno} examples/scattering/config.h Namespace VFP

First, we define the dimensionality of the problem. Since the solution does not
depend on either $\mathbf{x}$ nor $p$, we just use one space dimension.

@snippet{lineno} examples/scattering/config.h Var dim

Next, we define which terms of the VFP equation we use. To this end, we define a
`static constexpr` which the **compiler** can use to determine which terms are
activated. Selecting only the relevant terms here results in big performance
gains, even though setting the respective terms to 0 at runtime would produce
the same results. In this example, we only have a scattering term, and no $p$
dependence. We therefore only activate the `collision` term, while all other
terms (and the $p$ dependence are turned off by default.)

@snippet{lineno} examples/scattering/config.h Var vfp flags

To define the initial values, we use the class @dealref{Function} provided by
@dealii. We define our own class @vfpref{InitialValueFunction} that inherits all
properties of the parent class @dealref{Function}. Keeping the style of @dealii,
we keep the dimension as a template parameter (even though it is at the moment
defined by the compile time variable `dimension`).

As we have seen before, the initial condition depends on only one parameter,
`f0`. But in this example, we will extend the scope of this function, by using
it as the analytic solution we can compare to. Therefore, the function will also
depend on `nu` as a function of time. The function has to provide a value for
each component $f_{lms}$. Therefore, it is a **vector-valued** function with
$(l_{\rm max} +1)^2$ components. To convert between the system index $i$ and the
spherical harmonic indices $l, m, s$ we need a mapping given by `lms_indices`.

@snippet{lineno} examples/scattering/config.h Initial value constructor

After defining the constructor of the function, we have to define its value. For
this we override the
@dealref{vector_value(),classFunction,ae316ebc05d21989d573024f8a23c49cb}
function of the parent class. At a *point* `p` it has to provide a *value* `f`
for each component. In this example, value is given by

$$
  f_{lms}(t) = f_{lms, 0} \exp\left(-\nu \frac{l(l + 1)}{2} t\right)\,.
$$

@snippet{lineno} examples/scattering/config.h Initial value vector value

The definition of the scattering frequency works similar. The only difference
is, that the scattering frequency is a scalar function, therefore we use the
function @dealref{value_list(),classFunction,abe86ee7f7f12cf4041d1e714c0fb42f3}
to get the scattering frequency at multiple points in one function call.

@snippet{lineno} examples/scattering/config.h Scattering freq

Since @sapphire expects a definition of the function @vfpref{Source},
@vfpref{MagneticField} and @vfpref{BackgroundVelocityField} in `config.h`,
we have to implement them here. We can however leave the implementation empty,
since the functions won't be used in this example.

@snippet{lineno} examples/scattering/config.h Source

Empty implementation of the magnetic field:

@snippet{lineno} examples/scattering/config.h Magnetic field

Empty implementation of the background velocity field:

@snippet{lineno} examples/scattering/config.h Velocity field

Last, we have to close the namespaces and the include guard.

@snippet{lineno} examples/scattering/config.h Close namespaces


### scattering.cpp {#mainFunction}

In this file, we will define the `main` function for the scattering
example. First, we include the header files we need:

@snippet{lineno} examples/scattering/scattering.cpp Includes

Then we define the `main` function:

@snippet{lineno} examples/scattering/scattering.cpp main func

We will run the program inside a try-catch block to catch any exceptions and
output them.

@snippet{lineno} examples/scattering/scattering.cpp try catch block begin

First, we need to initialize MPI for parallel runs:

@snippet{lineno} examples/scattering/scattering.cpp MPI init

We will use the @ref Sapphire::Utils::SapphireLogStream "SapphireLogStream"
class to output information to the console. The
@dealref{depth_console(),classLogStream,a8028e970ad8388596d625ed463894e98}
function allows specifying how detailed the output should be. We will set it to
2, which results in one line of output per time-step. For parallel runs, we will
turn of the output of all but the first process.

@snippet{lineno} examples/scattering/scattering.cpp Saplog

The program allows to specify the parameter file as a command line argument. In
case no argument is given, it will default to the file `parameter.prm`. When
starting the program, we output the configuration specified by command line
arguments.

@snippet{lineno} examples/scattering/scattering.cpp Command line argument

Next, we create all objects that declare runtime parameter. We prefix these
lines with `main` in the log stream, to silence the initialization process in a
run with low verbosity.

@snippet{lineno} examples/scattering/scattering.cpp Run time params

We now declare all parameters that the ParameterHandler should expect.

@note All parameters have to be declared before we can parse the parameter file.

@snippet{lineno} examples/scattering/scattering.cpp Declare params

After declaring the parameters, we can parse the parameter file and the
parameters of the objects

@snippet{lineno} examples/scattering/scattering.cpp Parse params

Finally, we can create the VFP equation solver and run it.

@snippet{lineno} examples/scattering/scattering.cpp VFP Solver

After the simulation ended, we can compute the error by comparing to the
analytic solution (implemented in @vfpref{InitialValueFunction}).

@snippet{lineno} examples/scattering/scattering.cpp L2 error

In case an exception is thrown, we will catch it, output it to `std::err` and
return with an error code.

@snippet{lineno} examples/scattering/scattering.cpp try catch block end


@todo Add section "Compiling the example" including the "CMakeLists.txt"


## Results {#results}

After compiling the example program, we can run it with

```bash
./scattering parameter.prm
```

With the [parameter file](@ref parameter) given below, we expect the following
output:

```bash
Sapphire::Start example scattering with parameter file "parameter.prm" on 1 processor(s) [2023/10/30 13:17:39]
WARNING: spatial advection is deactivated, but dim_cs > 0
  If the spatial advection term is deactivated,
  the distribution function is assumed to be homogeneous
  i.e. the dimension of the configuration space should be zero.
Sapphire:VFP::Run VFP equation solver.          [13:17:39]
Sapphire:VFP::Time step    1 at t = 0.00000     [13:17:39]
Sapphire:VFP::Time step    2 at t = 0.100000    [13:17:39]
Sapphire:VFP::Time step    3 at t = 0.200000    [13:17:39]
Sapphire:VFP::Time step    4 at t = 0.300000    [13:17:39]
Sapphire:VFP::Time step    5 at t = 0.400000    [13:17:39]
Sapphire:VFP::Time step    6 at t = 0.500000    [13:17:39]
Sapphire:VFP::Time step    7 at t = 0.600000    [13:17:39]
Sapphire:VFP::Time step    8 at t = 0.700000    [13:17:39]
Sapphire:VFP::Time step    9 at t = 0.800000    [13:17:39]
Sapphire:VFP::Time step   10 at t = 0.900000    [13:17:39]
Sapphire:VFP::Time step   11 at t = 1.00000     [13:17:39]
Sapphire:VFP::Time step   12 at t = 1.10000     [13:17:39]
Sapphire:VFP::Time step   13 at t = 1.20000     [13:17:39]
Sapphire:VFP::Time step   14 at t = 1.30000     [13:17:39]
Sapphire:VFP::Time step   15 at t = 1.40000     [13:17:39]
Sapphire:VFP::Time step   16 at t = 1.50000     [13:17:39]
Sapphire:VFP::Time step   17 at t = 1.60000     [13:17:39]
Sapphire:VFP::Time step   18 at t = 1.70000     [13:17:39]
Sapphire:VFP::Time step   19 at t = 1.80000     [13:17:39]
Sapphire:VFP::Time step   20 at t = 1.90000     [13:17:39]
Sapphire:VFP::The simulation ended.             [13:17:39]

+------------------------------------------+------------------+------------+------------------+
| Total wallclock time elapsed             |    0.1548s     0 |    0.1548s |    0.1548s     0 |
|                                          |                  |                               |
| Section                      | no. calls |   min time  rank |   avg time |   max time  rank |
+------------------------------------------+------------------+------------+------------------+
| DG matrix                    |        21 |   0.07757s     0 |   0.07757s |   0.07757s     0 |
| FE system                    |         1 |   0.01056s     0 |   0.01056s |   0.01056s     0 |
| Grid setup                   |         1 |  0.001882s     0 |  0.001882s |  0.001882s     0 |
| Mass matrix                  |         1 |  0.002406s     0 |  0.002406s |  0.002406s     0 |
| Output                       |        21 |    0.0196s     0 |    0.0196s |    0.0196s     0 |
| Project f onto the FEM space |         1 |  0.004495s     0 |  0.004495s |  0.004495s     0 |
| Setup                        |         1 |   0.02981s     0 |   0.02981s |   0.02981s     0 |
| Theta method                 |        20 |    0.1114s     0 |    0.1114s |    0.1114s     0 |
+------------------------------------------+------------------+------------+------------------+
Sapphire:VFP::Compute the global error
Sapphire::L2 error: 0.000654138
```

@todo Explain console output

The results of the run will be saved in the `results/scattering` folder. In the
folder are two different kinds of files. First, the `log.prm` is a log of the
parameter file used for the run. This allows to reproduce/identify the results
of a simulation at a later stage. Second, the `solution_*.vtu` files with the
solution at the different time steps. These files can be visualized using
visualization software like [ParaView](https://www.paraview.org) or
[VisIt](https://visit-dav.github.io/visit-website/).

In our example the result is not very interesting. The values for $f_{lms}(x,t)$
are calculated for $x$ between $-1$ and $1$. As expected, the solution decays
with time as $f_{lms}(t) = \exp\left(- \frac{l(l + 1)}{2} t \right)$.

The following image shows the results at $t=2$ for $l = 0,1,2,3$:

![Results at t=2](/path/to/image)

@todo Use relative picture path, or path to git.


### Parameter file {#parameter}

The parameter file can either be a `.prm`, `.json` or `.xml` file. Here, we use
a `.prm` for easy readability. For more information on the different input
formats, we refer to the @dealref{ParameterHandler} documentation.

\include{lineno} scattering/parameter.prm


## The plain program {#plainCode}


### config.h {#plainConfig}

\include{lineno} scattering/config.h


### scattering.cpp {#plainMain}

\include{lineno} scattering/scattering.cpp

---

@author Florian Schulze (<florian.schulze@mpi-hd.mpg.de>)
@date 2023-10-31
