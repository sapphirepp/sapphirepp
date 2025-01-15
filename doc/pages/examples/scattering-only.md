# Advanced example: Scattering only {#scattering-only}

@tableofcontents

## Introduction {#introduction-scattering-only}

This example provides an in-depth exploration of @sapphire using a simple case
study. It's designed for advanced users and new developers of @sapphire who are
already familiar with the basic concepts presented in the [quick
start](#quick-start) and [parallel shock](#parallel-shock) examples. If you're
new to @sapphire, we recommend starting with these examples.

In this tutorial, we'll examine a simplified version of the VFP equation that
includes only a scattering term:

$$
  \frac{\partial f}{\partial t} =  \frac{\nu}{2} \Delta_{\theta, \varphi} f \,.
$$

We'll use a constant scattering frequency $\nu = \nu_0$. When we expand this
equation in spherical harmonics (as detailed in @cite Schween2025), we get:

$$
  \partial_{t} f_{lms} + \nu \, \frac{l(l + 1)}{2} f_{lms} = 0 \,.
$$

This shows that the expansion coefficients $f_{lms}$ are independent of each
other. Therefore, the solution to this equation is an exponential decay:

$$
  f_{lms}(t) = f_{lms, 0} \exp\left(-\nu \frac{l(l + 1)}{2} t \right) \,,
$$

where $f_{lms, 0} = f_{lms}(t=0)$.

## Implementation {#implementation-scattering-only}

This section provides a detailed walkthrough of the example's implementation in
@sapphire. You can locate the relevant files in the
`sapphirepp/examples/scattering-only` directory.

The example requires four files:

- [config.h](#config-scattering-only): Defines the physical setup, including
  the initial condition, scattering frequency, magnetic field, and more.
- [scattering-only.cpp](#main-scattering-only): Implements the `main` function
  and executes the simulation.
- [CMakeLists.txt](#compiling-scattering-only): The CMake file that outlines
  the compilation process for the example.
- [parameter.prm](#example-parameter-scattering-only): Specifies the runtime
  parameters.

We will go through each of these files line by line, highlighting the
conventions and best practices employed in @sapphire.

@note The `config.h` file is included in the @sapphire library. As such, it must
  retain its name and define specific variables and functions used by @sapphire.

### config.h {#config-scattering-only}

We start by investigating the `config.h` file line by line. The complete file
can be found at the bottom of this page and in the
`sapphirepp/examples/scattering-only` directory.

First, we set up the include guard and import the necessary dependencies:

- `cmath`: Provides access to mathematical functions like `std::exp()`.
- `vector`: Enables use of the `std::vector` class.

- `deal.II/base/exceptions.h`: Facilitates throwing exceptions.
- `deal.II/base/point.h`: Provides the @dealref{Point} class.
- `deal.II/base/function.h`: Allows implementation of functions as
  @dealref{Function} classes.
- `deal.II/base/parameter_handler.h`: Enables reading parameters from a
  parameter file.
- `deal.II/lac/vector.h`: Provides the @dealref{Vector} class.

- `pde-system.h`: Implements basic functionality of the VFP PDE system, like the
  @ref sapphirepp:VFP::PDESystem::create_lms_indices "create_lms_indices()"
  function.
- `vfp-flags.h`: Implements the enum @ref sapphirepp::VFP::VFPFlags "VFPFlags".
- `sapphirepp-logstream.h`: Facilitates output to the console using different
  levels of verbosity.

@snippet{lineno} examples/scattering-only/config.h Includes

To avoid naming conflicts between variables and classes, we encapsulate
everything related to @sapphire within the @ref sapphirepp namespace:

@snippet{lineno} examples/scattering-only/config.h Namespace sapphirepp

#### Custom Runtime Parameters {#parameter-scattering-only}

In this example, we define two custom runtime parameters using the @ref
sapphirepp::PhysicalParameters "PhysicalParameters" class, as introduced in the
[parallel shock](#parallel-shock) example. These parameters are the scattering
frequency $\nu$ and the initial value of the expansion coefficients $f_0$
(assuming the same initial value for all coefficients).

These parameters are defined as **public** variables at the start of the @ref
sapphirepp::PhysicalParameters "PhysicalParameters" class:

@snippet{lineno} examples/scattering-only/config.h PhysicalParameters

@note The `// !!!EDIT HERE!!!` comments are left as a reference where the user
  modifies the code as presented in the [parallel shock](#parallel-shock)
  example.

Next, we **declare** the parameters in a @dealref{ParameterHandler} object in
the @ref sapphirepp::PhysicalParameters::declare_parameters()
"declare_parameters()" function. We enter the subsection `Physical parameters`
at the start and leave it at the end. Notice the use of the @ref
sapphirepp::saplog "saplog" stream to output debug information, with
@dealref{LogStream::Prefix,classLogStream_1_1Prefix} controlling the verbosity
of the output.

@snippet{lineno} examples/scattering-only/config.h Declare parameters

Finally, **parsing** the parameters is straightforward in the @ref
sapphirepp::PhysicalParameters::parse_parameters() "parse_parameters()"
function. Again, we use the @ref sapphirepp::saplog "saplog" stream to output
debug information.

@snippet{lineno} examples/scattering-only/config.h Parse parameters

#### VFP equation {#dimension-scattering-only}

Next, we define static variables and functions related to the VFP equation.
These are grouped together under the @ref sapphirepp::VFP "VFP" namespace for
organization and ease of access.

@snippet{lineno} examples/scattering-only/config.h Namespace VFP

First, we specify the dimensionality of the problem. Given that the solution is
independent of both $\mathbf{x}$ and $p$, we set `dim = 1` (since `dim = 0` is
not allowed in the code). This corresponds to one space dimension `dim = dim_cs
= 1` as the momentum dependence is not activated. This variable is defined as a
`static constexpr` to ensure its value is known at compile time.

@snippet{lineno} examples/scattering-only/config.h Dimension

Next, we specify the terms of the VFP equation to be used. The `vfp_flags`
variable is also defined as `static constexpr`, enabling the compiler to
determine the active terms of the VFP equation at compile time. Selecting only
the relevant terms can significantly boost performance, even though setting the
respective terms to zero at runtime would yield the same results. In this
example, we only have a scattering term, @ref
sapphirepp::VFP::VFPFlags::collision "collision". All other terms (and the $p$
dependence) are deactivated by default.

@snippet{lineno} examples/scattering-only/config.h VFP Flags

#### Initial condition {#initial-condition-scattering-only}

To define the initial condition, we implement it as a @dealii @dealref{Function}
class. We create the class @ref sapphirepp::VFP::InitialValueFunction
"InitialValueFunction", which inherits all properties of the parent class
@dealref{Function}. In line with the @dealii style, we keep the dimension as a
template parameter, even though it is defined in the compile-time variable
`dimension`. This function needs to provide a value for each component
$f_{lms}$, making it a **vector-valued** function with `system_size =
(l_max+1)^2` components.

In the constructor, we pass the @ref sapphirepp::PhysicalParameters
"PhysicalParameters" and store it in the `private const` variable `prm`. This
allows us to access the runtime parameters within the function. Additionally, we
create a mapping between the system index $i$ and the spherical harmonic indices
$l, m, s$ using the @ref sapphirepp::VFP::PDESystem::create_lms_indices
"create_lms_indices()" function.

@snippet{lineno} examples/scattering-only/config.h InitialValueFunction constructor

We define the function value at a `point` in the @ref
sapphirepp::VFP::InitialValueFunction::vector_value "vector_value()" function,
overriding the
@dealref{vector_value(),classFunction,ae316ebc05d21989d573024f8a23c49cb}
function of the parent class. This function needs to provide a value for each
component $f_{lms}$ at the given point, which is stored in the results vector
`f`.

To catch implementation errors, we check that the size of the results vector `f`
matches the number of components `system_size`, using the
@dealref{Assert,group__Exceptions,ga70a0bb353656e704acf927945277bbc6} macro,
which is only activated in debug mode. Note that we use
`static_cast<void>(point)` to avoid a compiler warning about an unused variable.
We then loop over all indices $i$ and set the result vector `f` to $f_{0}$.

Finally, we close the function definition, after defining the private variable
`prm` and the mapping `lms_indices`.

@snippet{lineno} examples/scattering-only/config.h InitialValueFunction value

#### Scattering frequency {#scattering-frequency-scattering-only}

The scattering frequency is defined similarly to the initial condition. However,
the scattering frequency is a scalar function, so we use the @dealref{Function}
class with only one component.

@snippet{lineno} examples/scattering-only/config.h ScatteringFrequency constructor

To enhance performance by reducing the number of function calls, we employ the
@ref sapphirepp::VFP::ScatteringFrequency::value_list "value_list()" function.
This function allows us to obtain the function value at multiple points in a
single function call. Therefore, we override the
@dealref{value_list(),classFunction,abe86ee7f7f12cf4041d1e714c0fb42f3} function,
which takes a vector of `points` and the return vector `scattering_frequencies`
as arguments. The `component` argument can be ignored as we only have a scalar
function.

Initially, we verify that the size of the `scattering_frequencies` vector
matches the number of points. Then, we iterate over all points and set the
scattering frequency to $\nu$.

Again, we close the function definition, after defining the private variable
`prm`.

@snippet{lineno} examples/scattering-only/config.h ScatteringFrequency value

#### Source term {#source-term-scattering-only}

We need to provide definitions for the @ref sapphirepp::VFP::Source "Source",
@ref sapphirepp::VFP::MagneticField "MagneticField", and @ref
sapphirepp::VFP::BackgroundVelocityField "BackgroundVelocityField" functions in
`config.h` as expected by @sapphire. However, since these functions are not
utilized in this example, we can leave their implementations empty.

The @ref sapphirepp::VFP::Source "Source" function is a vector-valued function,
where each component corresponds to the expansion coefficients $S_{lms}$. Its
implementation closely mirrors that of the @ref
sapphirepp::VFP::InitialValueFunction "InitialValueFunction". Therefore, we skip
the detailed explanation and directly provide the code.

@snippet{lineno} examples/scattering-only/config.h Source

#### Magnetic field {#magnetic-field-scattering-only}

The @ref sapphirepp::VFP::MagneticField "MagneticField" function is a
vector-valued function with three components, representing the magnetic field
vector $\mathbf{B}=(B_x,B_y,B_z)^\top$. It is important to note that even if the
configuration space dimension is less than three, we still use three components.
This allows for a magnetic field that points out of the plane. Again, we skip
the detailed explanation as the implementation is similar to @ref
sapphirepp::VFP::InitialValueFunction "InitialValueFunction".

@snippet{lineno} examples/scattering-only/config.h MagneticField

#### Velocity field {#velocity-scattering-only}

The @ref sapphirepp::VFP::BackgroundVelocityField "BackgroundVelocityField"
function is a vector-valued function with three components,
$\mathbf{u}=(u_x,u_y,u_z)^\top$, similar to the @ref
sapphirepp::VFP::MagneticField "MagneticField". However, as demonstrated in the
[quick start](#quick-start) example, we need to provide additional functions for
the divergence, material derivative, and Jacobian of the velocity field. We
achieve this by using the @dealref{Function} class and defining custom
functions.

The constructor is similar to the @ref sapphirepp::VFP::MagneticField
"MagneticField".

@snippet{lineno} examples/scattering-only/config.h BackgroundVelocityField

Inside the velocity field class, we define four different functions:

1. Background velocity field value $\mathbf{u}(\mathbf{x})$:

   We use the @ref sapphirepp::VFP::BackgroundVelocityField::vector_value
   "vector_value()" function to define all components of the velocity field
   value at one `point`. This function overrides the
   @dealref{vector_value(),classFunction,ae316ebc05d21989d573024f8a23c49cb}
   function of the parent class.

   @snippet{lineno} examples/scattering-only/config.h BackgroundVelocityField value

2. Background velocity divergence $\nabla \cdot \mathbf{u}(\mathbf{x})$:

   We use the custom @ref
   sapphirepp::VFP::BackgroundVelocityField::divergence_list "divergence_list()"
   function to define the divergence of the velocity field at multiple `points`
   at once. This mirrors the
   @dealref{value_list(),classFunction,abe86ee7f7f12cf4041d1e714c0fb42f3}.

   @snippet{lineno} examples/scattering-only/config.h BackgroundVelocityField divergence

3. Background velocity material derivative $\frac{\mathrm{D}
   \mathbf{u}}{\mathrm{D} t}$:

   We use the custom @ref
   sapphirepp::VFP::BackgroundVelocityField::material_derivative_list
   "material_derivative_list()" function to define the material derivative of
   the velocity field at multiple `points` at once.

   @snippet{lineno} examples/scattering-only/config.h BackgroundVelocityField material derivative

4. Background velocity Jacobian $\frac{\partial u_{x}}{\partial x}$:

   We use the custom @ref
   sapphirepp::VFP::BackgroundVelocityField::jacobian_list "jacobian_list()"
   function to define the Jacobian of the velocity field at multiple `points` at
   once.

   After defining these functions, we define the private parameters and close the
   class definition.

   @snippet{lineno} examples/scattering-only/config.h BackgroundVelocityField Jacobian

Finally, we close the namespaces and the include guard.

@snippet{lineno} examples/scattering-only/config.h Close namespaces

### scattering-only.cpp {#main-scattering-only}

The second file, `scattering-only.cpp`, houses the `main` function. It closely
mirrors the main executable for @sapphire, `sapphirepp/sapphirepp.cpp`. Since
the setup is already implemented in the `config.h` file, this file usually
remains unchanged across different examples. However, in this case, we will
incorporate additional functionality to compare the numerical solution with the
analytic solution.

The file begins with the inclusion of several header files:

- `deal.II/base/mpi.h` and `mpi.h` enable MPI parallelization.
- `deal.II/base/parameter_handler.h` allows reading parameters from a parameter
  file.
- `config.h` is the configuration file we previously implemented.
- `vfp-solver.h` implements the VFP equation solver, which is the main
  interface to @sapphire.
- `vfp-parameters.h` declares the parameters related to VFP equation
  solver, such as the time step size and final time.
- `output-parameters.h` declares the parameters related to output, like the
  output directory and the output frequency.
- `sapphirepp-logstream.h` enables console and logfile output
  with varying levels of verbosity.

@snippet{lineno} examples/scattering-only/scattering-only.cpp Includes

#### Analytic solution {#analytic-solution-scattering-only}

We start by implementing a @dealref{Function} for the analytic solution.
Following the [coding conventions](#coding-conventions), we use an internal
namespace to avoid naming conflicts with other functions.

@snippet{lineno} examples/scattering-only/scattering-only.cpp AnalyticSolution namespace

The analytic solution is a vector valued function with each component
representing one expansion coefficient. The constructor is similar to the @ref
sapphirepp::VFP::InitialValueFunction "InitialValueFunction". One additional
parameter is the time at which the solution should be evaluated.

@snippet{lineno} examples/scattering-only/scattering-only.cpp AnalyticSolution constructor

The @dealref{vector_value(),classFunction,ae316ebc05d21989d573024f8a23c49cb}
function defines the value of the analytic solution at a given `point`. We can
convert the system index $i$ to $l$ using `l = lms_indices[i][0]`, and access
the time via `t = this->get_time()`. To recall, the analytic solution is given
by:

$$
  f_{lms}(t) = f_{0} \exp\left(-\nu \frac{l(l + 1)}{2} t \right) \,.
$$

@snippet{lineno} examples/scattering-only/scattering-only.cpp AnalyticSolution value

#### Main function {#main-function-scattering-only}

The `main` function is defined next, with `argc` and `argv` used to process
command-line arguments.

@snippet{lineno} examples/scattering-only/scattering-only.cpp Main function

The program is run within a try-catch block to catch any exceptions and output
them. The namespaces @ref sapphirepp and @ref sapphirepp::VFP are used to avoid
prefixing the respective classes.

@snippet{lineno} examples/scattering-only/scattering-only.cpp Try-Catch begin

Before proceeding, MPI needs to be initialized for parallel runs using the
@dealref{MPI_InitFinalize,classUtilities_1_1MPI_1_1MPI__InitFinalize} utility.
This utility initializes MPI, PETSc, and p4est at the start and shuts them down
at the end of the program. The last argument of the constructor specifies the
number of threads to be used per MPI process. We avoid TBB parallelization and
only use one thread per process.

@snippet{lineno} examples/scattering-only/scattering-only.cpp MPI initialization

The custom log stream @ref sapphirepp::saplog "saplog"
is used to output status and debug information to the console and logfiles.
The verbosity level and logfile can be specified via the command line,
see @ref sapphirepp::Utils::SapphireppLogStream::init "saplog.init()".
The default verbosity of `2` corresponds to progress information.

@snippet{lineno} examples/scattering-only/scattering-only.cpp Saplog

The parameter file is specified as a command-line argument. Therefore, the
command-line arguments are processed, defaulting to the file `parameter.prm` if
no argument is given.

@snippet{lineno} examples/scattering-only/scattering-only.cpp Command line argument

Next, all objects that declare runtime parameters are created. These include the
@ref sapphirepp::VFP::VFPParameters "VFPParameters", @ref
sapphirepp::Utils::OutputParameters "OutputParameters", and @ref
sapphirepp::PhysicalParameters "PhysicalParameters" classes. The
@dealref{ParameterHandler} is used to handle the parameter file.

@snippet{lineno} examples/scattering-only/scattering-only.cpp Run time parameters

Before reading the parameter file, all parameters that the
@dealref{ParameterHandler} should expect are declared.

@snippet{lineno} examples/scattering-only/scattering-only.cpp Declare parameters

After declaring the parameters, the parameter file is parsed and the parameters
of the objects are set.

@snippet{lineno} examples/scattering-only/scattering-only.cpp Parse parameters

Finally, we can create the VFP equation solver @ref sapphirepp::VFP::VFPSolver
"VFPSolver" and run the simulation.

@snippet{lineno} examples/scattering-only/scattering-only.cpp VFP Solver

#### Comparing with the analytic solution {#error-scattering-only}

After the simulation finishes, we can evaluate the error by comparing the
results with the analytic solution. To do this, we first need to create the
analytic solution function.

@snippet{lineno} examples/scattering-only/scattering-only.cpp Create AnalyticSolution

We can then compute the L2 error using the @ref
sapphirepp::VFP::VFPSolver::compute_global_error "compute_global_error()"
function. Additionally, we compute the norm of the solution to obtain the
relative error.

@snippet{lineno} examples/scattering-only/scattering-only.cpp Calculate L2-error

We also want to write the analytic solution to a file. To this end, we need to
create a @dealref{MPI::Vector,classPETScWrappers_1_1MPI_1_1Vector} object. We
store the finite element coefficients of the analytic solution in this vector
using the
@dealref{interpolate(),namespaceVectorTools,a32b68e11070496fcaf6c2389f088d09c}
function. Note that the interpolation ensures that the values at the nodes are
given by the function values. This is not the case for a projected function,
which the weak solution approximates. We refer to the discussion in the
[visualization](#visualization) section.

@snippet{lineno} examples/scattering-only/scattering-only.cpp Calculate analytic solution

Finally, we write the analytic solution to a file using the @dealref{DataOut}
object and the @ref sapphirepp::Utils::OutputParameters::write_results
"write_results()" function. We can make use of the @ref
sapphirepp::VFP::PDESystem::create_component_name_list
"create_component_name_list()" function to name the components of the
vector.

@snippet{lineno} examples/scattering-only/scattering-only.cpp Output analytic solution

In case an exception is thrown,
we catch it and use the
@ref sapphirepp::Utils::SapphireppLogStream::print_error "saplog.print_error()"
method to output it to both the log stream and `std::err`.
Then we return with an error code.

@snippet{lineno} examples/scattering-only/scattering-only.cpp Try-Catch end

### CMakeLists.txt {#compiling-scattering-only}

To compile the example, we need to create a `CMakeLists.txt` file.
This file is used by the [CMake](https://cmake.org) build system
to generate makefiles that compile the program.
Since the `config.h` file is included in the main library,
we must link against it.
This requirement means that @sapphire cannot be treated as an independent library.
As a consequence, all examples or applications of @sapphire
must be included alongside the main library.

To simplify the process,
we defined a macro `sapphirepp_setup_target_vfp`
which links all necessary libraries and files to the target.
As long as applications are located in the `sapphirepp/examples` folder,
there are only two steps to follow:

1. Add the application as a subdirectory to the
   `sapphirepp/examples/CMakeLists.txt` file:

    ```cmake
    add_subdirectory(scattering-only)
    ```

2. Add the following `CMakeLists.txt` in the application folder:

   @include{lineno} examples/scattering-only/CMakeLists.txt

After these steps, re-run `cmake` to update the build system:

```shell
cd sapphirepp
rm -rf build
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DEXAMPLES=ON
```

@note For development, we recommend using the Debug mode. You can switch between
  Debug and Release mode by changing the `CMAKE_BUILD_TYPE` variable.

This process will create a `sapphirepp/build/examples/scattering-only` folder.
To build the executable for the example, execute `make` in this folder:

```shell
cd sapphirepp/build/examples/scattering-only
make
```

This will create the executable `scattering-only`, which can be run with:

```shell
./scattering-only parameter.prm
```

### parameter.prm {#example-parameter-scattering-only}

The final file required is the parameter file. This file can be in `.prm`,
`.json`, or `.xml` format. For more details on the input formats, we refer to
the @dealref{ParameterHandler} documentation. For a reference on the syntax used
in the `.prm` file, see the [parallel shock](#parallel-shock) example.

In this example, we use a scattering frequency of $\nu = 0.1$ and an initial
condition of $f_0=1$. We set $l_{\rm max}=9$, which results in a CFL condition
of:

$$
  \Delta t \leq \frac{2}{\nu \, l_{\rm max} (l_{\rm max} + 1)} = 0.222 \,.
$$

Given this, we select a time step size of $\Delta t = 0.1$. We also opt to use a
fourth-order explicit Runge-Kutta method, which should provide adequate
precision for this simulation.

@include{lineno} examples/scattering-only/parameter.prm

## Executing the example {#execute-scattering-only}

Once the example program is compiled, it can be executed with the following
command:

```shell
cd sapphirepp/build/examples/scattering-only
./scattering-only parameter.prm
```

The console output should resemble the following:

```output
Sapphire::Start Sapphire++ v1.1.0-dev with 1 MPI process(es) [2023/12/22 15:34:41]
WARNING: spatial advection is deactivated, but dim_cs > 0
  If the spatial advection term is deactivated,
  the distribution function is assumed to be homogeneous
  i.e. the dimension of the configuration space should be zero.
Sapphire:VFP::Run VFP equation solver.          [15:34:41]
Sapphire:VFP::Time step    0 at t = 0.00000     [15:34:41]
...
Sapphire:VFP::Time step    9 at t = 0.900000    [15:34:44]
Sapphire:VFP::Simulation ended at t = 1.00000   [15:34:44]
Sapphire:VFP::Compute the global error
Sapphire:VFP::Compute the weighted norm
Sapphire::L2_error = 0.000199350, L2_norm = 4.54748, rel error = 4.38375e-05
Sapphire:OutputParameters::Writing results at time_step 0
```

The output begins with a startup message and a warning. The warning is due to
the lack of spatial dependence in this simple scattering-only example, despite
`dimension = 1` being used in the code. This warning can be ignored for this toy
problem, but in real applications, it might indicate a faulty setup.

The @ref sapphirepp::VFP::VFPSolver "VFPSolver" then starts and evolves the
solution in 10 time steps. The relative error is approximately $5 \times
10^{-5}$, which is expected for a 4th order Runge-Kutta scheme with a time step
size of $0.1$.

The program creates a `results` folder in the current directory, with a
`scattering-only` subdirectory. This folder contains a log file of the
parameters (`log.prm`), the analytic solution (`analytic_solution_0000.vtu`),
and the numerical solution at each time step (`solution_*.vtu`).

The `*.vtu` files can be opened in visualization software like @paraview or
[VisIt](https://visit-dav.github.io/visit-website/). These files contain all
expansion coefficients `f_lms` along with the MPI `subdomain`. Note that the
solution `f_lms` is a `dim = dim_ps` dimensional quantity in mixed coordinate
reduced phase space.

## Results {#results-scattering-only}

The $f_{l00}(t)$ solution yields the following result:

![Results GIF](https://sapphirepp.org/img/examples/scattering-only.gif)

When compared with the analytic solution at $t=1$, we expect the following
values:

$$
  f_{lms}^{\rm analytic} = f_0 \exp\left(-\nu \frac{l(l + 1)}{2} t \right) \,.
$$

| $l$ | $f_{lms}(t=1)$ | $f_{lms}^{\rm analytic}$ |
|----:|:---------------|:-------------------------|
| 0   | 1.0000         | 1.0000                   |
| 1   | 0.9048         | 0.9048                   |
| 2   | 0.7408         | 0.7408                   |
| 3   | 0.5488         | 0.5488                   |
| 4   | 0.3679         | 0.3679                   |
| 5   | 0.2231         | 0.2231                   |
| 6   | 0.1225         | 0.1225                   |
| 7   | 0.0608         | 0.0608                   |
| 8   | 0.0273         | 0.0273                   |
| 9   | 0.0111         | 0.0111                   |

The numerical and analytic solutions align within this precision, as illustrated
in the following plot:

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
