# Advanced example: Convergence study {#convergence-study}

@tableofcontents

## Introduction {#introduction-convergence-study}

This example presents a convergence study for @sapphire, utilizing a
one-dimensional test case with an analytical solution. It also provides an
overview of advanced features in @sapphire. It is recommended to familiarize
yourself with @sapphire through the previous examples before starting this one,
particularly the [quick start](#quick-start) and [scattering
only](#scattering-only) examples.

We consider a simple one one-dimensional test case with a static magnetic field
$\mathbf{B} = B_0 \hat{\mathbf{e}}_z$ in a static medium, $\mathbf{u} = 0$.
Assuming a collisionless plasma, the VFP equation simplifies to:

$$
  \frac{\partial f}{\partial t} + \mathbf{v} \cdot \nabla_{x} f +
  q \mathbf{v} \cdot \left( \mathbf{B} \times \nabla_{p} f \right) =
  0\, .
$$

Using a cut-off in the spherical harmonic expansion order $l_{\rm max}=1$, the
system of equations for the expansion coefficients $f_{lms}$ is given by:

\begin{align*}
   & \partial_{t} f_{000} +
  \frac{v}{\sqrt{3}} \partial_{x} f_{100}
  = 0                       \\
   & \partial_{t} f_{110} -
  \omega f_{100}
  = 0                       \\
   & \partial_{t} f_{100} +
  \frac{v}{\sqrt{3}} \partial_{x} f_{000} +
  \omega f_{110}
  = 0 \\
  & \partial_{t} f_{111} = 0
\end{align*}

Here, we introduced the gyrofrequency $\omega = \frac{B_0 q}{\gamma m}$, with
the magnetic field $B_0$, the charge $q$, the Lorentz factor $\gamma$, and the
mass $m$. In the following, we drop $f_{111}$ as it decouples and has a trivial
solution.

In @cite Schween2025, we derived an analytical solution for this system of
equations assuming a periodic box of length $L$. We obtained the following
solution for the expansion coefficients $f_{lms}$:

\begin{align*}
  f_{000}(t,x)  & = \sum_{n} \frac{c_n^2}{\sqrt{\omega^2 + c_n^2}}
  \left[
    \cos(\sqrt{\omega^2 + c_n^2} t) - 1 +
    \frac{\omega^2 + c_n^2}{c_n^2}
    \right]
  \left[ A_n \cos(k_n x) - B_n \sin(k_n x) \right]                      \\
  f_{110}(t,x)  & = \sum_{n} \frac{\omega c_n}{\sqrt{\omega^2 + c_n^2}}
  \left[ 1 - \cos(\sqrt{\omega^2 + c_n^2} t) \right]
  \left[ A_n \sin(k_n x) + B_n \cos(k_n x) \right]                      \\
  f_{100}(t, x) & = \sum_{n} c_n \sin(\sqrt{\omega^2 + c_n^2} t)
  \left[ A_n \sin(k_n x) + B_n \cos(k_n x) \right]
\end{align*}

$A_n$ and $B_n$ are the Fourier coefficients given by the initial condition for
$f_{000}$. The wave number $k_n = \frac{2 \pi n}{L}$ and the constant $c_n =
\frac{v}{\sqrt{3}} k_n$ are determined by the velocity $v$ and the box length $L$.

This analytic solution can be used to compare with the numerical solution at an
arbitrary time $t$ to assess the convergence of the numerical solution.

## Implementation {#implementation-convergence-study}

In the following, we will show the implementation of the convergence study
example in @sapphire. As always, the setup is divided up into two main files:

- [config.h](#config-convergence-study) which defines the setup of the VFP
  system and runtime parameters.
- [convergence-study.cpp](#main-convergence-study) which implements the `main`
  function to executes the simulation and compares the result to the analytic
  solution.

In this example, we will only highlight the main parts of the implementation.
For a comprehensive overview, please refer to the [scattering
only](#scattering-only) example.

### config.h {#config-convergence-study}

#### VFP equation {#dimension-convergence-study}

Given that the example considers a one-dimensional problem with decoupling of
the momenta, we can set the dimensionality of the problem to `dim = dim_ps =
dim_cs = 2`.

@snippet{lineno} examples/vfp/convergence-study/config.h Dimension

In the VFP equation, we only have the @ref sapphirepp::VFP::VFPFlags::rotation
"rotation", $q \mathbf{v} \cdot \left( \mathbf{B} \times \nabla_{p} f \right)$,
and @ref sapphirepp::VFP::VFPFlags::spatial_advection "spatial advection",
$(\mathbf{u} + \mathbf{v}) \cdot \nabla_{x} f$, terms. For the latter, we are in
the spatial case of a static background plasma, $\mathbf{u} = 0$. The magnetic
field is time independent, hence we can use
@ref sapphirepp::VFP::VFPFlags::time_independent_fields "time independent" flag.

@snippet{lineno} examples/vfp/convergence-study/config.h VFP Flags

#### Runtime Parameters {#parameter-convergence-study}

The main parameters describing out setup are the magnetic field strength $B_0$,
and the box length $L$. For simplicity, we hardcode the coefficients $A_n$ and
$B_n$ for the initial conditions. We leave the implementation of these as
runtime parameters as an exercise for the interested reader.

One peculiarity is, that the box length $L$ is already given by the domain size
defined in the @ref sapphirepp::VFP::VFPParameters "VFPParameters" class.
However, in order to access this parameter in the @ref
sapphirepp::VFP::InitialValueFunction "InitialValueFunction", we need to copy it
to the @ref sapphirepp::PhysicalParameters "PhysicalParameters" class. The same
goes for the velocity $v$, gamma factor $\gamma$, charge $q$, and mass $m$.

1. **Define**

   We start by defining the runtime parameters, including a copy of the @ref
   sapphirepp::VFP::VFPParameters "VFPParameters" $L$, $v$, $\gamma$, $q$, and
   $m$. We will copy their values in the [main
   function](#setup-convergence-study).

   @snippet{lineno} examples/vfp/convergence-study/config.h Define runtime parameter

2. **Declare**
  
   As in the [gyro-advection](#gyro-advection) example, we use the trick to use
   $\frac{B_0}{2 \pi}$ instead of $B_0$ as an input parameter. Since the main
   frequency in the problem is the gyrofrequency $\omega = \frac{B_0 q}{\gamma
    m}$, we want to ensure that a final time $T_{\rm final} = n \times T_g = n
    \times \frac{2 \pi}{B_0} \frac{\gamma m}{q}$ corresponds to exactly $n$ full
   gyroperiods.

   @snippet{lineno} examples/vfp/convergence-study/config.h Declare runtime parameter

3. **Parse**

   Last, we parse the runtime parameter $B_0$:

   @snippet{lineno} examples/vfp/convergence-study/config.h Parse runtime parameter

#### Analytic solution {#initial-condition-convergence-study}

We "misuse" the @ref sapphirepp::VFP::InitialValueFunction
"InitialValueFunction" to implement the analytic solution at all times. Setting
$t=0$ yields the initial condition, but setting $t>0$ yields the analytic
solution at a given time. To recall the analytic solution, we have:

\begin{align*}
  f_{000}(t,x)  & = \sum_{n} \frac{c_n^2}{\sqrt{\omega^2 + c_n^2}}
  \left[
    \cos(\sqrt{\omega^2 + c_n^2} t) - 1 +
    \frac{\omega^2 + c_n^2}{c_n^2}
    \right]
  \left[ A_n \cos(k_n x) - B_n \sin(k_n x) \right]                      \\
  f_{110}(t,x)  & = \sum_{n} \frac{\omega c_n}{\sqrt{\omega^2 + c_n^2}}
  \left[ 1 - \cos(\sqrt{\omega^2 + c_n^2} t) \right]
  \left[ A_n \sin(k_n x) + B_n \cos(k_n x) \right]                      \\
  f_{100}(t, x) & = \sum_{n} c_n \sin(\sqrt{\omega^2 + c_n^2} t)
  \left[ A_n \sin(k_n x) + B_n \cos(k_n x) \right]
\end{align*}

with $\omega = \frac{B_0 q}{\gamma m}$, $k_n = \frac{2 \pi n}{L}$ and $c_n =
\frac{v}{\sqrt{3}} k_n$. As mentioned before, we hardcode the coefficients $A_n$ and
$B_n$ for simplicity. We choose a simple sin wave $B_1 = 1$. To ensure
positivity of the analytic solution, we use the offset $A_0 = 2$.

We start the implementation by defining some useful constants. To access the
current time, we use the
@dealref{get_time(),classFunctionTime,ae7d37ddb04314b38cf67c6cba22923f6} method
of the @dealref{dealii::Function,classFunction} class. We define the
coefficients $A_n$ and $B_n$ as array, containing `A[0] = A_0`, `A[1] = A_1`,
etc.

@snippet{lineno} examples/vfp/convergence-study/config.h Initial value constants

Now, we implement the solution by looping over the modes $n=0$ until $n_{\rm
max}$. Note that we took special care in the $f_{000}$ term to avoid division by
$c_n$, which evaluates to zero for $n=0$.

@snippet{lineno} examples/vfp/convergence-study/config.h Initial value

#### Magnetic field {#magnetic-field-convergence-study}

The implementation of the magnetic field $\mathbf{B} = B_0 \hat{\mathbf{e}}_z$
is straightforward:f

@snippet{lineno} examples/vfp/convergence-study/config.h Magnetic field

### convergence-study.cpp {#main-convergence-study}

In the `main` function, we will implement the comparison of the numerical and
analytic solution. Because we want to compare the solution at all time steps, we
will not use the @ref sapphirepp::VFP::VFPSolver::run "run()" function, but
implement the time loop ourselves. We will also output add the analytic solution
to the output, showing a comparison between `projection` and `interpolation` of
the solution (see the [visualization](#visualization) section).

#### Setup {#setup-convergence-study}

The setup in the main function is fairly standard. We initialize MPI, and the
@ref sapphirepp::saplog "saplog" stream. We then parse the runtime parameters.

@snippet{lineno} examples/vfp/convergence-study/convergence-study.cpp Main function setup

One speciality is that we need to copy the runtime parameters $L$, $v$,
$\gamma$, $q$, and $m$ to the @ref sapphirepp::PhysicalParameters
"PhysicalParameters" class. In addition, we check if the assumptions of $l_{\rm
max} = 1$ and periodic boundary conditions are met.

@snippet{lineno} examples/vfp/convergence-study/convergence-study.cpp Copy VFP parameter

As we want to calculate the numerical error at every time step, we create a
`csv` file to store the results.

@snippet{lineno} examples/vfp/convergence-study/convergence-study.cpp Create error file

Now we create and set up the @ref sapphirepp::VFP::VFPSolver "VFPSolver" object.
Notice that we are not using the @ref sapphirepp::VFP::VFPSolver::run "run()"
function, but implement the time loop ourselves.

@snippet{lineno} examples/vfp/convergence-study/convergence-study.cpp Setup vfp_solver

Last, we have to prepare the analytic solution. To this end create an instance
of the @ref sapphirepp::VFP::InitialValueFunction "InitialValueFunction" and a
@dealref{Vector,classPETScWrappers_1_1MPI_1_1Vector}. In addition, we create a
@dealref{ComponentSelectFunction} in order to on;y compare the $f_{000}$ modes
when calculating the error.

@snippet{lineno} examples/vfp/convergence-study/convergence-study.cpp Setup analytic solution

#### Time loop {#time-loop-convergence-study}

Next comes the integral part of the simulation: the time loop. Every time step
we have to perform the following steps:

- Output the solution, numeric and analytic
- Calculate the error
- Advance the time, by performing an Euler or Runge-Kutta step

To keep track of the time steps, we use the @dealref{DiscreteTime} class.

@snippet{lineno} examples/vfp/convergence-study/convergence-study.cpp Time loop

We only want to output the solution every Nth time steps, where N is the
`output_frequency`. Ensuring this with an `if` statement, we add three
different vectors to the output:

- the numeric solution
- a projection of analytic solution onto the finite element DG space
- an interpolation of the analytic solution

@snippet{lineno} examples/vfp/convergence-study/convergence-study.cpp Output solution

Next we calculate the error and the norm of the solution. Notice that the
`weight` function allows us to only select the $f_{000}$ mode. Having both the
$L^2$ error, $\lVert f_{000, h}^{n} - f_{000}(t^n) \rVert_{L^2}$, and the $L^2$
norm, $\lVert f_{000, h}^{n} \rVert_{L^2}$, we can calculate the relative error
as

$$
  \text{rel. error} := \frac{\lVert f_{000, h}^{n} - f_{000}(t^n)
  \rVert_{L^2}}{\lVert f_{000, h}^{n} \rVert_{L^2}} \,.
$$

We output this results and save them to the `csv` file.

@snippet{lineno} examples/vfp/convergence-study/convergence-study.cpp Calculate error

Finally, we advance the time by performing an Euler or Runge-Kutta step. We
select the method according to the `time_stepping_method` runtime parameter, and
use the methods implemented in the @ref sapphirepp::VFP::VFPSolver "VFPSolver".

@snippet{lineno} examples/vfp/convergence-study/convergence-study.cpp Time step

#### End simulation {#end-simulation-convergence-study}

After finishing the time loop, we output the final results, and calculate the
final error.

@note In @sapphire, we output the results at the *start* of each time step.
  Therefore, we additionally need to output the final results after the last
  time step.

@snippet{lineno} examples/vfp/convergence-study/convergence-study.cpp Last time step

Finally, we can end the simulation and close the `csv` file.

@snippet{lineno} examples/vfp/convergence-study/convergence-study.cpp End simulation

## Results {#results-convergence-study}

Showing the results for a low resolution run, illustrates the difference between
*projection* and *interpolation* of the analytic solution. While the
interpolation gives the exact values at the nodes, the projection approximates
the solution in the finite element space. Using a discontinuous Galerkin method,
the projection is not continuous, but has jumps at the nodes. We can see that
the *weak* numerical solution is approximating the *projection* of the analytic
solution. Consequently, the numeric solution is not continuous and shows jumps
at the nodes. We discuss this in more detail in the
[visualization](#visualization) section.

![Results t0](https://sapphirepp.org/img/examples/convergence-study/convergence_study_t0.png)
![Results tEnd](https://sapphirepp.org/img/examples/convergence-study/convergence_study_tEnd.png)

Running the simulation for different time steps, $\Delta t$, we can assess the
convergence of the numerical solution. As expected, the `ERK4` method is fourth
order, while the forward Euler method is only first order. We can also study the
convergence of the numerical solution with spatial resolutions $\Delta
x$ for different polynomial degrees $k$. In this case, we expect the error to
follow a $(\Delta x)^{k+1}$ scaling, assuming the projection error dominates.

<p float="center">
  <img src="https://sapphirepp.org/img/examples/convergence-study/convergence_study_dt.png" alt="Convergence dt" width="49%">
  <img src="https://sapphirepp.org/img/examples/convergence-study/convergence_study_dx.png" alt="Convergence dx" width="49%">
</p>

## The plain program {#plain-convergence-study}

### config.h {#plain-config-convergence-study}

@include{lineno} examples/vfp/convergence-study/config.h

### convergence-study.cpp {#plain-main-convergence-study}

@include{lineno} examples/vfp/convergence-study/convergence-study.cpp

<div class="section_buttons">

| Previous              |
|:----------------------|
| [Examples](#examples) |

</div>

---

@author Florian Schulze (<florian.schulze@mpi-hd.mpg.de>)
@date 2024-03-07
