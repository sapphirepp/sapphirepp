# Quick-start {#quick-start}

@tableofcontents

## Introduction {#introduction-quick-start}

This quick-start guide will help you get started with @sapphire. While it does
not cover all the features of @sapphire, it provides the basics to write your
first application. We will demonstrate @sapphire using a simplified version of
CR acceleration at a parallel shock. Note that this example is primarily
intended as a technical demonstration. For a physically accurate model of CR
acceleration in a parallel shock, refer to the [parallel shock
example](#parallel-shock). We strongly recommend going through the [parallel
shock example](#parallel-shock) after finishing this
[quick-start](#quick-start), as it introduces additional concepts in @sapphire,
such as parameter files. Furthermore, we provide a [tutorial for visualizing
results using ParaView](#paraview-tutorial) as well as general [tips and tricks
for visualization](#visualization-tips). More examples that demonstrate advanced
features of @sapphire can be found in the [examples](#examples) section.

The example we are considering uses the following setup: We have a planar shock
that the particles will cross, to gain energy by the first-order Fermi mechanism
@cite Drury1983 @cite Kirk1994 @cite Schween2024b. We will prescribe the
background plasma in the shock frame, resulting in a stationary velocity field
$\mathbf{u}(\mathbf{x})$ with the shock compression ratio $r$ at $x = 0$:

\begin{align}
  \mathbf{u}(\mathbf{x}) = \hat{\mathbf{e}}_x
  \begin{cases}
      u_{\mathrm{sh}}           & \text{if} & x < 0 \\
      \frac{u_{\mathrm{sh}}}{r} & \text{if} & x > 0
  \end{cases}
\end{align}

The transport of energetic charged (test-)particles in this plasma is described
by the Vlasov-Fokker-Planck (VFP) equation:

$$
  \frac{\partial f}{\partial t} + (\mathbf{u} + \mathbf{v}) \cdot \nabla_{x} f -
  \gamma m \frac{\mathrm{D} \mathbf{u}}{\mathrm{D} t} \cdot \nabla_{p}f -
  \mathbf{p} \cdot\nabla_{x} \mathbf{u}\cdot \nabla_{p} f =
  \frac{\nu}{2} \Delta_{\theta, \varphi} f + S\, .
$$

Here, $f(\mathbf{x}, \mathbf{p})$ is the distribution function of the energetic
particles with momentum $\mathbf{p} = \gamma m \mathbf{v}$.
$\frac{\mathrm{D}}{\mathrm{D} t} = \frac{\partial}{\partial t} + \mathbf{u}
\cdot \frac{\partial}{\partial \mathbf{x}}$ denotes the material derivative,
$\nu(p)$ and $S(\mathbf{x}, \mathbf{p})$ the scattering frequency and source
respectively. (We neglect the magnetic field, as the gyro-motion is
perpendicular to the plane).

@sapphire efficiently solves this equation by exploiting that the distribution
function $f$ is in many physical environments almost isotropic. It thus
expresses $f$ as a series of spherical harmonics, i.e.

$$
 f(t, \mathbf{x}, \mathbf{p}) = \sum^{l_{\rm max}}_{l = 0} \sum^{l}_{m = 0}
 \sum^{1}_{s = 0} f_{lms}(t, \mathbf{x}, p) Y_{lms}(\theta,\varphi) \,,
$$

where $\mathbf{p} = p \, (\cos\theta, \sin\theta \cos\varphi, \sin\theta
\sin\varphi)^T$ is used. The zeroth order term is the isotropic part and higher
order terms represent the anisotropies of the distribution function. @sapphire
computes the expansion coefficients $f_{lms}(\mathbf{x}, p, t)$ up to a maximum
order $l_{\rm max}$. Hence, it reduces the dimensionality of the problem from
the full six dimensional phase space $(\mathbf{x}, \mathbf{p})$ to a four
dimensional **reduced phase space** $\mathbf{\xi} = (\mathbf{x}, p)^{T}$. We can
reduce the dimensionality even further, noticing that our problem only depends
on one spatial dimension. Our spatial domain for this problem $(\mathbf{x}) =
(x)$, called **configuration space**, has therefore the dimensionality `dim_cs =
1`. The **computational domain** for this problem therefore has the
dimensionality `dim = dim_ps = dim_cs + 1 = 2`, equalling the dimensionality of
the reduced phase space `dim_ps`.

@note The computational domain is always the reduced phase space. Therefore, the
      dimension of the grid and (most) quantities `dim`, equals the dimension of
      the reduced phase space, `dim = dim_ps`.

In practice this means, that if we refer to a `point` in @sapphire, the first
`dim_cs` components are the spatial coordinates $\mathbf{x} = (x,y,z)$ and the
last component corresponds to the momentum $p$. Since @sapphire uses a
logarithmic momentum variable $\ln(p)$ by default, we have the following in this
example:

```cpp
    const double x     = point[0];
    const double log_p = point[1];
```

@note Due to limitations of @dealii, the maximum computational dimension is
      currently limited to `dim = dim_ps = 3`. This means that one either has to
      assume an at most 2D spatial domain (`dim_cs = 2` $\implies$ `dim_ps =
      3`), or a three-dimensional but momentum independent problem (`dim_ps =
      dim_cs = 3`).

## Implementation {#implementation-quick-start}

Now that we have covered the introduction, we can implement the equation in
@sapphire. To do this, we need to modify the `sapphirepp/include/config.h` file.
This file contains the physical setup for the simulation. You can find the
complete source code on
[GitHub](https://github.com/sapphirepp/sapphirepp/blob/main/include/config.h).
Although the file is over 600 lines long, most of it is boilerplate code that
you can ignore. The lines you need to adapt for your specific problem are marked
with `// !!!EDIT HERE!!!`. In the following we go through these lines step by
step.

### Dimensionality and VFPFlags {#dimension-quick-start}

First, we have to specify the dimension of the **reduced phase space** for the
problem. As elaborated above, the configuration space for this problem is 1
dimensional, corresponding to a 2D reduced phase space. We therefore set the
`dimension` of the computational domain to 2 in `sapphirepp/include/config.h`
(see the corresponding line numbers below):

@snippet{lineno} include/config.h Dimension

Next, we have to specify which terms of the VFP equation are used. Because we do
not consider a magnetic field in this example, the VFP equation is,

$$
  \frac{\partial f}{\partial t} + (\mathbf{u} + \mathbf{v}) \cdot \nabla_{x} f -
  \gamma m \frac{\mathrm{D} \mathbf{u}}{\mathrm{D} t} \cdot \nabla_{p}f -
  \mathbf{p} \cdot\nabla_{x} \mathbf{u}\cdot \nabla_{p} f =
  \frac{\nu}{2} \Delta_{\theta, \varphi} f + S\, .
$$

In total, we have four effects at play. First, we have @ref
sapphirepp::VFP::VFPFlags::spatial_advection "spatial advection", $(\mathbf{u} +
\mathbf{v}) \cdot \nabla_{x} f$, due to the background plasma velocity
$\mathbf{u}(\mathbf{x},t)$ and the particle velocity $\mathbf{v}(\mathbf{p})$.
Next, we have a @ref sapphirepp::VFP::VFPFlags::momentum "momentum" term,
$\left( \gamma m \frac{\mathrm{D} \mathbf{u}}{\mathrm{D} t} + \mathbf{p}
\cdot\nabla_{x} \mathbf{u} \right) \cdot \nabla_{p} f$, describing the
(perceived) acceleration of particles in the mixed coordinate frame due to
changes in background velocity. The third effect are @ref
sapphirepp::VFP::VFPFlags::collision "collisions" (i.e. elastic binary
interactions) of the particles that we can model with an effective scattering
term, $\frac{\nu}{2} \Delta_{\theta, \varphi} f$, which only changes the
direction $\theta, \varphi$ of the particles but not their energy. Last, we have
the injection of particles at a @ref sapphirepp:VFP::VFPFlags::source "source"
$S$. In @sapphire, we activate these effects by listing them in the @ref
sapphirepp::VFP::VFPFlags "vfp_flags" variable (separated by a bitwise-or `|`
and prefixed with the namespace `VFPFlags::`):

@snippet{lineno} include/config.h VFP Flags

As you can see, we added two additional flags, @ref
sapphirepp::VFP::VFPFlags::time_independent_fields "time_independent_fields" and
@ref sapphirepp::VFP::VFPFlags::time_independent_source
"time_independent_source", specifying that both the velocity field and the
source are time independent. This results in noticeable performance gains, as
@sapphire only needs to calculate the corresponding DG-matrices once instead of
every time step. A full list of all possible flags can be found in the @ref
sapphirepp::VFP::VFPFlags documentation.

### Scattering frequency {#scattering-frequency-quick-start}

Now that we have stated the general problem, we need to specify the scattering
frequency $\nu(p)$, the source $\mathbf{S}(x,p)$, and the background velocity
field $\mathbf{u}(x)$. We start with the @ref
sapphirepp::VFP::ScatteringFrequency "scattering frequency". For simplicity, we
assume it is independent of momentum. Using @ref
sapphirepp::VFP::ReferenceValues "dimensionless units", we set

$$
  \nu = 100 \,.
$$

It is worth mentioning the units used in @sapphire. @sapphire employs
dimensionless units throughout the code, meaning both the input (in terms of
user-defined functions in `config.h`) and the output are dimensionless.
Generally, length scales are defined in terms of the gyro-radius of a reference
particle, $x^{*} = x/r_{g,0}$, and time is given in multiples of the reference
gyro-frequency, $t^{*} = t \omega_{g,0}$. Momentum is defined in terms of the
rest mass of a test particle, $p^{*} = p/m_{0} c$, while velocity is given as a
fraction of $c$, $v^{*} = v/c$. An overview of the reference values can be found
in @ref sapphirepp::VFP::ReferenceValues.

We note again, that this setup serves as a technical demonstration. A more
realistic setup with a momentum-depend scattering frequency is shown in the
[parallel shock example](#parallel-shock). In this example, the implementation
in @sapphire corresponds to only one line in `sapphirepp/include/config.h`:

@snippet{lineno} include/config.h Scattering frequency

@note You might notice the index `q_index` in the line above. For performance
      reasons, the @ref sapphirepp::VFP::ScatteringFrequency
      "ScatteringFrequency" function calculates the scattering frequencies for a
      list of `points` rather than for one individual `point`. Consequently, the
      function returns a list of scattering frequency values, each corresponding
      to a point in the list, indexed by `q_index`.

### Source term {#source-term-quick-start}

Next, we specify the injection of particles at the shock, $x = 0$. In analytic
models, particles are typically injected isotropically at a discrete momentum
$p_0$ @cite Drury1983 @cite Kirk1994 @cite Schween2024b,

$$
  S(x, \mathbf{p}) = Q \delta(x) \delta(p - p_0) \,.
$$

Here, $Q$ represents the injection rate, and $\delta$ denotes the delta
distribution. Numerically, we approximate these delta distributions as Gaussians
with standard deviations $\sigma_x$ and $\sigma_p$, respectively
@cite Schween2024b,

$$
  S(x, \mathbf{p}) = \frac{Q}{2 \pi \sigma_p \sigma_x}
  e^{-\frac{\left(p - p_0 \right)^2}{2 \sigma_p^2}}
  e^{-\frac{x^2}{2 \sigma_x^2}} \,.
$$

To implement this source function in @sapphire, we decompose it into its
coefficients $S_{lms}$ using spherical harmonics, similar to the distribution
function,

$$
  S(t, \mathbf{x}, \mathbf{p}) = \sum^{l_{\rm max}}_{l = 0} \sum^{l}_{m = 0}
  \sum^{1}_{s = 0} S_{lms}(t, \mathbf{x}, p) Y_{lms}(\theta,\varphi) \,.
$$

Since particles are injected isotropically in this example, only $S_{000}$
contributes:

\begin{align}
  & S(x, \mathbf{p}) = S_{000}(x, p) Y_{000}(\theta,\varphi)
    = \frac{1}{\sqrt{4 \pi}} S_{000}(x, p) \\
  \implies
  & S_{000}(x, p) = \sqrt{4 \pi} \frac{Q}{2 \pi \sigma_p \sigma_x}
    e^{-\frac{\left(p - p_0 \right)^2}{2 \sigma_p^2}}
    e^{-\frac{x^2}{2 \sigma_x^2}} \\
  & S_{lms}(x,p) = 0 \,, \quad \text{otherwise}
\end{align}

For this example, we use an injection momentum of $p_0 = 1$, an amplitude of $Q
= 0.1$, and standard deviations $\sigma_x = 0.1$ and $\sigma_p = 0.001$, again
using @ref sapphirepp::VFP::ReferenceValues "dimensionless units". In @sapphire,
this corresponds to the following implementation:

@snippet{lineno} include/config.h Source

Given the complexity of this implementation, additional comments may be helpful.
First, we define the parameters $p_0$, $Q$, $\sigma_p$, and $\sigma_x$ as `const
double` with their respective values. Next, we retrieve the values of $x$ and
$p$ at the current `point`. As noted earlier, `point[0]` corresponds to $x$,
while `point[1]` corresponds to $\ln(p)$.

In the next step, it is important to note that @sapphire collapses the three
indices $l$, $m$, and $s$ into a single index $i(l,m,s)$. This approach is
advantageous as it allows to store all coefficients $S_{lms} = S_{i(l,m,s)}$ in
a single @dealref{Vector}. To convert from the index `i` to the indices $l$,
$m$, and $s$, we can utilize the `lms_indices` map (see @ref
sapphirepp::VFP::PDESystem::create_lms_indices() "here" for more details).

Finally, we implement the coefficients of the source term $S_{lms}$ by
separating the isotropic part, $l = m = s = 0$, from the vanishing anisotropic
part.

### Velocity field {#velocity-quick-start}

Last, we need to implement the background velocity field
$\mathbf{u}(\mathbf{x})$. Since the velocity field of a shock, as given above,
is discontinuous, we instead parametrize it using a $\tanh$ profile for
numerical implementation @cite Schween2024b,

$$
  \mathbf{u}(x) = \frac{u_{\mathrm{sh}}}{2r} \left(
  (1-r) \tanh\left(\frac{x}{d_{\mathrm{sh}}}\right)
  + (1+r)\right) \hat{\mathbf{e}}_x\,.
$$

Here, $r$ is the compression ratio, $u_{\mathrm{sh}}$ is the shock velocity, and
$d_{\mathrm{sh}}$ is the shock width. This profile smooths the shock over a
width $d_{\mathrm{sh}}$, recovering the discontinuous profile in the limit
$d_{\mathrm{sh}} \to 0$,

\begin{align}
  \mathbf{u}(\mathbf{x}) \overset{d_{\mathrm{sh}} \to 0}{\to}
  \frac{u_{\mathrm{sh}}}{2r} \hat{\mathbf{e}}_x
  \begin{cases}
      -(1-r) + (1+r) & \text{if} & x < 0 \\
      (1-r) + (1+r)  & \text{if} & x > 0
  \end{cases}
  = \frac{u_{\mathrm{sh}}}{2r} \hat{\mathbf{e}}_x
  \begin{cases}
      2r & \text{if} & x < 0 \\
      2  & \text{if} & x > 0
  \end{cases} \,.
\end{align}

In this example, we choose the following set of parameters in @ref
sapphirepp::VFP::ReferenceValues "dimensionless units": $r = 4$, corresponding
to a strong shock, $u_{\mathrm{sh}} = 0.1$, and $d_{\mathrm{sh}} = 0.001$. The
implementation in @sapphire is given by:

@snippet{lineno} include/config.h Background velocity value

Note that even though we specified the problem as 1D in configuration space,
@sapphire requires us to specify $u_y$ and $u_z$ components. These correspond to
a homogeneous flow out of the plane. In this example, we set them to zero and
only specify the $u_x$ (i.e., `velocity[0]`) component. Again, we can access the
$x$ value via `point[0]`.

Investigating the VFP equation above, we notice that it not only depends on the
value of velocity but also on different combinations of its derivatives
$\frac{\partial u_{x_i}}{\partial x_j}$. Therefore, we also need to implement
these expressions in @sapphire. We start by prescribing the divergence, which
can be computed as,

$$
  \nabla \cdot \mathbf{u}(x) = \frac{u_{\mathrm{sh}}}{2 r} \frac{1-r}{d_{\mathrm{sh}}}
  \left(1 - \tanh\left(\frac{x}{d_{\mathrm{sh}}}\right)^2\right) \,.
$$

@snippet{lineno} include/config.h Background velocity divergence

@note The @ref sapphirepp::VFP::BackgroundVelocityField::divergence_list
      "divergence_list" function uses a list of `points` to calculate the
      divergence for multiple points at once. Therefore, we use the index
      `q_index` to access individual points and their corresponding divergence
      values. The $x$ value is consequently given by `points[q_index][0]`.

Similarly, we can calculate and implement the material derivative,

$$
  \frac{\mathrm{D} \mathbf{u}}{\mathrm{D} t}
  = \frac{\partial \mathbf{u}}{\partial t}
  + \mathbf{u} \cdot \frac{\partial \mathbf{u}}{\partial \mathbf{x}}
  = \left(\frac{u_{\mathrm{sh}}}{2r}\right)^2 \left( (1-r) \tanh\left(\frac{x}{d_{\mathrm{sh}}}\right)
  + (1+r)\right) \frac{1-r}{d_{\mathrm{sh}}} \left(1 -
  \tanh\left(\frac{x}{d_{\mathrm{sh}}}\right)^2\right) \,.
$$

@snippet{lineno} include/config.h Background velocity material derivative

Last, we also need the Jacobian $\frac{\partial u_{x_i}}{\partial x_j}$ of the
velocity field to compute the numerical flux. We only have a $x$-$x$-component
in the Jacobian which is given by:

$$
  \frac{\partial u_{x}}{\partial x} = \frac{u_{\mathrm{sh}}}{2 r} \frac{1-r}{d_{\mathrm{sh}}}
  \left(1 - \tanh\left(\frac{x}{d_{\mathrm{sh}}}\right)^2\right) \,.
$$

@snippet{lineno} include/config.h Background velocity Jacobian

<br>

This concludes the implementation of the example in @sapphire. The next step is
to compile the code to prepare it for running the simulation.

### Compile and run {#compile-quick-start}

After modifying `sapphirepp/include/config.h` to implement the physical setup,
we need to recompile @sapphire. First, we use [CMake](https://cmake.org/) to
locate and link all dependencies and prepare the compilation process. To keep
source code and executables separate, we create a `sapphirepp/build` directory
to store the compiled executable and related files. To prepare for compilation,
execute the following commands (this may be redundant if you have already
followed the instructions in [Compiling @sapphire](#compilation)):

```shell
cd sapphirepp
mkdir build
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
```

Upon successful execution, you should see output similar to:

```text
Sapphire++ version 1.1.0
-- The C compiler identification is GNU 11.4.0
-- The CXX compiler identification is GNU 11.4.0
...
-- Using the deal.II-9.5.2 installation found at /home/schulze/.local/lib/dealii
...
-- Preparing Sapphire++
-- Preparing VFP module
-- Configuring done (0.3s)
-- Generating done (0.0s)
-- Build files have been written to: /home/schulze/Documents/PhD/Code/sapphirepp/build
```

<br>

Next, we have to compile @sapphire. Change to the `sapphirepp/build` directory
and execute:

```shell
cd build
make
```

To speed up the compilation, you can use multiple processes in parallel by
running:

```shell
make -j N
```

where `N` is the number of processes to use. You should see output similar to
the following:

```text
[ 10%] Building CXX object src/utils/CMakeFiles/UtilsLib.dir/output-parameters.cpp.o
[ 20%] Building CXX object src/utils/CMakeFiles/UtilsLib.dir/sapphirepp-logstream.cpp.o
[ 30%] Linking CXX static library libUtilsLib.a
[ 30%] Built target UtilsLib
[ 40%] Building CXX object CMakeFiles/sapphirepp.dir/sapphirepp.cpp.o
[ 50%] Building CXX object CMakeFiles/sapphirepp.dir/src/vfp/particle-functions.cpp.o
[ 60%] Building CXX object CMakeFiles/sapphirepp.dir/src/vfp/pde-system.cpp.o
[ 70%] Building CXX object CMakeFiles/sapphirepp.dir/src/vfp/upwind-flux.cpp.o
[ 80%] Building CXX object CMakeFiles/sapphirepp.dir/src/vfp/vfp-parameters.cpp.o
[ 90%] Building CXX object CMakeFiles/sapphirepp.dir/src/vfp/vfp-solver.cpp.o
[100%] Linking CXX executable sapphirepp
[100%] Built target sapphirepp
```

<br>

Once compiled, we can run the `sapphirepp` executable inside the
`sapphirepp/build` folder using the following command:

```shell
./sapphirepp
```

A successful execution will produce output similar to the following:

```text
Sapphire::Start Sapphire++ v1.1.0 with 1 MPI process(es) [2024/7/31 10:43:50]
Sapphire:VFP::Setup VFP equation solver.        [10:43:50]
Sapphire:VFP::Time step      0 at t = 0.00000   [10:43:50]
Sapphire:VFP::Time step      1 at t = 1.00000   [10:43:51]
Sapphire:VFP::Time step      2 at t = 2.00000   [10:43:51]
Sapphire:VFP::Time step      3 at t = 3.00000   [10:43:51]
Sapphire:VFP::Time step      4 at t = 4.00000   [10:43:51]
Sapphire:VFP::Time step      5 at t = 5.00000   [10:43:52]
Sapphire:VFP::Time step      6 at t = 6.00000   [10:43:52]
Sapphire:VFP::Time step      7 at t = 7.00000   [10:43:52]
Sapphire:VFP::Time step      8 at t = 8.00000   [10:43:53]
Sapphire:VFP::Time step      9 at t = 9.00000   [10:43:53]
Sapphire:VFP::Time step     10 at t = 10.0000   [10:43:53]
...
Sapphire:VFP::Time step    195 at t = 195.000   [10:44:46]
Sapphire:VFP::Time step    196 at t = 196.000   [10:44:47]
Sapphire:VFP::Time step    197 at t = 197.000   [10:44:47]
Sapphire:VFP::Time step    198 at t = 198.000   [10:44:47]
Sapphire:VFP::Time step    199 at t = 199.000   [10:44:48]
Sapphire:VFP::Simulation ended at t = 200.000   [10:44:48]
```

<br>

@sapphire also supports high performance multiprocessing with
[Open-MPI](https://www.open-mpi.org). To activate parallel computation on `N`
processors, use the `mpirun` command:

```shell
mpirun -np N ./sapphirepp
```

You can verify if the parallelization is successful by checking the startup
message of @sapphire:

```text
Sapphire::Start Sapphire++ v1.1.0 with N MPI process(es) [2024/7/31 10:43:50]
...
```

## Results {#results-quick-start}

Running the simulation will create a `results` directory in `sapphirepp/build`.
Listing the contents of this folder with:

```shell
ls results
```

we see the following files:

```text
log.prm
solution_0000.pvtu  solution_0000.0.vtu solution_0000.1.vtu ... solution_0000.N-1.vtu
solution_0001.pvtu  solution_0001.0.vtu solution_0001.1.vtu ... solution_0001.N-1.vtu
solution_0002.pvtu  solution_0002.0.vtu solution_0002.1.vtu ... solution_0002.N-1.vtu
solution_0003.pvtu  solution_0003.0.vtu solution_0003.1.vtu ... solution_0003.N-1.vtu
...
solution_0199.pvtu  solution_0199.0.vtu solution_0199.1.vtu ... solution_0199.N-1.vtu
solution_0200.pvtu  solution_0200.0.vtu solution_0200.1.vtu ... solution_0200.N-1.vtu
```

The first file, `log.prm`, contains the parameters used for the run. An
introduction to parameter files is given in the [parallel
shock](#parameter-files-parallel-shock) example. The other files,
`solution_XXXX.*`, contain the simulation results. The `_XXXX` suffix indicates
the time step the result corresponds to. Each time step has multiple associated
files:

- `solution_XXXX.pvtu`: Contains metadata for the results at each time step,
  primarily collecting the distributed data from different processors in a
  multiprocessor run.
- `solution_XXXX.n.vtu`: Contains part of the solution handled by processor
  number `n`. In a single-threaded run, you will only have the
  `solution_XXXX.0.vtu` files. Otherwise, the index `n` runs up to `N-1`, with
  `N` being the number of processors.

@note @sapphire supports other output formats, like `hdf5`.

@note @sapphire overwrites previous results upon rerunning. We recommend using a
      unique identifier for each run. For more details, refer to the
      documentation on [Parameter files](#parameter-files-parallel-shock) and
      @ref sapphirepp::Utils::OutputParameters::simulation_id
      "OutputParameters".

<br>

To visualize the results, we recommend using dedicated visualization software
like @paraview or [VisIt](https://visit-dav.github.io/visit-website/). A
[detailed tutorial](#paraview-tutorial) on how to visualise the results using
@paraview is given in the next section on [Visualization](#visualization). Here,
we want to present the results to give a first impression.

When opening the files as a time series, the data is displayed on a 2D grid. The
x-coordinate of this grid corresponds to $x$, while the y-coordinate corresponds
to $\ln p$. The output includes a scalar field for each expansion coefficient
$f_{lms}(x, \ln p)$, labelled `f_lms`. The most significant component is the
zeroth-order coefficient, $f_{000}$, which represents the isotropic part of the
distribution function. In the animation below, we illustrate how $f_{000}(t, x,
\ln p)$ evolves over time using a logarithmic colour scale:

![2d time series](https://sapphirepp.org/img/quick-start/quick-start-2D.gif)

We can see, that the particles are injected at $\ln p_0 = 0$ and gain energy
through diffusive shock acceleration, populating the phase space $\ln p > 0$.
Since the particles do not suffer any losses, the low momentum phase space ($\ln
p < 0$) remains unpopulated. The particles are advected downstream ($x > 0$)
with the background flow, while their density in the upstream ($x < 0$) is
exponentially suppressed. As the injection rate is constant, the distribution
function reaches a steady state by the end of the simulation.

@note As mentioned above, the output of @sapphire is in @ref
sapphirepp::VFP::ReferenceValues "dimensionless units".

We again want to emphasize, that this example mainly serves as technical
demonstration and therefore only uses a low-fidelity simulation with coarse
resolution. A physically accurate, high resolution simulation can be found in
the [parallel shock example](#parallel-shock). There, we also provide a detailed
discussion of the steady state solution and comparison with analytic models.

To further investigate the $x$ dependence of the solution, we can examine a
cut-out along $x$ at constant $\ln p$:

![Steady state f(x) plot](https://sapphirepp.org/img/quick-start/quick-start-f-x.png)

Again, we can see that $f_{000}$ is constant for downstream ($x > 0$) but decays
exponentially upstream, $f_{000}(x) \propto e^{x}$ for $x < 0$. Additionally,
the figure illustrates the non-vanishing dipole anisotropy $f_{100}(x)$. It
exhibits a discontinuity at $x = 0$: For $x > 0$, the dipole anisotropy vanishes
($f_{100}(x) = 0$), indicating an isotropic particle distribution. For $x < 0$,
$f_{100}$ is negative ($f_{100}(x) < 0$), which corresponds to an excess of
particles moving in negative $x$ direction relative to the background flow.
However, since @sapphire employs mixed coordinates, this does not necessarily
mean that particles are preferably moving towards the upstream. Instead, it
suggests that particles are "swimming against the stream".

@note @sapphire uses **mixed coordinates**. This means that the momentum
      $\mathbf{p}$ is defined in the frame of the background plasma
      $\mathbf{u}(\mathbf{x})$. One can transform between the momentum in a
      laboratory frame $\tilde{\mathbf{p}}$ and the mixed coordinate momentum
      $\mathbf{p}$ using: $\mathbf{p} = \tilde{\mathbf{p}} - m
      \mathbf{u}(\mathbf{x})$.

When comparing with analytic models of parallel shock acceleration, we expect
the momentum ($p$) dependence of the distribution function to follow a power-law
for $p > p_0$, specifically $f_{000}(p) \propto p^{-4}$ @cite Drury1983
@cite Kirk1994 @cite Schween2024b. To demonstrate that this example accurately
captures this behaviour, we present a log-log plot of $p^4 f_{000}$, using a
cut-out at constant $x$:

![Steady state f(ln(p)) plot](https://sapphirepp.org/img/quick-start/quick-start-f-lnp.png)

<br>

To continue, we recommend reviewing the [ParaView tutorial](#paraview-tutorial)
and the [parallel shock example](#parallel-shock). In the [ParaView
tutorial](#paraview-tutorial), we provide a beginner-friendly guide on creating
the plots shown above using @paraview. The [parallel shock
example](#parallel-shock) discusses how to use and define runtime parameters,
among other topics. Additionally, it serves as a high-fidelity physical model of
acceleration in parallel shocks.

<div class="section_buttons">

| Previous       |                            Next |
|:---------------|--------------------------------:|
| [Home](#index) | [Visualization](#visualization) |

</div>

---

@author Florian Schulze (<florian.schulze@mpi-hd.mpg.de>)
@date 2024-07-31
