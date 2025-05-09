# Parallel shock {#parallel-shock}

@tableofcontents

## Introduction {#introduction-parallel-shock}

In this example we simulate the time-dependent transport and acceleration of
particles at a parallel shock.
The setup is similar to the [quick start](#quick-start), but we use an advanced
and more realistic configuration.
You can find a detailed description of the setup in @cite Schween2025.

The astrophysical shock is assumed to move along the $x$-axis.
In the rest frame of the shock wave, the background plasma's velocity is represented as
$\mathbf{U}_1 = U_{\rm sh} \hat{\mathbf{e}}_x$ upstream ($x<0$) and
$\mathbf{U}_1 = U_{\rm sh}/r \hat{\mathbf{e}}_x$ downstream ($x>0$).
Here, $r$ is the compression ratio of the shock and $U_{\rm sh}$ is its velocity.
As this represents a simple plain shock wave,
we only require a single dimension in configuration space for the simulation.
We also presuppose a constant magnetic field that is parallel to the shock normal,
i.e. $\mathbf{B} = B_0 \hat{\mathbf{e}}_x$.
This kind of setup is referred to as a parallel shock.

In astrophysical shocks, particles can be accelerated to high energies
through the process of diffusive shock acceleration @cite Drury1983 @cite Kirk1994.
This involves injecting particles at the shock with a momentum of $p_{\rm inj}$ at a rate $Q$.
These particles originate from the thermal background plasma and propagate around the shock.
These particles may be scattered by Alvén waves travelling through the background plasma,
a process we model with the help of scattering frequency $\nu$.
In @cite Schween2025 a constant scattering frequency is assumed,
i.e. it is independent of $p$ and $x$:

$$
  \nu = \nu_0  \, .
$$

For the case of a constant scattering frequency,
an analytic solution for the steady-state distribution function in mixed coordinates is presented in
@cite Drury1983 :

\f[
  f(x, p, \theta) = \frac{Q}{p_{\rm inj}U_{1}} \frac{3r}{r - 1}
  \left(\frac{p}{p_{\rm inj}}\right)^{-3r/(r - 1)} \times
  \begin{cases}
    e^{3U_{1} \nu / V^{2} x} \left[1 - 3 U_{1}/V \cos\theta\right]
    &\quad \text{ for } x < 0 \\
    1 &\quad \text{ for } x \geq 0
  \end{cases}
\f]

where $V = p/\gamma m$ is the particle velocity in mixed coordinates.
Decomposing this solution into spherical harmonics,
we find that the isotropic part of the distribution function is given by

\f[
  f_{000}(x, p) = \sqrt{4\pi} \frac{Q}{p_{\rm inj}U_{1}} \frac{3r}{r - 1}
  \left(\frac{p}{p_{\rm inj}}\right)^{-3r/(r - 1)} \times
  \begin{cases}
    e^{3U_{1} \nu / V^{2} x} &\quad \text{ for } x < 0 \\
    1 &\quad \text{ for } x \geq 0
  \end{cases}
\f]

Moreover, in the upstream we have a non-trivial anisotropy, namely

\f[
  f_{100}(x, p) = -3 U_{1}/V \sqrt{\frac{4\pi}{3}}
  \frac{Q}{p_{\rm inj}U_{1}} \frac{3r}{r - 1}
  \left(\frac{p}{p_{\rm inj}}\right)^{-3r/(r - 1)}
  e^{3U_{1} \nu / V^{2} x} \quad \text{ for } x < 0 \,,
\f]

while the downstream is isotropic, i.e.

$$
  f_{100} = 0 \quad \text{ for } x \geq 0 \,.
$$

An **approximate** analytic expression for the temporal evolution of the isotropic part of the distribution function $f$
has been given in @cite Drury1991 and @cite Forman1983 .
At the shock wave

\f[
  F(t, x = 0, p) = f_{000}(x = 0, p) \phi(t)
\f]

where $F$ is the time-dependent isotropic part of the distribution function and

\f[
  \phi(t) = \frac{1}{2}\left[\exp\left(\frac{2 c^{2}_{1}}{c_{2}}\right) \mathrm{erfc}\left(\sqrt{\frac{c^{3}_{1}}{2 t c_{2}}} + \sqrt{\frac{c^{\phantom{3}}_{1} t}{2 c_{2}}}\right)  +
\mathrm{erfc}\left(\sqrt{\frac{c^{3}_{1}}{2 t c_{2}}} - \sqrt{\frac{c^{\phantom{3}}_{1}t}{2 c_{2}}}\right)\right] \,.
\f]

The constants $c_1$ and $c_2$ are the mean and the variance of the acceleration time.
They are

\f[
  \begin{align}
  c_1 &= \frac{r}{2 U^{2}_{1} \nu} \frac{r + 1}{r - 1} \ln\left(\frac{1 + p^{2}}{1 + p^{2}_{\rm inj}}\right) \, \\
  c_2 &= \frac{1}{3\nu^{2}} \frac{r}{U^{4}_{1}} \frac{r^{3} + 1}{r - 1}
    \left[ \frac{1}{1 + p^{2}} - \frac{1}{1 + p^{2}_{\rm inj}} + \ln\left(\frac{1 + p^{2}}{1 + p^{2}_{\rm inj}}\right)\right] \,.
  \end{align}
\f]

More details can be found in @cite Schween2025 .

## Parameter files {#parameter-files-parallel-shock}

A feature of @sapphire, not covered in the [quick start](#quick-start) guide,
is the use of parameter files.
Its implementation uses the @dealref{ParameterHandler} class of @dealii.
Parameter files allow you to define simulation parameters that are set at runtime,
referred to hereafter as **runtime parameters**.
This is in contrast to compile-time parameters.

Runtime parameters are read from a parameter file,
such as `parameter.prm`, at the start of the simulation.
You can pass a parameter file to @sapphire as the **first** command-line argument, e.g.

```shell
./build/sapphirepp parameter-file.prm
```

This enables you to modify simulation parameters without recompiling the code,
which is particularly useful for parameter studies.
@sapphire includes several predefined runtime parameters,
such as the grid dimensions and the expansion order of the spherical harmonic expansion of $f$.

Parameters in @sapphire are organised into categories and subcategories,
primarily focusing on `Output` and `VFP` equation parameters.
Key parameters are:

| Parameter Name        | Category              | Description                                                    |
|:----------------------|:----------------------|:---------------------------------------------------------------|
| Results directory     | `Output`              | Specifies the directory for output files.                      |
| Format                | `Output`              | Defines the output format (e.g., `vtu/pvtu/hdf5`).             |
| Output frequency      | `Output`              | Determines the frequency of output.                            |
| Simulation identifier | `Output`              | Name of the simulation run, i.e. subfolder for the simulation. |
| Expansion order       | `VFP:Expansion`       | Sets the expansion order in spherical harmonics.               |
| Polynomial degree     | `VFP:Finite Elements` | Indicates the polynomial degree of the finite element method.  |
| Grid type             | `VFP:Mesh`            | Specifies the grid type (e.g., `Shock grid/Hypercube/File`).   |
| Final time            | `VFP:Time stepping`   | Sets the final time of the simulation.                         |
| Time step size        | `VFP:Time stepping`   | Determines the time step size of the simulation.               |
| Method                | `VFP:Time stepping`   | Defines the time stepping method (e.g., `CN/ERK4/...`).        |

A comprehensive and descriptive sample parameter file, `template-parameter.prm`,
is located in the `sapphirepp` folder.
This file is automatically generated when running @sapphire without any arguments.
Here's an example of the syntax used in `.prm` files:

```parameter
subsection Output
  set Results folder        = ./results
  set Format                = pvtu
  set Output frequency      = 1
  set Simulation identifier =
end
subsection VFP
  subsection Expansion
    set Expansion order = 1
  end
  subsection Finite element
    set Polynomial degree = 1
  end
  subsection Mesh
    set Grid type = Shock grid
  end
  subsection Time stepping
    set Method         = CN
    set Final time     = 200
    set Time step size = 1.0
  end
end
```

@note Note that section and parameter names are case-sensitive and can contain spaces.
Moreover, not all parameters listed in the template file are used.
It depends on the activated @ref sapphirepp::VFP::VFPFlags, which ones are actually relevant.

For more details on the predefined parameters in @sapphire
we also refer to the documentation of the @ref sapphirepp::Utils::OutputParameters
and the @ref sapphirepp::VFP::VFPParameters class.

@note The default behaviour of @sapphire is to overwrite results from previous runs.
We recommend using the `Simulation identifier` parameter to give each  run a unique ID,
e.g. `parallel-shock-1`.

Users can also define their own parameters to describe their physical setup.
An example of this can be found in the [Implementation](#parameter-parallel-shock) section.

## Implementation {#implementation-parallel-shock}

The implementation of this example can be found in the example directory `sapphirepp/examples/vfp/parallel-shock`.
The following sections explain it in detail.
As shown in the [quick start](#quick-start),
there are only a few lines in the `config.h` file that need to be adjusted.

### VFP equation {#dimension-parallel-shock}

As can be seen in the analytic solution presented above,
the distribution function $f$ depends only on $(x, p, \theta)$.
This means that we need a single dimension in configuration space.
Together with the $p$ dependence this results in a two-dimensional reduced phase space.
The $\theta$ dependence is handled by the spherical harmonics with $l>0$.
Because we don't have any dependence on $\varphi$, spherical harmonics with $m>0$ will be zero.

@snippet{lineno} examples/vfp/parallel-shock/config.h Dimension

Next, we need to specify which terms of the VFP equation are necessary to model the transport
and acceleration of particles at a parallel shock.
We actually have to solve the full VFP equation,

$$
  \frac{\partial f}{\partial t} + (\mathbf{U} + \mathbf{v}) \cdot \nabla_{x} f -
  \gamma m \frac{\mathrm{D} \mathbf{U}}{\mathrm{D} t} \cdot \nabla_{p}f -
  \mathbf{p} \cdot\nabla_{x} \mathbf{U}\cdot \nabla_{p} f +
  q \mathbf{v} \cdot \left( \mathbf{B} \times \nabla_{p} f \right) =
  \frac{\nu}{2} \Delta_{\theta, \varphi} f + S \,.
$$

We therefore need to activate all flags in the `vfp_flags` variable.
In addition, we know that the velocity field and magnetic field are time independent.
This also holds true for the source term.

@snippet{lineno} examples/vfp/parallel-shock/config.h VFP Flags

### User-defined runtime parameters {#parameter-parallel-shock}

This section outlines the process of implementing user-defined runtime parameters
that describe their physical setup.
The example parameters include:

| Name                 | Symbol        | Code                | Description                            |
|:---------------------|:--------------|:--------------------|:---------------------------------------|
| Shock velocity       | $U_{\rm sh}$  | `u_sh`              | Velocity of the shock.                 |
| Compression ratio    | $r$           | `compression_ratio` | Compression ratio of the shock.        |
| Injection momentum   | $p_{\rm inj}$ | `p_inj`             | Momentum of the injected particles.    |
| Injection rate       | $Q$           | `Q`                 | Rate of injected particles.            |
| Scattering frequency | $\nu_0$       | `nu0`               | Scattering frequency of the particles. |
| Magnetic field       | $B_0$         | `B0`                | Magnetic field strength.               |

Additionally, there are several numerical parameters that are not part of the physical setup:

| Name                   | Symbol        | Code          | Description                                    |
|:-----------------------|:--------------|:--------------|:-----------------------------------------------|
| Shock width            | $d_{\rm sh}$  | `shock_width` | Width of the shock.                            |
| Injection width in $p$ | $\sigma_p$    | `sig_p`       | Width of the injection in momentum.            |
| Injection width in $x$ | $\sigma_x$    | `sig_x`       | Width of the injection in configuration space. |
| Injection position     | $x_{\rm inj}$ | `x_inj`       | Position of the injection.                     |

Implementing these parameters in @sapphire involves a two-step process:

1. **Define the parameters**
   within the @ref sapphirepp::PhysicalParameters "PhysicalParameters" class
   in the `config.h` file.
   This step informs the compiler about their existence
   and sets default values to be used
   if the parameter is not specified in the parameter file.

   @snippet{lineno} examples/vfp/parallel-shock/config.h Define runtime parameter

2. **Add the parameters to the parameter file**
  to ensure that the parameter parser expects them
  and automatically sets their values.
  The @dealii class @dealref{ParameterHandler} `prm`
  and its @dealref{add_parameter,classParameterHandler,a04b75c02037d19fd7fd781785fcefc79} method
  are used for this purpose:

   ```cpp
   prm.add_parameter("entry", parameter, "Description of the parameter");
   ```

   The @ref sapphirepp::PhysicalParameters::declare_parameters() "declare_parameters()" method
   is edited to include the following code:

   @snippet{lineno} examples/vfp/parallel-shock/config.h Declare runtime parameter

@note An alternative approach for defining parameters,
      offering more fine-grained control over parsing,
      is introduced in the [gyro motion with advection](#parameter-gyro-advection) example.

### Scattering frequency {#scattering-frequency-parallel-shock}

We now implement the functions modelling the physical scenario.
We start with the scattering frequency.
As stated above, we assume a constant scattering frequency, i.e. $\nu = \nu_0$.
This translates to

@snippet{lineno} examples/vfp/parallel-shock/config.h Scattering frequency

The function calculates the scattering frequency at multiple `points` in the $(x,\ln p)$--domain.
We use the index `q_index` to refer to an individual point.
Notice that we access the runtime parameter $\nu_0$
via the (previously modified) @ref sapphirepp::PhysicalParameters "PhysicalParameters" class `prm`.

### Source term {#source-term-parallel-shock}

Furthermore, we would like to inject mono-energetic particles at the shock.
To model this we use a Gaussian profile centred at $(x_{\rm inj}, p_{\rm inj})$
with width $\sigma_x$ and $\sigma_p$,

$$
  S(x, p) = \frac{Q}{2\pi\sigma_x\sigma_p} \exp\left[
    -\frac{(x - x_{\rm inj})^2}{2\sigma_x^2} -
    \frac{(p - p_{\rm inj})^2}{2\sigma_p^2}
  \right] \, .
$$

As already explained in the [quick start](#quick-start),
the source term needs to be decomposed into spherical harmonics.
As we inject an isotropic distribution of particles, we can set all $l>0$ components to zero.
Notice that we use the index `i` to refer to the components of the spherical harmonic decomposition,
$i(l,m,s)$, at a single `point` in the $(x,\ln p)$--domain.

@snippet{lineno} examples/vfp/parallel-shock/config.h Source

We emphasise that the second component of point, i.e. `point[1]`,
corresponds to the momentum variable $\ln p$.

### Magnetic field {#magnetic-field-parallel-shock}

As described above, we use a constant magnetic field parallel to the shock,
$\mathbf{B} = B_0 \hat{\mathbf{e}}_x$.
It does not influence the steady-state solution, but nevertheless we include it for completeness.

@snippet{lineno} examples/vfp/parallel-shock/config.h Magnetic field

### Velocity field {#velocity-parallel-shock}

The velocity field is the same as in the [quick start](#quick-start) example.
We use a $\tanh$ profile to model the transition at the shock wave.
The only modification is, that we now use the user-defined runtime parameters $U_{\rm sh}$,
$r$ and $d_{\rm sh}$.
Otherwise, $\mathbf{U}(\mathbf{x})$, its divergence $\nabla \cdot \mathbf{U}(\mathbf{x})$,
the material derivative $\mathrm{D} \mathbf{U}/\mathrm{D} t$
and the Jacobian $\partial U_{x}/\partial x$ stay the same.

1. Background velocity field value $\mathbf{U}(\mathbf{x})$:

   @snippet{lineno} examples/vfp/parallel-shock/config.h Background velocity value

2. Background velocity divergence $\nabla \cdot \mathbf{U}(\mathbf{x})$:

   @snippet{lineno} examples/vfp/parallel-shock/config.h Background velocity divergence

3. Background velocity material derivative $\mathrm{D}\mathbf{U}/\mathrm{D} t$:

   @snippet{lineno} examples/vfp/parallel-shock/config.h Background velocity material derivative

4. Background velocity Jacobian $\partial U_{x}/ \partial x$:

   @snippet{lineno} examples/vfp/parallel-shock/config.h Background velocity Jacobian

### Compile and run {#compile-parallel-shock}

To run the simulation, implement the above functions in the `config.h` file and recompile @sapphire:

```shell
make --directory=build
```

Alternatively, you can use the implementation in the `examples/vfp/parallel-shock` folder.
This executable will be named `parallel-shock` instead of `sapphirepp`.
Note that for this to work, @sapphire must be configured with the `-DEXAMPLES=ON` option.

```shell
cmake -S . -B build -DEXAMPLES=ON
make --directory=build parallel-shock
```

We recommend running the simulation with the parameters given in
`examples/vfp/parallel-shock/parameter.prm`:

@include examples/vfp/parallel-shock/parameter.prm

Run the simulation with:

```shell
mpirun -n 4 ./build/examples/vfp/parallel-shock/parallel-shock examples/vfp/parallel-shock/parameter.prm
```

We provide a [Paraview Python](#paraview-python) script to generate the plots:

```shell
pvbatch examples/vfp/parallel-shock/pvplot.py results/parallel-shock
```

## Results {#results-parallel-shock}

Before we analyse the results,
we check that the solution has converged to a steady state.
To demonstrate this, we show the time series of the 2D data in a logarithmic scale.

![2D time series](https://sapphirepp.org/img/examples/parallel-shock/parallel-shock-2D.gif)
![Plot of shock region](https://sapphirepp.org/img/examples/parallel-shock/shock-region.png)

From the first figure we can conclude that the solution reached steady state even for large momenta.
The second figure shows a cut-out of the solution in the shock region.
This serves an illustrative purpose demonstrating the stretched `Shock grid` used in the simulation.

The analytic solution predicts that the energy spectrum of the particles at shock
follows a $p^{-3r/(r - 1)}$ power law for $p > p_{\rm inj}$.
Given a compression ratio $r=4$, we anticipate a $p^{-4}$ power law.
We plot $f_{000}(x = 0, \ln p)$ and compare it with our expectation.

![Steady state f(ln(p)) plot](https://sapphirepp.org/img/examples/parallel-shock/particle-spectrum.png)

In the upstream region ($x<0$), we expect an exponential cut-off,
which is dependent on the scattering frequency $\nu$.
As discussed in the introduction, we also observe an upstream anisotropy, $f_{100}$.
Note that the anisotropic part disappears in the downstream region.
When comparing our simulation results to the analytic solution, we find they are in good agreement.
The figure below shows $f_{000}(x, p = 10)$ and $f_{100}(x, p = 10)$:

![Steady state f(x) plot](https://sapphirepp.org/img/examples/parallel-shock/spatial-distribution.png)

In the last figure we present the temporal evolution of the energy spectrum at the shock
for a specific $p$, namely $p=10$.
The simulation results are compared to the **approximate** analytical solution
mentioned in the introduction.

![f(t) plot](https://sapphirepp.org/img/examples/parallel-shock/temporal-evolution.png)

We note that there exists an **exact** analytical solution
for a specific choice of the scattering frequency in the up- and downstream,
see, for example, eq. 16 and 21 in @cite Drury1991 .
In @cite Dann2025 it is demonstrated that @sapphire reproduces it accurately.

<div class="section_buttons">

| Previous              |
|:----------------------|
| [Examples](#examples) |

</div>

---

@author Florian Schulze (<florian.schulze@mpi-hd.mpg.de>)
@author Nils Schween (<nils.schween@mpi-hd.mpg.de>)
@date 2025-03-11
