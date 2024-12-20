# Parallel shock {#parallel-shock}

@tableofcontents

## Introduction {#introduction-parallel-shock}

In this example we simulate the transport and acceleration of particles in a
parallel shock. The setup is similar to the [quick start](#quick-start), but we
use an advanced and more realistic configuration. You can find a detailed
description of the setup and an analytic solution in @cite Schween2025.

We're examining an astrophysical shock moving along the $x$-axis. The background
plasma's velocity is represented as $\mathbf{U}_1 = U_{\rm sh}
\hat{\mathbf{e}}_x$ upstream ($x<0$) and $\mathbf{U}_1 = U_{\rm sh}/r
\hat{\mathbf{e}}_x$ downstream ($x>0$). Here, $r$ stands for the compression
ratio and $U_{\rm sh}$ for the shock velocity. As this represents a simple plain
shock, we only require a single dimension in configuration space for the
simulation. We also employ a constant magnetic field, parallel to the shock
normal, denoted as $\mathbf{B} = B_0 \hat{\mathbf{e}}_x$. This kind of setup is
referred to as a parallel shock.

In astrophysical shocks, particles can be accelerated to high energies through
the process of diffusive shock acceleration @cite Drury1983 @cite Kirk1994. This
involves injecting particles at the shock with a momentum of $p_{\rm inj}$ at a
rate $Q$. These particles originate from the thermal background plasma and
propagate around the shock. Particles scatter at Alvén waves in the magnetic
field, a process we describe by the scattering frequency $\nu$. In our scenario,
the scattering frequency exhibits a minor momentum dependency, which can be
represented as:

$$
  \nu = \nu_0 \left(\frac{p}{p_0}\right)^{1/3} \, .
$$

@todo We can also just use $\nu = \nu_0 \left(\frac{p}{p_0}\right)^{\alpha}$
  with $\alpha$ a free parameter.

In @cite Schween2025 we derive an analytic solution for the distribution
function in mixed coordinates:

\f[
  f(x, p, \theta) = \frac{Q}{p_{\rm inj}U_{1}} \frac{3r}{r - 1}
  \left(\frac{p}{p_{\rm inj}}\right)^{-3r/(r - 1)} \times
  \begin{cases}
    e^{3U_{1} \nu / V^{2} x} \left[1 - 3 U_{1}/V \cos\theta\right]
    \quad \text{ for } x < 0 \\
    1 \quad \text{ for } x \geq 0
  \end{cases}
\f]

where $V = \frac{p}{\gamma m}$ is the particle velocity in mixed coordinates.
Decomposing this solution into spherical harmonics, we find that the isotropic
part of the distribution function is given by

\f[
  f_{000}(x, p) = \sqrt{4\pi} \frac{Q}{p_{\rm inj}U_{1}} \frac{3r}{r - 1}
  \left(\frac{p}{p_{\rm inj}}\right)^{-3r/(r - 1)} \times
  \begin{cases}
    e^{3U_{1} \nu / V^{2} x} \quad \text{ for } x < 0 \\
    1 \quad \text{ for } x \geq 0
  \end{cases}
\f]

In the upstream we have an anisotropy,

$$
  f_{100}(x, p) = -3 U_{1}/V \sqrt{\frac{4\pi}{3}}
  \frac{Q}{p_{\rm inj}U_{1}} \frac{3r}{r - 1}
  \left(\frac{p}{p_{\rm inj}}\right)^{-3r/(r - 1)}
  e^{3U_{1} \nu / V^{2} x} \quad \text{ for } x < 0 \,,
$$

while the downstream is isotropic,

$$
  f_{100} = 0 \quad \text{ for } x \geq 0 \,.
$$

## Parameter files {#parameter-files-parallel-shock}

A key feature of @sapphire, not covered in the [quick start](#quick-start)
guide, is the use of parameter files. These files allow you to define the
simulation parameters that are set at runtime, referred to hereafter as
**runtime parameters**. This is in contrast to compile-time parameters.

Runtime parameters are read from a parameter file, such as `parameter.prm`, at
the start of the simulation. You can pass a parameter file to @sapphire using
the following syntax:

```shell
  ./sapphirepp parameter-file.prm
```

This feature enables you to modify the parameters without needing to recompile
the code, which is particularly useful for running parameter studies. @sapphire
includes several predefined runtime parameters, such as the grid dimensions and
the expansion order in spherical harmonics.

Parameters in @sapphire are organized into categories and subcategories,
primarily focusing on `Output` and `VFP` equation parameters. Key parameters
include:

| Parameter Name | Category | Description |
|:--------------|:---------|:------------|
| Results directory | `Output` | Specifies the directory for output files. |
| Format | `Output` | Defines the output format (e.g., `vtu/pvtu/hdf5`). |
| Output frequency | `Output` | Determines the frequency of output. |
| Simulation identifier | `Output` | Name of the simulation run, i.e. subfolder for the simulation. |
| Expansion order | `VFP:Expansion` | Sets the expansion order in spherical harmonics. |
| Polynomial degree | `VFP:Finite Elements` | Indicates the polynomial degree of the finite element method. |
| Grid type | `VFP:Mesh` | Specifies the grid type (e.g., `Shock grid/Hypercube/File`). |
| Final time | `VFP:Time stepping` | Sets the final time of the simulation. |
| Time step size | `VFP:Time stepping` | Determines the time step size of the simulation. |
| Method | `VFP:Time stepping` | Defines the time stepping method (e.g., `CN/ERK4/...`). |

A comprehensive and descriptive sample parameter file, `template-parameter.prm`,
is located in the `sapphirepp` folder. This file is automatically generated when
running @sapphire without any arguments. Here's an example of the syntax used in
`.prm` files:

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

@note Note that section and parameter names are case-sensitive and can contain
  spaces.

For more details on the predefined parameters in @sapphire we also refer to the
documentation of the @ref sapphirepp::Utils::OutputParameters and the @ref
sapphirepp::VFP::VFPParameters class.

@note The default behaviour of @sapphire is to overwrite results from previous
  runs. We recommend using the `Simulation identifier` parameter to give each
  run a unique ID, e.g. `parallel-shock-1`.

Users can also define their own parameters to describe the physical setup. An
example of this can be found in the [implementation](#parameter-parallel-shock)
section.

## Implementation {#implementation-parallel-shock}

An implementation of this example can be found in the example directory
`sapphirepp/examples/parallel-shock`. The following sections will explain the
implementation in detail. As shown in the [quick start](#quick-start), there are
only a few lines in the `config.h` file that need to be adjusted.

### VFP equation {#dimension-parallel-shock}

As we see in the analytic solution, the distribution function $f$ depends only
on $(x, p, \theta)$. This means that we only need a single dimension in
configuration space. Together with the $p$ dependence this results in a
two-dimensional reduced phase space. The $\theta$ dependence is handled by the
spherical harmonics with $l>0$. Because we don't have any
dependence on $\varphi$, spherical harmonics with $m>0$ will be zero.

@snippet{lineno} examples/parallel-shock/config.h Dimension

Next, we need to specify which parts of the VFP equation are present in our
problem. Here, we have the full VFP equation,

$$
  \frac{\partial f}{\partial t} + (\mathbf{u} + \mathbf{v}) \cdot \nabla_{x} f -
  \gamma m \frac{\mathrm{D} \mathbf{u}}{\mathrm{D} t} \cdot \nabla_{p}f -
  \mathbf{p} \cdot\nabla_{x} \mathbf{u}\cdot \nabla_{p} f +
  q \mathbf{v} \cdot \left( \mathbf{B} \times \nabla_{p} f \right) =
  \frac{\nu}{2} \Delta_{\theta, \varphi} f + S \,.
$$

We therefore need to list all flags in the `vfp_flags` variable. In addition, we
know that the velocity field and magnetic field time independent, as well as the
source.

@snippet{lineno} examples/parallel-shock/config.h VFP Flags

### Defining Custom Runtime Parameters {#parameter-parallel-shock}

This section outlines the process of defining custom runtime parameters that
describe the physical setup. The parameters include:

| Name | Symbol | Code | Description |
|:-----|:-------|:-----|:------------|
| Shock velocity | $U_{\rm sh}$ | `u_sh` | Velocity of the shock. |
| Compression ratio | $r$ | `compression_ratio` | Compression ratio of the shock. |
| Injection momentum | $p_{\rm inj}$ | `p_inj` | Momentum of the injected particles. |
| Injection rate | $Q$ | `Q` | Rate of injected particles. |
| Scattering frequency | $\nu_0$ | `nu0` | Scattering frequency of the particles. |
| Magnetic field | $B_0$ | `B0` | Magnetic field strength. |

Additionally, there are several numerical parameters that are not part of the
physical setup:

| Name | Symbol |  Code | Description |
|:-----|:-------|:------|:------------|
| Shock width | $d_{\rm sh}$ | `shock_width` | Width of the shock. |
| Injection width in $p$ | $\sigma_p$ | `sig_p` | Width of the injection in momentum. |
| Injection width in $x$ | $\sigma_x$ | `sig_x` | Width of the injection in configuration space. |
| Injection position | $x_{\rm inj}$ | `x_inj` | Position of the injection. |

Implementing these parameters in @sapphire involves a three-step process:

1. **Define** the parameters within the @ref sapphirepp::PhysicalParameters
   "PhysicalParameters" class in the `config.h` file. This step informs the
   compiler about these parameters.

   @snippet{lineno} examples/parallel-shock/config.h Define runtime parameter

2. **Declare** the parameters in the parameter file. This step ensures that the
   parameter parser expects these parameters. The @dealii concept of a
   @dealref{ParameterHandler} `prm` and the
   @dealref{declare_entry,classParameterHandler,a6d65f458be69e23a348221cb67fc411d}
   function are used for this purpose:

   ```cpp
   prm.declare_entry("name", "default value",    dealii::Patterns::Type(),
                     "Description of the parameter");
   ```

   All parameters related to the source are collected in a subsection. This is
   achieved by entering the subsection before declaring the parameters and
   leaving it afterwards. The @ref
   sapphirepp::PhysicalParameters::declare_parameters() "declare_parameters()"
   function is edited to include the following code:

   @snippet{lineno} examples/parallel-shock/config.h Declare runtime parameter

3. **Parse** the values defined in the parameter file. This step sets the values
   of the parameters according to the parameter file. The
   @dealref{ParameterHandler} getter function
   @dealref{get_double(),classParameterHandler,aeaf3c7846747695b1f327677e3716ec5}
   is utilized inside the @ref
   sapphirepp::PhysicalParameters::parse_parameters() "parse_parameters()"
   function and the following code is included:

   @snippet{lineno} examples/parallel-shock/config.h Parse runtime parameter

### Scattering frequency {#scattering-frequency-parallel-shock}

Now we implement the functions related to the physical setup, starting with the
scattering frequency. As stated above, we use a scattering with a weak momentum
dependency,

$$
  \nu = \nu_0 \left(\frac{p}{p_0}\right)^{1/3} \, .
$$

To get the momentum variable $\ln p$ we use the second component of the `point`
in reduced phase space that is given to function. One peculiarity here is, that
the function calculates the scattering frequency at multiple `points`. We use
the index `q_index` to refer to an individual point:

@snippet{lineno} examples/parallel-shock/config.h Scattering frequency

Notice that we use access the runtime parameter $\nu_0$ via the (previously
modified) @ref sapphirepp::PhysicalParameters "PhysicalParameters" class `prm`.

### Source term {#source-term-parallel-shock}

We want to inject mono-energetic particles at the shock. To model this we use a
Gaussian centred at $(x_{\rm inj}, p_{\rm inj})$ with width $\sigma_x$ and
$\sigma_p$,

$$
  S(x, p) = \frac{Q}{2\pi\sigma_x\sigma_p} \exp\left[
    -\frac{(x - x_{\rm inj})^2}{2\sigma_x^2} -
    \frac{(p - p_{\rm inj})^2}{2\sigma_p^2}
  \right] \, .
$$

As already explained in the [quick start](#quick-start), the source term needs
to be decomposed into spherical harmonics. As want to inject an isotropic
distribution of particles, we can set all $l>0$ components to zero. Notice that
here we use the index `i` to refer to the components of the spherical harmonics,
$i(l,m,s)$, at one `point`.

@snippet{lineno} examples/parallel-shock/config.h Source

### Magnetic field {#magnetic-field-parallel-shock}

As described above, we use a constant magnetic field parallel to the shock,
$\mathbf{B} = B_0 \hat{\mathbf{e}}_x$. It does not influence the steady state
solution, but nevertheless we include it for completeness.

@snippet{lineno} examples/parallel-shock/config.h Magnetic field

### Velocity field {#velocity-parallel-shock}

The velocity field is the same as in the [quick start](#quick-start) example. We
use a $\tanh$ profile to model the transition at the shock. The only
modification is, that we now use the custom runtime parameters $U_{\rm sh}$, $r$
and $d_{\rm sh}$. Otherwise, of the value $\mathbf{u}(\mathbf{x})$, the
divergence $\nabla \cdot \mathbf{u}(\mathbf{x})$, the material derivative
$\frac{\mathrm{D} \mathbf{u}}{\mathrm{D} t}$ and the Jacobian $\frac{\partial
u_{x}}{\partial x}$ stays the same.

1. Background velocity field value $\mathbf{u}(\mathbf{x})$:

   @snippet{lineno} examples/parallel-shock/config.h Background velocity value

2. Background velocity divergence $\nabla \cdot \mathbf{u}(\mathbf{x})$:

   @snippet{lineno} examples/parallel-shock/config.h Background velocity divergence

3. Background velocity material derivative $\frac{\mathrm{D}
   \mathbf{u}}{\mathrm{D} t}$:

   @snippet{lineno} examples/parallel-shock/config.h Background velocity material derivative

4. Background velocity Jacobian $\frac{\partial u_{x}}{\partial x}$:

   @snippet{lineno} examples/parallel-shock/config.h Background velocity Jacobian

### Compile and run {#compile-parallel-shock}

To run the simulation, implement the above functions in the `config.h` file and
recompile @sapphire:

```shell
  cd sapphirepp/build
  make
```

Alternatively, you can use the implementation in the `examples/parallel-shock`
folder. This executable will be named `parallel-shock` instead of `sapphirepp`.
Note that for this to work, @sapphire must be configured with the
`-DEXAMPLES=ON` option.

```shell
  cd sapphirepp/build/examples/parallel-shock
  make parallel-shock
```

We recommend running the simulation with the parameters given in
`examples/parallel-shock/parameter.prm`:

@include examples/parallel-shock/parameter.prm

Run the simulation with:

```shell
  cd sapphirepp/build/examples/parallel-shock
  ./parallel-shock parameter.prm
```

## Results {#results-parallel-shock}

The analytic solution predicts that the distribution function will follow a
$p^{-3r/(r - 1)}$ power law for $p > p_{\rm inj}$. Given a compression ratio
$r=4$, we anticipate a $p^{-4}$ power law. To visualize these results, we plot
$p^4 f(x,p)$, which should yield an approximately constant value in the
downstream region ($x>0$). In the upstream region ($x<0$), we expect an
exponential cut-off, which is dependent on the scattering frequency $\nu$.

![2D time series](https://sapphirepp.org/img/examples/parallel-shock-2d.gif)

As discussed in the introduction, we also observe an upstream anisotropy,
$f_{100}$. When comparing our simulation results to the analytic solution, we
find they are in good agreement. The figure below showcases $p^4 f(\ln(p))$:

![Steady state f(ln(p)) plot](https://sapphirepp.org/img/examples/parallel-shock-f-lnp.png)

In the last figure we present $p^4 f(x)$. Note that the anisotropic part
disappears in the downstream region.

![Steady state f(x) plot](https://sapphirepp.org/img/examples/parallel-shock-f-x.png)

<div class="section_buttons">

| Previous              |
|:----------------------|
| [Examples](#examples) |

</div>

---

@author Florian Schulze (<florian.schulze@mpi-hd.mpg.de>)
@date 2023-12-20
