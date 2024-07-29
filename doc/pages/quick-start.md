# Quick-start {#quick-start}

@tableofcontents

## Introduction {#introduction-quick-start}

This quick-start guide will help you to get started with @sapphire. Note that it
will not show all features of @sapphire, but covers the basics you need to write
your own application in @sapphire. We will demonstrate @sapphire on a simplified
version of the CR acceleration at a parallel shock. For a more detailed
implementation including advanced features of sapphire we refer to the
[examples](#examples).

Let's get started by considering the following setup: We have a planar shock that
the particles will cross, to gain energy by the first-order Fermi mechanism
@cite Drury1983 @cite Kirk1994 @cite Schween2024. We will prescribe the
background plasma in the shock frame, resulting in a stationary velocity field
$\mathbf{u}(\mathbf{x})$ with the shock compression ratio $r$ at $0$:

\begin{align}
  \mathbf{u}(\mathbf{x}) = u_1 \hat{\mathbf{e}}_x \text{, for } x < 0 \\
  \mathbf{u}(\mathbf{x}) = r \cdot u_1 \hat{\mathbf{e}}_x \text{, for } x > 0
\end{align}

The transport of energetic charged (test-)particles in this plasma is described
by the Vlasov-Fokker-Planck (VFP) equation:

$$
  \frac{\partial f}{\partial t} + (\mathbf{u} + \mathbf{v}) \cdot \nabla_{x} f -
  \gamma m \frac{\mathrm{D} \mathbf{u}}{\mathrm{D} t} \cdot \nabla_{p}f -
  \mathbf{p} \cdot\nabla_{x} \mathbf{u}\cdot \nabla_{p} f =
  \frac{\nu}{2} \Delta_{\theta, \varphi} f + S\, .
$$

Here, $f(\mathbf{x}, \mathbf{p})$ is the distribution function of the
energetic particles with momentum $\mathbf{p} = \gamma m \mathbf{v}$.
$\frac{\mathrm{D}}{\mathrm{D} t} = \frac{\partial}{\partial t} + \mathbf{u}
\cdot \frac{\partial}{\partial \mathbf{x}}$ denotes the material derivative,
$\nu(p)$ and $S(\mathbf{x}, \mathbf{p})$ the scattering frequency and source
respectively.

@sapphire efficiently solves this equation by exploiting that the distribution
function $f$ is in many physical environments almost isotropic. It thus
expresses $f$ as a series of spherical harmonics, i.e.

$$
 f(t, \mathbf{x}, \mathbf{p}) = \sum^{l_{\rm max}}_{l = 0} \sum^{l}_{m = 0}
 \sum^{1}_{s = 0} f_{lms}(t, \mathbf{x}, p) Y_{lms}(\theta,\varphi) \,,
$$

whose zeroth order term is the isotropic part and higher order terms represent
the anisotropies of the distribution function. @sapphire computes the expansion
coefficients $f_{lms}(\mathbf{x}, p, t)$ up to a maximum order $l_{\rm max}$.
Hence, it reduces the dimensionality of the problem from the full six
dimensional phase space $(\mathbf{x}, \mathbf{p})$ to a four dimensional
**reduced phase space** $\mathbf{\xi} = (\mathbf{x}, p)^{T}$. We can reduce the
dimensionality even further, noticing that our problem only depends on one
spatial dimension. Our spatial domain for this problem $(\mathbf{x}) = (x)$,
called **configuration space**, has therefore the dimensionality `dim_cs = 1`.
The **computational domain** for this problem therefore has the dimensionality
`dim = dim_ps = dim_cs + 1 = 2`, equalling the dimensionality of the reduced
phase space `dim_ps`.

@note The computational domain is always the reduced phase space. Therefore, the
      dimension of the grid and (most) quantities `dim`, equals the dimension of
      the reduced phase space, `dim = dim_ps`.

In practice this means, that if we refer to a `point` in @sapphire, the first
`dim_cs` components are the spatial coordinates $\mathbf{x} = (x,y,z)$ and the
last component is the momentum $p$. Hence, we have in this example:

```cpp
    const double x = point[0];
    const double p = point[1];
```

@note Due to limitations of @dealii, the maximum computational dimension is
      currently limited to `dim = dim_ps = 3`. This means that one either has to
      assume an at most 2d spatial domain (`dim_cs = 2` $\implies$ `dim_ps =
      3`), or a three-dimensional but momentum independent problem (`dim_ps =
      dim_cs = 3`).

## Implementation {#implementation-quick-start}

After this somewhat lengthy introduction, let's start to implement the above
stated equation in @sapphire. To this end, we have to modify the
`sapphirepp/include/config.h` file. This file contains the physical setup for
the simulation. In total this file is over 600 lines long, but as a user you can
ignore most of it as boilerplate code. There are in practise very few lines you have to
adapt to your problem, marked with a `// !!!EDIT HERE!!!` comment. Let's start
going through these lines step by step.

### Dimensionality and VFPFlags {#dimension-quick-start}

First, we have to specify the dimension of the **reduced phase space** for the
problem. As elaborated above, the configuration space for this problem is $1$
dimensional, corresponding two a 2d reduced phase space. We therefore set the
`dimension` of the computational domain to $2$ in `sapphirepp/include/config.h`
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

In total, we have four effects at play that are important. First, we have @ref
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
direction $\mathbf{p} = p \left( \sin \theta \cos \varphi, \sin \theta \sin
\varphi, \cos \theta \right)$ of the particles but not their energy. Last, we
have the injection of particles at a @ref sapphirepp:VFP::VFPFlags::source
"source" $S$. In @sapphire, we activate these effects by listing them in the @ref
sapphirepp::VFP::VFPFlags "vfp_flags" variable (separated by a bitwise or `|`):

@snippet{lineno} include/config.h VFP Flags

As you can see, we added two additional flags, @ref
sapphirepp::VFP::VFPFlags::time_independent_fields and @ref
sapphirepp::VFP::VFPFlags::time_independent_source, specifying that both the
velocity field and the source are time independent. This results in noticeable
performance gains, as @sapphire only needs to calculate the corresponding
DG-matrices once.

### Scattering frequency {#scattering-frequency-quick-start}

Now that the general problem is stated, we have to specify the scattering
frequency $\nu(p)$, the source $\mathbf{S}(x,p)$ and the background velocity
field $\mathbf{u}(x)$. Let's start with the @ref
sapphirepp::VFP::ScatteringFrequency "scattering frequency". For simplicity, we
assume that its independent of momentum, and set it to $\nu = 100$. For a more
realistic setup with a momentum dependence we refer to the advanced [parallel
shock example](#parallel-shock). The corresponding implementation in @sapphire
is only line in `sapphirepp/include/config.h`:

@snippet{lineno} include/config.h Scattering frequency

### Source term {#source-term-quick-start}

Next, we specify the injection of particles at the shock ($x = 0$). We use a
Gaussian distribution to inject particles at a momentum $p_0 = 1$:

$$
  S(x,p) = \frac{Q}{2 \pi \sigma_p \sigma_x}
           e^{-\frac{\left(p - p_0 \right)^2}{2 \sigma_p^2}}
           e^{-\frac{x^2}{2 \sigma_x^2}} \,.
$$

We choose an amplitude of $Q=0.1$, and standard deviations $\sigma_x = 0.1$ and
$\sigma_p = 0.001$. In @sapphire, we have to specify the coefficient
$S_{i(l,m,s)}(t, \mathbf{x}, p)$ in reduced phase space of a spherical harmonic
decomposition of the source in full phase space $S(t, \mathbf{x}, \mathbf{p})$.
Since we inject particles isotropically, we just have

$$
 S(t, \mathbf{x}, \mathbf{p}) = \sum^{l_{\rm max}}_{l = 0} \sum^{l}_{m = 0}
 \sum^{1}_{s = 0} S_{lms}(t, \mathbf{x}, p) Y_{lms}(\theta,\varphi)
 = S_{000}(x, p) Y_{lms}(\theta,\varphi)
 = \frac{1}{\sqrt{4 \pi}} S_{000}(x, p) \,.
$$

Notice the factor $\frac{1}{\sqrt{4 \pi}}$ we have to account for.

@snippet{lineno} include/config.h Source

### Velocity field {#velocity-quick-start}

We use a $\tanh$ profile to prescribe the transition of the background velocity
in the shock region:

$$
  \mathbf{u}(x) = \frac{u_s}{2r} \left( (1-r) \tanh\left(\frac{x}{d_s}\right)
  + (1+r)\right) \hat{\mathbf{e}}_x\,.
$$

Here, $r = 4$ is the compression ratio, $u_s = 0.1$ the shock velocity and $d_s
= 0.001$ the shock width. Even though we specified the problem is 1d in
configuration space, @sapphire allows prescribing flows outside the plane by
setting $u_y$ and $u_z$-components. In our case, we only have a flow in
$x$-direction, so we can set the other components to zero:

@snippet{lineno} include/config.h Background velocity value

As we can see above, the VFP equation not only depends on the value of velocity,
but also on different combinations of its derivatives $\frac{\partial
u_{x_i}}{\partial x_j}$. Therefore, we also have to implement these expressions
in @sapphire. We start by prescribing the divergence, which we can compute to
be,

$$
  \nabla \cdot \mathbf{u}(x) = \frac{u_s}{2 r} \frac{1-r}{d_s}
  \left(1 - \tanh\left(\frac{x}{d_s}\right)^2\right) \,.
$$

@snippet{lineno} include/config.h Background velocity divergence

The material derivative is given by,

$$
  \frac{\mathrm{D} \mathbf{u}}{\mathrm{D} t}
  = \frac{\partial \mathbf{u}}{\partial t}
  + \mathbf{u} \cdot \frac{\partial \mathbf{u}}{\partial \mathbf{x}}
  = \left(\frac{u_s}{2r}\right)^2 \left( (1-r) \tanh\left(\frac{x}{d_s}\right)
  + (1+r)\right) \frac{1-r}{d_s} \left(1 -
  \tanh\left(\frac{x}{d_s}\right)^2\right) \,.
$$

@snippet{lineno} include/config.h Background velocity material derivative

Last, we also need the Jacobian $\frac{\partial u_{x_i}}{\partial x_j}$ of the
velocity field to compute the numerical flux. We only have a $x$-$x$-component
in the Jacobian which is given by:

$$
  \frac{\partial u_{x}}{\partial x} = \frac{u_s}{2 r} \frac{1-r}{d_s}
  \left(1 - \tanh\left(\frac{x}{d_s}\right)^2\right) \,.
$$

@snippet{lineno} include/config.h Background velocity Jacobian

### Compile and run {#compile-quick-start}

After making modifications to the physical setup, you'll need to recompile
@sapphire. We utilize [CMake](https://cmake.org/) for this purpose, adhering to
the principle of keeping source code and executables separate. As part of this
process, we'll create a `sapphirepp/build` directory to house the compiled
executable. Execute the following commands to compile @sapphire:

```shell
  cd sapphirepp
  mkdir build
  cmake -S . -B build
  cd build
  make
```

Once compiled, you can run the executable inside the `build` folder using the
command below:

```shell
  ./sapphirepp
```

To use multiple processors, you can use `mpirun`:

```shell
  mpirun -np N ./sapphirepp
```
  
where `N` is the number of processors you want to use.

## Results {#results-quick-start}

Running the simulation will produce a `results` directory containing the output
as `solution_*.pvtu` files. These files can be visualized using visualization
software like [ParaView](https://www.paraview.org) or
[VisIt](https://visit-dav.github.io/visit-website/).

Opening the files as a time series, we can see that they contain the expansion
coefficients `f_lms`. Let's have a look at the zeroth order coefficient
$f_{l=0,m=0,s=0}(x, \ln p)$/`f_000`. For better visibility, we show a 2d
log-plot of $p^4 f_{000}(x,p)$:

![2d time series](https://sapphirepp.org/img/examples/quick-start-2d.gif)

We can see that the particles are injected at $\ln p_0 = 0$ and gain energy by
diffusive shock acceleration. The particles are advected with the background
flow, and the shock is visible as a discontinuity in the distribution function
at $x = 0$. The particles are injected at a constant rate, and the distribution
function reaches a steady state after some time.

A detailed discussion of the steady state solution is found in the [parallel
shock example](#parallel-shock). Here, we can notice that the distribution
function follows a $p^4$ power law for $p > p_0$:

![Steady state f(ln(p)) plot](https://sapphirepp.org/img/examples/quick-start-f-lnp.png)

In the $f(x)$ figure blow we see that, the distribution function is constant in
the downstream ($x > 0$). In the upstream ($x < 0$) it has an exponential cut
off:

![Steady state f(x) plot](https://sapphirepp.org/img/examples/quick-start-f-x.png)

To continue, we recommend having a look at the [parallel shock
example](#parallel-shock). There, we will discuss how to use and define runtime
parameters, among other topics.

<div class="section_buttons">

| Previous       |                            Next |
|:---------------|--------------------------------:|
| [Home](#index) | [Visualization](#visualization) |

</div>

---

@author Florian Schulze (<florian.schulze@mpi-hd.mpg.de>)
@date 2023-12-14
