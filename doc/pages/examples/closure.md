# Closure {#closure}

@tableofcontents

## Introduction {#introduction-closure}

This example is designed to give a feeling
for the effect of truncating the expansion at $l_{\rm max}$.
Furthermore, we will introduce the post-processing module
to reconstruct the full distribution function
$f(\mathbf{x}, p, \theta, \varphi)$.
To this end, we follow the example introduced in @cite Schween2024, Sec. 4.3.

We start off with an isotropic distribution of mono-energetic particles,
homogeneous in $y$--$z$--plane with a Gaussian profile in $x$ direction.
This distribution will evolve according to the following differential equation
$$
  \frac{\mathrm{d} f}{\mathrm{d} t} =
  \frac{\partial f}{\partial t} + V_{x}\frac{\partial f}{\partial x} =
  \frac{\nu}{2} \Delta_{\theta,\varphi} f \,
$$
with the initial condition
$$
  f(t=0,x) = \frac{1}{\sqrt{2 \pi} \sigma}
    \exp\left(-\frac{x^{2}}{2\sigma^{2}}\right) \,.
$$

Since we explore the effects of scattering only at the end of this example,
we start with $\nu = 0$.
Because all particles have the same energy,
we have to fix their Lorentz factor.
We choose $\gamma = 2$,
which is equivalent to setting the magnitude of their velocities to
$V = \sqrt{3}/2$.
We note that the projected velocity along the $x$--axis
depends on $\theta$, $V_{x} = V \cos\theta$.

The solution to the evolution equation,
keeping in mind that $\nu = 0$,
can be computed with the method of characteristics and is
$$
  f(t,x) = \frac{1}{\sqrt{2 \pi} \sigma} \exp\left(-\frac{(x - V \cos\theta t)^{2}}{2\sigma^{2}} \right) \,.
$$
Particles with values of $\theta$ close to zero or $\pi$
move faster along the $x$-axis
than particles with $\theta$ values around $\pi/2$.
This leads to a separation of particles with different directions of motion,
encoded by $\theta$.
This is demonstrated in the figure below.

@todo Figure from paper
![Analytical velocity distribution](https://sapphirepp.org/img/examples/parallel-shock-f-lnp.png)
Analytical velocity distribution of the particles at three different space and time points.

We plot a heatmap of the average particle number density at $(x=0,t=0)$,
at $(x=3 \sigma, t = 3 \sigma/ V = 6/\sqrt{3})$
and $(x=6 \sigma, t= 6 \sigma /V = 12/\sqrt{3})$.
We choose the points such that a particle with velocity $V_{x} = \sqrt{3}/2$
starting at $x = 0$ has reached three or six standard deviations $\sigma$
respectively.
In the second and third plot the separation of particles with different $\theta$
shows up as a hot spot at the pole of the sphere.
The spot gets narrower
the further out, and later we look at the distribution of particles.

We will show how to implement this example in @sapphire,
and compare the results for different $l_{\text{max}}$.

## Implementation {#implementation-closure}

The implementation of this example in the `examples/closure` directory.
In the following we discuss details,
assuming the reader already familiarised himself with @sapphire.

### VFP equation {#dimension-closure}

The example is one dimensional in space,
and uses mono-energetic particles $\mathbf{p} = \gamma m \mathbf{v}$.
The reduces phase space is therefore only 1 dimensional,
`dim = dim_ps = dim_cs = 1`.

@snippet{lineno} examples/closure/config.h Dimension

In this simplified example, we have the following VFP equation,

$$
  \frac{\mathrm{d} f}{\mathrm{d} t} =
  \frac{\partial f}{\partial t} + V_{x}\frac{\partial f}{\partial x} =
  \frac{\nu}{2} \Delta_{\theta,\varphi} f \,
$$

This results in the following terms we have to activate in @sapphire,
the @ref sapphirepp::VFP::VFPFlags::spatial_advection "spatial advection",
$(\mathbf{u} + \mathbf{v}) \cdot \nabla_{x} f$,
and @ref sapphirepp::VFP::VFPFlags::collision "collisions",
$\frac{\nu}{2} \Delta_{\theta, \varphi} f$.
Furthermore, we can activate the
@ref sapphirepp::VFP::VFPFlags::time_independent_fields "time independent fields".

@snippet{lineno} examples/closure/config.h VFP Flags

### Runtime Parameters {#parameter-closure}

Our setup is characterized by only two parameters,
the standard deviation $\sigma$ of the initial distribution,
and the scattering frequency $\nu$.
As demonstrated in the [parallel shock](#parallel-shock) example,
we will define these as runtime parameters.

1. **Define**

   We start by defining our runtime parameters:

   @snippet{lineno} examples/closure/config.h Define runtime parameter

2. **Declare**
  
   Next, we declare the parameters in the @dealref{ParameterHandler}:

   @snippet{lineno} examples/closure/config.h Declare runtime parameter

3. **Parse**

   Finally, we parse the runtime parameters:

   @snippet{lineno} examples/closure/config.h Parse runtime parameter

### Initial Condition {#initial-condition-closure}

As previously described is the initial condition given by
an isotropic Gaussian distribution,

$$
  f(t=0,x) = \frac{1}{\sqrt{2 \pi} \sigma}
    \exp\left(-\frac{x^{2}}{2\sigma^{2}}\right) \,.
$$

Consequently, we set all expansion coefficients $f_{l>0}$ to zero,
except for the isotropic component $f_{i(0,0,0)} = f_0$.

@snippet{lineno} examples/closure/config.h Initial value

### Scattering frequency {#scattering-frequency-closure}

For the scattering frequency we simply use the constant
given as runtime parameter.

@snippet{lineno} examples/closure/config.h Scattering frequency

## Phase space reconstruction {#phase-space-reconstruction-closure}

To perform the phase-space reconstruction,
we can use the post-processor module provided in
@ref sapphirepp::VFP::PhaseSpaceReconstruction "PhaseSpaceReconstruction".
To activate it,
we just have to give a list of points to the parameter file.
The points are defined in the `Phase space reconstruction` section
as a list of semicolon separated point,
called`reconstruction points`.
Furthermore, we need to give the number of $\theta$ and $\varphi$ values
to use for the reconstruction on the momentum unit-sphere,
$(n_{p_x}, n_{p_y}, n_{p_z})$.

```parameter
subsection Phase space reconstruction
  set reconstruction points = -3; 3
  set n_phi                 = 75
  set n_theta               = 75
end
```

@note For higher dimensional runs,clo
      the reconstruction points are given as a list of points
      in reduced phase space,
      e.g. `x1, y1, p1; x2, y2, p2; x3, y3, p3; ...`

The post-processor will now output additional files in the `results` directory,
named `surface_plot_distribution_function_point_XX_t_XXXX.dat`
and `spherical_density_map_point_XX_t_XXXX.dat`.
The files save the values of $f(\theta, \varphi)$
and $f(n_{p_x}, n_{p_y}, n_{p_z})$
for each point `XX` and time step `XXXX` respectively.

We provide a gnuplot script,
[`examples/closure/heatmap.gp`](https://github.com/sapphirepp/sapphirepp/tree/main/examples/closure/heatmap.gp),
to create heatmaps from these files:

```shell
gnuplot -e "input_file='results/closure/spherical_density_map_point_00_t_0000.dat'; output_file='results/closure/heatmap.png'" heatmap.gp
```

## Results {#results-closure}

@todo Results

### Example parameter file {#example-parameter-closure}

The full parameter file to reproduce the results is given below.

@include{lineno} examples/closure/parameter.prm

<div class="section_buttons">

| Previous              |
|:----------------------|
| [Examples](#examples) |

</div>

---

@author Florian Schulze (<florian.schulze@mpi-hd.mpg.de>)
@date 2024-12-14
