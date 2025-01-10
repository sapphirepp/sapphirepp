# Closure {#closure}

@tableofcontents

## Introduction {#introduction-closure}

This example is designed to give a feeling
for the effect of truncating the expansion at $l_{\rm max}$.
Furthermore, we will introduce the post-processing module
to reconstruct the full distribution function
$f(\mathbf{x}, p, \theta, \varphi)$.
To this end, we follow the example introduced in @cite Schween2025, Sec. 4.3.

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

![Analytical velocity distribution](https://sapphirepp.org/img/examples/closure/analytical_heatmaps.png)
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

The implementation of this example can be found in the `examples/closure` directory.
In the following we discuss the details,
assuming that the reader already familiarised himself with @sapphire.

### VFP equation {#dimension-closure}

The example is one dimensional in space,
and uses mono-energetic particles  with momentum $\mathbf{p} = \gamma m \mathbf{v}$.
The reduced phase space is therefore only 1 dimensional,
`dim = dim_ps = dim_cs = 1`.

@snippet{lineno} examples/closure/config.h Dimension

As stated above, in this example, we use the subsequent VFP equation,

$$
  \frac{\mathrm{d} f}{\mathrm{d} t} =
  \frac{\partial f}{\partial t} + V_{x}\frac{\partial f}{\partial x} =
  \frac{\nu}{2} \Delta_{\theta,\varphi} f \,
$$

This results in the following terms that we have to activate in @sapphire,
the @ref sapphirepp::VFP::VFPFlags::spatial_advection "spatial advection",
$(\mathbf{u} + \mathbf{v}) \cdot \nabla_{x} f$,
and @ref sapphirepp::VFP::VFPFlags::collision "collisions",
$\frac{\nu}{2} \Delta_{\theta, \varphi} f$.
Furthermore, to decrease the computation time we include 
@ref sapphirepp::VFP::VFPFlags::time_independent_fields "time independent fields" flag.

@snippet{lineno} examples/closure/config.h VFP Flags

### Runtime Parameters {#parameter-closure}

Our setup is characterized by only two parameters,
the standard deviation $\sigma$ of the initial distribution,
and the scattering frequency $\nu$.
As demonstrated in the [parallel shock](#parallel-shock) example,
we will these as runtime parameters.

1. **Define**

   We start by declaring our runtime parameters:

   @snippet{lineno} examples/closure/config.h Define runtime parameter

2. **Declare**
  
   Next, we inform the @dealref{ParameterHandler} that the declared parameters exist:

   @snippet{lineno} examples/closure/config.h Declare runtime parameter

3. **Parse**

   Finally, we parse them:

   @snippet{lineno} examples/closure/config.h Parse runtime parameter

### Initial Condition {#initial-condition-closure}

As previously stated, the initial condition is given by
an isotropic Gaussian distribution,

$$
  f(t=0,x) = \frac{1}{\sqrt{2 \pi} \sigma}
    \exp\left(-\frac{x^{2}}{2\sigma^{2}}\right) \,.
$$

Consequently, we set all expansion coefficients $f_{l>0}$ to zero,
except for the isotropic component $f_{i(0,0,0)} = \sqrt{4\pi} f(t=0, x)$.

@snippet{lineno} examples/closure/config.h Initial value

### Scattering frequency {#scattering-frequency-closure}

For the scattering frequency we simply use the constant
given as a runtime parameter.

@snippet{lineno} examples/closure/config.h Scattering frequency

## Phase space reconstruction {#phase-space-reconstruction-closure}

To perform the phase-space reconstruction at speficic points in reduced phase
space $\xi = (x, y, p)^{T}$, 
we can use the post-processor module provided in
@ref sapphirepp::VFP::PhaseSpaceReconstruction "PhaseSpaceReconstruction".
To activate it,
we just have to give a list of points to the parameter file.
The points are defined in the `Phase space reconstruction` section
as a semicolon separated list,
called`reconstruction points`.
Furthermore, we need to set the number of $\theta$ and $\varphi$ values
used for the reconstruction on the momentum unit-sphere,
$(n_{p_x}, n_{p_y}, n_{p_z})$.

```parameter
subsection Phase space reconstruction
  set reconstruction points = -6; -3; 0; 3; 6
  set n_phi                 = 75
  set n_theta               = 75
end
```

@note For higher dimensional runs,
      the reconstruction points are given as a list of the following form,
      e.g. `x1, y1, p1; x2, y2, p2; x3, y3, p3; ...`

The post-processor will now output additional files in the `results` directory,
named `surface_plot_distribution_function_point_XX_t_XXXX.dat`
and `spherical_density_map_point_XX_t_XXXX.dat`.
The files save the values of $f(t,\xi,\theta, \varphi)$, where $\xi$ is one of
the points in the list given in the parameter file. Moreover,
$f(t,\xi,n_{p_x},n_{p_y}, n_{p_z})$ 
is also stored for each point `XX` and time step `XXXX` respectively.

We provide a gnuplot script,
[`examples/closure/heatmap.gp`](https://github.com/sapphirepp/sapphirepp/tree/main/examples/closure/heatmap.gp),
to create heatmaps from these files:

```shell
gnuplot -e "input_file='results/closure/spherical_density_map_point_00_t_0000.dat'; output_file='results/closure/heatmap.png'" heatmap.gp
```

## Results {#results-closure}

As we want to compare the solution from @sapphire
to the analytic solution given above at $(x=0,t=0)$,
$(x=3 \sigma, t = 3 \sigma/ V = 6/\sqrt{3})$
and $(x=6 \sigma, t= 6 \sigma /V = 12/\sqrt{3})$,
we choose a time-step
$\Delta t = \frac{1}{\sqrt{3} \cdot 100}$.
The corresponding time steps are then $N=0$, $N=600$ and $N=1200$ respectively.

In a first run, we use $l_{\text{max}} = 3$:

```parameter
subsection Expansion
  set Expansion order = 3
end
```

After compiling and running the example,
we can visualize the results at the predefined points
using the provided gnuplot script, i.e. we execute

```shell
gnuplot -e "input_file='results/closure/spherical_density_map_point_02_t_0000.dat'; output_file='results/closure/heatmap_l3_x0_t0.png'" heatmap.gp
gnuplot -e "input_file='results/closure/spherical_density_map_point_03_t_0600.dat'; output_file='results/closure/heatmap_l3_x3_t3.png'" heatmap.gp
gnuplot -e "input_file='results/closure/spherical_density_map_point_04_t_1200.dat'; output_file='results/closure/heatmap_l3_x6_t6.png'" heatmap.gp
```

@note The index for the point follows the same order
      as the definition in the parameter file.
      The used position and time can be checked
      in the header of the output files.

<p float="center">
  <img src="https://sapphirepp.org/img/examples/closure/heatmap_l3_x0_t0.png" alt="Numerical solution for l_max=3, x=0, t=0" width="30%">
  <img src="https://sapphirepp.org/img/examples/closure/heatmap_l3_x3_t3.png" alt="Numerical solution for l_max=3, x=3, t=3/V" width="30%">
  <img src="https://sapphirepp.org/img/examples/closure/heatmap_l3_x6_t6.png" alt="Numerical solution for l_max=3, x=6, t=6/V" width="30%">
</p>

Taking a closer look at the results, in particular the colorbar,
we see that the numerical solution does not match
the analytical solution presented above.
The distribution is too smooth and,
especially for the point $(x=6 \sigma, t=6/\sqrt{3})$,
it underestimates the peak at the pole of the sphere. Furthermore,
the numerically distribution function is negative.

We can improve the results by increasing the expansion order
to $l_{\text{max}} = 4$:

```parameter
subsection Expansion
  set Expansion order = 4
end
```

```shell
gnuplot -e "input_file='results/closure/spherical_density_map_point_02_t_0000.dat'; output_file='results/closure/heatmap_l4_x0_t0.png'" heatmap.gp
gnuplot -e "input_file='results/closure/spherical_density_map_point_03_t_0600.dat'; output_file='results/closure/heatmap_l4_x3_t3.png'" heatmap.gp
gnuplot -e "input_file='results/closure/spherical_density_map_point_04_t_1200.dat'; output_file='results/closure/heatmap_l4_x6_t6.png'" heatmap.gp
```

Repeating the steps above,
running the simulation and visualizing the results,
we obtain the following heatmaps.

<p float="center">
  <img src="https://sapphirepp.org/img/examples/closure/heatmap_l4_x0_t0.png" alt="Numerical solution for l_max=4, x=0, t=0" width="30%">
  <img src="https://sapphirepp.org/img/examples/closure/heatmap_l4_x3_t3.png" alt="Numerical solution for l_max=4, x=3, t=3/V" width="30%">
  <img src="https://sapphirepp.org/img/examples/closure/heatmap_l4_x6_t6.png" alt="Numerical solution for l_max=4, x=6, t=6/V" width="30%">
</p>

We see a slight improvement.

For a more direct comparison of the effects of different choices for
$l_{\text{max}}$, we plot the heatmaps at a single point, namely at $(x=6 \sigma,
t=6/\sqrt{3})$, for the different $l_{\text{max}}$. The left plot shows the
the particle distribution with $l_{\text{max}} = 3$, the plot in the middle
shows the $l_{\text{max}} = 4$ case and for the right heatmap, we ran a simulation with
$l_{\text{max}} = 11$. 

<p float="center">
  <img src="https://sapphirepp.org/img/examples/closure/heatmap_l3_x6_t6.png" alt="Numerical solution for l_max=3, x=6, t=6/V" width="30%">
  <img src="https://sapphirepp.org/img/examples/closure/heatmap_l4_x6_t6.png" alt="Numerical solution for l_max=4, x=6, t=6/V" width="30%">
  <img src="https://sapphirepp.org/img/examples/closure/heatmap_l11_x6_t6.png" alt="Numerical solution for l_max=11, x=6, t=6/V" width="30%">
</p>

Finally, we want to emphasize that using a scattering frequency $\nu > 0$,
suppresses the higher order terms.
Therefore, a lower $l_{\text{max}}$ is sufficient to obtain accurate results.

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
@author Nils Schween (<nils.schween@mpi-hd.mpg.de>)
@date 2024-12-14
