# Gyro motion with advection {#gyro-advection}

@tableofcontents


## Introduction {#introduction-gyro-advection}

This example provides an overview of how magnetic fields function within
@sapphire. We will explore a simple scenario of an isotropic particle
distribution in a constant magnetic field, advected in a constant background
plasma flow. This example closely follows section 4.2 in @cite Schween2024.

We use a constant magnetic field oriented in the $z$-direction,

$$
  \mathbf{B} = B_0 \hat{\mathbf{e}}_z \,,
$$

and a constant background plasma flow

$$
  \mathbf{u} = u_0 \hat{\mathbf{e}}_x + u_0 \hat{\mathbf{e}}_y \,.
$$

In this example, we will not take into account scattering, a background velocity
field, or a source term. This leads us to the following VFP equation:

$$
  \frac{\partial f}{\partial t} + (\mathbf{u} + \mathbf{v}) \cdot \nabla_{x} f +
  q \mathbf{v} \cdot \left( \mathbf{B} \times \nabla_{p} f \right) = 0 \,.
$$

It's evident that different momenta decouple, so we will examine a discrete
momentum $\mathbf{p} = \gamma m \mathbf{v}$. A single particle in this scenario
will gyrate in the magnetic field with a gyroradius

$$
  r_g = \frac{p_{\perp} c}{q B} = \frac{\gamma m c v_{\perp}}{q B} \,,
$$

and a gyrofrequency

$$
  \omega_g = \frac{q B}{\gamma m c} \,.
$$

Given that the $B$-field is oriented in the $z$-direction, the motion will occur
solely in the $x$-$y$-plane. This allows us to reduce the problem's
dimensionality to 2D. As initial condition, we use a Gaussian isotropic
distribution of particles concentrated at the origin:

\begin{align*}
  f_{l=0}(\mathbf{r}, t=0) &= e^{- \frac{r^2}{2 \sigma^2}} \,, \\
  f_{l>0}(\mathbf{r}, t=0) &= 0 \,.
\end{align*}

Here, $\sigma$ represents the spread of the Gaussian, and $r = \sqrt{x^2 + y^2}$
is the radial distance from the origin. After one gyroperiod,

$$
  T_g = \frac{2 \pi}{\omega_g} = \frac{2 \pi \gamma m c}{q B_0} \,,
$$

each particle completes a full gyro cycle and returns to its initial position.
Therefore, after one gyroperiod the particle distribution reverts to its initial
condition, except being advected by $\mathbf{u} \times T_g$. Using a periodic
box of size $L = n \times u_0 T_g$ with $n \in \mathbb{N}$, the particle
distribution will return to its initial condition after $n$ gyroperiods:

$$
  f_{l=0}(\mathbf{r}, t = n \times T_g) = f_{l=0}(\mathbf{r}, t = 0)
  = e^{- \frac{r^2}{2 \sigma^2}} \,.
$$

It's important to note that truncating the expansion of the distribution
function at $l_{\rm max}$ results in deviations. Therefore, $f_{l>0}(\mathbf{r},
t = n \times T_g) \neq 0$ even though the initial condition is isotropic.


## Implementation {#implementation-gyro-advection}

You can find the implementation of this example in the `examples/gyro-advection`
directory. The subsequent sections provide a concise guide on how to implement
this example in @sapphire. This guide assumes that you have already gone through
the [quick start](#quick-start) and [parallel-shock](#parallel-shock) examples,
and are familiar with the basic concepts of @sapphire.


### VFP equation {#dimension-gyro-advection}

As previously discussed, the motion is confined to the $x$-$y$-plane. Therefore,
we use a 2d configuration space, $(x,y)$. Given that we're examining a
discrete momentum $\mathbf{p} = \gamma m \mathbf{v}$, we can equate the reduced
phase space to the configuration space, setting `dim = dim_ps = dim_cs = 2`.

@snippet{lineno} examples/gyro-advection/config.h Dimension

Revisiting the VFP equation,

$$
  \frac{\partial f}{\partial t} + (\mathbf{u} + \mathbf{v}) \cdot \nabla_{x} f +
  q \mathbf{v} \cdot \left( \mathbf{B} \times \nabla_{p} f \right) = 0 \,,
$$

we have the time derivative and two additional terms. The first term represents
@ref sapphirepp::VFP::VFPFlags::spatial_advection "spatial advection",
$(\mathbf{u} + \mathbf{v}) \cdot \nabla_{x} f$, while the second term denotes
@ref sapphirepp::VFP::VFPFlags::rotation "rotation" due to magnetic
interactions, $q \mathbf{v} \cdot \left( \mathbf{B} \times \nabla_{p} f
\right)$. Moreover, the magnetic field is
@ref sapphirepp::VFP::VFPFlags::time_independent_fields "time independent".
Consequently, we activate these three @ref sapphirepp::VFP::VFPFlags
"VFP flags":

@snippet{lineno} examples/gyro-advection/config.h VFP Flags


### Runtime Parameters {#parameter-gyro-advection}

Our setup is characterized by three parameters: the magnetic field strength
$B_0$, the background plasma flow $u_0$, and the initial spread of the
distribution function $\sigma$. The additional parameters $m$, $q$, and $\gamma$
are already implemented in @sapphire. We will employ the
@ref sapphirepp::PhysicalParameters "PhysicalParameters" class to define,
declare, and parse our custom runtime parameters, as demonstrated in the
[parallel shock](#parallel-shock) example.

1. **Define**

   We start by defining our runtime parameters:

   @snippet{lineno} examples/gyro-advection/config.h Define runtime parameter

2. **Declare**
  
   Note that we use $\frac{B_0}{2 \pi}$ instead of $B_0$ as an input parameter.
   This ensures that $T_{\rm final} = n \times \gamma$ corresponds to exactly
   $n$ full gyroperiods $T_g$ (assuming $m = q = c =1$).

   @snippet{lineno} examples/gyro-advection/config.h Declare runtime parameter

3. **Parse**

   Finally, we parse the runtime parameters:

   @snippet{lineno} examples/gyro-advection/config.h Parse runtime parameter


### Initial Condition {#initial-condition-gyro-advection}

The initial condition is an isotropic Gaussian distribution, as previously
described. Consequently, we set all expansion coefficients $f_{l>0}$ to zero,
except for the isotropic component $f_{i(0,0,0)} = f_0$. To calculate $r^2 =
x^2+y^2$, we employ the
@dealref{point.norm_square(),classTensor,ac9829e8f74544262c6b0ec1f39429029}
function.

@snippet{lineno} examples/gyro-advection/config.h Initial value


### Magnetic field {#magnetic-field-gyro-advection}

The implementation of the magnetic field is straightforward. We assign
$\mathbf{B} = B_0 \hat{\mathbf{e}}_z$ within the @ref
sapphirepp::VFP::MagneticField "MagneticField" class:

@snippet{lineno} examples/gyro-advection/config.h Magnetic field


### Velocity field {#velocity-gyro-advection}

In our VFP equation, we have a constant velocity field $\mathbf{u} (\mathbf{x})
= u_0 \hat{\mathbf{e}}_x + u_0 \hat{\mathbf{e}}_y$. As demonstrated in the
[quick start](#quick-start) example, we not only need to implement the value,
but also the divergence, material derivative, and Jacobian of the velocity
field.

1. Background velocity field value $\mathbf{u}(\mathbf{x}) = u_0
   \hat{\mathbf{e}}_x + u_0 \hat{\mathbf{e}}_y$:

   @snippet{lineno} examples/gyro-advection/config.h Background velocity value

2. Background velocity divergence $\nabla \cdot \mathbf{u}(\mathbf{x}) = 0$:

   @snippet{lineno} examples/gyro-advection/config.h Background velocity divergence

3. Background velocity material derivative $\frac{\mathrm{D}
   \mathbf{u}}{\mathrm{D} t} = 0$:

   @snippet{lineno} examples/gyro-advection/config.h Background velocity material derivative

4. Background velocity Jacobian $\frac{\partial u_{x}}{\partial x} = 0$:

   @snippet{lineno} examples/gyro-advection/config.h Background velocity Jacobian


## Results {#results-gyro-advection}

We execute the simulation with parameters set to $B_0 = 2\pi$, $u_0 = 0.1$,
$\gamma = 2$, $l_{\rm max} = 3$, $L = 6$, and $T_{\rm final} = 200$, which
corresponds to $100$ gyroperiods. The comprehensive parameter file is provided
below. The results are visualized in the following figures, showing the
isotropic part of the distribution function $f_{000}$ at different time steps.

<p float="center">
  <img src="https://sapphirepp.org/img/examples/gyro-advection/advective_t000.png" alt="gyro-advection t=000" width="30%">
  <img src="https://sapphirepp.org/img/examples/gyro-advection/advective_t040.png" alt="gyro-advection t=040" width="30%">
  <img src="https://sapphirepp.org/img/examples/gyro-advection/advective_t080.png" alt="gyro-advection t=080" width="30%">
</p>
<p float="center">
  <img src="https://sapphirepp.org/img/examples/gyro-advection/advective_t120.png" alt="gyro-advection t=120" width="30%">
  <img src="https://sapphirepp.org/img/examples/gyro-advection/advective_t160.png" alt="gyro-advection t=160" width="30%">
  <img src="https://sapphirepp.org/img/examples/gyro-advection/advective_t200.png" alt="gyro-advection t=200" width="30%">
</p>

The distribution function exhibits a pulsation with a period of $T_g$ and is
advected with $\mathbf{0}$. The impact of the cut-off of in $l_{\rm max}$ can be
observed by showing the residual between $t=0$ and $t=200$. Below we show the
results with and without advection, for different $l_{\rm max}$. It can be seen,
that a higher $l_{\rm max}$ leads to a smaller residual.

<p float="center">
  <img src="https://sapphirepp.org/img/examples/gyro-advection/residual_advective_l3.png" alt="residual advective l_max=3" width="30%">
  <img src="https://sapphirepp.org/img/examples/gyro-advection/residual_static_l3.png" alt="residual static l_max=3" width="30%">
  <img src="https://sapphirepp.org/img/examples/gyro-advection/residual_static_l5.png" alt="residual static l_max=5" width="30%">
</p>


### Example parameter file {#example-parameter-gyro-advection}

The parameter file below can be used to run the simulation. Be aware that this
is a high resolution simulation and may take a long time to execute (~200 core
hours). Templates for lower resolution simulations are provided in
`tests/gyro-advection`.

@include{lineno} examples/gyro-advection/parameter.prm


<div class="section_buttons">

| Previous              |
|:----------------------|
| [Examples](#examples) |

</div>


---

@author Florian Schulze (<florian.schulze@mpi-hd.mpg.de>)
@date 2024-03-06
