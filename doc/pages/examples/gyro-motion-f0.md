# Gyro motion of an isotropic distribution {#gyro-motion-f0}

@tableofcontents


## Introduction {#introduction-gyro-motion-f0}

This example provides an overview of how magnetic fields function within
@sapphire. We will explore a simple scenario of an isotropic particle
distribution in a constant magnetic field, focusing primarily on the isotropic
part of the distribution function, $f_0$.

Consider a scenario with particles in a constant external magnetic field,

$$
  \mathbf{B} = B_0 \hat{\mathbf{e}}_z \,.
$$

In this example, we will not take into account scattering, a background velocity
field, or a source term. This leads us to the following VFP equation:

$$
  \frac{\partial f}{\partial t} + \mathbf{v} \cdot \nabla_{x} f -
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
Therefore, after $n \times T_g$ with $n \in \mathbb{N}$, the particle
distribution reverts to its initial condition:

$$
  f_{l=0}(\mathbf{r}, t = n \times T_g) = f_{l=0}(\mathbf{r}, t = 0)
  = e^{- \frac{r^2}{2 \sigma^2}} \,.
$$

It's important to note that truncating the expansion of the distribution
function at $l_{\rm max}$ results in deviations. Therefore, $f_{l>0}(\mathbf{r},
t = n \times T_g) \neq 0$ even though the initial condition is isotropic.


## Implementation {#implementation-gyro-motion-f0}

You can find the implementation of this example in the `examples/gyro-motion-f0`
directory. The subsequent sections provide a concise guide on how to implement
this example in @sapphire. This guide assumes that you have already gone through
the [quick start](#quick-start) and [parallel-shock](#parallel-shock) examples,
and are familiar with the basic concepts of @sapphire.


### VFP equation {#dimension-gyro-motion-f0}

As previously discussed, the motion is confined to the $x$-$y$-plane. Therefore,
we use a 2d configuration space, $(x,y)$. Given that we're examining a
discrete momentum $\mathbf{p} = \gamma m \mathbf{v}$, we can equate the reduced
phase space to the configuration space, setting `dim = dim_ps = dim_cs = 2`.

@snippet{lineno} examples/gyro-motion-f0/config.h Dimension

Revisiting the VFP equation,

$$
  \frac{\partial f}{\partial t} + \mathbf{v} \cdot \nabla_{x} f -
  q \mathbf{v} \cdot \left( \mathbf{B} \times \nabla_{p} f \right) = 0 \,,
$$

we have the time derivative and two additional terms. The first term represents
@ref sapphirepp::VFP::VFPFlags::spatial_advection "spatial advection",
$\mathbf{v} \cdot \nabla_{x} f$, while the second term denotes @ref
sapphirepp::VFP::VFPFlags::magnetic "magnetic" interactions, $q \mathbf{v} \cdot
\left( \mathbf{B} \times \nabla_{p} f \right)$. Moreover, the magnetic field is
@ref sapphirepp::VFP::VFPFlags::time_independent_fields "time independent".
Consequently, we activate these three @ref sapphirepp::VFP::VFPFlags
"VFP flags":

@snippet{lineno} examples/gyro-motion-f0/config.h VFP Flags


### Runtime Parameters {#parameter-gyro-motion-f0}

Our setup is characterized by two parameters: the magnetic field strength $B_0$
and the initial spread of the distribution function $\sigma$. The additional
parameters $m$, $q$, and $\gamma$ are already implemented in @sapphire. We will
employ the @ref sapphirepp::PhysicalParameters "PhysicalParameters" class to
define, declare, and parse our custom runtime parameters, as demonstrated in the
[parallel shock](#parallel-shock) example.

1. **Define**

   We start by defining our runtime parameters:

   @snippet{lineno} examples/gyro-motion-f0/config.h Define runtime parameter

2. **Declare**
  
   Note that we use $\frac{B_0}{2 \pi}$ instead of $B_0$ as an input parameter.
   This ensures that `Final time = gamma` corresponds to one full gyroperiod
   $T_g$ (assuming $m = q = c =1$).

   @snippet{lineno} examples/gyro-motion-f0/config.h Declare runtime parameter

3. **Parse**

   Finally, we parse the runtime parameters:

   @snippet{lineno} examples/gyro-motion-f0/config.h Parse runtime parameter


### Initial Condition {#initial-condition-gyro-motion-f0}

The initial condition is an isotropic Gaussian distribution, as previously
described. Consequently, we set all expansion coefficients $f_{l>0}$ to zero,
except for the isotropic component $f_{i(0,0,0)} = f_0$. To calculate $r^2 =
x^2+×y^2$, we employ the
@dealref{point.norm_square(),classTensor,ac9829e8f74544262c6b0ec1f39429029}
function.

@snippet{lineno} examples/gyro-motion-f0/config.h Initial value


### Magnetic field {#magnetic-field-gyro-motion-f0}

The implementation of the magnetic field is straightforward. We assign
$\mathbf{B} = B_0 \hat{\mathbf{e}}_z$ within the @ref
sapphirepp::VFP::MagneticField "MagneticField" class:

@snippet{lineno} examples/gyro-motion-f0/config.h Magnetic field


### Velocity field {#velocity-gyro-motion-f0}

In our VFP equation, we don't have a velocity field $\mathbf{u} (\mathbf{x})$.
However, the @ref sapphirepp::VFP::VFPFlags::spatial_advection
"VFPFlags::spatial_advection" flag activates the full $(\mathbf{u} + \mathbf{v})
\cdot \nabla_x f$ term. Therefore, we need to explicitly set $\mathbf{u} =
\mathbf{0}$. As demonstrated in the [quick start](#quick-start) example, we need
to define the value, divergence, material derivative, and Jacobian of the
velocity field.

1. Background velocity field value $\mathbf{u}(\mathbf{x})$:

   @snippet{lineno} examples/gyro-motion-f0/config.h Background velocity value

2. Background velocity divergence $\nabla \cdot \mathbf{u}(\mathbf{x})$:

   @snippet{lineno} examples/gyro-motion-f0/config.h Background velocity divergence

3. Background velocity material derivative $\frac{\mathrm{D}
   \mathbf{u}}{\mathrm{D} t}$:

   @snippet{lineno} examples/gyro-motion-f0/config.h Background velocity material derivative

4. Background velocity Jacobian $\frac{\partial u_{x}}{\partial x}$:

   @snippet{lineno} examples/gyro-motion-f0/config.h Background velocity Jacobian


## Results {#results-gyro-motion-f0}

We execute the simulation with parameters set to $B_0 = 2\pi$, $\gamma = 2$,
$l_{\rm max} = 3$, and $T_{\rm final} = 20$, which corresponds to ten gyro
periods. The comprehensive parameter file is provided below. The resulting
figures from the simulation are displayed as follows:

![2D time series](https://sapphirepp.org/img/examples/gyro-motion-f0-2d.gif)

The distribution function exhibits a pulsation with a period of $T_g$. If we
reduce the expansion order to $l_{\rm max} = 1$, a deviation from the initial
condition is observed, as anticipated. The subsequent figure illustrates the
results for $l_{\rm max} = 1$.

![2D time series l_max=1](https://sapphirepp.org/img/examples/gyro-motion-f0-2d-lmax1.gif)


### Example parameter file {#example-parameter-gyro-motion-f0}

The parameter file below can be used to run the simulation:

@include{lineno} examples/gyro-motion-f0/parameter.prm


<div class="section_buttons">

| Previous              |
|:----------------------|
| [Examples](#examples) |

</div>


---

@author Florian Schulze (<florian.schulze@mpi-hd.mpg.de>)
@date 2023-12-21
