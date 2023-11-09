# Gyro motion part 1: Isotropic distribution {#gyro-motion-f0}

@tableofcontents


## Introduction {#introduction-gyro-motion-f0}

In this example, we want to give an overview on how magnetic fields work in
@sapphire. Here, we will focus on only the isotropic part of the distribution,
$f_0$. We have a second example focusing on the anisotropies $f_1$:
@todo Link to gyro-motion-f1 example.

Let us consider a set-up, with particles in a constant external magnetic field,

$$
  \mathbf{B} = B_0 \hat{\mathbf{e}}_z \,.
$$

Since we don't consider scattering, background velocity field nor a source term,
we arrive at the following VFP equation:

$$
  \frac{\partial f}{\partial t} + \mathbf{v} \cdot \nabla_{x} f -
  q \mathbf{v} \cdot \left( \mathbf{B} \times \nabla_{p} f \right) = 0 \,.
$$

It can be seen, that different momenta decouple, so we will have a look at a
discrete momentum $\mathbf{p} = \gamma m \mathbf{v}$.  
Considering a single particle, it will gyrate in the magnetic field with a
gyroradius

$$
  r_g = \frac{p_{\perp}}{q B} = \frac{\gamma m v_{\perp}}{q B}
$$

and a gyrofrequency

$$
  \omega_g = \frac{q B}{\gamma m} \,.
$$

With the $B$-field orientated in $z$-direction, the motion will only take place
in the $x$-$y$-plane. This allows us to reduce the dimensionality of the problem
to 2D. As initial condition, we choose an isotropic distribution of particles
concentrated at the origin:

\begin{align*}
  %Gaussian
  f_0(\mathbf{r}, t=0) &= e^{- \frac{r^2}{2 \sigma^2}} \,, \\
  f_{l>0}(\mathbf{r}, t=0) &= 0 \,.
\end{align*}

Here $\sigma$ is the spread of the initial Gaussian, and $r = \sqrt{x^2 + y^2}$
is the radial distance from the origin. After one gyroperiod,

$$
  T_g = \frac{2 \pi}{\omega_g} = \frac{2 \pi q B_0}{\gamma m} \,,
$$

every particle performed a full gyro cycle and returned to the initial position,
undistinguishable from its initial state. Therefore, the distribution function
after $n \times T_g$ with $n \in \mathbb{N}$ is the same as the initial
condition:

$$
  f_l(\mathbf{r}, t = n \times T_g) = f_l(\mathbf{r}, t = 0) \,
$$

i.e.

\begin{align*}
  %Gaussian
  f_0(\mathbf{r}, t = n \times T_g) &= e^{- \frac{r^2}{2 \sigma^2}} \,, \\
  f_{l>0}(\mathbf{r}, t = n \times T_g) &= 0 \,.
\end{align*}

@todo Is it correct, that $f_{l>0}$ returns to $0$?


## Implementation {#implementation-gyro-motion-f0}


### config.h {#config-gyro-motion-f0}


### gyro-motion-f0.cpp {#main-gyro-motion-f0}


## Results {#results-gyro-motion-f0}


### Parameter file {#parameter-gyro-motion-f0}

@include{lineno} examples/gyro-motion-f0/parameter.prm


## The plain program {#plain-gyro-motion-f0}


### config.h {#plain-config-gyro-motion-f0}

@include{lineno} examples/gyro-motion-f0/config.h


### gyro-motion-f0.cpp {#plain-main-gyro-motion-f0}

@include{lineno} examples/gyro-motion-f0/gyro-motion-f0.cpp

---

@author Florian Schulze (<florian.schulze@mpi-hd.mpg.de>)
@date 2023-11-09
