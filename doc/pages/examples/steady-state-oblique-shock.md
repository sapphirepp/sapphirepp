# Steady-state oblique shock {#steady-state-oblique-shock}

@tableofcontents

In @cite Shirin2025 we suggest that curved or broken emission spectra of X- and gamma ray sources
could be explained with a combination of a magnetic field, 
which is inclined with respect to the shock normal, 
and a $p$-dependence of the scattering frequency, 
that is in agreement with Iroshnikov–Kraichnan MHD turbulence,
or an enhanced scattering zone in the precursor of the shock respectively.

In this example, we provide the @sapphire setup 
that we used to obtain the results presented in @cite Shirin2025 . 
The example focuses on the relevant aspects only. 
For a deeper understanding, a read of the [parallel shock](#parallel-shock) example
and its [steady-state version](#steady-state-parallel-shock) may be useful.

## Introduction {#introduction-steady-state-oblique-shock}

We study the acceleration of energetic particles at a collisionless shock with
compression ratio $r$ that propagates with velocity $\mathbf{U}$ through an astrophysical plasma,
e.g. the interstellar medium. 
The particles are injected at a constant rate $Q$ at the position of the shock wave.
The injected particle distribution is mono-energetic and isotropic.
The astrophysical plasma is permeated with a mean magnetic field $\mathbf{B}$.
Furthermore it is turbulent, 
which causes the energetic particles to frequently change their direction of motion;
a process that we model with the scattering frequency $\nu$.

We summarise the physical scenario in the following table:
<div style="text-align:center">
| Physical scenario:   |                                                                                                                                                                                                               |
|----------------------|---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| Velocity profile     | $$\mathbf{U} = \begin{cases}U_{\text{sh}} \mathbf{e}_x  &\text{ for } x < 0\\  \frac{U_{\text{sh}}}{r} \mathbf{e}_x & \text{ for } x > 0\end{cases}$$                                                         |
| Magnetic field       | $$\mathbf{B}(x) = B_0 \begin{cases} \cos\theta_n \mathbf{e}_x + \sin\theta_n \mathbf{e}_z & \text{ for } x < 0 \\  \cos\theta_n \mathbf{e}_x + r \sin\theta_n \mathbf{e}_z & \text{ for } x > 0 \end{cases}$$ |
| Scattering frequency | $$\nu(x) = \nu_0 B p^{\alpha} $$                                                                                                                                                                              |
| Source               | $$ S(x,p) = \frac{Q}{4\pi p^2} \frac{1}{2\pi \sigma_x \sigma_p}\exp\left[-\left(\frac{(x - x_{\text{inj}})^2}{2 \sigma^2_x} + \frac{(p -p_\text{inj})^2}{2\sigma^2_p}\right)\right]$$                         |
</div>

The shock velocity is $U_{\text{sh}}$. 
The angle $\theta_n = \mathbf{n} \cdot \mathbf{B}/ B$ is the angle between 
the magnetic field $\mathbf{B}$ and the normal $\mathbf{n}$ of the shock surface.
The scattering frequency is normalised with $\nu_{0}$
and its $p$-dependence is characterised by the parameter $\alpha$.

## Implementation {#implementation-steady-state-oblique-shock}

The `config.h` and the `parameter.prm` file of this example are stored in the
directory `sapphirepp/examples/vfp/steady-state-oblique-shock`
and can, for example, be used to explore specific X- or gamma ray sources.

Considering that we are interested in the effects of an oblique magnetic field 
in combination with a $p$- and $x$- dependent scattering frequency, 
we restrict our discussion of the implementation to these functions. 

### VFP equation {#vfp-steady-state-oblique-shock}

Though, before we look at the source code of the magnetic field and scattering frequency function,
we state the VFP equation that we solve, i.e. 

$$
  (\mathbf{U} + \mathbf{v}) \cdot \nabla_{x} f -
  \gamma m \frac{\mathrm{D} \mathbf{U}}{\mathrm{D} t} \cdot \nabla_{p}f -
  \mathbf{p} \cdot\nabla_{x} \mathbf{U}\cdot \nabla_{p} f +
  q \mathbf{v} \cdot \left( \mathbf{B} \times \nabla_{p} f \right) =
  \frac{\nu}{2} \Delta_{\theta, \varphi} f + S \,.
$$

Therefore, the list of `vfp_flags` is

@snippet{lineno} examples/vfp/steady-state-oblique-shock/config.h VFP Flags

Note, we included the
@ref sapphirepp::VFP::VFPFlags::scaled_distribution_function "scaled_distribution_function" flag.
This means that we compute $g = p^{3} f$.
Moreover, we use $\ln p$ instead of $p$, i.e. $g = g(x, \ln p)$. 

### Magnetic field {#magnetic-field-steady-state-oblique-shock}

We orient the magnetic field in the $x$-$z$ plane. The shock normal points into the $x$-direction.  
For a strong shock, namely a shock whose Alfvén Mach number ($M_{A}$) is much greater than one,
the $\mathbf{B}$-field component tangential to the shock surface $B_{z}$,
is enhanced by the compression ratio $r = 4$.

Note, the transition from its up- to its downstream value is _not_ discontinuous:
The Maxwell-Faraday equation implies for a time-independent $\mathbf{B}$-field that
$$
 0 = \nabla \times \mathbf{E} = \begin{pmatrix}0 \\ 0 \\ - \partial_{x} (U_x B_z)\end{pmatrix}\,,
$$
where we used the ideal MHD approximation $\mathbf{E} = - \mathbf{U} \times \mathbf{B}$ 
and the assumption that $U_x$ and $B_z$ do not change along $z$-direction.
Additionally, it is necessary to continuously vary the velocity $\mathbf{U}$ through the shock;
the energy change of the particles due to the change of reference frames on shock crossings
is reflected in the VFP equation in terms containing the derivative of $\mathbf{U}$, 
which is not defined for a discontinuous velocity profile. We chose
$$
\mathbf{U} = \frac{U_{\text{sh}}}{2r}\left[ (r+1) - (r-1) \tanh\left(\frac{x}{L_s}\right) \right]\mathbf{e}_x \,.
$$

Because $U_x B_z = \text{const.}$, a possible choice for $B_z$ is
$$
B_z = B_{0} \sin\theta_{n} 2r \left[ (r+1) - (r-1) \tanh\left(\frac{x}{L_s}\right) \right]^{-1} \, .
$$

The implementation of the `MagneticField` function looks like 

@snippet{lineno} examples/vfp/steady-state-oblique-shock/config.h Magnetic field

The parameter `obliqueness` is the angle between the magnetic field and the shock normal, i.e.  $\theta_n$ in degrees, 
the `compression ratio` is $r$ and the `shock_width` is $L_s$. 
The parameters can be changed in the `parameter.prm` file. 

### Scattering frequency {#scattering-frequency-steady-state-oblique-shock}

It is the interplay between the oblique magnetic field 
and the $p$ and $x$ dependence of the scattering frequency $\nu$
that produces a curved or broken spectrum. 

We start with a $p$-dependent scattering frequency, namely
$$
  \nu(p) = \nu_0 B p^{\alpha} \, ,
$$

where $\nu_0$ is a parameter which can be freely set and $B$ is the magnitude of the magnetic field.
The case $\alpha = -1$ corresponds to a scattering frequency 
that is proportional to the gyro frequency of the particles, i.e. $\nu \propto \omega_g$. 
This scaling of the scattering frequency with momentum is called "Bohm scaling".
$\nu = \omega_g$ is the Bohm limit and corresponds to $\nu_0 = 1$. 
We note that a plasma in which particles are scattered more often is not magnetised anymore.
The Bohm scaling case in the context of oblique shock was, for example, studied in @cite Bell2011 .
No deviations from a power-law spectrum are expected.

In this example, we choose $\alpha$ to equal $-1/2$. 
This corresponds to Iroshnikov–Kraichnan MHD turbulence of the background astrophysical plasma. 

In a second step, we investigate a $x$- and $p$-dependent scattering frequency. 
We choose the $x$-dependence such that it models an enhanced scattering zone in the precursor of the shock wave.
The idea is to mimic the self-excitation of magnetic field perturbations ahead of the shock. 
How these perturbations change when being advected through the shock is not trivial, 
we therefore experiment with two cases:

__Case (a)__:  
The scattering frequency is 

$$
	\nu(x, p) = \nu_0 B \left[\frac{\zeta_1 + \zeta_2}{2} - \frac{\zeta_1 - \zeta_2}{2} 
	\tanh\left(\frac{x + x_{\text{pre}}}{L_T}\right)\right] p^{-1}
$$

This choice implies that $\eta = \nu/\omega_g$ is constant across the shock.
Note that we returned to $\alpha = -1$ to isolate the effects of the enhanced scattering zone.
The spatial profile of the scattering frequency and of $\eta$
are plotted in the vicinity of the shock (grey dashed line).

<div style="text-align:center;">
<img alt="Case A: The ratio of the scattering frequency to the gyro frequency, called eta, is constant across the shock."
src="https://sapphirepp.org/img/examples/steady-state-oblique-shock/scattering_regime_case_a.png">
</div>


__Case (b)__:  

The scattering frequency is 
$$
	\nu(x, p) = \nu_0 B_0 \left[\frac{\zeta_1 + \zeta_2}{2} - \frac{\zeta_1 - \zeta_2}{2} 
	\tanh\left(\frac{x + x_{\text{pre}}}{L_T}\right)\right] p^{-1}
$$

In contrast to case (a), it is the scattering frequency $\nu$ that does not jump at the shock. 

<div style="text-align:center;">
<img alt="Case B: The scattering frequency is constant across the shock."
src="https://sapphirepp.org/img/examples/steady-state-oblique-shock/scattering_regime_case_b.png">
</div>

The implementation looks like

@snippet{lineno} examples/vfp/steady-state-oblique-shock/config.h Scattering frequency

Note, we use $\ln p$ instead of $p$ and, thus, $\nu(\ln p) = \nu_0 B \exp(\alpha \ln p)$.
Furthermore, `points[q_index][1]` is $\ln p$.

We introduced a set of parameters, namely `enhanced scattering zone`, `nu0`, `alpha`, 
`case identifier`, `zeta one`, `zeta two`, `transition point` and `transition length`.

If the `enhanced scattering zone` parameter is set to `true`,
the scattering frequency depends on $x$ and $p$.
The `case identifier` distinguishes between case (a) and case (b).
The parameters `zeta one` and `zeta two` determine the final values of
the scattering frequency (or $\eta$) in the up- and downstream respectively.
The `transition point` determines where the transition between $\zeta_1$ and $zeta_2$ happens.
It could be interpreted as the size of the enhanced scattering zone.
The `transition length` determines the width of the used tanh-profile.
If the `enhanced scattering zone` parameter is set to `false`,
the scattering frequency only depends on $p$.


### Compile and run {#compile-steady-state-oblique-shock}

If you configured @sapphire with the `-DEXAMPLES=ON` option,
see [Compilation](#compilation),
the oblique shock example is compiled
and the executable can be found in the
`build/examples/vfp/steady-state-oblique-shock` folder.
The executable is called `steady-state-oblique-shock`

We supply the user with a physically reasonable set of parameters that are
listed in the `examples/vfp/steady-state-oblique-shock/parameter.prm` file:

@include examples/vfp/steady-state-oblique-shock/parameter.prm

Run the simulation with:

```shell
mpirun -n 6 ./build/examples/vfp/steady-state-oblique-shock/steady-state-oblique-shock examples/vfp/steady-state-oblique-shock/parameter.prm
```

The plots can be created with the command:

```shell
pvbatch examples/vfp/steady-state-parallel-shock/pvplot.py results/steady-state-parallel-shock
```

## Results {#results-steady-state-oblique-shock}


<div class="section_buttons">

| Previous              |
|:----------------------|
| [Examples](#examples) |

</div>

---

@author Nils Schween (<nils.schween@mpi-hd.mpg.de>)
@date 2025-02-05
