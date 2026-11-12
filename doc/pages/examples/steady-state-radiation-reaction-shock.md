# Steady-state radiation reaction shock {#steady-state-radiation-reaction-shock}

@tableofcontents

In @cite Grebe2026 we compare the results of analytical estimates of the
cut-off momentum of electrons accelerated at astrophysical shocks
and damped by synchrotron radiation,
obtained in the diffusion approximation,
to the results of simulations including higher-order anisotropies.

In this example, we provide the @sapphire setup
that we used to obtain the results presented in @cite Grebe2026.

## Introduction {#introduction-steady-state-radiation-reaction-shock}

The setup closely follows the [oblique shock](#steady-state-oblique-shock) example.

Compared to the scenario presented there,
here we apply Bohm-scaling scattering,
i.e. $\alpha = -1$,
and do not need the parameters `enhanced scattering zone`,
`case identifier`, `zeta one`, `zeta two`, `transition point` and `transition length`.
We thus removed them from the example.

Furthermore, we change the injection from Gaussian injection at the shock to
&ldquo; lower-$p$-boundary &rdquo; injection.
This means that we prescribe an isotropic particle distribution at the lower $p$ boundary
along the total spatial domain.
As the acceleration term in the transport equation is proportional to the
gradient of the fluid velocity,
this only injects particles in the shock zone.
The resulting spectrum is comparable to the Gaussian injection at the shock
used in [oblique shock](#steady-state-oblique-shock) example.

We detail the implementation in the following section.

## Implementation {#implementation-steady-state-radiation-reaction-shock}

The `config.h` and the `parameter.prm` file of this example are stored in the
directory `sapphirepp/examples/vfp/steady-state-radiation-reaction-shock`.

### VFP equation {#vfp-steady-state-radiation-reaction-shock}

We solve the VFP equation including the radiation reaction force
$\mathbf{F}_R$, i.e.

$$
\frac{\partial f}{\partial t} + (\mathbf{U} + \mathbf{v}) \cdot \nabla_x f
-\gamma m \frac{\mathrm{D} \mathbf{U}}{\mathrm{D} t} \cdot \nabla_p f
- (\mathbf{p} \cdot \nabla_x) \mathbf{U} \cdot \nabla_p f
+\nabla_p \cdot (\mathbf{F}_R f)
+ q (\mathbf{v} \times \mathbf{B} ) \cdot \nabla_p f
= \frac{\nu}{2}\,\Delta_{\theta \varphi} f \, .
$$

The list of activated `vfp_flags` is

@snippet{lineno} examples/vfp/steady-state-radiation-reaction-shock/config.h VFP Flags

We include again the
@ref sapphirepp::VFP::VFPFlags::scaled_distribution_function "scaled_distribution_function" flag, note that the
@ref sapphirepp::VFP::VFPFlags::source "source" flag
was removed,
and the new flag
@ref sapphirepp::VFP::VFPFlags::radiation_reaction "radiation_reaction".

### Injection {#injection-steady-state-radiation-reaction-shock}

We prescribe an isotropic particle distribution along the lower $p$ boundary across the total
spatial domain and normalize the distribution with `prm.Q`.
The lower $p$ boundary sets the injection momentum.
Here, the `prm.ln_p_inj` is instead employed to scale the momentum variable,
see below.

The implementation of the `BoundaryValueFunction` function is

@snippet{lineno} examples/vfp/steady-state-radiation-reaction-shock/config.h Boundary value

The variable `boundary_id` is used to distinguish the different boundaries of the computational domain,
cf. "Boundary indicator" and "Colorization" in the [deal.II Glossary](https://dealii.org/current/doxygen/deal.II/DEALGlossary.html).
The id two corresponds to the lower $p$ boundary.

The factor `std::sqrt(4*std::numbers::pi)` accounts for the normalization of the
zeroth order spherical harmonic in the expansion of the distribution function
and the factor $(p/p_{inj})^3$ reflects that we work with a scaled distribution function.

In the parameter file, we set furthermore
`set lower p = inflow`.

### Radiation Reaction Force {#radiation-reaction-force-steady-state-radiation-reaction-shock}

A detailed derivation of the radiation reaction term as implemented
in @sapphire,
as well as validation thereof,
can be found at [radiation reaction](#radiation-reaction)
and @cite Grebe2026.

As we use electrons in this example,
we use a negative charge $-q^*$,
and set the @ref sapphirepp::VFP::ReferenceValues "reference units"
for the simulation to
the electron mass $\underline{m} = 9.1093837 \times 10^{-31} \, \mathrm{kg}$
and elementary charge $\underline{q} = 1.602176634 \times 10^{−19} \, \mathrm{C}$.

```prm
subsection Particle properties
  set Charge = -1
  set Mass   = 1.0
end
subsection Reference units
  set Charge                  = 1.602176634e-19
  set Magnetic field strength = 1.e-10
  set Mass                    = 9.1093837e-31
end
```

We recommend that the spatial grid extends at least 7 diffusion lengths
of the maximum momentum upstream
and 3 diffusion lengths downstream.

### Compile and run {#compile-steady-state-radiation-reaction-shock}

If you configured @sapphire with the `-DEXAMPLES=ON` option,
see [Compilation](#compilation),
the radiation reaction shock example is compiled
and the executable can be found in the
`build/examples/vfp/steady-state-radiation-reaction-shock` folder.
The executable is called `steady-state-radiation-reaction-shock`.

The parameter file is
`examples/vfp/steady-state-radiation-reaction-shock/parameter.prm`.

@include examples/vfp/steady-state-radiation-reaction-shock/parameter.prm

Each run can be started with

```shell
mpirun -n 32 ./build/examples/vfp/steady-state-radiation-reaction-shock/steady-state-radiation-reaction-shock examples/vfp/steady-state-radiation-reaction-shock/parameter.prm
```

@note The size of the simulation is such that it requires the usage of many cores
and a lot of memory,
typically only available on clusters.

The following @sapplot script creates a $2D$-plot of the distribution function,
the scaled, normalized distribution $p^4 f_{000}$ at the shock,
as well as the spatial distribution at fixed momentum value.

It furthermore provides .csv files of the momentum and spatial distribution.

```shell
conda activate sapplot
export SAPPHIREPP_RESULTS=$(pwd)/results
python examples/vfp/steady-state-radiation-reaction-shock/plot_steady_state_radiation_reaction_shock.py
```

## Results {#results-steady-state-radiation-reaction-shock}

The position and steepness of the cut-off is sensitive to the obliquity,
shock velocity, gyrofrequency-averaged scattering rate and the magnetic field strength.
This is illustrated in the following plot.

The fiducial shock parameters are $U/c = 0.05$, $\nu_0 = 0.05$,
$\theta_\mathrm{B} = 30^\circ$ and $B_1 = 10 \, \mu \mathrm{G}$.

The first panel shows how the cut-off moves to higher momenta when comparing
a magnetic field parallel to the shock normal or at
${\theta_\mathrm{B} = 30 ^\circ}$.
The spectral index decreases due to the obliquity,
and we observe the expected pile up near the shock region.

The cut-off also increases for higher shock velocity,
as shown in the second panel, and for increased scattering rate at the chosen
obliquity,
as shown in the third panel.

The last panel shows how an enhanced magnetic field increases radiation damping
and thereby decreases the cut-off momentum.

<div style="text-align:center;">
<img alt="Dependence of particle distribution at the shock with and without radiation damping on the shock parameters."
src="https://sapphirepp.org/img/examples/steady-state-radiation-reaction-shock/shockparameter_influence.png">
</div>

<div class="section_buttons">

| Previous              |
|:----------------------|
| [Examples](#examples) |

</div>

---

@author Marie-Ann Grebe (<magrebe@mpi-hd.mpg.de>)
@author Nils Schween (<nils.schween@mpi-hd.mpg.de>)
@date 2026-11-11
