# Steady-state acceleration with synchrotron damping {#steady-state-radiation-reaction-shock}

@tableofcontents

In @cite Grebe2026 we compare the results of analytical estimates in the 
diffusion approximation of the cut-off momentum of electrons accelerated at 
astrophysical shocks and damped by synchrotron radiation,
 to the results of simulations including higher-order ansiotropies. 

In this example, we provide the @sapphire setup
that we used to obtain the results presented in @cite Grebe2026 .
The example focuses on the relevant aspects only.

## Introduction {#introduction-steady-state-radiation-reaction-shock}

The setup closely follows the [oblique shock](#oblique-shock) example. 

Compared to the physical scenario presented there, 
we apply Bohm-scaling scattering, 
i.e. $\alpha = 1$, 
and do not require the parameters `enhanced scattering zone`, `alpha`, 
`case identifier`, `zeta one`, `zeta two`, `transition point` and `transition length`.

Furthermore, we change the injection from Gaussian injection at the shock to 
lower $p$ boundary injection. 
This means that we inject the istropic distribution at the lower $p$ boundary 
along the total spatial domain. 
As the acceleration term in the transport equation is proportional to the 
gradient of the fluid velocity, 
this produces the same spectrum at the shock as Gaussian injection would. 

We detail the implementation in the following section.


## Implementation {#implementation-steady-state-radiation-reaction-shock}

The `config.h` and the `parameter.prm` file of this example are stored in the
directory `sapphirepp/examples/vfp/steady-state-radiation-reaction-shock`.

### VFP equation {#vfp-steady-state-radiation-reaction-shock}

We solve the VFP equation including the radiation reaction force 
$\mathbf{F}_R$, i.e. 

\begin{multline}
\frac{\partial f}{\partial t} + (\mathbf{U} + \mathbf{v}) \cdot \nabla_x f
-\gamma m \frac{\mathrm{D} \mathbf{U}}{\mathrm{D} t} \cdot \nabla_p f
- (\mathbf{p} \cdot \nabla_x) \mathbf{U} \cdot \nabla_p f \\
+\nabla_p \cdot (\mathbf{F}_R f)
+ q (\mathbf{v} \times \mathbf{B} ) \cdot \nabla_p f
= \frac{\nu}{2}\,\Delta_{\theta \varphi} f + S \, .
\end{multline}

The list of activated `vfp_flags` is 

@snippet{lineno} examples/vfp/steady-state-radiation-reaction-shock/config.h VFP Flags

We include again the 
@ref sapphirepp::VFP::VFPFlags::scaled_distribution_function "scaled_distribution_function" flag.
and work in logarithmic momenta.

Note that the 
@ref sapphirepp::VFP::VFPFlags::source "source" flag
is disabled.

### Injection {#injection-steady-state-radiation-reaction-shock}

We choose to inject particles along the lower $p$ boundary across the total 
spatial domain. 

The implementation of the `BoundaryValueFunction` function including lower $p$ 
injection looks like

@snippet{lineno} examples/vfp/steady-state-radiation-reaction-shock/config.h Boundary Value Function

We normalize the distribution to the lower $p$ boundary,
which we set equal to the injection momentum `prm.ln_p_inj` specified in the 
parameter file, see below. 

The factor `std::sqrt(4*std::numbers::pi)` accounts for the normalization. 

In the parameter file, we set for the lower $p$ boundary 
`set lower p = inflow`.

### Radiation Reaction Force {#radiation-reaction-force-steady-state-radiation-reaction-shock}

The details of the derivation of the radiation reaction term as implemented 
in @sapphire
as well as validation thereof can be found at [radiation reaction](#radiation-reaction)
and @cite Grebe2026.

As we use dimensionless units for the radiation reaction term, 
we introduce a new subsection to the parameter file, 

```  subsection Particle properties
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
`examples/vfp/steady-state-radiation-reaction-shock/parameter.prm` .

@include examples/vfp/steady-state-radiation-reaction-shock/parameter.prm

Each run can be started with

```shell
mpirun -n 32 ./build/examples/vfp/steady-state-radiation-reaction-shock/steady-state-radiation-reaction-shock examples/vfp/steady-state-radiation-reaction-shock/parameter.prm
```

@note The size of the simulation is such that it requires the usage of many cores
and a lot of memory,
typically only available on clusters.

The plots can be created with a @sapplot script:

```shell
conda activate sapplot
export SAPPHIREPP_RESULTS=$(pwd)/results
python examples/vfp/steady-state-radiation-reaction-shock/plot_steady_state_radiation_reaction_shock.py
```

## Results {#results-steady-state-radiation-reaction-shock}

The position and steepness of the cut-off is sensitive to the obliquity,
 shock velocity, gyrofrequency-averaged scattering rate and the magnetic field 
 strength. This is illustrated in the following plot. 

The fiducial shock parameters are $U/c = 0.05$, $\nu_0 = 0.05$, 
$\theta_\mathrm{B} = 30^\circ$ and $B_1 = 10 \, \mu \mathrm{G}$.

The first panel shows how the cut-off moves to higher momenta when the magentic
field has an angle of $\theta_\mathrm{B} = 30 ^\circ$ to the shock normal. 
The spectral index decreases as well, 
and we observe the expected pile up near the shock region. 

The cut-off momentum also increases for higher shock velocity, 
as shown in the second panel, and for increased scattering rate at the chosen
obliquity, as shown in the third panel.

The last panel shows how radiation damping decreases the cut-off momentum 
when the magnetic field is enhanced. 


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

@author Nils Schween (<nils.schween@mpi-hd.mpg.de>)
@author Marie-Ann Grebe (<magrebe@mpi-hd.mpg.de>)
@date 2025-10-19
