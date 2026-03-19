# Radiation Reaction {#radiation-reaction}

@tableofcontents

@sapphire simulates the evolution of charged particles in a background plasma. 
The particles are modelled with a phase-space density called the distribution function, 
which is denoted with $f(t, \mathbf{x}, \mathbf{p})$.
Its evolution is computed by solving a Vlasov-Fokker-Planck (VFP) equation.

This page documents what term we add to the VFP equation to include the **radiation reaction force**.
An accelerated charged particle emits electromagnetic radiation and, as a consequence,
experiences a reaction force due to the loss of energy and momentum carried away by the emitted radiation.
This force is called the radiation reaction force.
We use the Landau-Lifshitz (LL) formulation of the Lorentz-Abraham-Dirac radiation–reaction force.

Additionally, we write down how the radiation reaction force modifies the discontinuous Galerkin (dG) weak formulation used in @sapphire
and we validate its implementation against analytical solutions.

## From the LL radiation reaction force to the VFP equation

The back-reaction on the particles because of the emission of photons
is described by the radiation reaction equation as discussed by Landau and Lifshitz
@cite landau1975classical.
It supplements the Lorentz force law with an additional damping force, namely

$$ 
\mathbf{F}_R = \frac{\mu_0q^{3}\gamma}{6\pi mc} 
 \Bigg\{ \left( \frac{\partial}{\partial t} + \mathbf{v} \cdot \nabla \right)(\mathbf{E}+\mathbf{v} \times \mathbf{B}) \Bigg\} + \frac{\mu_0q^{4}}{6\pi m^{2}c} \Bigg\{\frac{1}{c^2} \mathbf{E} (\mathbf{v} \cdot \mathbf{E}) - \mathbf{B} \times (\mathbf{E} +\mathbf{v} \times \mathbf{B}) \Bigg\} - \frac{\mu_0q^{4} \gamma^2}{6\pi m^{2}c^{3}} \Bigg\{\left( \mathbf{E} + \mathbf{v} \times \mathbf{B} \right)^{2} - \frac{1}{c^{2}} (\mathbf{E} \cdot \mathbf{v})^{2} \Bigg\} \mathbf{v} .
$$

Here, $\mu_0$ denotes the magnetic vacuum permeability, 
$\mathbf{E}$ and $\mathbf{B}$ are the electromagnetic fields of the background plasma, 
$q$ is the charge of the particle and $\mathbf{v}$ its velocity.
Moreover, $\gamma = (1 - \beta^{2})^{-1/2}$ is its Lorentz factor 
and $\beta = \| \boldsymbol{\beta} \| = \| \mathbf{v}\|/c$.

Following @cite tamburini2010radiation , the first curly-bracketed term,
containing derivatives of the electromagnetic fields, is dropped,
whereas the second term of the radiation reaction force $\mathbf{F}_R$ is kept.
We note that for ultra-relativistic electrons, 
namely electrons with Lorentz factors $\gamma \gg 1$, 
it is the third term that captures the dominant contribution the effects from the radiation reaction, in particular the energy loss.

Before including the radiation reaction force in the VFP equation, we further simplify it.
In @sapphire we use mixed-coordinates, 
i.e we use a different coordinate system for the configuration and the momentum space. 
For example, our $\mathbf{x}$-coordinates may refer to the rest frame of a shock wave that propagates through the background plasma 
whereas the $\mathbf{p}$-coordinates always refer to particle momenta as given in the rest frame of the background plasma. 
Furthermore, the background plasma is modelled using the **ideal magnetohydrodynamic (MHD) approximation**, 
which assumes infinite plasma conductivity. The approximation implies that 

$$
\mathbf{E} = - \mathbf{U} \times \mathbf{B} \, ,
$$

where $\mathbf{U}$ is the plasma velocity. 

If we use an apostrophe to denote quantities in the rest frame of the plasma 
(the fluid rest frame as we sometimes say), then $\mathbf{E}' = 0$ 
and the radiation reaction force reduces to 
(keeping in mind that we dropped the term in $\mathbf{F}_R$
that contains the derivatives of the electromagnetic fields)

$$
\mathbf{F}'_R = \sigma \Big( \mathbf{B}' (\mathbf{B}' \cdot \boldsymbol{\beta}') + \gamma'^2 \boldsymbol{\beta}' ((\boldsymbol{\beta}'\cdot \mathbf{B}')^2 - B'^2) \Big),
$$

with the coefficient

$$
\sigma =  \frac{\mu_0 q^4}{6\pi m^2} .
$$

Note that the first term in the expression for $\mathbf{F}'_R$ contains the projection of the velocity along the magnetic field, 
i.e. it gives the rate at which the momentum component parallel to the $\mathbf{B}$-field changes, 
while the second term accounts for the change of perpendicular momentum component
responsible for energy losses.

We, finally, add the radiation reaction force $\mathbf{F}'_R$ to the VFP equation used in @sapphire, i.e. the modified VFP equation is

$$
\frac{\partial f}{\partial t} + (\mathbf{U} + \mathbf{v}) \cdot \nabla_x f
-\gamma m \frac{\mathrm{D} \mathbf{U}}{\mathrm{D} t} \cdot \nabla_p f
- (\mathbf{p} \cdot \nabla_x) \mathbf{U} \cdot \nabla_p f
+\nabla_p \cdot (\mathbf{F}_R f)
+ q (\mathbf{v} \times \mathbf{B} ) \cdot \nabla_p f
= \frac{\nu}{2}\,\Delta_{\theta \varphi} f + S \, .
$$

We emphasise that we dropped the apostrophe again; in @sapphire it is implicit that all quantities related to momentum are given in the fluid rest frame. For a more detailed discussion of the above equation see @cite Schween2024a , in particular Sec. 2. 

## Spherical harmonic representation of the LL radiation reaction force

In @sapphire we use a truncated (real) spherical harmonic expansion of the distribution function $f$
, i.e.

$$
f(t,\mathbf{x}, p, \theta, \varphi) = \sum^{l_{\mathrm{max}}}_{l = 0} \sum^{1}_{s = 0} \sum^{l}_{m = s} f_{lms}(t, \mathbf{x}, p) Y_{lms}(\theta, \varphi) \,.
$$

Note that we use a spherical coordinate system for the momentum space coordinates.
Our polar direction is the $x$-direction.

This implies that it is necessary to derive partial differential equations (PDEs) for the expansion coefficients $f_{lms}$. This can be done with the operator-based method introduced in @cite Schween2024a . 
The method boils to down to replacing the operators appearing in the VFP equation with their matrix representations in the space of spherical harmonics.

For the LL radiation reaction force this results in the following replacement:

$$
\nabla_p \cdot (\mathbf{F}_R f)
\longrightarrow \frac{\sigma}{p^2} \mathbf{M}_1 \partial_p \left(p^2 \gamma^2 \beta  \mathbf{f}\right)
-\frac{\sigma \beta}{p} \mathbf{M}_2 \mathbf{f} \,,
$$

where $\mathbf{f} = \left(f_{000}, f_{110}, f_{100}, f_{111}, \dots \right)^{T}$ is a vector
whose components are the expansion coefficients.
Additionally, we introduced the matrices 

$$ \mathbf{M}_1=\mathbf{A}^a \mathbf{A}^b B_a B_b - B^2\mathbf{1} 
\quad \text{and} \quad 
\mathbf{M}_2=\frac{3}{2}\mathbf{A}^a\mathbf{A}^bB_aB_b-\frac{1}{2}B^2\mathbf{1} \,,
$$ 

where $\mathbf{1}$ is the identity matrix and where we implicitly sum over the indices $a$ and $b$. 
Moreover, $\mathbf{A}^{i}$ are the (real) representation matrices of the direction operators
and $\boldsymbol{\Omega}^{i}$ are the ones of the angular momentum operator. 

The following relations have been used to arrive at the above form

$$
\varepsilon_{abc} \mathbf{A}^b \boldsymbol{\Omega}^c \mathbf{A}^a  = -2 \mathbf{1} \quad \text{and} \quad \varepsilon_{abc} B_{a}  B_{d} \mathbf{A}^{b} \boldsymbol{\Omega}^{c}\mathbf{A}^{d} = -\frac{1}{2}\mathbf{A}^a \mathbf{A}^b B_a B_b - \frac{1}{2}B^2 \mathbf{1}\, .
$$

The details of the derivation and a proof that the two identities hold can be found in @cite Harsh2026 .

In the next section, we state which form the radiation reaction force term has in the weak formulation of the VFP equation. 
Though before doing this, we rewrite it using the dimensionless units of @sapphire, see @ref sapphirepp::VFP::ReferenceValues "Reference Values". 
This yields

$$
\frac{3}{2 \underline{\tau}^*_R} \frac{q^{*4}}{m^{*3}} \left[ \frac{1}{p^{*2}} \mathbf{M}^{*}_1 \partial_p \left(p^{*3} \gamma  \mathbf{f}\right)
-\frac{1}{\gamma} \mathbf{M}^{*}_2 \mathbf{f} \right] \,,
$$

where the asterisk indicates that the physical quantities are dimensionless.
The charge and the mass are expressed in terms of a reference charge $\underline{q}$ and a reference mass $\underline{m}$ respectively. 
Furthermore, The reference momentum is $\underline{m}c$ 
and the dimensionless _energy independent_ radiation-reaction timescale is defined by 
$$
\underline{\tau}_R = \frac{9\pi \underline{m}^3 c \underline{\omega}_g}{\mu_0 \underline{q}^4\underline{B}^2} \,,
$$ 
where the inverse of  $\underline{\omega}_g = \underline{q} \underline{B}/ \underline{m}$ is used as the reference time. For the case of an electron with the default reference magnetic field of $1\mu G$, the radiation-reaction timescale is $\underline{\tau}_R \approx 1.3619 \times 10^{22}$.

Note that we drop the asterisk in the following equations;
it is implicit that everything is done @sapphire units.

## Weak formulation

@sapphire uses the dG method for the spatial discretisation of the system of PDEs used to compute the expansion coefficients $\mathbf{f}$. 
The heart of a dG method is the weak formulation of the PDEs. 
Here we merely state it. For more details we refer the reader to @cite Schween2025 .

Because we allow users of @sapphire to choose between a linear and logarithmic momentum variable ($p$ and $\ln p$ respectively),
 there are variants of the weak formulation. 
Additionally, it is possible to scale the distribution function, i.e. $\mathbf{g} = p^{\alpha} \mathbf{f}$, where $\alpha$ is the scaling exponent.

### Linear momentum variable

$$
\begin{split}
& \frac{3}{2 \underline{\tau}_R} \frac{q^4}{m^3} \Bigg\{ \int_{\partial T} \mathbf{h} \cdot
\left[\left( \mathbf{A}^{a}\mathbf{A}^{b}B_{a}B_{b}-B^2\mathbf{1}\right) p \gamma
\right]
\mathbf{f} n_{p} \\
&{}- \int_{T} \frac{\partial \mathbf{h}}{\partial p}\cdot
\left[
  \left(\mathbf{A}^{a}\mathbf{A}^{b}B_{a}B_{b}-B^2\mathbf{1}\right) p \gamma
\right]
\mathbf{f}  \\
&{}+ \int_{T} \mathbf{h}\cdot
\left[
  \left((2-\alpha) \gamma-\frac{3}{2\gamma}\right)\mathbf{A}^{a}\mathbf{A}^{b}B_{a}B_{b}
  -\left((2-\alpha)\gamma - \frac{1}{2\gamma}\right)B^2\mathbf{1}
 \right] 
\mathbf{f}  \Bigg\}
\end{split}
$$


### Logarithmic momentum variable

$$
\begin{split}
& \frac{3}{2 \underline{\tau}_R} \frac{q^4}{m^3} \Bigg\{ \int_{\partial T} \mathbf{h} \cdot
\left[\left( \mathbf{A}^{a}\mathbf{A}^{b}B_{a}B_{b}-B^2\mathbf{1}\right) \gamma
 \right]
\mathbf{f} n_{p} \\
&{}- \int_{T} \frac{\partial \mathbf{h}}{\partial \ln p} \cdot
\left[
  \left(\mathbf{A}^{a}\mathbf{A}^{b}B_{a}B_{b}-B^2\mathbf{1}\right) \gamma
\right] 
\mathbf{f}  \\
&{}+ \int_{T} \mathbf{h}\cdot
\left[
  \left((3-\alpha) \gamma-\frac{3}{2\gamma}\right)\mathbf{A}^{a}\mathbf{A}^{b}B_{a}B_{b}
  - \left((3-\alpha)\gamma - \frac{1}{2\gamma}\right)B^2\mathbf{1}
\right] 
\mathbf{f}  \Bigg\}
\end{split}
$$

## Validation

To validate our implementation of the radiation reaction force, 
we restrict ourselves to the synchrotron radiation of ultra-relativistic electrons, 
namely electrons with $\gamma \in  \mathcal{O}(10^6)$.

We distinguish two cases: 
1. __Pure cooling__, i.e. we instruct @sapphire to solve
$$
\frac{\partial \mathbf{f}}{\partial t} + \frac{3}{2 \underline{\tau}_R} \frac{q^{4}}{m^{3}} \left[ \frac{1}{p^{2}} \mathbf{M}_1 \partial_p \left(p^{3} \gamma  \mathbf{f}\right)
-\frac{1}{\gamma} \mathbf{M}_2 \mathbf{f} \right]  = \mathbf{S} \,,
$$
which models the evolution of a homogeneous ($\mathbf{x}$-independent) particle distribution that changes only because of the radiation reaction force.
2. __Cooling and gyromotion__, i.e. we numerically solve
$$
\frac{\partial \mathbf{f}}{\partial t} - \omega_{a} \boldsymbol{\Omega}^{a} + \frac{3}{2 \underline{\tau}_R} \frac{q^{4}}{m^{3}} \left[ \frac{1}{p^{2}} \mathbf{M}_1 \partial_p \left(p^{3} \gamma  \mathbf{f}\right)
-\frac{1}{\gamma} \mathbf{M}_2 \mathbf{f} \right]  = \mathbf{0} \,,
$$
where $\omega_a = q B_a / \gamma m $ is a gyro frequency 
and the $\boldsymbol{\Omega}^{a}$ matrices introduce the magnetic force that causes the particles to gyrate about the $\mathbf{B}$-field.

In both cases we choose the magnetic field of the background plasma to be aligned with the $x$-axis,
i.e. $\mathbf{B} = B \mathbf{e}_x$ 
and we truncate the spherical harmonic expansion of $f$ at $l_{\mathrm{max}} = 1$,
i.e. we only include a dipole anisotropy. 
The matrices $\mathbf{M}_1$ and $\mathbf{M}_2$ reduce to

$$
\mathbf{M}_1 = B^2
\begin{pmatrix}
				-2/3		&  0   &  0				&  0   \\
				0			& -4/5  &  0				&  0   \\
				0			&  0   & -2/5				&  0   \\
				0			&  0   &  0				& -4/5
\end{pmatrix}
\quad \text{and} \quad

\mathbf{M}_2 = B^2
\begin{pmatrix}
				0		&  0   &  0				&  0   \\
				0			& -1/5  &  0				&  0   \\
				0			&  0   & 2/5				&  0   \\
				0			&  0   &  0				& -1/5
\end{pmatrix} \,.
$$

The explicit form of the rotation generation matrix is 

$$
\boldsymbol{\Omega}_x
=
\begin{pmatrix}
				0		&  0   &  0				&  0   \\
				0			& 0  &  0				&  1   \\
				0			&  0   & 0				&  0   \\
				0			&  -1   &  0				& 0 
\end{pmatrix} \,.
$$

### Pure cooling

If the electron mass $m_e$ and the elementary change $e$ are used as reference quantities, we can set $q = m = 1$. 
Furthermore, for ultra-relativistic electrons $v \approx 1$ and $p \approx \gamma$.

In this approximation, the equation for the isotropic part of the distribution function is 

$$
\partial_t f_{000} - \frac{B^2}{\underline{\tau}_R p^2}
\partial_p \left( p^4 f_{000} \right) = S_{000} \,.
$$

This equation can be solved with the [methods of characteristics](https://en.wikipedia.org/wiki/Method_of_characteristics): Assuming no source, i.e. $S_{000} = 0$, and the initial condition $f_{000}(t = 0, p) = k_{000}(p)$ with $k_{000}(p)$ being an arbitrary function, the solution is 

$$
g_{000}(t, p) = h_{000}\left(\frac{p}{1 - p B^2 t/\underline{\tau}_R}\right) \,,
$$

where $g_{000} = p^4 f_{000}$ and $h_{000} = p^4 k_{000}$, i.e. we scaled the distribution function 
($\alpha = 4$).
Note in @sapphire the scaling exponent $\alpha$ is set to three. 

This shows that the cooling of the particles, namely the energy loss due to the radiation reaction force, happens at a different rate for particles with different energies.
This is, as expected for the LL radiation reaction force, in agreement with [Larmor's formula](https://en.wikipedia.org/wiki/Larmor_formula);
the radiated power is proportional to $\gamma^2$.

This can be highlighted with the introduction of a cooling time scale. We define
$$
 \tau_{\mathrm{c}}(p) = \frac{\underline{\tau}_R}{p B^2} \,.
$$

Note that if we had not included the factor $2/3$ in the definition of $\underline{\tau}_R$,
it would have been explicit in $\tau_c$. 

We proceed analogously with the equation for the dipole component $f_{100}$ 
and only state that the procedure is the same for $f_{110}$ and $f_{111}$. 
It is

$$
\partial_t f_{100} - \frac{3 B^2}{5 \underline{\tau}_R p^2}
\partial_p \left( p^4 f_{100} \right) - \frac{3 B^2}{5 \underline{\tau}_R p} f_{100} = 0 \,.
$$

It differs structurally from the equation for the isotropic part $f_{000}$ in having an additional reaction term.
Though, the methods of characteristics can be still be used to solve it.
The solution is 

$$
	g_{100}(t, p) = h_{100}\left(\frac{p}{1 - 3 t/5 \tau_c}\right)
	\exp\left(\frac{3}{5 p^2} \frac{t}{\tau_c} \left[1 - \frac{3}{10} \frac{t}{\tau_c}\right] \right) \,.
$$

We emphasise that the exponential correction factor is strongly suppressed.

We now show the @sapphire and the analytic solutions for the initial conditions $k_{000}(p) = N p^{-4} \exp(-p/p_{\mathrm{max}})$, 
i.e. a power-law with exponential cut-off and an arbitrary normalisation $N$ 
and 
$k_{100}(p) = \frac{N}{2 \sqrt{3}} p^{-4} \exp(-p/p_{\mathrm{max}})$ .

<CENTER>
<img src="https://sapphirepp.org/img/implementation/synchrotron/isotropic-cooling.gif"
alt="Time-dependent isotropic synchrotron cooling vs analytic." width="60%"/>
</CENTER>

The numerical curve follows the analytic cooling trajectory exactly.

At the end of the pure-cooling test case, we demonstrate that a continuous injection of mono-energetic particles, that are subsequently cooled, 
results in a power-law energy spectrum with spectral index minus four. 

Therefore, we again solve the equation for the isotropic part $f_{000}$.
This time there are no particles at $t = 0$, i.e. $k_{000} = 0$. 
However, there is a non-trivial source term, namely

$$
S_{000}(t, p) = \frac{Q}{\sqrt{4 \pi} p^2_{\mathrm{inj}}} \delta\left(p - p_{\mathrm{inj}}\right) H(t)\, ,
$$

where $p_{\mathrm{inj}}$ is the injection energy of the particles, $Q$ is the injection rate 
and $H(t)$ is a [Heaviside step function](https://en.wikipedia.org/wiki/Heaviside_step_function). 

Keeping in mind that the (weak) derivative of a step function is a Delta distribution, 
and that $\delta(f(x)) = \sum_i \delta(x - x_i)/|f'(x_i)|$, where $x_i$ are the roots of $f(x)$, 
the method of characteristics gives

$$
g_{000}(t, p)
= \frac{Q}{\sqrt{4\pi}}  p_{\mathrm{inj}} \tau_{c,\mathrm{inj}}
H\left(p - p_{\mathrm{inj}}\right)
H\left(t - (\tau_c - \tau_{c, \mathrm{inj}})\right)\,.
$$

$\tau_{c,\mathrm{inj}}$ is the cooling time corresponding to the injection energy.
The difference $\tau_c - \tau_{c, \mathrm{inj}}$ is the time it takes to decrease the energy of a particle from its injection energy to the energy $p$.
The second Heaviside step function ensures that there are no particles with momenta smaller than $p$ at time $t$: If a particle has a momentum smaller $p$, 
then the corresponding difference between its cooling time and the cooling time of the injection momentum is larger than the one for $p$ 
and $t$ smaller than this difference 
and, hence, the step function evaluates to zero.

Picture 

### Cooling and gyromotion

In this test case, we include the magnetic force in our computations. 
This leads to the inclusion of the matrix $\boldsymbol{\Omega}_x$ in our system of PDEs
that determines the expansion coefficients $\mathbf{f}$. 
Physically, the particles now rotate about the magnetic field with angular frequency $\omega_x = B/\gamma$ 
_and_ loose energy due to the back reaction of the synchrotron radiation that they emit.

The explicit form of $\boldsymbol{\Omega}_x$ shows that only the equations for the dipole components $f_{110}$ and $f_{111}$ change. 
The evolution of isotropic part and the dipole component pointing into the $x$-direction is unaffected. 

To derive an analytic expression for the dipole components perpendicular to the $\mathbf{B}$-field,
we define 

$$d_\perp = p^4 (f_{110}+ \mathrm{i} f_{111})$$

and add up the equations for $f_{110}$ and $f_{111}$, 
multiplied in agreement with the above definition.
This yields

$$
\partial_t d_{\perp} - \frac{ 6 B^2 p^2}{ 5 \underline{\tau}_R} \partial_p d_{\perp}
= \left(\frac{3 B^2}{10 \underline{\tau}_R p} -\mathrm{i} \omega_x\right) d_{\perp} \,.
$$

Its solution is 

$$
d_{\perp}(t, p) = d_{\perp,0}\left( \frac{p}{1 - 6 t/5 \tau_c}\right)
\exp(-\mathrm{i} \omega_x t) 
\exp\left(\frac{3}{10 p^2} \frac{t}{\tau_c} \left[\frac{3}{5} \frac{t}{\tau_c} - 1 \right]\right)
\,,
$$

where $d_{\perp, 0} = h_{110}(p) + \mathrm{i} h_{111}(p)$ are the initial conditions.

<CENTER>
<img src="https://sapphirepp.org/img/implementation/synchrotron/anisotropy-rotation.gif"
alt="Anisotropic l=1 cooling with rotation: precession and damping." width="60%"/>
</CENTER>

Turning on rotation couples $(f_{110},f_{111})$ and produces oscillatory exchange with amplitude damping at the cooling rate where, the phase follows $\omega_g$. @sapphire matches the analytic solution.

<div class="section_buttons">

| Previous                          |
| :-------------------------------- |
| [Implementation](#implementation) |

</div>

---

@author Harsh Goyal (<harsh.goyal@mpi-hd.mpg.de>)
@author Nils Schween (<nils.schween@mpi-hd.mpg.de>)
@date 2026-03-19
