# Radiation Reaction {#radiation-reaction}

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
Additionally, we introduced the matrices $\mathbf{M}_1=\mathbf{A}^a \mathbf{A}^b B_a B_b - B^2\mathbf{1}$ and $\mathbf{M}_2=\frac{3}{2}\mathbf{A}^a\mathbf{A}^bB_aB_b-\frac{1}{2}B^2\mathbf{1}$, 
where $\mathbf{1}$ is the identity matrix. 
Moreover, $\mathbf{A}^{i}$ are the representation matrices of the direction operators
and $\boldsymbol{\Omega}^{i}$ are the ones of the angular momentum operator. 

The following relations have been used to arrive at the above form

$$
\varepsilon_{abc} \mathbf{A}^b \boldsymbol{\Omega}^c \mathbf{A}^a  = -2 \mathbf{1} \quad \text{and} \quad \varepsilon_{abc} B_{a}  B_{d} \mathbf{A}^{b} \boldsymbol{\Omega}^{c}\mathbf{A}^{d} = -\frac{1}{2}\mathbf{A}^a \mathbf{A}^b B_a B_b - \frac{1}{2}B^2 \mathbf{1}\, .
$$

The details of the derivation and a proof that the two identities hold can be found in @cite Harsh2026 .

In the next section, we state which form the radiation reaction force term has in the weak formulation of the VFP equation. 
Though before doing this, we rewrite it using the dimensionless units of @sapphire . This yields

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

@sapphire uses the dG method for the spatial discretisation of the system of PDEs used to compute the expansion coefficients $\mathbf{h}$. 
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

## Analytic Solutions

We add the radiation reaction force straightforwardly to the momentum part of the VFP equation, i.e.

$$
\frac{\partial f}{\partial t}  + \dots = -\nabla_p \cdot (\mathbf{F}_R f).
\quad (6)
$$

### Isotropic and anisotropic solutions (ultra-relativistic)

To expose anisotropy at first order, we use the linear expansion
$$
f = f_0 + \boldsymbol{f_1} \cdot \boldsymbol{n}, \quad (17)
$$

where $\boldsymbol{n}$ denotes the unit vector in the direction of the particle velocity (i.e.\ $\boldsymbol{n} \equiv \boldsymbol{v}/|\boldsymbol{v}| = \boldsymbol{\beta}/|\boldsymbol{\beta}|$).

Projecting the radiation reaction term onto the isotropic($f_0$) and dipole($f_1$) sectors gives

$$
\left(\frac{\partial f_0}{\partial t}\right)_{\!\rm sync}
=\frac{1}{p^2}\frac{\partial}{\partial p}\!\left[\frac{2\sigma B^2}{3}\,p^2\,\gamma^2\beta\,f_0\right].
\quad (18)
$$

$$
\left(\frac{\partial f_{1i}}{\partial t}\right)_{\!\rm sync}
= -\,\frac{\sigma}{5p^2}\frac{\partial}{\partial p}\!\left[ -4B^2p^2(\gamma^2\beta f_{1i})
+2p^2 B_i B_j(\gamma^2\beta f_{1j})\right]
-\frac{\sigma}{5p}\!\left[-3B_iB_j(\beta f_{1j})+B^2(\beta f_{1i})\right].
\quad (19)
$$

Using $\gamma^2\beta \simeq (p/mc)^2$, Eq. (18) reduces to an advection–reaction
equation in $p$ with solution (characteristics)

$$
p^3 f_0(t,p)=k\!\left(\frac{p}{1-t/t_{\rm cool}}\right)\,\frac{1}{1-t/t_{\rm cool}},\qquad
t_{\rm cool}=\frac{3(m c)^2}{2\sigma B^2 p}
=\frac{9\pi m^3 c}{\mu_0 e^4 B^2\,\gamma}\ \text{(SI)}.
\quad (20)
$$

Similarly, writing $\mathbf{f}_1=(f_{1x},f_{1y},f_{1z})^\top$ and diagonalising $B\otimes B=\mathbf{V}\Lambda\mathbf{V}^\top$,
the rotated components $\tilde{\mathbf{f}}_1=\mathbf{V}^\top\mathbf{f}_1$ decouple.

### Including rotation (gyromotion)

Rotation couples the real pair $(f_{110},f_{111})$ through $\omega_g=eB/(\gamma m)$.
Defining $f_\perp=f_{110}+i f_{111}$ yields a single complex advection reaction
equation with a $-i\omega_g f_\perp$ term.

## Examples

### 1) Pure cooling, isotropic and anisotropic test

We compare numerical $f_{000}(t,p)$ and  $f_{100}(t,p)$ against analytical results for uniform $B$, $p\in[p_{\min},p_{\max}]$,
and initial $k(p)=p^{-1}\exp(-p/p_{\max})$ (i.e. $f_0\propto p^{-4}$ with cutoff; scaling index $\alpha=3$).
The numerical curve follows the analytic cooling trajectory exactly.

<CENTER>
<img src="https://sapphirepp.org/img/implementation/synchrotron/isotropic-cooling.gif"
alt="Time-dependent isotropic synchrotron cooling vs analytic." width="60%"/>
</CENTER>

### 2) l = 1 anisotropy with rotation enabled

Turning on rotation couples $(f_{110},f_{111})$ and produces oscillatory exchange with amplitude damping at the cooling rate where, the phase follows $\omega_g$. @sapphire matches the analytic solution.

<CENTER>
<img src="https://sapphirepp.org/img/implementation/synchrotron/anisotropy-rotation.gif"
alt="Anisotropic l=1 cooling with rotation: precession and damping." width="60%"/>
</CENTER>

<div class="section_buttons">

| Previous                          |
| :-------------------------------- |
| [Implementation](#implementation) |

</div>

---

@author Harsh Goyal (<harsh.goyal@mpi-hd.mpg.de>)
@date 2025-11-13
