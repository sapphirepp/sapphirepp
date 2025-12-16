# Radiation Reaction {#radiation-reaction}

@sapphire solves the Vlasov–Fokker–Planck (VFP) system in a spherical–harmonic basis.
This page documents how **radiation reaction** is derived, inserted
into the DG weak formulation, and validated against analytic solutions. The implementation
uses the reduced Landau–Lifshitz (LL) radiation–reaction force and supports *linear/log*
momentum coordinates and *unscaled/scaled* distributions.

## From Landau–Lifshitz to VFP

An accelerated charged particle emits electromagnetic radiation and, as a consequence,
experiences a recoil force due to the loss of energy and momentum carried away by the 
emitted radiation. This back-reaction is described by the Landau–Lifshitz radiation-reaction equation,
which supplements the Lorentz force law with an additional radiation-reaction term.

The following equation for the radiation damping force has been adapted from @cite landau1975classical .

$$
\mathbf{F_R} = \frac{\mu_0q^{3}\gamma}{6\pi mc} 
 \Bigg\{ \left( \frac{\partial}{\partial t} + \mathbf{v} \cdot \nabla \right)(\mathbf{E}+\mathbf{v} \times \mathbf{B}) \Bigg\} + \frac{\mu_0q^{4}}{6\pi m^{2}c} \Bigg\{\frac{1}{c^2} \mathbf{E} (\mathbf{v} \cdot \mathbf{E}) - \mathbf{B} \times (\mathbf{E} +\mathbf{v} \times \mathbf{B}) \Bigg\} - \frac{\mu_0q^{4} \gamma^2}{6\pi m^{2}c^{3}} \Bigg[\left( \mathbf{E} + \mathbf{v} \times \mathbf{B} \right)^{2} - \frac{1}{c^{2}} (\mathbf{E} \cdot \mathbf{v})^{2} \Bigg] \mathbf{v} .
\quad (1)
$$

Following @cite tamburini2010radiation , the first curly–bracketed term containing derivatives of the fields is dropped because, for ultra-relativistic electrons $$(\gamma \gg 1) \quad (2),$$ its contribution is orders of magnitude smaller than the dominant quadratic terms. The remaining part of the Landau–Lifshitz force therefore captures the principal radiation–reaction effects responsible for synchrotron energy loss.

We now derive the expression of the synchrotron cooling term relevant to the evolution of the particle distribution. The analysis is carried out under the **ideal magnetohydrodynamic (MHD) approximation**, which assumes infinite plasma conductivity. In this limit, the electric field in the fluid rest frame vanishes, and in the laboratory frame it satisfies
$$
\mathbf{E} = -\,\mathbf{u}\!\times\!\mathbf{B}.
\quad (3)
$$

where, u denotes the bulk plasma (fluid) velocity. This approximation ensures that the electromagnetic fields evolve in a perfectly conducting medium, a valid assumption for many astrophysical plasmas where the magnetic Reynolds number is extremely large.

Substituting the ideal-MHD condition and neglecting the derivative term as justified above, the Landau–Lifshitz force reduces to

$$
\mathbf{F_R} = \sigma \Big( \mathbf{B} (\mathbf{B} \cdot \boldsymbol{\beta}) + \gamma^2 \boldsymbol{\beta} ((\boldsymbol{\beta}\cdot \mathbf{B})^2 - B^2) \Big),
\quad (4)
$$

with the coupling
$$
\sigma =  \frac{\mu_0q^4}{6\pi m^2} .
\quad (5)
$$
Here, $$\gamma = (1 - \beta^{2})^{-1/2} \quad (6) $$ is the Lorentz factor, and $$\boldsymbol{\beta} = \mathbf{v}/c \quad (7) $$ is the dimensionless velocity vector of the charged particle. The first term in Eq.(4) represents the projection of the velocity along the magnetic field, while the second term accounts for the perpendicular component responsible for synchrotron energy losses.

Finally, the loss terms due to the radiation raection are represented as

$$
\left(\frac{\partial f}{\partial t}\right)
=-\nabla_p \cdot (\mathbf{F_R}f).
\quad (8)
$$

To expose anisotropy at first order, we use the linear expansion
$$
f = f_0 + \boldsymbol{f_1} \cdot \boldsymbol{n}, \quad (9)
$$

where $\boldsymbol{n}$ denotes the unit vector in the direction of the particle velocity (i.e.\ $\boldsymbol{n} \equiv \boldsymbol{v}/|\boldsymbol{v}| = \boldsymbol{\beta}/|\boldsymbol{\beta}|$).

Projecting the Landau–Lifshitz radiation–reaction term onto the isotropic($f_0$) and dipole($f_1$) sectors gives

$$
\left(\frac{\partial f_0}{\partial t}\right)_{\!\rm sync}
=\frac{1}{p^2}\frac{\partial}{\partial p}\!\left[\frac{2\sigma B^2}{3}\,p^2\,\gamma^2\beta\,f_0\right].
\quad (10)
$$

$$
\left(\frac{\partial f_{1i}}{\partial t}\right)_{\!\rm sync}
= -\,\frac{\sigma}{5p^2}\frac{\partial}{\partial p}\!\left[ -4B^2p^2(\gamma^2\beta f_{1i})
+2p^2 B_i B_j(\gamma^2\beta f_{1j})\right]
-\frac{\sigma}{5p}\!\left[-3B_iB_j(\beta f_{1j})+B^2(\beta f_{1i})\right].
\quad (11)
$$

The same physics in the spherical–harmonic operator form gives an equivalent expression which reads (for each $(\ell,m)$)

$$
\left(\partial_t f_{\ell}^{m}\right)_{\rm sync}
= -\frac{\sigma}{p^2}\,\mathbf{M}_1\,\partial_p\!\left(p^2\gamma^2\beta\, f_{\ell}^{m}\right)
+\frac{\sigma}{p}\,\mathbf{M}_2\,(\beta\,f_{\ell}^{m}),
\quad (12)
$$

with $\mathbf{M}_1=A^aA^bB_aB_b-B^2\mathbf{1}$, $\mathbf{M}_2=\frac{3}{2}A^aA^bB_aB_b-\frac{1}{2}B^2\mathbf{1}. \quad (13)$  

Finally, @sapphire uses the **modified VFP equation** with the synchrotron term moved to the LHS, i.e.

$$
\frac{\partial f}{\partial t} + (u+v)\!\cdot\!\nabla_x f
-\gamma m\,\frac{D u}{D t}\!\cdot\nabla_p f
- p.(\nabla_x u)\,\nabla_p f
+\nabla_p \cdot (\mathbf{F_R}f)
+ q\,\boldsymbol{v}\!\cdot\!(\boldsymbol{B}\times\nabla_p f)
= \frac{\nu}{2}\,\Delta_{\theta,\varphi} f + S.
\quad (14)
$$

This is the theoretical form implemented.

## Weak formulation (dimensionless units)

The following weak forms show how the radiation reaction term enters the DG formulation for
different combinations of **momentum coordinates** (linear or logarithmic)
and **distribution scaling** non-scaled (f) or scaled ($g=p^{\alpha}f$).

### Linear, Non-Scaled

$$
\begin{split}
& \text{coeff} \times \Bigg\{ \int_{\partial T} \boldsymbol{\Phi}_{i} \cdot
\left[\left(\!\mathbf{A}^{a}\mathbf{A}^{b}B_{a}B_{b}-B^2\mathbf{1}\right) p \gamma
\!\right]\!
\zeta_{j}\boldsymbol{\Phi}_{j} n_{p} \\
&{}- \int_{T} \partial_{p}\boldsymbol{\Phi}_{i}\cdot
\left[\!
  \left(\mathbf{A}^{a}\mathbf{A}^{b}B_{a}B_{b}-B^2\mathbf{1}\right) p \gamma
\!\right]\!
\zeta_{j}\boldsymbol{\Phi}_{j} \\
&{}+ \int_{T} \boldsymbol{\Phi}_{i}\cdot
\left[\!
  \left(4 \gamma-\frac{5}{2 \gamma}\right)\mathbf{A}^{a}\mathbf{A}^{b}B_{a}B_{b}
  +\left(\frac{3}{2 \gamma}-4\gamma\right)B^2\mathbf{1}
\!\right]
\zeta_{j}\boldsymbol{\Phi}_{j} \Bigg\}.
\end{split}
\quad (15)
$$

### Logarithmic, Non-Scaled

$$
\begin{split}
& \text{coeff} \times \Bigg\{ \int_{\partial T} \boldsymbol{\Phi}_{i} \cdot
\left[\left(\!\mathbf{A}^{a}\mathbf{A}^{b}B_{a}B_{b}-B^2\mathbf{1}\right) \gamma
\!\right]\!
\zeta_{j}\boldsymbol{\Phi}_{j} n_{p} \\
&{}- \int_{T} \frac{\partial \boldsymbol{\Phi}_{i}}{\partial \ln p}\boldsymbol{\Phi}_{i}\cdot
\left[\!
  \left(\mathbf{A}^{a}\mathbf{A}^{b}B_{a}B_{b}-B^2\mathbf{1}\right) \gamma
\!\right]\!
\zeta_{j}\boldsymbol{\Phi}_{j} \\
&{}+ \int_{T} \boldsymbol{\Phi}_{i}\cdot
\left[\!
  \left(4 \gamma-\frac{5}{2 \gamma}\right)\mathbf{A}^{a}\mathbf{A}^{b}B_{a}B_{b}
  +\left(\frac{3}{2 \gamma}-4\gamma\right)B^2\mathbf{1}
\!\right]
\zeta_{j}\boldsymbol{\Phi}_{j} \Bigg\}.
\end{split}
\quad (16)
$$

### Linear, Scaled

$$
\begin{split}
& \text{coeff} \times \Bigg\{ \int_{\partial T} \boldsymbol{\Phi}_{i} \cdot
\left[\left(\!\mathbf{A}^{a}\mathbf{A}^{b}B_{a}B_{b}-B^2\mathbf{1}\right) p \gamma
\!\right]\!
\zeta_{j}\boldsymbol{\Phi}_{j} n_{p} \\
&{}- \int_{T} \partial_{p}\boldsymbol{\Phi}_{i}\cdot
\left[\!
  \left(\mathbf{A}^{a}\mathbf{A}^{b}B_{a}B_{b}-B^2\mathbf{1}\right) p \gamma
\!\right]\!
\zeta_{j}\boldsymbol{\Phi}_{j} \\
&{}+ \int_{T} \boldsymbol{\Phi}_{i}\cdot
\left[\!
  \left((4-\alpha) \gamma-\frac{5}{2 \gamma}\right)\mathbf{A}^{a}\mathbf{A}^{b}B_{a}B_{b}
  +\left(\frac{3}{2 \gamma}-(4-\alpha)\gamma\right)B^2\mathbf{1}
\!\right]
\zeta_{j}\boldsymbol{\Phi}_{j} \Bigg\}.
\end{split}
\quad (17)
$$

### Logarithmic, Scaled

$$
\begin{split}
& \text{coeff} \times \Bigg\{ \int_{\partial T} \boldsymbol{\Phi}_{i} \cdot
\left[\left(\!\mathbf{A}^{a}\mathbf{A}^{b}B_{a}B_{b}-B^2\mathbf{1}\right) \gamma
\!\right]\!
\zeta_{j}\boldsymbol{\Phi}_{j} n_{p} \\
&{}- \int_{T} \frac{\partial \boldsymbol{\Phi}_{i}}{\partial \ln p}\boldsymbol{\Phi}_{i}\cdot
\left[\!
  \left(\mathbf{A}^{a}\mathbf{A}^{b}B_{a}B_{b}-B^2\mathbf{1}\right) \gamma
\!\right]\!
\zeta_{j}\boldsymbol{\Phi}_{j} \\
&{}+ \int_{T} \boldsymbol{\Phi}_{i}\cdot
\left[\!
  \left((4-\alpha) \gamma-\frac{5}{2 \gamma}\right)\mathbf{A}^{a}\mathbf{A}^{b}B_{a}B_{b}
  +\left(\frac{3}{2 \gamma}-(4-\alpha)\gamma\right)B^2\mathbf{1}
\!\right]
\zeta_{j}\boldsymbol{\Phi}_{j} \Bigg\}.
\end{split}
\quad (18)
$$

With, $$\text{coeff} = \frac{\mu_0 q^4\underline{B}^2}{6\pi m^3 c \underline{\omega_g}} \quad (19)$$ 

here,
$$ \mu_0 \: \text{is the permeability of free space,}\\ q \: \text{is the particle charge,}\\ m \:  \text{is the particle mass,}\\ c \:  \text{is the speed of light,}\\ \underline{B} \: \text{is the reference magnetic field,}\\ \underline{\omega_g} \:  \text{is the reference gyrofrequency.}
$$


## Analytic checks

### Isotropic solution (ultra-relativistic)

Let $g_0=p^3 f_0$. Using $\gamma^2\beta \simeq (p/mc)^2$, Eq. (10) reduces to an advection–reaction
equation in $p$ with solution (characteristics)

$$
g_0(t,p)=k\!\left(\frac{p}{1-t/t_{\rm cool}}\right)\,\frac{1}{1-t/t_{\rm cool}},\qquad
t_{\rm cool}=\frac{3(m c)^2}{2\sigma B^2 p}
=\frac{9\pi m^3 c}{\mu_0 e^4 B^2\,\gamma}\ \text{(SI)}.
\quad (5)
$$

### l=1 anisotropy (no rotation)

Writing $\mathbf{f}_1=(f_{1x},f_{1y},f_{1z})^\top$ and diagonalising $B\otimes B=\mathbf{V}\Lambda\mathbf{V}^\top$,
the rotated components $\tilde{\mathbf{f}}_1=\mathbf{V}^\top\mathbf{f}_1$ decouple. Each component’s scaled
variable $g_{1i}=p^3\tilde f_{1i}$ obeys a scalar advection–reaction equation with coefficients set by
the eigenvalue $\lambda_i$.

### Including rotation (gyromotion) and complexification

Rotation couples the real pair $(f_{110},f_{111})$ through $\omega_g=eB/(\gamma m)$.
Defining $f_\perp=f_{110}+i f_{111}$ and $g_\perp=p^3 f_\perp$ yields a single complex advection–reaction
equation with a $-i\omega_g f_\perp$ term; the cooled spectrum is modulated by phase rotation at
$\omega_g$.

## Examples

### 1) Pure cooling, isotropic test

We compare numerical $g_0(t,p)$ against Eq. (5) for uniform $B$, $p\in[p_{\min},p_{\max}]$,
and initial $k(p)=p^{-1}\exp(-p/p_{\max})$ (i.e. $f_0\propto p^{-4}$ with cutoff; scaling index $\alpha=3$).
The numerical curve follows the analytic cooling trajectory exactly.

<CENTER>
<img src="https://sapphirepp.org/img/implementation/synchrotron/isotropic-cooling.gif"
alt="Time-dependent isotropic synchrotron cooling vs analytic." width="60%"/>
</CENTER>

### 2) l = 1 anisotropy (no rotation)

Initialise $(k_{1x},k_{1y},k_{1z})=\big(\tfrac{0.5}{\sqrt{3}}p^{-1}e^{-p/p_{\max}},0,0\big)$; the parallel mode and the
two transverse modes evolve with different damping rates governed by the eigenvalues of
$B\otimes B$. Numerics reproduce the analytic prediction.

<CENTER>
<img src="https://sapphirepp.org/img/implementation/synchrotron/anisotropy-l1.png"
alt="l=1 anisotropic cooling: numerical vs analytic." width="60%"/>
</CENTER>

### 3) l = 1 anisotropy with rotation enabled

Turning on rotation couples $(f_{110},f_{111})$ and produces oscillatory exchange with amplitude
damping at the cooling rate; the phase follows $\omega_g$. @sapphire matches the analytic complex
solution.

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
