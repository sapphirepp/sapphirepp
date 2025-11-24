# Synchrotron cooling {#synchrotron-radiation-implementation}

@sapphire solves the Vlasov–Fokker–Planck (VFP) system in a spherical–harmonic basis.
This page documents how **synchrotron radiation reaction (cooling)** is derived, inserted
into the DG weak formulation, and validated against analytic solutions. The implementation
uses the reduced Landau–Lifshitz (LL) radiation–reaction force and supports *linear/log*
momentum coordinates and *unscaled/scaled* distributions (with $g=p^{3}f$).

## From Landau–Lifshitz to VFP

Projecting the LL radiation–reaction term onto the isotropic and $l=1$ dipole sectors gives
(for electrons; SI units via $\sigma=\mu_0 q^4/(6\pi m^2)$)

$$
\left(\frac{\partial f_0}{\partial t}\right)_{\!\rm sync}
=\frac{1}{p^2}\frac{\partial}{\partial p}\!\left[\frac{2\sigma B^2}{3}\,p^2\,\gamma^2\beta\,f_0\right].
\quad (1)
$$

$$
\left(\frac{\partial f_{1i}}{\partial t}\right)_{\!\rm sync}
= -\,\frac{\sigma}{5p^2}\frac{\partial}{\partial p}\!\left[ -4B^2p^2(\gamma^2\beta f_{1i})
+2p^2 B_i B_j(\gamma^2\beta f_{1j})\right]
-\frac{\sigma}{5p}\!\left[-3B_iB_j(\beta f_{1j})+B^2(\beta f_{1i})\right].
\quad (2)
$$

The same physics in the spherical–harmonic operator form reads (for each $(\ell,m)$)

$$
\left(\partial_t f_{\ell}^{m}\right)_{\rm sync}
= -\frac{\sigma}{p^2}\,\mathbf{M}_1\,\partial_p\!\left(p^2\gamma^2\beta\, f_{\ell}^{m}\right)
+\frac{\sigma}{p}\,(\mathbf{M}_2+2\mathbf{S})\,(\beta\,f_{\ell}^{m}),
\quad (3)
$$

with $\mathbf{M}_1=A^aA^bB_aB_b-B^2\mathbf{1}$, $\mathbf{M}_2=i\varepsilon_{abc}B^aA^b\Omega^c\,B^dA_d$,
and $\mathbf{S}=A^aA^bB_aB_b$. The Cartesian and spherical forms are equivalent via the
basis change $\mathbf{C}$: $\mathbf{C}^{-1}\mathbf{M}_1\mathbf{C}=\mathbf{M}_3$,
$\mathbf{C}^{-1}(\mathbf{M}_2+2\mathbf{S})\mathbf{C}=\mathbf{M}_4$.

Finally, @sapphire uses the **modified VFP equation** with the synchrotron term moved to the LHS, i.e.

$$
\frac{\partial f}{\partial t} + (u+v)\!\cdot\!\nabla_x f
-\gamma m\,\frac{D u}{D t}\!\cdot\nabla_p f
- p\!:\!(\nabla_x u)\,\nabla_p f
-\left(\frac{\partial f}{\partial t}\right)_{\!\rm sync}
+ q\,\boldsymbol{v}\!\cdot\!(\boldsymbol{B}\times\nabla_p f)
= \frac{\nu}{2}\,\Delta_{\theta,\varphi} f + S.
\quad (4)
$$

This is the theoretical form implemented.

## Weak formulation (dimensionless units)

The following weak forms show how the synchrotron term enters the DG formulation for
different combinations of **momentum coordinates** (linear or logarithmic)
and **distribution scaling** (non-scaled or scaled).

### Linear, Non-Scaled

$$
\begin{split}
&+ \int_{\partial T} \boldsymbol{\Phi}*{i} \cdot
\left[\!
  \left(\mathbf{A}^{a}\mathbf{A}^{b}B*{a*}B_{b*}-B_*^2\mathbf{1}\right)
  \frac{3\gamma p_* q_*^4}{2\bar{\tau}_{s*} m_*^2}
\!\right]\!
\zeta_{j}\boldsymbol{\Phi}*{j} n*{p*} \\
&{}- \int_{T} \partial_{p*}\boldsymbol{\Phi}*{i}\cdot
\left[\!
  \left(\mathbf{A}^{a}\mathbf{A}^{b}B*{a*}B_{b*}-B_*^2\mathbf{1}\right)
  \frac{3\gamma p_* q_*^4}{2\bar{\tau}_{s*} m_*^2}
\!\right]\!
\zeta_{j}\boldsymbol{\Phi}*{j} \\
&{}+ \int*{T} \boldsymbol{\Phi}*{i}\cdot
\frac{3 q**^4}{2\bar{\tau}*{s*} m**^2}
\left[\!
  \left(2\gamma-\frac{3}{\gamma}\right)\mathbf{A}^{a}\mathbf{A}^{b}B_{a*}B_{b*}
  +\left(\frac{1}{\gamma}-2\gamma\right)B_*^2\mathbf{1}
\!\right]
\zeta_{j}\boldsymbol{\Phi}_{j}.
\end{split}
$$

### Logarithmic, Non-Scaled

$$
\begin{split}
&+ \int_{\partial T} \boldsymbol{\Phi}*{i}\cdot
\frac{3 q**^4 \gamma}{2\bar{\tau}_{s*} m_*^2}
\left(\mathbf{A}^{a}\mathbf{A}^{b}B_{a*}B_{b*}-B_*^2\mathbf{1}\right)
\zeta_{j}\boldsymbol{\Phi}*{j} n*{\ln p*} \\
&{}- \int_{T} \frac{\partial \boldsymbol{\Phi}*{i}}{\partial \ln p**}\cdot
\frac{3 q_*^4 \gamma}{2\bar{\tau}*{s*} m**^2}
\left(\mathbf{A}^{a}\mathbf{A}^{b}B_{a*}B_{b*}-B_*^2\mathbf{1}\right)
\zeta_{j}\boldsymbol{\Phi}*{j} \\
&{}+ \int*{T} \boldsymbol{\Phi}*{i}\cdot
\frac{3 q**^4}{2\bar{\tau}*{s*} m**^2}
\left[\!
  \left(3\gamma-\frac{3}{\gamma}\right)\mathbf{A}^{a}\mathbf{A}^{b}B_{a*}B_{b*}
  +\left(\frac{1}{\gamma}-3\gamma\right)B_*^2\mathbf{1}
\!\right]
\zeta_{j}\boldsymbol{\Phi}_{j}.
\end{split}
$$

### Linear, Scaled

$$
\begin{split}
&+ \int_{\partial T} \boldsymbol{\Phi}*{i}\cdot
\left[\!
  \left(\mathbf{A}^{a}\mathbf{A}^{b}B*{a*}B_{b*}-B_*^2\mathbf{1}\right)
  \frac{3 q_*^4 \gamma p_*}{2\bar{\tau}_{s*} m_*^2}
\!\right]\!
\zeta_{j}\boldsymbol{\Phi}*{j} n*{p*} \\
&{}- \int_{T} \partial_{p*}\boldsymbol{\Phi}*{i}\cdot
\left[\!
  \left(\mathbf{A}^{a}\mathbf{A}^{b}B*{a*}B_{b*}-B_*^2\mathbf{1}\right)
  \frac{3 q_*^4 \gamma p_*}{2\bar{\tau}_{s*} m_*^2}
\!\right]\!
\zeta_{j}\boldsymbol{\Phi}*{j} \\
&{}- \int*{T} \boldsymbol{\Phi}*{i}\cdot
\frac{3 q**^4}{2\bar{\tau}*{s*} m**^2}
\left[\!
  \left((\alpha-2)\gamma+\frac{3}{\gamma}\right)\mathbf{A}^{a}\mathbf{A}^{b}B_{a*}B_{b*}
  -\left((\alpha-2)\gamma+\frac{1}{\gamma}\right)B_*^2\mathbf{1}
\!\right]
\zeta_{j}\boldsymbol{\Phi}_{j}.
\end{split}
$$

### Logarithmic, Scaled

$$
\begin{split}
&+ \int_{\partial T} \boldsymbol{\Phi}*{i}\cdot
\frac{3 q**^4 \gamma}{2\bar{\tau}_{s*} m_*^2}
\left(\mathbf{A}^{a}\mathbf{A}^{b}B_{a*}B_{b*}-B_*^2\mathbf{1}\right)
\zeta_{j}\boldsymbol{\Phi}*{j} n*{\ln p*} \\
&{}- \int_{T} \frac{\partial \boldsymbol{\Phi}*{i}}{\partial \ln p**}\cdot
\frac{3 q_*^4 \gamma}{2\bar{\tau}*{s*} m**^2}
\left(\mathbf{A}^{a}\mathbf{A}^{b}B_{a*}B_{b*}-B_*^2\mathbf{1}\right)
\zeta_{j}\boldsymbol{\Phi}*{j} \\
&{}- \int*{T} \boldsymbol{\Phi}*{i}\cdot
\frac{3 q**^4}{2\bar{\tau}*{s*} m**^2}
\left[\!
  \left((\alpha-3)\gamma+\frac{3}{\gamma}\right)\mathbf{A}^{a}\mathbf{A}^{b}B_{a*}B_{b*}
  -\left((\alpha-3)\gamma+\frac{1}{\gamma}\right)B_*^2\mathbf{1}
\!\right]
\zeta_{j}\boldsymbol{\Phi}_{j}.
\end{split}
$$

## Momentum coordinate and scaling

The solver supports either **linear** $p$ or **logarithmic** $\ln p$ grids and can evolve either the
**unscaled** distribution $f$ or the **scaled** distribution $g=p^3 f$. Evolving $g$ improves conditioning
on wide momentum ranges and aligns directly with the analytic solutions below.

## Analytic checks

### Isotropic solution (ultra-relativistic)

Let $g_0=p^3 f_0$. Using $\gamma^2\beta \simeq (p/mc)^2$, Eq. (1) reduces to an advection–reaction
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
