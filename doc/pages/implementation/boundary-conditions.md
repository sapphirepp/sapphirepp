# Boundary conditions {#boundary-conditions}

@sapphire numerically solves a system of advection-reaction equations
to obtain a solution to the Vlasov--Fokker--Planck equation.
It applies the discontinuous Galerkin method. This means that it computes a weak solution.
Consequently, the prescribed boundary conditions also only hold weakly,
i.e. almost everywhere on the boundary of the computational domain.

For example, the weak formulation of a steady-state advection-reaction equation, i.e.

$$
\begin{split}
	\boldsymbol{\beta}(\mathbf{x}) \cdot \nabla f  + r(\mathbf{x}) f &= s(\mathbf{x}) \, ,\\
	f & = h \text{ on } \partial D^{-} \, ,
\end{split}
$$

where $D^{-} \equiv \{x \in D \mid \boldsymbol{\beta} \cdot \nabla \mathbf{n} < 0 \}$ is the inflow boundary, is

$$
	\text{Find f} \in V \text{ subject to } a(f,w) = \int_D s w \, \mathrm{d}^3 x + \int_{D^-} \left(\boldsymbol{\beta} \cdot \mathbf{n} \right) h w \, \mathrm{d}A
	\quad \text{ for all } w \in V
$$

with

$$
	a(f,w) \equiv \int_D \left( \boldsymbol{\beta} \cdot \nabla f \right) w \, \mathrm{d}^3 x + \int_D r \nabla f w \, \mathrm{d}^3 x + \int_{D^-} \left(\boldsymbol{\beta} \cdot \mathbf{n} \right) f w \, \mathrm{d}A \,.
$$

Note the boundary condition enters the weak formulation, i.e. it enters as an integral,
which implies that it holds almost everywhere,
see equations (2.13) -- (2.16) and eq. (2.21) in @cite Pietro_MathematicalAspectsDG .
Moreover, the boundary condition prescribe the inflow,
namely $(\boldsymbol{\beta} \cdot \mathbf{n} )f$, only.
The outflow is a consequence of the problem itself, e.g. of the flow field $\boldsymbol{\beta}$,
and cannot be prescribed.

@sapphire solves a __system__ of advection-reaction equations
resulting from a spherical harmonic expansion of the distribution function $f$, 
see equations (14) and (15) in @cite Schween2025 .
However, the inclusion of the boundary conditions into the (discrete) weak formulation
can be generalised to the case of a system of partial differential equations (PDE)s.
In that case, the flux (or the flow) through the boundary of the domain is not a scalar but a vector. The corresponding term in the weak formulation of the system of equations is
$$
	\int_{\partial D} \mathbf{J}_{B}(\mathbf{f}, \mathbf{h}) \, \mathrm{d}A \, ,
$$
where $\mathbf{f} = \left(f_{000}, f_{110}, f_{100}, f_{111}, \dots\right)$ is a vector
whose components are the expansion coefficients $f_{lms}$,
and $\mathbf{h}$ are their __user-defined__ values on the boundary.

In the _discrete_ weak formulation, the explicit form of $\mathbf{J}_{B}$ depends on the choice of the numerical flux used. An intuitive choice is the upwind flux, i.e.

$$
	\mathbf{J}_{B}(\mathbf{f}, \mathbf{h}) = \mathbf{W}\left(\boldsymbol{\Lambda}_{+}\mathbf{W}^{T} \mathbf{f} + \boldsymbol{\Lambda}_{-} \mathbf{W}^{T} \mathbf{h}\right) \quad \text{with } (\mathbf{n} \cdot \boldsymbol{\beta}) \mathbf{W} = \mathbf{W} \boldsymbol{\Lambda} \text{ and } \boldsymbol{\Lambda} = \boldsymbol{\Lambda}^{+} + \boldsymbol{\Lambda}^{-} \,,
$$
see eq. (29) in @cite Schween2025 for details.
The salient point is that in the case of a system of equations, there will be inflowing ($\boldsymbol{\Lambda}^{-}$)
_and_ outflowing ($\boldsymbol{\Lambda}^{+}$) components,
depending on the normal $\mathbf{n}$ of the boundary and the advection matrices $\boldsymbol{\beta}$. The outflow is determined by the values of the expansion coefficients $\mathbf{f}$ on the "inner" side of the boundary surface and cannot be prescribed. The inflow is determined by the values of $\mathbf{h}$ on the "outer" side of the boundary surface and can be user defined. Different boundary conditions result from the different choices for $\mathbf{h}$ . For example, if $\mathbf{h} = 0$, we speak of a _zero inflow_ or _outflow_ boundary conditions.

In the next sections, we go through possible choices for $\mathbf{h}$ and illustrate their consequences by means of examples.

## Inflow and zero inflow (outflow) boundary conditions

Inflow boundary conditions results from setting $\mathbf{h}$ independently from $\mathbf{f}$.
For example, we may choose the zeroth component of $\mathbf{h}$ to be a constant in space and time;
setting $h_0 = \sqrt{4 \pi}$ at $x = -L$, 
enforces an inflow that stems from an isotropic distribution $f$ of value $1$ 
at the left $x$ boundary located at $-L$. 
Since the distribution function $f$ is isotropic at the boundary,
we expect that half of the particles, namely all particles with $\theta < \pi/2$, 
enter the computational domain.
As a side remark, we note that representing such a phase-space distribution of particles
requires the inclusion of many spherical harmonics in the expansion of $f$.

To illustrate the inflow produced by an isotropic distribution $f$, 
we use the @ref sapphirepp::VFP::VFPFlags "VFPFlags" to only include the temporal evolution
and the advection term. As stated above, we also set $h_0 = \sqrt{4\pi}$ at $x = -L$ for all $t$.
This results in the following set of equations

$$
	\frac{\partial \mathbf{f}}{\partial t} + v \mathbf{A}_x \partial_x \mathbf{f} = \mathbf{0}
	\quad \text{with } \mathbf{f}^{-}(t, x = -L) = \mathbf{h}^{-} \text{ and } \mathbf{f}_{0}(x) \equiv \mathbf{f}(t = 0, x ) = 0
$$
where $\mathbf{f}^{-}$ denotes the inflowing components of the expansion coefficients. 

To determine the inflow , namely the flux through the boundary at $x = -L$, we integrate
over a small control volume, e.g. $[-L, -L + \epsilon]$, close to the boundary:

$$
	\frac{\partial }{\partial t} \int^{-L + \epsilon}_{-L} \mathbf{f} \, \mathrm{d} x
	+ v \mathbf{A}_x \mathbf{f} \Big |^{-L + \epsilon}_{-L} \,.
$$

At $x = -L$ , we have $- v \mathbf{A}_x \mathbf{f}(t, x = -L)$.
We emphasise that, in the case of a multidimensional problem, 
the minus is consequence of the outward pointing normal $\mathbf{n}$ of the boundary surface.
If we diagonalize $- v \mathbf{A}_x$, i.e. 

$$
	-v \mathbf{A}_x = \mathbf{V}_x \boldsymbol{\Lambda}_x \mathbf{V}^{T}_{x} 
	= \mathbf{V}_x \left(\boldsymbol{\Lambda}^{+}_{x} + \boldsymbol{\Lambda}^{-}_x \right) \mathbf{V}^{T}_x \,,
$$
we can distinguish between inflowing components, 
i.e. the flow is in the opposite direction of $\mathbf{n}$ 
and they correspond to the negative eigenvalues $\boldsymbol{\Lambda}^{-}_{x}$,
and the outflowing component belonging to the positive eigenvalues $\boldsymbol{\Lambda}^{+}_{x}$. 
This distinction may become clearer 
when looking at the following transformation of the above system of PDEs

$$
	\frac{\partial \tilde{\mathbf{f}}}{\partial t}
	- \boldsymbol{\Lambda}_x \frac{\partial \tilde{\mathbf{f}}}{\partial x} = \mathbf{0} \,,
$$

where $\tilde{\mathbf{f}} \equiv \mathbf{V}^{T}_{x} \mathbf{f}$ are the _characteristic variables_. Note that $\mathbf{V}^{-1} = \mathbf{V}^{T}$, because $\mathbf{A}^{T}_{x} = \mathbf{A}_x$.
The analytical solution of the above equations is 

$$
	\tilde{f}_i(x,t) = \tilde{f}_{0,i}(x + \lambda_{i} t) \, .
$$

$\mathbf{f}_{0}(x)$ are the initial conditions (or the boundary conditions) 
of the expansion coefficients.
This solution implies that negative eigenvalues, $\lambda_{i} < 0$, 
move the transformed initial (or boundary) conditions, $\tilde{f}_{0,i}$, 
into the positive $x$-direction, 
i.e. into the opposite direction of the outward normal of the boundary surface at $x = -L$.

We now formalise the concept of inflowing components $\mathbf{f}^{-}$ with the following definition

$$
	\mathbf{f}^{-} = \mathbf{V}_{x}\boldsymbol{\mathbb{1}}^{-} \mathbf{V}^{T}_x \mathbf{f} \, ,
$$

where $\boldsymbol{\mathbb{1}}^{-}$ is a matrix with ones on the diagonal
where the diagonal elements of $\boldsymbol{\Lambda}^{-}_{x}$ are non-zero. 

We wrap up and use the boundary condition $h_0 = \sqrt{4 \pi}$ at $x = -L$ in @sapphire 
and compare it to the analytical solution. 
We use an expansion order $l_{\mathrm{max}} = 9$ and set the velocity $v = \sqrt{8}/3$.
The result is shown in the following animation:

<img src="https://sapphirepp.org/img/implementation/boundary-conditions/inflow-bc.gif" alt="drawing" width="90%"/>

The dashed lines show the derived analytical solution. 
Because we chose $f$ to be isotropic at the boundary, 
particles with velocity component $v_x = v \cos\theta \in [0, v]$ enter the computational domain.
The ones with velocities close to $v$ cross the domain fastest.
This can be seen during the first time steps.
Though, the steady inflow of particles will eventually lead to a state
in which particles cover the complete velocity range $[0, v]$.
This state is steady and reached at the end of the animation.

During the animation's first moments, the plotted expansion coefficients are step functions. 
This is an effect of a finite expansion order.
If we solved the equation for $f$ directly, i.e. 

$$
	\frac{\partial f}{\partial t} + v \cos\theta \frac{\partial f}{\partial x} \quad 
	\text{with } f(0,x) = 0 \text{ and } f(t,-L) = 1 \,,
$$

we would get 

$$
	f(t,x,\theta) =
	\begin{cases}
		1 &\text{for } x < -L + v \cos\theta t \\
		0 &\text{for } x \geq -L + v \cos\theta t
	\end{cases}\; .
$$

This implies that $f_{000}$ of an infinite order expansion is 

$$
	f_{000}(t,x) = \int Y_{000} f(t,x,\theta) \mathrm{d}\Omega
	= 
	\begin{cases}
	\sqrt{\pi}\left(1 - \frac{x + L}{vt}\right) &\text{for } x < -L + vt \\
	0 &\text{for } x \geq -L + vt
	\end{cases}	\; ,
$$
where we used that $f = 1$ when $\cos\theta > (x + L)/vt$.
We also plotted this function; it's label is `f_000_no_expansion`.
It is a linear function; that is approximated by a step function and the higher the expansion order,
the better the approximation. 
This means that we would need an infinite angular resolution to capture it exactly. 

In applications, the last statement is crucial.
If there is no scattering, the expansion order must be very high.
A fact that becomes clear when we reconstruct phase space at some point $x$ 
after steady state has been reached:

The plot compares the phase space reconstruction for $l_{\mathrm{max}} = 9$ 
and $l_{\mathrm{max}} = 63$. 

<p float="left">
  <img src="https://sapphirepp.org/img/implementation/boundary-conditions/l_max9.png" width="48%" />
  <img src="https://sapphirepp.org/img/implementation/boundary-conditions/l_max63.png" width="48%" />
</p>

A low expansion order leads to a distribution function with negative values. 
Moreover, if there is no scattering 
a difference between odd and even $l_{\mathrm{max}}$ can be observed,
for details see Sec. 2 in  @cite Garret2016.

To conclude we look at the same example, 
but we include scattering, we directly compute the steady-state solution and use the zero inflow boundary condition. 
We solve

$$
	 v \mathbf{A}_x \partial_x \mathbf{f} = \nu \mathbf{C} \quad \text{with } \mathbf{f}^{-}(x = -L) = \mathbf{h}^{-} \text{ and } \mathbf{f}^{-}(x = L) = \mathbf{0}  \,.
$$

For $l_{\mathrm{max}} = 1$, the explicit system of equations is 

$$
	\frac{v}{\sqrt{3}} 
	\begin{pmatrix}
	0 & 0 & 1 & 0 \\
	0 & 0 & 0 & 0 \\ 
	1 & 0 & 0 & 0 \\
	0 & 0 & 0 & 0 
	\end{pmatrix} \partial_{x} \mathbf{f}
	= -\nu 
\begin{pmatrix}
	0 & 0 & 0 & 0 \\
	0 & 1 & 0 & 0 \\ 
	0 & 0 & 1 & 0 \\
	0 & 0 & 0 & 1 
   	\end{pmatrix} \mathbf{f} \, .
$$
Hence, 

$$
\begin{split}
f_{000} &= -\frac{\sqrt{3}\nu}{v} f_{100} x + \text{const.} \\
f_{110} &= \text{const.} \\ 
f_{100} &=\text{const.} \\
f_{111}  &= \text{const.} \\ 
\end{split}
$$

For the boundary conditions, we again set $h_{0} = \sqrt{4 \pi}$ and use the eigenvectors $\mathbf{V}_{x}$ to see that

$$
\mathbf{f}^{-}(x = -L) = \frac{1}{2}
\begin{pmatrix}
	f_{000} + f_{100}\\
	0 \\
	f_{000} + f_{100} \\
	0 
\end{pmatrix}
= \sqrt{\pi}
\begin{pmatrix}
1 \\
0 \\
1 \\
0 \\
\end{pmatrix}
\quad
\text{and}
\quad
\mathbf{f}^{-}(x = L) = \frac{1}{2}
\begin{pmatrix}
	f_{000} - f_{100}\\
	0 \\
	-f_{000} + f_{100} \\
	0 
\end{pmatrix}
= 
\begin{pmatrix}
0 \\
0 \\
0 \\
0 \\
\end{pmatrix}
$$

This can be condensed in the two boundary conditions $ f_{000} + f_{100} = 2\sqrt{\pi}$ at $x = -L$ and $f_{000} = f_{100}$ at $x = L$. 

## Reflecting boundary conditions

This implies that implementing reflexive boundary conditions, requires us to find a vector $\mathbf{g}$ that "reflects" $\mathbf{f}_{h}$. For the sake of simplicity, we restrict our attention to the boundaries at $x = 0$ and at $x = L $, where $L$ is the size of our computational domain in $x$-direction. A reflection at $y-z$-plane is a transformation that flips the sign of the $x$-component of, say, a vector $\mathbf{w} \in \mathbb{R}^{3}$. Formally,

$$
\mathbf{w}'=
\begin{pmatrix}
	-1 & 0 & 0 \\
	0  & 1 & 0 \\
    0  & 0 & 1
\end{pmatrix} \mathbf{w} \,.
$$

To apply this transformation to our spherical harmonic expansion of the distribution function $f$, we need quantities that transform like vectors. These are the moments of $f$. For example,

$$
\begin{split}
  j^{i}(\mathbf{x}) &= \int f v^{i} \, \mathrm{d}^{3} p \\
  T^{ij}(\mathbf{x}) & = \int f v^{i}v^{j} \, \mathrm{d}^{3} p \, .
\end{split}
$$

The moments of $f$ can be computed straightforwardly, $f$ is expanded in the Cartesian tensors, i.e.
$$
  f(t, \mathbf{x}, \mathbf{p}) = \sum^{\infty}_{l_{\mathrm{max}}} F^{(l)}_{i_{1} \cdots i_{l}} \frac{p^{i_{1}} \cdots p^{i_{l}}}{p^{l}} \,.
$$

The second moment is then a sum of the zeroth order expansion coefficient and the second order expansion coefficient. The third moment is a sum of the first order expansion coefficient and the third order expansion coefficient. In mathematics

$$
\begin{split}
	n         & = c_{1} \int F^{(0)} p^{2} \mathrm{d}^{3}p  \\
        j^{i} & = c_{2} \int F^{(1)}_{i} p \mathrm{d}^{3} p \\
    T^{ij}    & = c_{3} \int F^{(0)} p^{2} \mathrm{d}^{3}p  + \int F^{(2)}_{ij} h(p) \mathrm{d}^{3} p
\end{split}
$$

```math
\alpha + \beta = \gamma
```

The Cartesian tensor expansion coefficients can be expressed in spherical harmonics. For example,

$$
\begin{split}
  j^{x} &= c f_{100} \\
  j^{y} &= -c f_{110} \\
  j^{z} &= -c f_{111}
\end{split}
$$

```math
\begin{split}
 \alpha &= \gamma \\
 \beta &= \delta
\end{split}
```
or for the second order expansion coefficient, $F_{yz} = \frac{1}{4} \sqrt{\frac{15}{\pi}} f_{221}$, $F_{xy} = -\frac{1}{4} \sqrt{\frac{15}{\pi}} f_{210}$, $F_{xz} = -\frac{1}{4} \sqrt{\frac{15}{\pi}} f_{211}$, $F_{zz} = -\frac{1}{4} \sqrt{\frac{15}{\pi}} f_{220} - \frac{1}{4} \sqrt{\frac{5}{\pi}} f_{200}$ and $F_{xx} = \frac{1}{2} \sqrt{\frac{5}{\pi}} f_{000}$.

We can now apply our transformation to this tensors, to get our vector $\mathbf{g}$.
For example, $\mathbf{j}'$ gets a minus sign in its $x$-component. For a vector $\mathbf{g}$, we can take the value $\mathbf{f}_{h,1}$ on the inner side of our computational domain, and flip the sign of the $f_{100}$ component. Repeating this procedure with the second moment, we see that the sign of the $F_{xx}$ component is unchanged, because the reflection was applied twice. However, the sign of $F_{xy}$ and $F_{xz}$ changes. Hence, we have to flip the sign of the component $f_{210}$ and $f_{211}$. This can be done with higher order moments as well. There will always be a sign flip whenever the is an odd number of $x$s, e.g. $F_{xxx}$, will flip the signs of the spherical harmonic expansion coefficients that contain it.

In the end, we can work out which components of $\mathbf{f}_{h,1}$ have to change their signs. We formalise this with a signature $\mathbf{g} = \mathbf{\eta} \mathbf{f}_{h,1}$.

Using the above procedure to determine this signature matrix is complicated. Luckily there is an easier, if we take a look at momentum space. In this space we use spherical coordinates, and we reflect  our velocity vector $\mathbf{v}$, by changing $\theta \rightarrow \pi - \theta$. Now keeping in mind that all the moments are computed with the velocities in it,
we can apply our transformation, before integrating them out.
(This needs a clearer argument.) We can directly get the signature by using working out the sing changes corresponding to $Y^{m}_{l}(\pi - \theta, \varphi)$.

The reflection at the other boundaries should work analogously.

## Continuous boundary conditions
