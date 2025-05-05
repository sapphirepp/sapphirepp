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
resulting from a spherical harmonic expansion of the distribution function $f$, s
ee equations (14) and (15) in @cite Schween2025 .
However, the inclusion of the boundary conditions into the (discrete) weak formulation
can be generalised to the case of a system of partial differential equations (PDE)s.
In that case, the flux (or the flow) through the boundary of the domain is not a scalar but a vector.The corresponding term in the weak formulation of the system of equations is
$$
	\int_{\partial D} \mathbf{J}_{B}(\mathbf{f}, \mathbf{h}) \, \mathrm{d}A \, ,
$$
where $\mathbf{f} = \left(f_{000}, f_{110}, f_{100}, f_{111}, \dots\right)$ is a vector
whose components are the expansion coefficients $f_{lms}$,
and $\mathbf{h}$ are their __user-defined__ values on the boundary.

In the _discrete_ weak formulation, the explicit form of $\mathbf{J}_{B}$ depends on the choice of the numerical flux used. An intuitive choice is the upwind flux, i.e.

$$
	\mathbf{J}_{B}(\mathbf{f}, \mathbf{h}) = \mathbf{W}\left(\boldsymbol{\Lambda}_{+}\mathbf{W}^{T} \mathbf{f} + \boldsymbol{\Lambda}_{-} \mathbf{W}^{T} \mathbf{h}\right) \quad \text{with } (\mathbf{n} \cdot \boldsymbol{\beta}) \mathbf{W} = \mathbf{W} \boldsymbol{\Lambda} \text{ and } \boldsymbol{\Lambda} = \boldsymbol{\Lambda}_{+} + \boldsymbol{\Lambda}_{-} \,,
$$
see eq. (29) in @cite Schween2025 for details.
The salient point is that in the case of a system of equations, there will be inflowing ($\boldsymbol{\Lambda}_{-}$)
_and_ outflowing ($\boldsymbol{\Lambda}_{+}$) components,
depending on the normal $\mathbf{n}$ of the boundary and the advection matrices $\boldsymbol{\beta}$. The outflow is determined by the values of the expansion coefficients $\mathbf{f}$ on the "inner" side of the boundary surface. The inflow is determined by the values of $\mathbf{h}$ on the "outer" side of the boundary surface. Different boundary conditions result from the different choices for $\mathbf{h}$ . For example, if $\mathbf{h} = 0$, we speak of a _zero inflow_ or _outflow_ boundary conditions.

In next sections, we go through possible choices for $\mathbf{h}$ and illustrate their consequences by means of examples.

## Inflow and zero inflow (outflow) boundary conditions

Inflow boundary conditions results from setting $\mathbf{h}$ independently from $\mathbf{f}$. For example, we may choose the zeroth component of $\mathbf{h}$ to be a constant; setting $h_0 = \sqrt{4 \pi}$ at $x = -L$, enforces an isotropic distribution of value $1$ at the left $x$ boundary located at $-L$. If we use the @ref sapphirepp::VFP::VFPFlags "VFPFlags" to only include the temporal evolution and the advection term, this boundary condition results in the following Riemann problem

$$
	\frac{\partial \mathbf{f}}{\partial t} + v \mathbf{A}_x \partial_x \mathbf{f} = 0
	\quad \text{with } f_{0} =
	\begin{cases}
		\sqrt{4 \pi} &\text{for } x \leq -L \\
		 0 &\text{for } x > -L
	\end{cases}
$$

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
