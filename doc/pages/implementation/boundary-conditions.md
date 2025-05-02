# Boundary conditions {#boundary-conditions}

@sapphire numerically solves a system of advection-reaction equations
to obtain a solution to the Vlasov--Fokker--Planck equation.
It applies the discontinuous Galerkin method. This means that it computes a weak solution.
Consequently, the prescribed boundary conditions also only hold weakly,
i.e. almost everywhere on the boundary of the computational domain.

For example, the weak formulation of an advection-reaction equations, i.e.

$$
\begin{split}
	\frac{\partial f}{\partial t} + \boldsymbol{\beta}(\mathbf{x}) \cdot \nabla f  + r(\mathbf{x}) f &= s(\mathbf{x}) \, ,\\
	f & = h \text{ on } \partial D^{-}
\end{split}
$$

where $D^{-} \equiv \{x \in D \mid \boldsymbol{\beta} \cdot \nabla \mathbf{n} < 0 \}$ is the inflow boundary,

$$
	\text{Find f} \in V \text{ subject to } a(f,w) = \int_D s w \, \mathrm{d}^3 x + \int_{D^-} \left(\boldsymbol{\beta} \cdot \mathbf{n} \right) h w
	\text{ for all } w \in V
$$

with

$$
a(f,w) \equiv \int_D \left( \boldsymbol{\beta} \cdot \nabla f \right) w \, \mathrm{d}^3 x + \int_D r \nabla f w \, \mathrm{d}^3 x + \int_{D^-} \left(\boldsymbol{\beta} \cdot \mathbf{n} \right) f w \, \mathrm{d}A \,.
$$

Note the boundary condition enters the weak formulation, i.e. it enters as an integral,
which implies that it holds almost everywhere,
see equations (2.13) -- (2.16) and eq. (2.21) in @cite Pietro_MathematicalAspectsDG .
Moreover, the boundary condition prescribe the inflow,
i.e. $(\boldsymbol{\beta} \cdot \mathbf{n} f)$, only.
The outflow is a consequence of the problem itself, e.g. of the flow field $\boldsymbol{\beta}$,
and cannot be prescribed.

@sapphire solves a __system__ of advection-reaction equations, see equations (14) and (15) in
@cite Schween2025 .
However, the inclusion of the boundary conditions into the (discrete) weak formulation
can be generalised to the system of partial differential equations (PDE)s.

In the case of the system, the flux (or the flow) through the boundary of the
domain is not a scalar but a vector.  There will be Additionally, @sapphire computes an
approximation to the solution of the advection-reaction equations.  Hence, the
flux is a numerical flux $\mathbf{J}_{B}$. At boundary cell interfaces the numerical approximation to the solution
has two different values on each side of the interface, say $\mathbf{f}_{h,1}$
and $\mathbf{f}_{h,2}$. The numerical flux is a function of both of them,
i.e. $\mathbf{J}_{F}(\mathbf{f}_{h,1}, \mathbf{f}_{h,2})$.

At the boundary faces of our computational domain, the value for the approximate solution outside the domain is unknown. Choosing it determines the kind of boundary condition. For example, setting it to zero constitutes the "zero inflow" (or "outflow") boundary condition. In formulas

$$
  \mathbf{J}_{B}(\mathbf{f}_{h,1}, \mathbf{g}) \equiv \mathbf{J}_{F}(\mathbf{f}_{h,1}, \mathbf{g})
$$
with $\mathbf{g} = 0$.

## Inflow and zero inflow boundary conditions

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
