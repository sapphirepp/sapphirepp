# Advanced example: Convergence study {#convergence-study}


## Analytical solution

These are only notes for internal use, we should remove them later.

We derive an analytic solution for the static `gyro` setup with a truncation at
$l=1$ using a Fourier transform (FT).

The VFP equation is given by:

$$
  \frac{\partial f}{\partial t} + \mathbf{v} \cdot \nabla_{x} f -
  q \mathbf{v} \cdot \left( \mathbf{B} \times \nabla_{p} f \right) = 0 \,.
$$

$$
  \mathbf{B} = B_0 \hat{\mathbf{e}}_z \,.
$$

We only work in 1D.

We have two different notations, the @sapphire notation and Brian's.

The expansion of the distribution function is given by:

$$
 f(t, \mathbf{x}, \mathbf{p}) = \sum^{l_{\rm max}}_{l = 0} \sum^{l}_{m = 0}
 \sum^{1}_{s = 0} f_{lms}(t, \mathbf{x}, p) Y_{lms}(\theta,\varphi) \,,
$$

We get a system of equations:

$$
  \partial_{t}\pmb{f} +
  \left(
    U^{a}\mathbf{1} +
    V \mathbf{A}^{a}
  \right) \partial_{x_{a}}\pmb{f} -
  \left(
    \gamma m \frac{\mathrm{d} U_{a}}{\mathrm{d} t} \mathbf{A}^{a} -
    p \frac{\partial U_{b}}{\partial x^{a}} \mathbf{A}^{a}\mathbf{A}^{b}
  \right) \partial_{p} \pmb{f} +
  \left(
    \frac{1}{V} \epsilon_{abc} \frac{\mathrm{d} U^{a}}{\mathrm{d} t}
      \mathbf{A}^{b}\mathbf{\Omega}^{c} +
    \epsilon_{bcd} \frac{\partial U_{b}}{\partial x^{a}}\mathbf{A}^{a}
      \mathbf{A}^{c}\mathbf{\Omega}^{d}
  \right) \pmb{f} -
  \omega_{a}\mathbf{\Omega}^{a}\pmb{f} +
  \nu \mathbf{C}\pmb{f}
  = 0
$$

In our case and with $l_{\rm max}=1$ this simplifies to:

$$
  \partial_{t}\pmb{f} +
  V \mathbf{A}^{a} \partial_{x_{a}}\pmb{f} -
  \omega_{a}\mathbf{\Omega}^{a}\pmb{f}
  = 0
$$

Where $\omega_{a} = \omega \delta_{a,z}$, $\omega = \frac{B_0 q}{\gamma m}$.
In 1D we have $\partial_{x}\pmb{f} = \partial_{z}\pmb{f} = 0$. Therefore:

$$
  \partial_{t}\pmb{f} +
  V \mathbf{A}^{x} \partial_{x}\pmb{f} -
  \omega \mathbf{\Omega}^{z}\pmb{f}
  = 0
$$

Putting in the components of $\mathbf{A}$ and $\mathbf{\Omega}$, we get:

$$
  \partial_{t} f_{000} +
  \frac{V}{\sqrt{3}} \partial_{x} f_{100}
  = 0
$$

$$
  \partial_{t} f_{110} -
  \omega f_{100}
  = 0
$$

$$
  \partial_{t} f_{100} +
  \frac{V}{\sqrt{3}} \partial_{x} f_{000} +
  \omega f_{110}
  = 0
$$

$$
  \partial_{t} f_{111} = 0
$$

Here, we used
$A^{x}_{20}=A^{x}_{02}=\sqrt{\frac{(l+m)(l-m)}{(2l+1)(2l-1)}}=1/\sqrt{3}$,
and $\Omega^{z}_{12} = -\Omega^{z}_{21} = 1$ with all other components
vanishing.
Be aware that the ordering of components in the vector is
$\pmb{f} = (f_{000}, f_{110}, f_{100}, f_{111})$.

$f_{111} = \rm const.$ decouples from the system of equations, and we choose the
initial conditions so that $f_{111} = 0$.

Acting with $\partial_{t}$ on the third equation and inserting the other
equations, we arrive at:

$$
  \partial_{t}^2 f_{100} -
  \frac{V^2}{3} \partial_{x}^2 f_{100} +
  \omega^2 f_{100}
  = 0
$$

We apply the ansatz of factorization for $f_{110}$:

$$
  f_{100}(t, x) = \phi(t) \psi(x)
$$

Inserting this into the equation above and denoting
$\partial_{t} \phi(t) = \dot{\phi}$ and $\partial_{x} \psi(x) = \psi'(x)$,
we get:

$$
  \ddot{\phi(t)} \psi(x) -
  \frac{V^2}{3} \phi(t) \psi''(x) +
  \omega^2 \phi(t) \psi(x)
  = 0
$$

$$
  \implies
  \frac{\ddot{\phi(t)} + \omega^2 \phi(t)}{\phi(t)}
  = \frac{V^2}{3} \frac{\psi''(x)}{\psi(x)}
  = - c_n^2 = \rm const.
$$

These equations are respectively solved by:

$$
  \phi(t) = \phi_{0,n} \sin(\sqrt{\omega^2 + c_n^2} t) +
    \phi_{1,n} \cos(\sqrt{\omega^2 + c_n^2} t)
$$

We use initial conditions, such that $f_{100}(t=0) = 0 \implies \phi_{1,n} = 0$.

$$
  \psi(x) = \psi_{0,n} \sin\left(\frac{\sqrt{3} c_n}{V} x\right) +
    \psi_{1,n} \cos\left(\frac{\sqrt{3} c_n}{V} x\right)
$$

The full solution is given by a superposition of these Fourier modes:

$$
  f_{100}(t, x) = \sum_{n} \phi_{0,n} \sin(\sqrt{\omega^2 + c_n^2} t)
    \left[ \psi_{0,n} \sin(k_n x) + \psi_{1,n} \cos(k_n x) \right]
$$

where we introduced the wave number $k_n = \frac{\sqrt{3} c_n}{V}$.

One can notice, that $\phi_{0,n}$, $\psi_{0,n}$ and $\psi_{1,n}$ are degenerate
once. In the following, we choose $\phi_{0,n} = 1$.

Enforcing periodic boundary conditions in a box with length $L$, we get:

$$
  k_n L = \frac{\sqrt{3} c_n}{V} L \overset{!}{=} 2 \pi n
  \implies c_n = \frac{V}{\sqrt{3} L} 2 \pi n \,, k_n = \frac{2 \pi n}{L}
$$

with $n \in \mathbb{N}$.

Inserting this solution in the equation for $f_{000}$, we get:

$$
  f_{000}(t,x) - f_{000}(t=0,x)
  = -\int_{0}^{t} \frac{V}{\sqrt{3}} \partial_{x} f_{100}(\tilde{t}, x)
    \mathrm{d}\tilde{t}
  = \sum_{n} \frac{V}{\sqrt{3}} \frac{1}{\sqrt{\omega^2 + c_n^2}}
    \left[ \cos(\sqrt{\omega^2 + c_n^2} t) - 1 \right] \frac{\sqrt{3} c_n}{V}
    \left[ \psi_{0,n} \cos(k_n x) - \psi_{1,n} \sin(k_n x) \right]
$$

We choose the following initial condition for $f_{000}$ (wisely choosing the
prefactor):

$$
  f_{000}(t=0, x) = \sum_{n} \sqrt{\omega^2 + c_n^2}
    \left[ A_n \cos(k_n x) - B_n \sin(k_n x) \right]
$$

Using this gives us an initial condition for the derivative of $f_{100}$:

$$
  \partial_{t} f_{100} \rvert_{t=0} +
  \frac{V}{\sqrt{3}} \partial_{x} f_{000} \rvert_{t=0} +
  \omega f_{110} \rvert_{t=0}
  = 0
$$

Using $f_{110}(t=0,x) = 0$, this results in:

$$
  \sum_{n} \sqrt{\omega^2 + c_n^2}
    \left[ \psi_{0,n} \sin(k_n x) + \psi_{1,n} \cos(k_n x) \right]
  = -\frac{V}{\sqrt{3}} \sum_{n} \sqrt{\omega^2 + c_n^2} k_n
    \left[ - A_n \sin(k_n x) - B_n \cos(k_n x) \right]
$$

Comparing the coefficients we get:

$$
  \psi_{n,0} = \frac{V}{\sqrt{3}} k_n A_n = c_n A_n
$$

$$
  \psi_{n,1} = c_n B_n
$$

We can get $f_{110}$ by integrating $f_{100}$ in time. Using the initial
condition $f_{110}(t=0) = 0$, we get:

$$
  f_{110}(t,x) = \int_{0}^{t} \omega f_{100}(\tilde{t}, x) \mathrm{d}\tilde{t}
  = \sum_{n} \frac{\omega c_n}{\sqrt{\omega^2 + c_n^2}}
    \left[ 1 - \cos(\sqrt{\omega^2 + c_n^2} t) \right]
    \left[ A_n \sin(k_n x) + B_n \cos(k_n x) \right]
$$

To summarize, we have:

$$
  f_{000}(t,x) = \sum_{n} \frac{c_n^2}{\sqrt{\omega^2 + c_n^2}}
    \left[
      \cos(\sqrt{\omega^2 + c_n^2} t) - 1 +
      \frac{\omega^2 + c_n^2}{c_n^2}
    \right]
    \left[ A_n \cos(k_n x) - B_n \sin(k_n x) \right]
$$

$$
  f_{110}(t,x) = \sum_{n} \frac{\omega c_n}{\sqrt{\omega^2 + c_n^2}}
    \left[ 1 - \cos(\sqrt{\omega^2 + c_n^2} t) \right]
    \left[ A_n \sin(k_n x) + B_n \cos(k_n x) \right]
$$

$$
  f_{100}(t, x) = \sum_{n} c_n \sin(\sqrt{\omega^2 + c_n^2} t)
    \left[ A_n \sin(k_n x) + B_n \cos(k_n x) \right]
$$

$$
  f_{111}(t, x) = 0
$$


---

Last, we check again, that this solution is correct by inserting it into the
third equation of the system of equations:

$$
  \partial_{t} f_{100} +
  \frac{V}{\sqrt{3}} \partial_{x} f_{000} +
  \omega f_{110}
  = 0
$$

Inserting all functions we arrive at:

$$
  \sum_{n} c_n \sqrt{\omega^2 + c_n^2} \cos(\sqrt{\omega^2 + c_n^2} t)
    \left[ A_n \sin(k_n x) + B_n \cos(k_n x) \right] +
  \sum_{n} \frac{V}{\sqrt{3}}\frac{c_n^2}{\sqrt{\omega^2 + c_n^2}}
    \left[
      \cos(\sqrt{\omega^2 + c_n^2} t) - 1 +
      \frac{\omega^2 + c_n^2}{c_n^2}
    \right]
    k_n \left[ -A_n \sin(k_n x) - B_n \cos(k_n x) \right] +
  \sum_{n} \frac{\omega^2 c_n}{\sqrt{\omega^2 + c_n^2}}
    \left[ 1 - \cos(\sqrt{\omega^2 + c_n^2} t) \right]
    \left[ A_n \sin(k_n x) + B_n \cos(k_n x) \right]
  = 0
$$

Comparing the coefficients we have:

$$
  c_n \sqrt{\omega^2 + c_n^2} \cos(\sqrt{\omega^2 + c_n^2} t) -
  \frac{c_n^3}{\sqrt{\omega^2 + c_n^2}}
    \left[
      \cos(\sqrt{\omega^2 + c_n^2} t) - 1 +
      \frac{\omega^2 + c_n^2}{c_n^2}
    \right] +
  \frac{\omega^2 c_n}{\sqrt{\omega^2 + c_n^2}}
    \left[ 1 - \cos(\sqrt{\omega^2 + c_n^2} t) \right]
  = 0
$$

$$
  \implies
  \left[
    \omega^2 + c_n^2 - c_n^2 - \omega^2
  \right] \cos(\sqrt{\omega^2 + c_n^2} t) +
  c_n^2 - (\omega^2 + c_n^2) + \omega^2
  = 0 \quad \rm (correct)
$$
