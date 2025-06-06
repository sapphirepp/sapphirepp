# Steady-state parallel shock {#steady-state-parallel-shock}

@tableofcontents

This example builds on the [parallel shock](#parallel-shock) example.

## Introduction {#introduction-steady-state-parallel-shock}

In the [parallel shock](#parallel-shock) example, we showed how to simulate the
_time-dependent_ acceleration and propagation of particles at a parallel shock
front. It could be seen that after some time a steady state was reached, i.e.
the phase-space distribution of the particles did not change anymore with time.
The example, that you are currently reading, demonstrates how steady-state
solutions of the VFP equation can be obtained directly. Though various physical
scenarios lead to a steady state, we once more opted to compute the particle
distribution for the case of a parallel shock.

As in the [quick start](#quick-start) and in the [parallel
shock](#parallel-shock) example, we examine a collisionless shock with
compression ratio $r$ that propagates through an astrophysical plasma,
e.g. the interstellar medium. Particles are injected at a constant rate $Q$ with
an energy determined by $p_{\text{inj}}$ at the position of the shock wave. They
are scattered in the up- and downstream at a frequency $\nu$ and, thus, undergo
diffusive shock acceleration. Since a detailed explanation of the described
scenario is given in the two mentioned examples, we summarise the physical
setup in the following table:
<div style="text-align:center">
| Physical scenario:   |                                                                                                                                                                                       |
|----------------------|---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| Velocity profile     | $$\mathbf{U} = U_{\text{sh}} \hat{\mathbf{e}}_x  \text{ for } x < 0 \text{ and } \mathbf{U} = U_{\text{sh}}/r \hat{\mathbf{e}}_x \text{ for } x > 0 $$                                |
| Magnetic field       | $$\mathbf{B}(x) = B_0 \hat{\mathbf{e}}_x$$                                                                                                                                            |
| Scattering frequency | $$\nu(x) = \nu_0 B_0 p^{-1} $$                                                                                                                                                        |
| Source               | $$ S(x,p) = \frac{Q}{4\pi p^2} \frac{1}{2\pi \sigma_x \sigma_p}\exp\left[-\left(\frac{(x - x_{\text{inj}})^2}{2 \sigma^2_x} + \frac{(p -p_\text{inj})^2}{2\sigma^2_p}\right)\right]$$ |
</div>

## Implementation {#implementation-steady-state-parallel-shock}

The implementation of the sketched physical scenario can be found in the
directory `sapphirepp/examples/vfp/steady-state-parallel-shock`. Since the
implementation is essentially the same as in the [parallel
shock](#parallel-shock) example, we restrain our discussion to the changes
necessary to directly compute the steady-state solution with @sapphire.

### VFP equation {#dimension-steady-state-parallel-shock}

Not very surprising, the main difference to the time-dependent parallel shock
example is the exclusion of the time derivative of the distribution function in
the VFP equation. We, namely, solve

$$
  (\mathbf{u} + \mathbf{v}) \cdot \nabla_{x} f -
  \gamma m \frac{\mathrm{D} \mathbf{u}}{\mathrm{D} t} \cdot \nabla_{p}f -
  \mathbf{p} \cdot\nabla_{x} \mathbf{u}\cdot \nabla_{p} f +
  q \mathbf{v} \cdot \left( \mathbf{B} \times \nabla_{p} f \right) =
  \frac{\nu}{2} \Delta_{\theta, \varphi} f + S \,.
$$

Therefore, the list of `vfp_flags` is

@snippet{lineno} examples/vfp/steady-state-parallel-shock/config.h VFP Flags

We note that we did _not_ include the
@ref sapphirepp::VFP::VFPFlags::time_evolution "time_evolution" flag. Also, it
was not necessary to explicitly state that the velocity field, the magnetic
field and the source term are time-independent.

### Scattering frequency {#scattering-frequency-steady-state-parallel-shock}

The steady-state parallel shock scenario also differs from the time-dependent
case in the momentum dependence of the scattering frequency. We set

$$
  \nu(p) = \nu_0 B_0 p^{-1} \, ,
$$

where $\nu_0$ is a parameter which can be freely set. Because $h = \omega_g/\nu$, where $\omega_g = q B/\gamma m$ is the gyro frequency, is the hall parameter, $1/nu_0$  gives approximately the number of gyrations a particle performs before being scattered.
$\nu_0 = 1$ is called the "Bohm limit" or "Bohm scaling".
We note that a plasma in which particles are scattered more often is not magnetised anymore.
Furthermore, for the above choice of $\nu(p)$ the Hall parameter becomes independent of $p$.

The implementation looks like

@snippet{lineno} examples/vfp/steady-state-parallel-shock/config.h Scattering frequency

Note that we chose to work with $\ln p$ instead of $p$ and, hence,
$\nu(\ln p) = \nu_0 * B_0 \exp(-1 \ln p)$.
We remember the user that the last component of `point` in reduced phase space
is the magnitude of the momentum variable.
Since the reduced phase space is $(x, \ln p)$, `points[q_index][1]` is $\ln p$.

Notice that $\nu_0$ is a runtime parameter that is set in the supplied parameter
file.

### Source term {#source-term-steady-state-parallel-shock}

We also slightly changed the source term. We, namely, included the factor
$1/4\pi p^2$. Thus we still inject isotropic and almost mono-energetic particles
at the shock. However, the additional factor gives the user-defined rate $Q$ a
clearer physical interpretation. The reason is that the _configuration space_
number density of injected particles is given by

$$
  n_{\text{inj}} = \int S \, \mathrm{d}^3 p = 4\pi \int S p^2 \, \mathrm{d} p = \frac{Q}{\sqrt{2
  \pi} \sigma_x} \exp\left(-\frac{(x - x_{\text{inj}})^2}{2 \sigma^2_x}\right) \,,
$$
which is now directly set by the parameter $Q$.

Keeping in mind that we actually solve a system of PDEs that determines the
expansion coefficients of the spherical harmonic expansion of the distribution
function $f$, we need to also expand the source term in spherical harmonics.
This is a straightforward computation, because the source term is
independent of $\theta$ and $\varphi$, i.e. we inject an isotropic particle
distribution,

$$
  S_{000} = \int Y_{000} S \, \mathrm{d}\Omega = \sqrt{4\pi} S \,.
$$

All $S_{lms}$ with $l > 0$ are zero.

@snippet{lineno} examples/vfp/steady-state-parallel-shock/config.h Source

We again note that we use the index `i` to refer to the components of the
spherical harmonic decomposition of the source term. $i = 0$ corresponds $l = 0
, m = 0$ and $s = 0$.

### Compile and run {#compile-steady-state-parallel-shock}

If you configured @sapphire with the `-DEXAMPLES=ON` option,
see [Compilation](#compilation),
the steady-state shock example is compiled
and the executable can be found in the
`build/examples/vfp/steady-state-parallel-shock` folder.
The executable is called `steady-state-parallel-shock`

We supply the user with a physically reasonable set of parameters that are
listed in the `examples/vfp/steady-state-parallel-shock/parameter.prm` file:

@include examples/vfp/steady-state-parallel-shock/parameter.prm

Run the simulation with:

```shell
mpirun -n 6 ./build/examples/vfp/steady-state-parallel-shock/steady-state-parallel-shock examples/vfp/steady-state-parallel-shock/parameter.prm
```

The plots can be created with the command:

```shell
pvbatch examples/vfp/steady-state-parallel-shock/pvplot.py results/steady-state-parallel-shock
```

## Results {#results-steady-state-parallel-shock}

In @cite Drury1983 eq. 3.24, Drury derives an analytic solution for the
isotropic part $f_{000}$ of the distribution function at the shock, i.e. at $x = 0$.
He assumes a constant scattering frequency. However, in the case of a parallel shock the $p$-dependence of the $\nu$ does not change the solution. Moreover, his source term is a
delta distribution in $x$ and $p$ (point injection) and his velocity profile has
a sharp discontinuity at the shock. In contrast to Drury, we model the velocity
discontinuity with a tanh-function and the delta distribution with a narrow
Gaussian profile. We, thus, do not expect that the analytic solution in @cite Drury1983
and our simulation results match exactly, but they should get closer
to each other the smaller the shock width and the narrower the injection
profile. Adapting his solution to our source term, we expect the simulation
result to be approximately

$$
 f_{000}(x = 0, p) = \frac{3}{\sqrt{4 \pi} p^{3}_{\text{inj}}} \frac{Q}{U_1 - U_2}
 \left(\frac{p}{p_{\text{inj}}}\right)^{-3 U_1/(U_1 - U_2)}
 = \frac{3Q}{\sqrt{4 \pi} U_1 p^{3}_{\text{inj}}} \frac{r}{r - 1} \left(\frac{p}{p_{\text{inj}}}\right)^{-3 r/(r - 1)}
$$

For the additional factor of $\sqrt{4\pi}$ see @cite Schween2025 eq. 61. Because
we set the compression ratio to $r=4$, we expect a $p^{-4}$ power law.

Drury also gives a formula for the spatial dependence of the isotropic part
of the distribution function $f_{000}$, see eq. 2.34 in @cite Drury1983, namely

\begin{equation}
  f_{000}(x, p = \hat{p}) =
  \begin{cases}
  f_{000}(x = 0, \hat{p})  \exp\left(3 U_{1} \nu(\hat{p})/\hat{v}^2 x \right) &\text{for } x < 0 \\
  f_{000}(x = 0, \hat{p}) &\text{for } x > 0
  \end{cases} \,
\end{equation}
where $\hat{p}$ is an arbitrary but specified value of the magnitude of the
momentum and $\hat{v}$ the corresponding velocity.

The following two plots show a comparison of the simulation results with the
expected analytical solution:

<div style="text-align:center;">
  <div style="display:inline-block;">
    <img alt="Particle spectrum"
  src="https://sapphirepp.org/img/examples/steady-state-parallel-shock/particle-spectrum.png"
   height=450>
  </div>
<div style="display:inline-block;">
    <img alt="Downstream particle distribution"
  src="https://sapphirepp.org/img/examples/steady-state-parallel-shock/spatial-distribution.png"
    height=450>
  </div>
</div>

Despite the differences in the simulation setup and the derivation of the
analytical solution, the analytical prediction and the simulation results agree
well.

## Discussion {#discussion-steady-state-parallel-shock}

Since we use the dG method in combination with a spherical harmonic expansion of
the distribution function $f$, we need to solve a large system of equations to
obtain the steady state solution. We apply an iterative method, namely the
[generalized minimal residual
method](https://en.wikipedia.org/wiki/Generalized_minimal_residual_method)
(GMRES) in combination with a Block-Jacobi preconditioner. We have not
implemented this algorithm, but we use the @dealii interface to the
[PETSc](https://petsc.org/release/) library.

This has multiple consequences for users of the steady-state solver. Firstly, an
iterative method requires termination criteria. These are the number of
iterations, the magnitude of the relative (or absolute) preconditioned residual,
i.e. $\| \mathbf{B}(\mathbf{b} - \mathbf{A}\mathbf{x})\|$ and the relative
increase in the residual, see [PETSc
documentation](https://petsc.org/main/manual/ksp/#convergence-tests) for more
details. The number of iterations and the relative tolerance are currently
hard-coded, i.e. $5000$ and $1 \times 10^{-8} * \|\mathbf{b}\|_2$ respectively.
This implies that if a user computes a distribution functions whose values are smaller than
the tolerance, the solver will not convergence. This happens, for example, if
the $p$-range covers many orders of magnitudes. This issue can be addressed with
a scaled distribution function, see Sec.
[Extenstion](#scaled-f-steady-state-parallel-shock).

Secondly, the Block-Jacobi preconditioner is not particularly designed for the
system of equations that we solve. This means that we need many iterations.
Moreover, the number of used blocks in the preconditioner is determined by the
number of processes, i.e. MPI ranks, used to solve the system of equations. Its
quality, thus, is probably better for a low number of ranks. However, a single
compute node does not have infinite memory, and it might be necessary to use
many cores to actually solve a problem.

@note If you, by chance, have experience with preconditioners for advection-reaction
  systems (or Friedrichs' systems), please contact me. We are interested in
  developing a more robust preconditioner

If a user is familiar with PETSc, she can try to change the iterative method and
the used preconditioner via the command line using PETSc commands. For example,
each Jacobi block is solved with an
[incomplete LU
factorization](https://en.wikipedia.org/wiki/Incomplete_LU_factorization) (ILU),
if a complete LU factorization of the blocks is wanted, the program can be run
with

```shell
mpirun -n 4 ./build/examples/vfp/steady-state-parallel-shock/steady-state-parallel-shock examples/vfp/steady-state-parallel-shock/parameter.prm -ksp_monitor -sub_pc_type lu -sub_ksp_type preonly
```

In general, the option `-ksp_monitor` allows a user to monitor the convergence of
the iterative method. More information about the iterative method and the
preconditioner can be obtained with the option `-ksp_view`.

@note Unfortunately, it is not possible to overwrite the tolerance for the
  residual `-ksp_rtol`. The reason is that @dealii uses its own routine, namely
  `SolverBase::convergence_test()`, to check for convergence.

## Extension: A scaled distribution function {#scaled-f-steady-state-parallel-shock}

As explained at the beginning of the [Disucssion
section](#discussion-steady-state-parallel-shock), it might be that the value of
the distribution function $f$ is below the termination tolerance of the
iterative method, i.e. $f < 1 \times 10^{-8} * \|\mathbf{b}\|$.
If this is the case, it is possible to scale the distribution function with the factor $p^{s}$. We
decided to call $s$ the `scaling_spectral_index` and it is set to three, i.e. $s = 3$.
Moreover, we define $g = p^3 f$.
To scale the distribution function it is necessary to add an additional `VFPFlag`. We
included a subfolder `scaled` in the `steady-state-parallel-shock` example
folder.

It is only the `config.h` file that changes. The additional VFPFlag
is called `scaled_distribution_function`. Hence the variable `vfp_flags` becomes

@snippet{lineno} examples/vfp/steady-state-parallel-shock/scaled/config.h VFP Flags

If we scale the distribution function, we also have to scale the source $S$,
i.e. we have to multiply it with $p^{3}$. This gives

@snippet{lineno} examples/vfp/steady-state-parallel-shock/scaled/config.h Source

The steady-state shock example with a scaled distribution function can be run
with

```shell
mpirun -n 6 build/examples/vfp/steady-state-parallel-shock/scaled/steady-state-parallel-shock-scaled examples/vfp/steady-state-parallel-shock/scaled/parameter.prm &&
pvbatch examples/vfp/steady-state-parallel-shock/pvplot.py results/steady-state-parallel-shock
```

@warning It seems that scaling the distribution function increases the condition
number of the system of equations that we need to solve. This implies that the
GMRES method may **not** converge or take very long to converge. It might be necessary to increase the number of allowed iterations.

A result of the above simulation run is depicted in the following plot:

<div style="text-align:center;">
<img alt="Scaled particle spectrum at the shock"
src="https://sapphirepp.org/img/examples/steady-state-parallel-shock/scaled-particle-spectrum.png"
height=450>
</div>

The plot shows the scaled particle spectrum $g(0,p)$ multiplied with $p$ at the shock.
The slight deviation of the numerical solution's spectral index is a consequence of the fact that we
cannot use a discontinuous velocity profile. Instead we use a tanh-function and,
thus, have finite shock width. For a detailed discussion of the effect of a
finite shock width on the spectral index, we refer our users to @cite Achterberg2011.
The results in the cited paper were reproduced with @sapphire in @cite Dann2025.

<div class="section_buttons">

| Previous              |
|:----------------------|
| [Examples](#examples) |

</div>

---

@author Nils Schween (<nils.schween@mpi-hd.mpg.de>)
@date 2025-02-05
