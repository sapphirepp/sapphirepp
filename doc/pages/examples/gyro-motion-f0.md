# Gyro motion part 1: Isotropic distribution {#gyro-motion-f0}

@tableofcontents


## Introduction {#introduction-gyro-motion-f0}

In this example, we want to give an overview on how magnetic fields work in
@sapphire. Here, we will focus on only the isotropic part of the distribution,
$f_0$. We have a second example focusing on the anisotropies $f_1$:
@todo Link to gyro-motion-f1 example.

Let us consider a set-up, with particles in a constant external magnetic field,

$$
  \mathbf{B} = B_0 \hat{\mathbf{e}}_z \,.
$$

Since we don't consider scattering, background velocity field nor a source term,
we arrive at the following VFP equation:

$$
  \frac{\partial f}{\partial t} + \mathbf{v} \cdot \nabla_{x} f -
  q \mathbf{v} \cdot \left( \mathbf{B} \times \nabla_{p} f \right) = 0 \,.
$$

It can be seen, that different momenta decouple, so we will have a look at a
discrete momentum $\mathbf{p} = \gamma m \mathbf{v}$.  
Considering a single particle, it will gyrate in the magnetic field with a
gyroradius

$$
  r_g = \frac{p_{\perp} c}{q B} = \frac{\gamma m c v_{\perp}}{q B}
$$

and a gyrofrequency

$$
  \omega_g = \frac{q B}{\gamma m c} \,.
$$

With the $B$-field orientated in $z$-direction, the motion will only take place
in the $x$-$y$-plane. This allows us to reduce the dimensionality of the problem
to 2D. As initial condition, we choose an isotropic distribution of particles
concentrated at the origin:

\begin{align*}
  %Gaussian
  f_0(\mathbf{r}, t=0) &= e^{- \frac{r^2}{2 \sigma^2}} \,, \\
  f_{l>0}(\mathbf{r}, t=0) &= 0 \,.
\end{align*}

Here $\sigma$ is the spread of the initial Gaussian, and $r = \sqrt{x^2 + y^2}$
is the radial distance from the origin. After one gyroperiod,

$$
  T_g = \frac{2 \pi}{\omega_g} = \frac{2 \pi \gamma m c}{q B_0} \,,
$$

every particle performed a full gyro cycle and returned to the initial position.
Therefore, after $n \times T_g$ with $n \in \mathbb{N}$ the particle
distribution returns to its initial condition:

$$
  f_0(\mathbf{r}, t = n \times T_g) = f_0(\mathbf{r}, t = 0)
  = e^{- \frac{r^2}{2 \sigma^2}} \,.
$$

Note that truncating the expansion of the distribution function at $l_{\rm max}$
leads to deviations. Therefore, $f_{l>0}(\mathbf{r}, t = n \times T_g) \neq 0$
even though the initial condition is isotropic.


## Implementation {#implementation-gyro-motion-f0}


### config.h {#config-gyro-motion-f0}

We again start with some includes and open the namespace @ref Sapphire.

@snippet{lineno} examples/gyro-motion-f0/config.h Includes

Next, we define the runtime parameters in the @ref Sapphire::PhysicalProperties
"PhysicalProperties" class. We need the magnetic field strength $B_0$, and the
initial spread of the distribution function $\sigma$. The mass $m$, charge $q$
and velocity/gamma factor $\gamma$ of the particles are already implemented as
runtime parameters in the class @vfpref{VFPParameters}. Notice, that instead
of $B_0$ we use $\frac{B_0}{2 \pi}$ as an input parameter. This is, so that
`Final time = gamma` corresponds to one full gyroperiod $T_g$ (for $m = q = c
=1$).


@snippet{lineno} examples/gyro-motion-f0/config.h PhysicalProperties

The rest of the set-up is again collected in the namespace @ref Sapphire::VFP
"VFP".

@snippet{lineno} examples/gyro-motion-f0/config.h Namespace VFP

We set the dimensionality of the problem to 2D, since the motion is only in the
$x$-$y$-plane.

@snippet{lineno} examples/gyro-motion-f0/config.h Dimension

In this example, we only have the magnetic and spatial advection terms in the
VFP equation. Furthermore, the magnetic field is constant, so we can use the
@vfpref{VFPFlags::time_independent_fields} flag:

@snippet{lineno} examples/gyro-motion-f0/config.h VFP flags

The initial condition is an isotropic Gaussian distribution as described above.
In @sapphire, this is implemented in the @vfpref{InitialValueFunction} function.
We first set the system size to $(l+1)^2$ and the value of $\sigma$ in the
constructor.

@snippet{lineno} examples/gyro-motion-f0/config.h InitialValueFunction constructor

Then we calculate the value for each component $f_{i(l,m,s)}$ in the
@vfpref{InitialValueFunction::vector_value()} function. We first check if the
provided output vector `f` has the correct size using the
@dealref{AssertDimension(),group__Exceptions,ga9442b63275c9ef3fab29bc222831c49c}
macro provided by @dealii. Next, we set all entries in the vector to zero,
before we set only the isotropic component $f_0$ to a Gaussian distribution.

@snippet{lineno} examples/gyro-motion-f0/config.h InitialValueFunction value

The $B$-field is defined in the @vfpref{MagneticField} function. It has three
components, $x$, $y$, $z$, even if the dimension of the problem is only 2D.
Again, the run time parameter `B0` is provided by the @ref
Sapphire::PhysicalProperties "PhysicalProperties" class and set in the
constructor.

@snippet{lineno} examples/gyro-motion-f0/config.h MagneticField constructor

The value of the magnetic field is set in the @vfpref{MagneticField::value()}
function. We again first check if the provided output vector `magnetic_field`
has the correct size, before setting the components according to
$B(\mathbf{r},t) = B_0 \hat{\mathbf{e}}_z$.

@snippet{lineno} examples/gyro-motion-f0/config.h MagneticField value

The last term that @sapphire will use in the computation is the background
velocity field $\mathbf{u}$. Even though the background velocity field is zero
in our case, using the @vfpref{VFPFlags::spatial_advection} flag will activate
the full $(\mathbf{u} + \mathbf{v}) \cdot \nabla_x f$ term. We can however just
set $\mathbf{u} = 0$ explicitly. This is what we do in the
@vfpref{BackgroundVelocityField} function, we set its value $\mathbf{u}$
@vfpref{BackgroundVelocityField::vector_value()}, divergence $\nabla \cdot
\mathbf{u}$ @vfpref{BackgroundVelocityField::divergence_list()}, material
derivative $\frac{\mathrm{D} \mathbf{u}}{\mathrm{D} t}$
@vfpref{BackgroundVelocityField::material_derivative_list()} and Jacobi matrix
$\frac{\partial u_i}{\partial x_j}$
@vfpref{BackgroundVelocityField::jacobian_list()} to zero.

@snippet{lineno} examples/gyro-motion-f0/config.h BackgroundVelocityField

For the scattering frequency and the source term we can provide empty
implementations, since @sapphire will not use them in this example.

@snippet{lineno} examples/gyro-motion-f0/config.h ScatteringFrequency

@snippet{lineno} examples/gyro-motion-f0/config.h Source

Finally, we close the namespaces and terminate the include guard.

@snippet{lineno} examples/gyro-motion-f0/config.h Close namespaces


### gyro-motion-f0.cpp {#main-gyro-motion-f0}

In the `gyro-motion-f0.cpp` we define the main function for the executable. We
start by including the header files for MPI, @dealii and @sapphire.

@snippet{lineno} examples/gyro-motion-f0/gyro-motion-f0.cpp Includes

Next, we define the main function and open a `try-catch` block to abort the
program on exceptions with a meaningful error message.

@snippet{lineno} examples/gyro-motion-f0/gyro-motion-f0.cpp Main function

We initialize MPI using the @dealii wrapper
@dealref{MPI_InitFinalize,classUtilities_1_1MPI_1_1MPI__InitFinalize}
to initialize the MPI and PETSc and p4est. We use one thread per process, so all
parallelisation is done via MPI.

@note It is also possible, to use MPI for the communication between nodes, and
      use TBB to parallelise the computation on a single node.

@snippet{lineno} examples/gyro-motion-f0/gyro-motion-f0.cpp MPI initialization

Next, we create set the @ref Sapphire::Utils:SapphireLogStream "saplog" to only
show the progress of the simulation, and not the debug messages.

@snippet{lineno} examples/gyro-motion-f0/gyro-motion-f0.cpp Saplog

The parameter file to be used is given via a command line argument.

@snippet{lineno} examples/gyro-motion-f0/gyro-motion-f0.cpp Command line argument

As it is @dealii convention, we first need to create all `Parameters` objects
that handle runtime parameters.

@snippet{lineno} examples/gyro-motion-f0/gyro-motion-f0.cpp Run time parameters

We then declare the parameters in the @dealref{ParameterHandler},

@snippet{lineno} examples/gyro-motion-f0/gyro-motion-f0.cpp Declare parameters

and parse them from the parameter file.

@snippet{lineno} examples/gyro-motion-f0/gyro-motion-f0.cpp Parse parameters

Now we can create the @vfpref{VFPEquationSolver} and run the simulation.

@snippet{lineno} examples/gyro-motion-f0/gyro-motion-f0.cpp VFP solver

After the simulation finished, we can compare the result with our expectation.
Assuming the end time of the `Final time` of the simulation corresponds to
multiple of the gyroperiod, we expect the particle distribution $f_0$ to return
to its initial condition. But we only want to compare the $f_0$ component of the
solution with the initial condition, and neglect higher order components. To
this end mask the solution with the @dealref{ComponentSelectFunction}.

@snippet{lineno} examples/gyro-motion-f0/gyro-motion-f0.cpp L2 error

Finally, we end the `try-catch` block and terminate the program.

@snippet{lineno} examples/gyro-motion-f0/gyro-motion-f0.cpp Try-Catch end


## Results {#results-gyro-motion-f0}

@todo Add results section


### Parameter file {#parameter-gyro-motion-f0}

@include{lineno} examples/gyro-motion-f0/parameter.prm


## The plain program {#plain-gyro-motion-f0}


### config.h {#plain-config-gyro-motion-f0}

@include{lineno} examples/gyro-motion-f0/config.h


### gyro-motion-f0.cpp {#plain-main-gyro-motion-f0}

@include{lineno} examples/gyro-motion-f0/gyro-motion-f0.cpp

---

@author Florian Schulze (<florian.schulze@mpi-hd.mpg.de>)
@date 2023-11-09
