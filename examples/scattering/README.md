# Introduction and Basic Usage: Scattering {#scattering}

Table of contents:
\tableofcontents
\todo Table of contents is not working


## Introduction {#introduction}

In this example, we want to give an introduction to `Sapphire++` and show the
basic usage.
All examples are structured in the following way:
First, we give an introduction, showing which parts of the VFP equation we want
to solve, giving a physical motivation and show how `Sapphire++` translates this
problem to finite elements.
Next, we show how to implement the example in `Sapphire++`. We start with a
detailed description of the `config.h`, going through the file line by line giving
some explanatory comments. Similarly, we present the `scattering.cpp`
which implements the `main` function.
In the results section, we show the expected output of the program and explain
how to interpret it.

As an introductory example, we want to study the simple case with only a
scattering term in the VFP equation,
\todo Check sign in front of scattering operator
$$
  \frac{\partial f}{\partial t} =  \frac{\nu}{2} \Delta_{\theta, \varphi} f \,.
$$
Here $f(\mathbf{x}, \mathbf{p}, t)$ is the distribution function and $\nu$ is
the scattering frequency. The Laplacian $\Delta_{\theta, \varphi}$ only acts on
the directional part of the momentum $\mathbf{p} = (p \cos\theta, p \cos\varphi
\sin\theta, p \sin\varphi \cos\theta)^{\top}$.

`Sapphire++` uses an expansion real spherical harmonics, $Y_{lms}(\theta,\varphi)$,
to solve the VFP equation. The distribution function can therefore be written as
$$
  f(\mathbf{x}, \mathbf{p}) = \sum^{\infty}_{l = 0} \sum^{l}_{m = 0} \sum^{1}_{s
  = 0} f_{lms}(\mathbf{x}, p ,t) Y_{lms}(\theta,\varphi) \,.
$$
This expansion uses, that the distribution function is nearly isotropic, with
anisotropies described by the higher moments $f_{l>0}$. The VFP equation now
turns into a system of equations for the expansion coefficients
$f_{lms}(\mathbf{x}, p, t)$. Notice, that $f_{lms}$ depends only on the absolute
value of the momentum $p$. This reduces the dimensionality of the problem from
$(\mathrm{dim\_cs})^2$ to $\mathrm{dim\_ps} = \mathrm{dim\_cs}+1$. Here, we
introduce call the nomenclature configuration space $\mathrm{cs}$ for the
dimension of $\mathbf{x}$ and (reduced) phase space $\mathrm{ps}$ for
$(\mathbf{x}, p)$.

Collecting all expansion coefficients into a vector $\boldsymbol{f} = (f_{000},
f_{110}, f_{100}, f_{111}, f_{220}, f_{210}, f_{200}, f_{211}, f_{221} \dots)^
{\top}$ the VFP equation now becomes
$$
  \partial_{t} \boldsymbol{f} + \nu \, \mathbf{C} \boldsymbol{f} = 0 \,.
$$
The matrix $\mathbf{C}$ represents the scattering operator and is given by
\todo Insert citation to DG paper
$$
  (\mathbf{C})_{i(l'm's')j(l,m,s)} = \frac{l(l + 1)}{2}\delta_{l'l}\delta_{m'm}
  \delta_{s's} \,.
$$

Focusing on the equation for a single expansion coefficient, it becomes obvious
that in this example the different expansion coefficients are independent,
$$
  \partial_{t} f_{lms} + \nu \, \frac{l(l + 1)}{2} f_{lms} = 0 \,.
$$
This differential equation describes and exponential decay,
$$
  f_{lms}(t) = f_{lms, 0} \exp\left(-\nu \frac{l(l + 1)}{2} t \right) \,,
$$
where $f_{lms, 0} = f_{lms}(t=0)$. We will use this solution to check the
correctness of the numerical solution.


## Implementation {#code}

To implement the example in `Sapphire++`, we need to create two.
The first one is the header file `config.h` which implements compile time
flags and the physical setup. The later one consist of the definition of the
initial condition $f_{lms, 0}$ "InitialValueFunction", the scattering
frequency $\nu$ "ScatteringFrequency", the source $S$ "Source",
the magnetic field $\mathbf{B}$ "MagneticField" and the background
velocity field $\mathbf{u}$ "BackgroundVelocityField".
The second file is a `.cpp` file that implements the main function. The
since the physical setup is defined in the `config.h` file, the `main.cpp`
will be nearly identical for most use-cases of `Sapphire++`.
The last file that has to be created is the `Parameter.prm` which defines
the run time parameter of Sapphire.
