<p align="center">
<img
src="https://git.mpi-hd.mpg.de/schween/sapphire/media/branch/main/logo/sapphire-logo-text.png"
width=563 height=300 > 
</p>

# About

Sapphire is an acronym and stands for "**S**imulating **a**strophysical
**p**lasmas and **p**articles with **hi**gly **r**elativistic **e**nergies". It
is made to simulate the interaction of particles with a background plasma. To
this end it solves a Vlasov-Fokker-Planck equation in mixed coordinates, namely

$$
  \frac{\partial f}{\partial t} + (\mathbf{u} + \mathbf{v}) \cdot \nabla_{x}
  f - \gamma m \frac{\mathrm{D} \mathbf{u}}{\mathrm{D} t} \cdot \nabla_{p}f
  - \mathbf{p} \cdot\nabla_{x} \mathbf{u}\cdot \nabla_{p} f
  + q \mathbf{v} \cdot \left( \mathbf{B} \times \nabla_{p} f \right) =
  \frac{\nu}{2} \Delta_{\theta, \varphi} f\, .
$$

$f$ is a number density of particles in phase space. The above equation models how this
distribution evolves in time assuming that the particles interact with the a
background plasma, which is determined by means of its velocity $\mathbf{u}$ and
its magnetic field $\mathbf{B}$. The interaction is determined by the right-hand
side of the equation: It models particles whose direction of diffuses.

To solve the above equation, we expand $f$ using real spherical harmonics, namely
$$
 f(\mathbf{x}, \mathbf{p}) = \sum^{\infty}_{l = 0} \sum^{l}_{m = 0} \sum^{1}_{s
 = 0} f_{lms}(\mathbf{x}, p ,t) Y_{lms}(\theta,\varphi) \,.
$$

This turns the original equation into a system of equations determining the
coefficients $f_{lms}$ of the expansions. And this is the output of the
simulation.

The system of equations is solved applying a finite element method. We use the
discontinuous Galerkin method. Sapphire is implemented with the help of the
finite element library (deal.ii)[https://dealii.org/]. 

# Compiling Sapphire

# Quick start guide
