# Synchrotron cooling {#synchrotron-cooling}

@tableofcontents

This example adds synchrotron radiation reaction to the parallel-shock setup used in the @ref parallel-shock and @ref steady-state-parallel-shock examples.
It demonstrates how the Landau–Lifshitz (LL) radiation–reaction term modifies the particle spectrum in a one-dimensional parallel shock, enabling a direct comparison between radiative and non-radiative cases.

## Relation to the parallel-shock examples {#relation-parallel-shock}

The synchrotron-radiation example is built on top of the parallel-shock configuration described in @ref parallel-shock and in the steady-state formulation @ref steady-state-parallel-shock.
It reuses the same spatial shock setup and scattering physics, but adds an energy-loss channel via synchrotron cooling.
This allows one to quantify radiative effects by comparing with the base examples using identical runtime parameters wherever possible.

A detailed derivation and discussion of the weak formulation can be found on the implementation page @ref radiation-reaction.

## VFP equation with synchrotron cooling {#vfp-synchrotron-radiation}

In @sapphire the synchrotron term is added directly to the Vlasov–Fokker–Planck equation.
In compact vector notation the equation solved in this example reads

$$
\frac{\partial f}{\partial t} + (\mathbf{u} + \mathbf{v}) \cdot \nabla_{x} f - \gamma m \frac{\mathrm{D} \mathbf{u}}{\mathrm{D} t} \cdot \nabla_{p} f - \mathbf{p} : (\nabla_{x} \mathbf{u}) \, \nabla_{p} f - \left(\frac{\partial f}{\partial t}\right)_{\mathrm{sync}} + q \, \mathbf{v} \cdot \left( \mathbf{B} \times \nabla_{p} f \right) = \frac{\nu}{2} \, \Delta_{\theta,\varphi} f + S \, .
$$

Here $f(t,\mathbf{x},\mathbf{p})$ is the phase-space distribution, $\mathbf{u}$ is the background plasma velocity, $\mathbf{v}$ the particle velocity, and $\mathbf{B}$ the magnetic field.
The operator $\Delta_{\theta,\varphi}$ denotes the angular part of the Laplacian on the momentum sphere, and $S$ collects any additional source terms.
The synchrotron contribution $(\partial f / \partial t)_{\mathrm{sync}}$ is implemented using the reduced Landau–Lifshitz (LL) radiation–reaction force and is treated as part of the left-hand side of the equation.
This form is equivalent to the more detailed spherical-harmonic representation discussed in @ref radiation-reaction page.

## Implementation overview {#implementation-synchrotron-radiation}

The synchrotron-radiation example lives in the directory
`examples/vfp/synchrotron-radiation`
and consists of the following files:

- `config.h` – defines the geometry, background velocity field, magnetic field, and scattering model used for the parallel shock, as well as enabling synchrotron cooling.
- `synchrotron-radiation.cpp` – driver executable that sets up the VFP solver, reads the parameter file, and advances the solution in time.
- `parameter.prm` – runtime parameter file (in deal.II `ParameterHandler` format) controlling mesh, time-stepping, output, and physical reference units.

The actual synchrotron terms are implemented in the core VFP classes in the library and are activated via the VFP flags passed to the solver.
In particular, the example uses a combination of flags defined in
@ref sapphirepp::VFP::VFPFlags "VFP flags"
to include advection, scattering, and synchrotron cooling in the momentum equation.

Implementing parameters in @sapphire involves a two-step process:

1. **Define the parameters**
   within the @ref sapphirepp::PhysicalParameters "PhysicalParameters" class
   in the `config.h` file.
   This step informs the compiler about their existence
   and sets default values to be used
   if the parameter is not specified in the parameter file.

   @snippet{lineno} examples/vfp/synchrotron-radiation/config.h Define runtime parameter

2. **Add the parameters to the parameter file**
  to ensure that the parameter parser expects them
  and automatically sets their values.
  The @dealii class @dealref{ParameterHandler} `prm`
  and its @dealref{add_parameter,classParameterHandler,a04b75c02037d19fd7fd781785fcefc79} method
  are used for this purpose:

   ```cpp
   prm.add_parameter("entry", parameter, "Description of the parameter");
   ```

   The @ref sapphirepp::PhysicalParameters::declare_parameters() "declare_parameters()" method
   is edited to include the following code:

   @snippet{lineno} examples/vfp/synchrotron-radiation/config.h Declare runtime parameter

## Example parameter file {#example-parameter-synchrotron-radiation}

A typical parameter file for this example is distributed with @sapphire and can be found at
`examples/vfp/synchrotron-radiation/parameter.prm`.
For convenience, we reproduce it here:

@include{lineno} examples/vfp/synchrotron-radiation/parameter.prm

This configuration uses a one-dimensional grid in logarithmic momentum with 64 cells, electron-based reference units, and a fourth-order explicit Runge–Kutta time integrator with a fixed time step.
Output is written to the `./results` directory with base file name `solution` and simulation identifier `synchrotron-test-run-`.
Boundary conditions are chosen to mimic an effectively isolated system in momentum space over the simulated time interval.

## Running the example {#execute-synchrotron-radiation}

After configuring and building @sapphire with examples enabled, the synchrotron-radiation example can be executed from the build directory with

```shell
mpirun -n 1 ./build/examples/vfp/synchrotron-radiation/synchrotron-radiation ./examples/vfp/synchrotron-radiation/parameter.prm
```

The simulation writes a sequence of solution files to the results folder specified in the parameter file.
These can be inspected with @paraview; a simple workflow is described in @ref paraview-tutorial and @ref paraview-python.

## Results and comparison to non-radiative runs {#results-synchrotron-radiation}

Physically, synchrotron cooling introduces an energy-dependent loss term that steepens the high-momentum tail of the shock-accelerated distribution.
In the absence of cooling, the parallel-shock example produces a power-law in momentum over a wide range, limited only by the finite simulation time and momentum grid.

With synchrotron radiation enabled, particles at the highest momenta lose energy on a timescale comparable to, or shorter than, the acceleration time.
The resulting spectrum develops a cooling break and an exponential-like cutoff whose position depends on the magnetic-field strength, the scattering rate, and the total runtime.

By running the synchrotron-radiation example alongside the standard parallel-shock examples with otherwise identical parameters, one can directly compare:

- The evolution of the momentum spectrum at fixed positions relative to the shock.
- The spatial structure of the distribution upstream and downstream of the shock.
- The dependence of the cutoff momentum on magnetic-field strength and runtime.

A typical comparison of the isotropic spectrum at the shock, for simulations with and without synchrotron cooling, is shown below.

![Steady-state spectra at the shock with and without synchrotron cooling](https://sapphirepp.org/img/examples/synchrotron-radiation/synchrotron-spectrum-comparison.png)

These comparisons provide a convenient test of the synchrotron implementation and a starting point for more realistic applications, for example modelling synchrotron-cooled cosmic-ray populations in supernova remnants.

<div class="section_buttons">

| Previous              |
|:----------------------|
| [Examples](#examples) |

</div>

---

@author Harsh Goyal (<harsh.goyal@mpi-hd.mpg.de>)
@date 2025-11-13
