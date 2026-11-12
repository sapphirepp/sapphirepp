# Release Notes {#release-notes}

@tableofcontents

## Changes: 1.4.0

- **Breaking change:** The minimum required deal.II version is now 9.7.1.
  The minimum required Doxygen version is now 1.18.0.

### VFP

- Include the Landau--Lifschitz radiation-reaction (RR) force in the VFP equation  
  The RR force can, for example, be used to compute the energy loss
  of relativistic electrons due to their synchrotron emission.
  Its implementation is described [here](https://sapphirepp.org/latest/radiation-reaction.html)
  and an [example](https://sapphirepp.org/latest/steady-state-radiation-reaction-shock.html) is provided.
- Add the ability to checkpoint simulations  
  A simulation can be restarted after a failure.
  An explanation on how to checkpoint and how to resume simulation
  is given in the subsection "Resume simulation"
  on the [parallel-shock example](https://sapphirepp.org/latest/parallel-shock.html) page.
- **Breaking change:** Parameter `Output frequency` was replaced with `Output time step`  
  Instead of specifying the frequency at which simulation results are output,
  the user sets the length of time interval after which output is produced.
- Bug: Weak formulation of the spatial advection part  
  Simulations with `dim_cs > 1` and $\nabla \cdot U \neq 0$ were affected.  
  See [Commit 70484e3](https://github.com/sapphirepp/sapphirepp/commit/70484e3752c63c0fc2a717c9507d30cb883457c4).
- Bug: Weak formulation of the reaction part of the fictitious force term  
  Simulations with `VFPFlags::momentum` activated
  and a Jacobian velocity field with $(\mathbf{J}_{\mathbf{U}})_{ij} \neq =0$ for $i,j > 1$.  
  See [Commit c4ed003](https://github.com/sapphirepp/sapphirepp/commit/c4ed003ce16e819dc90338ad9ecd3a7ce16248cd)
- Improve the `install-dealii.sh` script
- Improve sapplot (link to changelog of sapplot)

**Full Changelog**: https://github.com/sapphirepp/sapphirepp/compare/v1.3.0...v1.4.0

## Changes: 1.3.0

- **Breaking change:** Change C++ standard from 17 to 20.  
  A recompilation of Sapphire++ and deal.II may be required.
- @sapphire can now be installed with Docker  
  Find the installation instructions [here](https://sapphirepp.org/latest/index.html#docker).

### VFP

- Add the local Lax-Friedrichs flux  
  Up to now users were restricted to the upwind flux.
  They can use the more diffusive but faster local Lax-Friedrichs flux
- Change the plotting scripts of the examples to use the [Sapphire++ — Plot library](https://plot.sapphirepp.org)
- Add a grid option that allows users to first create the grid and then refine it  
  This is much faster than creating a refined grid directly.
  In particular, if periodic boundary conditions are used.
- Add a steady-state oblique shock example
- Add examples demonstrating the boundary conditions
- Add VFP parameters that allow users to define the reference units  
  Currently @sapphire uses, for example, the proton mass as a reference mass.
  Now it is possible to change this to the electron mass or any other mass.
  It is also possible to set the reference magnetic field strength, the velocity and the charge.
- Add @ref numerical-parameters that control the solver tolerance and max. number of iterations
- Add convenience functions that simplify the read in of physical functions from files
- Fix: The matrices of the momentum term are now computed with $l_{\mathrm{max}} + 1$

### Utils

- Add a new output parameter `Debug input functions`  
  If set to `true` the velocity field $\mathbf{U}$, its derivatives,
  the magnetic field $\mathbf{B}$, the scattering frequency $\nu$
  and the source $\mathbf{s}$ are included in the ParaView output

**Full Changelog**: https://github.com/sapphirepp/sapphirepp/compare/v1.2.0...v1.3.0

## Changes: 1.2.0

### VFP

- Add the ability to directly solve for steady-state solutions
  - **Breaking change:**
    A new VFP flag called `VFPFlag::time_evolution` was added
    to indicate the inclusion (or exclusion) of the time derivative of $f$.
    Users are asked to adapt their codes: please look at the examples.
- Add the ability to scale the distribution function
  - We added `VFPFlag::scaled_distribution_function` to allow users to compute
    $p^3f$ instead of $f$.
    This is particularly useful if a large $p$-range is needed
    and there are not many particles with high energies.
- Add an example
  [Steady state parallel shock](https://sapphirepp.org/latest/steady-state-parallel-shock.html)
  to demonstrate the above two features
- Add two new boundary conditions and
  improve the documentation of the boundary conditions
  - **Breaking change:**
    Users can now set inflow and reflective boundary conditions.
    All boundary conditions are now documented on an extra
    [webpage](https://sapphirepp.org/latest/boundary-conditions.html).
- The phase space reconstruction of $f$
  now uses (more correctly) $\cos\theta$ instead of $\theta$
  - **Breaking change:**
    Renamed `PhaseSpaceReconstruction` to `ProbeLocation`.
- Improve the parameter handling of user-defined input functions
- Improve the [Parallel shock](https://sapphirepp.org/latest/parallel-shock.html) example

### Utils

- Add plotting routines
  - We use ParaView Python scripts to provide reproducible plotting routines.
    An introduction is given on an extra
    [webpage](https://sapphirepp.org/latest/paraview-python.html).
- Add logfiles
  - The log verbosity and logfile can be specified via command line options.
    See `./build/sapphirepp -h` for all options.
- Simplify `cmake` building routine
  - **Breaking change:**
    Add `cmake` macro `sapphirepp_setup_target`.
- Change folder structure
  - **Breaking change:**
    Examples are now in `examples/vfp/`.
    Parameters files are not copied to `build` folder.
    We recommended to always execute from the `sapphirepp` folder:
    `./build/sapphirepp parameter-template.prm`
- Add function to read in `.csv` files

<div class="section_buttons">

| Previous    |                                  Next |
| :---------- | ------------------------------------: |
| [IDE](#ide) | [Acknowledgements](#acknowledgements) |

</div>
