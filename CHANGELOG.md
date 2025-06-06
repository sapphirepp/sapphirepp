# Release Notes {#release-notes}

@tableofcontents

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
