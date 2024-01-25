# Visualization {#visualization}

@tableofcontents

@todo So far these are only notes. They need to be reworked.

Here are some tips and tricks for visualizing the output of @sapphire. We also
give some advice how to interpret the output.


## Finite Elements {#fe-visualization}

FE always does projections. Therefore, the FE presentation does not converge
point wise, but only in the integral sense. This means that the FE solution can
mismatch the analytical solution at every point, but still be an accurate
presentation in the FE sense.

The FE representation has polynomials of degree $k$ on each cell. But most
visualisation programme expect a point wise representation of the solutions.
Hence, when saving the results we "convert" it to point wise data, by giving the
values on each node of the FE mesh. If we have a polynomial of degree $k>1$ we
add additional nodes to the mesh, representing the polynomials inside the cells.
@note The nodes in the output mesh therefore do not correspond to the mesh used
      in the simulation, if $k>1$.

We use a DG method. This means that the solution is **discontinues** between
the cells. Even the projection of a continuous function will be discontinues.
Furthermore, each node is shared by multiple cells and therefore has multiple
values.

- `hdf5` format can only be used for `dim > 1`
- `vtu` format allows representing higher order polynomials on a cell, if we
  activate the `write_higher_order_cells` flag. We are not using it so far!


## VFP modules specialities {#vfp-visualization}

- @sapphire uses an expansion into spherical harmonics. It output these
  expansion coefficients, $f_{lms}$ labelled `f_lms`. It does *not* the give the
  distribution function $f$.
- @sapphire uses *mixed coordinates*. The corresponding laboratory frame
  momentum does depend on the position *and* direction of the particles (the
  background velocity is a vector field $\mathbf{u}(\mathbf{x}$). Because we
  expand the momentum dependence into spherical harmonics, there is no simple
  conversion between the absolute value of the momentum in laboratory frame and
  mixed coordinates. If you want to compare your result, always make sure that
  both solutions are given in the same coordinate system (since there is no easy
  way to convert @sapphire output to laboratory frame coordinates, we recommend
  to use mixed coordinates for the comparison).


## ParaView {#paraview-visualization}

- It uses the values given at the nodes and linearly interpolates between them
  inside each cell.
- ParaView has problems reading the mesh data (especially in the 1D case).
  Therefore, one should *not* use the options `Sample At Cell Boundaries` or
  `Sample At Segment Centers` when using the `PlotOverLine` feature. Instead,
  one should use the option `Sample Uniformly`. It is recommended to choose a
  `Resolution` equal to the number of cells in the mesh, to sample one point per
  cell and get a smooth representation of the solution without jumps of the DG
  method.
- When using the `PlotOverLine` feature in 2D or 3D, one should make sure to
  *not* sample along a cell boundary. Otherwise, the solution will have 4/8 or
  more values at a node (one from each neighbouring cell). ParaView is not able
  to interpret this and will show a wrong representation of the solution. We
  recommend using a small offset from the cell boundary when using this feature.
  Combined with using the `Sample Uniformly` option, this will ensure the
  solution is only sampled inside a cell, where the solution is well-defined.


<div class="section_buttons">

| Previous                    |                  Next |
|:----------------------------|----------------------:|
| [Quick-start](#quick-start) | [Examples](#examples) |

</div>


---

@author Florian Schulze (<florian.schulze@mpi-hd.mpg.de>)
@date 2024-01-25
