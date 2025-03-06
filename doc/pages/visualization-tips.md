# Tips and Tricks for Visualization {#visualization-tips}

@tableofcontents

Working with discontinuous Galerkin (DG) finite element methods can be counterintuitive
when it comes to visualizing and interpreting the results.
As the name suggests,
the solution is multivalued and discontinuous at cell edges,
a feature which contrasts the nature of Finite Difference (FD)
or Finite Volume (FV) methods.
These methods by definition only have one value per node/cell.
In this section, we want to provide you some tips, tricks
and lessons we had to learn on the way,
that might help with visualizing and interpreting the output of @sapphire.

## Discontinuous Galerkin Finite Elements: Projection vs. Interpolation {#dg-fe-visualization}

In finite element (FE) methods,
differential equation are solved by finding *weak* solutions to the problem.
This means, that the numeric solution $f_h(x)$
does not need to converge to the full solution $f$ point wise,
but only in the integral sense. We
express this in the following way:

$$
  (f, v) = (f_h, v) \quad \forall v \in V_h \,,
$$

where $v$ is a test function
and $(f, v) = \int f(x) v(x) dx$ is the inner product
of the function space $V_h$.
Using that $V_h$ is spanned by a set of basis functions $\phi_i(x)$,
Galerkin methods use them to project the solution $f_h$
into the function space $V_h$:

$$
  f_h(x) = \sum_i f_i \phi_i(x) \,.
$$

The numeric solution is now represented by the coefficients $f_i$.
We say that $f_h(x)$ is the *projection*
of $f(x)$ onto the FE space.
In the picture below
we show the difference between a *projection* $f_h$ and an *interpolation*.
In the interpolation,
the function has the exact value at the nodes of the mesh,
using a linear interpolation between the nodes.
This is not the case for the projection.
In fact, using a *discontinuous* Galerkin method,
the basis functions are only defined on one cell each.
Therefore, the projection is discontinuous at the cell boundaries.

![Projection vs interpolation](https://sapphirepp.org/img/examples/visualization/visualization_projection_interpolation.png)

Running @sapphire,
we observed that at later time steps $t^n$,
these discontinuities might be more pronounced in the numerical solution $f_h^n$,
than those of the projection at the given time step $f_h(t^n)$.
This becomes especially noticeable,
if the analytic solution would be zero at some time,
$f(t^n, x) = 0$.
An example is the $f_{110}$ mode
in the [convergence-study](#convergence-study) example.

![Discontinuities](https://sapphirepp.org/img/examples/visualization/visualization_discontinuity_k1.gif)

These jumps are a feature of the DG method and are not a sign of instability.
In fact, error estimates for DG methods include estimated for the jump terms,
see @cite Pietro_MathematicalAspectsDG.
A typical example is the $L^2$ error estimate of a DG method is,

$$
  \|u^N - u_h^N\|_{L^2}
  + \left( \sum_m^{N-1} \delta t \sum_F \int_F \frac{1}{2}
    |\beta \cdot n_F| [[u^m - u_h^m]]^2 \right)^{1/2}
  \lesssim C (\Delta x)^a + D (\Delta t)^b \,.
$$

The second term on the left-hand side
provides an estimate of the jumps at the cell boundaries,
$[[u]] = \sum_j (u(x^+_{j+\frac{1}{2}}) - u(x^-_{j+\frac{1}{2}}))$.
$a$, $b$, $C$, and $D$ are constants that bound the error.
Concluding, we can see that the discontinuities are a feature of the DG method,
that we can expect to be on the order of numerical errors.

## Polynomial basis functions and deal.II output {#dealii-visualization}

In @sapphire,
we use polynomial basis functions $\phi_i$ for the FE representation.
More precisely, we Lagrange polynomials of degree $k$ on each cell.
Representing these in a visualization programme,
like ParaView and VisIt,
turn out to be challenging.
Most programmes expect a point wise representation of the solution,
while the FE representation is given by the coefficients $f_i$.
Hence, when saving the results
@dealii "converts" it to point wise data,
by giving the values on each node of the FE mesh.
In the special case of first order polynomials $k=1$,
this results in an exact representation of the numerical solution,
as most visualisation programmes will linearly interpolate between the nodes.
For higher order polynomials,
@dealii provides the following workaround:
it adds additional nodes to the mesh,
representing the polynomials inside the cells.
There are some important points to note:

- Each (cell-edge) node has multiple values,
  one for each cell it is part of.
  (In 1D this corresponds to left- and right-sided values.)
- For higher polynomial degree $k>1$,
  the nodes in the output mesh do not correspond to the mesh used in the simulation,
  as @dealii introduces additional cell-internal nodes.
- The `hdf5` format can only be used for `dim > 1`.
- The `vtu` format allows representing higher order polynomials on a cell,
  if we activate the `write_higher_order_cells` flag.
  We are not using it so far!
- For `pvtu` files, the time variable is not saved correctly.
  Therefore, only the time step can be as a proxy for the simulation time.
- If one wants "smoother" results,
  we recommend using higher order polynomials.
  This is often more effective in reducing discontinuities
  than increasing the mesh resolution.
  See for example the figure below.

![Discontinuities for $k=2$](https://sapphirepp.org/img/examples/visualization/visualization_discontinuity_k2.png)

This figure shows again the $f_{110}$ mode
of the [convergence-study](#convergence-study) example,
but with a higher polynomial degree of $k=2$.
One can see that the discontinuities much less pronounced
than in the $k=1$ case
(figure in previous section).

## VFP modules output {#vfp-visualization}

In the VFP module in @sapphire,
we have to consider some additional points
when visualizing and interpreting the results:

- @sapphire uses an expansion into spherical harmonics.
  It outputs these expansion coefficients,
  $f_{lms}$ labelled `f_lms`.
  It does *not* give the distribution function $f$.
- @sapphire uses *mixed coordinates*.
  The corresponding laboratory frame momentum
  does depend on the position *and* direction of the particles
  (the background velocity is a vector field $\mathbf{u}(\mathbf{x})$).
  Because we expand the momentum dependence into spherical harmonics,
  there is no simple conversion between the absolute value of the momentum in laboratory frame
  and mixed coordinates.
  If you want to compare your result,
  always make sure that both solutions are given in the same coordinate system
  (since there is no easy way to convert @sapphire output to laboratory frame coordinates,
  we recommend to use mixed coordinates for the comparison).

## ParaView {#paraview-visualization}

Here we want to share some tips and tricks
we learned when using @paraview to visualize the results of @sapphire:

- @paraview uses the values given at the nodes
  and linearly interpolates between them inside each cell.
- To extract a 1D line from 2D/3D data,
  use the
  [`PlotOverLine` feature](https://www.paraview.org/paraview-docs/latest/python/paraview.simple.PlotOverLine.html).
  The best way to use this feature
  is to select the `Sample At Cell Boundaries` option
  to get the full discontinuous representation of the data.
  If you prefer a "smooth" result,
  the `Sample At Segment Centers` option is a good choice.
  (For $k=1$, this equals the cell average.)

  Unfortunately,
  the `PlotOverLine` feature has a few quirks we want to mention here:
  - Using the `Sample Uniformly` option
    can lead to sub-sampling and aliasing effects.
    We therefore discourage using this option.
  - In 1D,
    the `Sample At Cell Boundaries` and `Sample At Segment Centers` options
    do not work.
    For the 1D case,
    we recommend showing the original data in a `Line Chart View`.
    After opening a new `Line Chart View`,
    you can display the loaded data
    by clicking to the left of it in the `Pipeline Browser`
    (where the eye icon normally is,
    this is also possible here,
    even though the eye icon is not shown).
    Then deactivate the `Use Index For X Axis` option,
    and instead use `X Array Name = Points_X`.
  - The `Sample At Segment Centers` option includes two additional points,
    one at the start and one at the end of the defined line.
    This is normally not an issue,
    but it means that not all values are sampled at cell centres.
  - In 2D or 3D, make sure to *not* sample along a cell boundary.
    Otherwise, the solution will have 4/8 values at a node
    (one from each neighbouring cell).
    @paraview is not able to interpret this
    and will arbitrarily switch cell sides.
    We recommend using a small offset from the cell boundary.
  - For a large
    @ref sapphirepp::VFP::GridType::shock "shock grids",
    @paraview has some trouble automatically finding the cells
    for `Sample At Cell Boundaries/Segment Centers`.
    In this case,
    deactivating the `Compute Tolerance` option
    and manually setting the `Tolerance` to an appropriate value
    has proven effective.
- We often encountered problems when exporting an animation as `avi` or `ogv`.
  A workaround is to export the animation as a series of `png` files,
  and convert them with an external tool instead.
- As @paraview is in active development,
  you might occasionally encounter bugs
  (e.g., in the `PlotOverLine` feature).
  Sometimes updating to a newer version can fix the issue.
  Otherwise, we recommend checking the
  [ParaView Discourse](https://discourse.paraview.org/)
  and [ParaView Issues](https://gitlab.kitware.com/paraview/paraview/-/issues).

<div class="section_buttons">

| Previous                                         |                  Next |
|:-------------------------------------------------|----------------------:|
| [ParaView Python Introduction](#paraview-python) | [Examples](#examples) |

</div>

---

@author Florian Schulze (<florian.schulze@mpi-hd.mpg.de>)
@date 2024-03-08
