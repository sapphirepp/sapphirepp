# Visualization {#visualization}

@tableofcontents

@sapphire uses standardized output formats like `vtu`, `pvtu`, and `hdf5` to
save results. We recommend using dedicated tools like @paraview or
[VisIt](https://visit-dav.github.io/visit-website/) to visualize these results,
as they handle much of the heavy lifting involved in visualizing Finite Element
simulation data. We provide a [tutorial for ParaView](#paraview-tutorial) to
help new users become familiar with the basic concepts. If you are already
familiar with @paraview or [VisIt](https://visit-dav.github.io/visit-website/),
you can directly jump to the [tips and tricks for
visualization](#visualization-tips). Here, we explain some caveats of results
produced by the discontinuous Galerkin method utilized by @sapphire.

@subpage paraview-tutorial

@subpage visualization-tips

<div class="section_buttons">

| Previous                    |                                    Next |
|:----------------------------|----------------------------------------:|
| [Quick-start](#quick-start) | [ParaView Tutorial](#paraview-tutorial) |

</div>
