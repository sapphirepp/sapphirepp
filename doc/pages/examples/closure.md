# Closure {#closure}

@tableofcontents

## Introduction {#introduction-closure}

This example is designed to give a feeling
for the effect of truncating the expansion at $l_{\rm max}$.
Furthermore, we will introduce the post-processing module
to reconstruct the full distribution function
$f(\mathbf{x}, p, \theta, \varphi)$.

@todo Motivation and introduction

## Implementation {#implementation-closure}

The implementation of this example in the `examples/closure` directory.
In the following we discuss details,
assuming the reader already familiarised himself with @sapphire.

### VFP equation {#dimension-closure}

The example is one dimensional in space,
and uses mono-energetic particles $\mathbf{p} = \gamma m \mathbf{v}$.
The reduces phase space is therefore only 1 dimensional,
`dim = dim_ps = dim_cs = 1`.

@snippet{lineno} examples/closure/config.h Dimension

In this simplified example, we have the following VFP equation,

$$
  \frac{{\rm d} f}{{\rm d} t} =
  \frac{\partial f}{\partial t} + \mathbf{v} \cdot \nabla_{x} f =
  \frac{\nu}{2} \Delta_{\theta, \varphi} f \,.
$$

This results in the following terms we have to activate in @sapphire,
the @ref sapphirepp::VFP::VFPFlags::spatial_advection "spatial advection"
($(\mathbf{u} + \mathbf{v}) \cdot \nabla_{x} f$),
and @ref sapphirepp::VFP::VFPFlags::collision "collisions"
($\frac{\nu}{2} \Delta_{\theta, \varphi} f$).
Furthermore, we can activate the
@ref sapphirepp::VFP::VFPFlags::time_independent_fields "time independent fields".

@snippet{lineno} examples/closure/config.h VFP Flags

### Runtime Parameters {#parameter-closure}

Our setup is characterized by only two parameters,
the standard deviation $\simga$ of the initial distribution,
and the scattering frequency $\nu$.
As demonstrated in the [parallel shock](#parallel-shock) example,
we will define these as runtime parameters.

1. **Define**

   We start by defining our runtime parameters:

   @snippet{lineno} examples/closure/config.h Define runtime parameter

2. **Declare**
  
   Next, we declare the parameters in the @dealref{ParameterHandler}:

   @snippet{lineno} examples/closure/config.h Declare runtime parameter

3. **Parse**

   Finally, we parse the runtime parameters:

   @snippet{lineno} examples/closure/config.h Parse runtime parameter

### Initial Condition {#initial-condition-closure}

As previously described is the initial condition given by
an isotropic Gaussian distribution.
Consequently, we set all expansion coefficients $f_{l>0}$ to zero,
except for the isotropic component $f_{i(0,0,0)} = f_0$.

@snippet{lineno} examples/closure/config.h Initial value

### Scattering frequency {#scattering-frequency-closure}

For the scattering frequency we simply use the constant
given as runtime parameter.

@snippet{lineno} examples/closure/config.h Scattering frequency

### Magnetic field {#magnetic-field-closure}

@todo Remove magnetic field?

@snippet{lineno} examples/closure/config.h Magnetic field

### Phase space reconstruction {#phase-space-reconstruction-closure}

To perform the phase-space reconstruction,
we can use the post-processor module provided in
@ref sapphirepp::VFP::PhaseSpaceReconstruction "PhaseSpaceReconstruction".
To activate it,
we just have to give a list of points to the parameter file.

@todo Example

## Results {#results-closure}

@todo Results

### Example parameter file {#example-parameter-closure}

The full parameter file to reproduce the results is given below.

@include{lineno} examples/closure/parameter.prm

<div class="section_buttons">

| Previous              |
|:----------------------|
| [Examples](#examples) |

</div>

---

@author Florian Schulze (<florian.schulze@mpi-hd.mpg.de>)
@date 2024-12-14
