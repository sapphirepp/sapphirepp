# Licence and citation {#licence}

Sapphire++ is distributed under the [LGPL 3.0 license](https://github.com/sapphirepp/sapphirepp/blob/main/LICENSE).

If you use Sapphire++ in your research, please cite the following paper and the software:

> Schween, N. W. and Schulze, F. and Reville, B., Sapphire++: A Particle Transport Code Combining a Spherical Harmonic Expansion and the Discontinuous Galerkin Method, 2025, DOI: https://doi.org/10.1016/j.jcp.2024.113690

```bibtex
@article{Sapphirepp2025,
title = {Sapphire++: A particle transport code combining a spherical harmonic expansion and the discontinuous Galerkin method},
journal = {Journal of Computational Physics},
volume = {523},
pages = {113690},
year = {2025},
issn = {0021-9991},
doi = {https://doi.org/10.1016/j.jcp.2024.113690},
url = {https://www.sciencedirect.com/science/article/pii/S0021999124009380},
author = {Nils W. Schween and Florian Schulze and Brian Reville},
keywords = {Numerical methods, Vlasov-Fokker-Planck, Cosmic rays, Discontinuous Galerkin method, Spherical harmonics, Particle acceleration},
abstract = {We present Sapphire++, an open-source code designed to numerically solve the Vlasov–Fokker–Planck equation for astrophysical applications. Sapphire++ employs a numerical algorithm based on a spherical harmonic expansion of the distribution function, expressing the Vlasov–Fokker–Planck equation as a system of partial differential equations governing the evolution of the expansion coefficients. The code utilises the discontinuous Galerkin method in conjunction with implicit and explicit time stepping methods to compute these coefficients, providing significant flexibility in its choice of spatial and temporal accuracy. We showcase the code's validity using examples. In particular, we simulate the acceleration of test particles at a parallel shock and compare the results to analytical predictions. The Sapphire++ code Image 1 is available as a free and open-source tool for the community.}
}
```

The reference to the Sapphire++ software
can be found on [Zenodo](https://doi.org/10.5281/zenodo.14506754).

Sapphire++ builds on top of the [deal.II](https://www.dealii.org) library.
We therefore strongly encourage you to [cite the deal.II library](https://www.dealii.org/community/publications/) as well.

<div class="section_buttons">

| Previous                              |
|:--------------------------------------|
| [Acknowledgements](#acknowledgements) |

</div>
