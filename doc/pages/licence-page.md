# Licence {#licence}

Sapphire++ is distributed under the [LGPL 3.0 license](https://github.com/sapphirepp/sapphirepp/blob/main/LICENSE).

If you use this software in your research, please cite the following paper:

> Schween, N. W. and Schulze, F. and Reville, B., Sapphire++: A Particle Transport Code Combining a Spherical Harmonic Expansion and the Discontinuous Galerkin Method, 2025, DOI: https://doi.org/10.1016/j.jcp.2024.113690

Here's the BibTeX entry for the paper:

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

Sapphire++ builds on top of the [deal.II](https://www.dealii.org) library.
We therefore encourage you to [cite the deal.II library](https://www.dealii.org/community/publications/) as well:

> Pasquale Claudio Africa, Daniel Arndt, Wolfgang Bangerth, Bruno Blais, Marc Fehling, Rene Gassmöller, Timo Heister, Luca Heltai, Martin Kronbichler, Matthias Maier, Peter Munch, Magdalena Schreter-Fleischhacker, Jan Philipp Thiele, Bruno Turcksin, David Wells, Vladimir Yushutin  
> The deal.II Library, Version 9.6  
> Journal of Numerical Mathematics, vol. 32, pp. 369-380, 2024.  
> DOI: 10.1515/jnma-2024-0137; preprint; bibtex

```bibtex
@Article{2024:africa.arndt.ea:deal,
  author  = {Pasquale C. Africa and Daniel Arndt and Wolfgang Bangerth and Bruno Blais and
             Marc Fehling and Rene Gassm{\"o}ller and Timo Heister and Luca Heltai and
             Sebastian Kinnewig and Martin Kronbichler and Matthias Maier and Peter Munch and
             Magdalena Schreter-Fleischhacker and Jan P. Thiele and Bruno Turcksin and
             David Wells and Vladimir Yushutin},
  title   = {The deal.II library, Version 9.6},
  journal = {Journal of Numerical Mathematics},
  year    = 2024,
  volume  = 32,
  number  = 4,
  pages   = {369--380},
  doi     = {10.1515/jnma-2024-0137}
}
```

> D. Arndt, W. Bangerth, D. Davydov, T. Heister, L. Heltai, M. Kronbichler, M. Maier, J.-P. Pelteret, B. Turcksin, D. Wells  
> The deal.II finite element library: design, features, and insights  
> Computers & Mathematics with Applications, vol. 81, pages 407-422, 2021.  
> DOI: 10.1016/j.camwa.2020.02.022; preprint

```bibtex
@Article{dealii2019design,
  title   = {The {deal.II} finite element library: Design, features, and insights},
  author  = {Daniel Arndt and Wolfgang Bangerth and Denis Davydov and
             Timo Heister and Luca Heltai and Martin Kronbichler and
             Matthias Maier and Jean-Paul Pelteret and Bruno Turcksin and
             David Wells},
  journal = {Computers \& Mathematics with Applications},
  year    = {2021},
  DOI     = {10.1016/j.camwa.2020.02.022},
  pages   = {407-422},
  volume  = {81},
  issn    = {0898-1221},
  url     = {https://arxiv.org/abs/1910.13247}
}
```

<div class="section_buttons">

| Previous                              |
|:--------------------------------------|
| [Acknowledgements](#acknowledgements) |

</div>
