[![GitHub CI](https://github.com/sapphirepp/sapphirepp/actions/workflows/tests.yml/badge.svg)](https://github.com/sapphirepp/sapphirepp/actions/workflows/tests.yml)
[![clang-format Check](https://github.com/sapphirepp/sapphirepp/actions/workflows/clang-format-check.yml/badge.svg)](https://github.com/sapphirepp/sapphirepp/actions/workflows/clang-format-check.yml)
[![Monthly tests](https://github.com/sapphirepp/sapphirepp/actions/workflows/monthly-tests.yml/badge.svg)](https://github.com/sapphirepp/sapphirepp/actions/workflows/monthly-tests.yml)

<p align="center">
<img src="logo/sapphire-logo-text.png" width=563 height=300>
</p>


# About Sapphire++

Sapphire++ (<strong>S</strong>imulating
<strong>a</strong>strophysical <strong>p</strong>lasmas and
<strong>p</strong>articles with <strong>hi</strong>ghly
<strong>r</strong>elativistic <strong>e</strong>nergies in
C<strong>++</strong>) is a code to solve that solves the Vlasov-Fokker-Planck
equation for particle transport in plasma:

$$
  \frac{\partial f}{\partial t} + (\mathbf{u} + \mathbf{v}) \cdot \nabla_{x} f -
  \gamma m \frac{\mathrm{D} \mathbf{u}}{\mathrm{D} t} \cdot \nabla_{p}f -
  \mathbf{p} \cdot\nabla_{x} \mathbf{u}\cdot \nabla_{p} f +
  q \mathbf{v} \cdot \left( \mathbf{B} \times \nabla_{p} f \right) =
  \frac{\nu}{2} \Delta_{\theta, \varphi} f + S(\mathbf{x}, \mathbf{p}, t)\, .
$$


## Installation

Sapphire++ builds on top of the [deal.II](https://www.dealii.org) finite element
library. To install [deal.II](https://www.dealii.org) and other prerequisites
for Sapphire++ we provide an [installation script](scripts/install-dealii.sh).
Run the script by executing:

```shell
chmod u+x install-dealii.sh
./install-dealii.sh
```

To download Sapphire++ itself, you can either get the latest release as a
tarball on the [release page](https://github.com/sapphirepp/sapphirepp/releases)
or clone the repository:

```shell
git clone https://github.com/sapphirepp/sapphirepp
```

Afterwards you can compile Sapphire++, with `cmake` and `make`:

```shell
cd sapphirepp
export DEAL_II_DIR="path/to/deal.II"
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DEXAMPLES=ON
cd build && make
```

A detailed description of the installation process can be found in the
[documentation](https://sapphirepp.org).


## Getting started

To run Sapphire++ you need to provide a configuration file:
  
```shell
./sapphirepp parameter.prm
```

You enable parallel execution by using `mpirun`:

```shell
mpirun ./sapphirepp parameter.prm
```

We provide a documentation on our website [sapphire.org](https://sapphirepp.org)
together with [tutorials](https://sapphirepp.org/latest/examples.html) and
[examples](examples).


## Licence

Sapphire++ is distributed under the [LGPL 3.0 license](LICENSE).

If you use this software in your research, please cite the following paper:

> [Your Name], [Paper Title], [Journal], [Year], DOI: [DOI]

Here's the BibTeX entry for the paper:

```bibtex
@article{YourPaper,
  title={Paper Title},
  author={Your Name},
  journal={Journal},
  year={Year},
  doi={DOI}
}
```


## Contributing

We welcome contributions to Sapphire++. If you want to contribute code, please
follow the [contribution guidelines](CONTRIBUTING.md) and our
[code of conduct](CODE_OF_CONDUCT.md).
