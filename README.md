# cif2x

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![Python](https://img.shields.io/badge/python-3.9%2B-blue.svg)](https://www.python.org/)
![Version](https://img.shields.io/badge/version-1.1.0-blue.svg)
[![Documentation](https://img.shields.io/badge/docs-online-brightgreen.svg)](https://issp-center-dev.github.io/cif2x/manual/main/en/html/)

In recent years, the use of machine learning for predicting material properties and designing substances (known as materials informatics) has gained considerable attention.
The accuracy of machine learning depends heavily on the preparation of appropriate training data.
Therefore, the development of tools and environments for the rapid generation of training data is expected to contribute significantly to the advancement of research in materials informatics.

Cif2x is a tool to generate input files for first-principles calculation software. It takes crystal structure data in CIF format and input parameters as a template, and constructs the parts that vary depending on the type of material and conditions. It is capable of generating multiple input files tailored to specific computational conditions. Currently, it supports [VASP](https://www.vasp.at), [Quantum ESPRESSO](https://www.quantum-espresso.org), [OpenMX](http://www.openmx-square.org), and [AkaiKKR](http://kkr.issp.u-tokyo.ac.jp).

getcif is a tool to retrieve crystallographic information and other properties of materials from databases. The latest version of getcif provides access to Materials Project database. Users can search database and obtain information by specifying symmetry, composition, or physical properties of materials.


## Target applications

Quantum ESPRESSO, VASP, OpenMX, and AkaiKKR.

## Requirement

Python3 with pymatgen, qe-tools, AkaiKKRPythonUtil, and other library packages.

getcif requires mp-api, pymatgen, and other library packages.

## Install

- From source

``` bash
python3 -m pip install DIRECTORY_OF_THE_REPOSITORY
```

Installing the package provides the `cif2x` and `getcif` command-line tools.

## Usage

### cif2x

Generate input files from a CIF (or molecular structure) file and a YAML
parameter file. The target code is selected with `-t`
(`quantum_espresso`/`qe`, `vasp`, `openmx`, `akaikkr`; case-insensitive):

``` bash
cif2x -t quantum_espresso input.yaml structure.cif
```

### getcif

Query the Materials Project database and download structures and properties.
A Materials Project API key is required (supplied via a key file, the
environment, or pymatgen settings):

``` bash
getcif input.yaml
```

Example inputs are available under the `sample/` directory and in
`docs/tutorial/`. See the documentation below for the YAML input format.

## Documentation

- [User Guide (English)](https://issp-center-dev.github.io/cif2x/manual/main/en/html/)
- [ユーザーガイド (日本語)](https://issp-center-dev.github.io/cif2x/manual/main/ja/html/)

## Bug reports & questions

Please report bugs or ask questions on the
[GitHub issue tracker](https://github.com/issp-center-dev/cif2x/issues).

## License

The distribution of the program package and the source codes for cif2x follow
GNU General Public License version 3
([GPL v3](https://www.gnu.org/licenses/gpl-3.0.en.html)).

Copyright (c) <2023-> The University of Tokyo. All rights reserved.

This software was developed with the support of
"Project for Advancement of Software Usability in Materials Science"
of The Institute for Solid State Physics, The University of Tokyo.

## Official page

- [HTP-tools project site](https://www.pasums.issp.u-tokyo.ac.jp/htp-tools/)
- [Software repository](https://github.com/issp-center-dev/cif2x)

## Authors

Kazuyoshi Yoshimi (ISSP, Univ. of Tokyo), 
Tatsumi Aoyama (ISSP, Univ. of Tokyo), 
Yuichi Motoyama (ISSP, Univ. of Tokyo), 
Masahiro Fukuda (ISSP, Univ. of Tokyo), 
Kota Ido (ISSP, Univ. of Tokyo), 
Tetsuya Fukushima (AIST), 
Shusuke Kasamatsu (Yamagata University), 
Takashi Koretsune (Tohoku University), 
Taisuke Ozaki (ISSP, Univ. of Tokyo)
