****************************************************************
Introduction
****************************************************************

What is cif2x?
----------------------------------------------------------------

In recent years, the use of machine learning for predicting material properties and designing substances (known as materials informatics) has gained considerable attention.
The accuracy of machine learning depends heavily on the preparation of appropriate training data.
Therefore, the development of tools and environments for the rapid generation of training data is expected to contribute significantly to the advancement of research in materials informatics.

Cif2x is a tool that generates input files for first-principles calculations from cif files.
It constructs parts that vary depending on the type of material and computational conditions from crystal structure data, using input parameters as a template.
It is capable of generating multiple input files tailored to specific computational conditions.
Currently, it supports `VASP <https://www.vasp.at>`, `Quantum ESPRESSO <https://www.quantum-espresso.org>`, and `OpenMX <http://www.openmx-square.org>`,
with plans to support `AkaiKKR <http://kkr.issp.u-tokyo.ac.jp>` in the future.

License
----------------------------------------------------------------

The distribution of the program package and the source codes for cif2x follow GNU General Public License version 3 (GPL v3) or later.

Contributors
----------------------------------------------------------------

This software was developed by the following contributors.

-  ver.1.0-alpha (Released on 2023/12/28)

   -  Developers

      -  Kazuyoshi Yoshimi (The Instutite for Solid State Physics, The University of Tokyo)

      -  Tatsumi Aoyama (The Instutite for Solid State Physics, The University of Tokyo)

      -  Yuichi Motoyama (The Instutite for Solid State Physics, The University of Tokyo)

      -  Masahiro Fukuda (The Instutite for Solid State Physics, The University of Tokyo)

      -  Kota Ido (The Instutite for Solid State Physics, The University of Tokyo)

      -  Tetsuya Fukushima (The National Institute of Advanced Industrial Science and Technology (AIST))

      -  Shusuke Kasamatsu (Yamagata University)

      -  Takashi Koretsune (Tohoku University)

   -  Project Corrdinator

      -  Taisuke Ozaki (The Instutite for Solid State Physics, The University of Tokyo)


Copyright
----------------------------------------------------------------

.. only:: html

  |copy| *2023- The University of Tokyo. All rights reserved.*

  .. |copy| unicode:: 0xA9 .. copyright sign

.. only:: latex

  :math:`\copyright` *2023- The University of Tokyo. All rights reserved.*


This software was developed with the support of "Project for advancement of software usability in materials science" of The Institute for Solid State Physics, The University of Tokyo.

Operating environment
----------------------------------------------------------------

This tool was tested on the following platforms:

- Ubuntu Linux + python3

