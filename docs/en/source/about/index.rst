****************************************************************
Introduction
****************************************************************

What is HTP-tools?
----------------------------------------------------------------

In recent years, the use of machine learning for predicting material properties and designing substances (known as materials informatics) has gained considerable attention.
The accuracy of machine learning depends heavily on the preparation of appropriate training data.
Therefore, the development of tools and environments for the rapid generation of training data is expected to contribute significantly to the advancement of research in materials informatics.

HTP-Tools is a package which is specifically designed to support high throughput computations.
It includes 'moller', a scripting tool for managing large-scale computations on supercomputers and clusters, and 'cif2x', which generates input files for first-principles calculations from cif files.
cif2x is compatible with `VASP <https://www.vasp.at>`_ , `Quantum ESPRESSO  <https://www.quantum-espresso.org>`_ , and `OpenMX <http://www.openmx-square.org>`_
with future support for `AkaiKKR <http://kkr.issp.u-tokyo.ac.jp>`_ planned.
In addition, cif2x allows the generation of multiple input files tailored to specific computational conditions.

Using HTP tools, researchers can easily perform a large number of first-principles calculations on high-performance computers.
In the future, there are plans to develop tools that are compatible with various simulation software,
such as solvers for quantum lattice models, to further enrich the environment for generating training data for machine learning.

License
----------------------------------------------------------------

The distribution of the program package and the source codes for HTP-tools follow GNU General Public License version 3 (GPL v3) or later.

Contributors
----------------------------------------------------------------

This software was developed by the following contributors.

-  ver.1.0-beta (Released on 2023/12/08)

   -  Developers

      -  Kazuyoshi Yoshimi (The Instutite for Solid State Physics, The University of Tokyo)

      -  Tatsumi Aoyama (The Instutite for Solid State Physics, The University of Tokyo)

      -  Yuichi Motoyama (The Instutite for Solid State Physics, The University of Tokyo)

      -  Masahiro Fukuda (The Instutite for Solid State Physics, The University of Tokyo)

      -  Tetsuya Fukushima (The National Institute of Advanced Industrial Science and Technology (AIST))

      -  Kota Ido (The Instutite for Solid State Physics, The University of Tokyo)

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

HTP-tools was tested on the following platforms:

- Ubuntu Linux + python3

