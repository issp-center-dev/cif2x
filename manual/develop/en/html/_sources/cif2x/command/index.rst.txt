Command reference
================================================================

cif2x
----------------------------------------------------------------

  Generate input files for first-principles calculation software

SYNOPSIS:

  .. code-block:: bash

    cif2x [-v][-q] -t target input_yaml material.cif
    cif2x -h
    cif2x --version

DESCRIPTION:

  This program reads an input parameter file specified by ``input_yaml`` and a crystal data file specified by ``material.cif``, and generates a set of input files for first-principles calculation software. In the current version, the supported software includes Quantum ESPRESSO, VASP, and OpenMX.
  It takes the following command line options.

  - ``-v``

    increases verbosity of the runtime messages. When specified multiple times, the program becomes more verbose.
    
  - ``-q``

    decreases verbosity of the runtime messages. It cancels the effect of ``-v`` option, and when specified multiple times, the program becomes more quiet.

  - ``-t`` *target*

    specifies the target first-principles calculation software. The supported software for *target* is listed as follows:

    - ``QE``, ``espresso``, ``quantum_espresso``: generates input files for Quantum ESPRESSO.

    - ``VASP``: generates input files for VASP.

    - ``OpenMX``: generates input files for OpenMX

  - ``input_yaml``

    specifies an input parameter file in YAML format.

  - ``material.cif``

    specifies crystal structure data file. It is in CIF (Crystallographic Information Framework) format, or other format supported by pymatgen.

  - ``-h``

    displays help and exits.

  - ``--version``

    displays version information.
