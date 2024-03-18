Installation and basic usage
================================================================

**Prerequisite**

  Input file generator for first-principles calculation ``cif2x`` included in HTP-tools requires the following programs and libraries:

  - python 3.x
  - pymatgen module
  - ruamel.yaml module
  - f90nml module
  - qe-tools module
  - numpy module
  - pandas module
  - monty module
  - OpenBabel module (optional)
  - AkaiKKRPythonUtil module

**Official pages**

  - `GitHub repository <https://github.com/issp-center-dev/cif2x>`_

**Downloads**

  cif2x can be downloaded by the following command with git:

  .. code-block:: bash

     $ git clone https://github.com/issp-center-dev/cif2x.git

**Installation**

  Once the source files are obtained, you can install cif2x by running the following command. The required libraries will also be installed automatically at the same time. 

  .. code-block:: bash

     $ cd ./cif2x
     $ python3 -m pip install .

  The executable file ``cif2x`` will be installed.
  You may need to add ``--user`` option next to ``install`` keyword above in case you are not allowed to install packages system-wide.

  AkaiKKRPythonUtil module need to be installed separately. The source package is available from `the repository <https://github.com/AkaiKKRteam/AkaiKKRPythonUtil>`_. Then follow the steps below to install the module along with the required seaborn module:

  .. code-block:: bash

     $ git clone https://github.com/AkaiKKRteam/AkaiKKRPythonUtil.git
     $ cd ./AkaiKKRPythonUtil/library/PyAkaiKKR
     $ python3 -m pip install .
     $ python3 -m pip install seaborn


**Directory structure**

  ::

     .
     |-- LICENSE
     |-- README.md
     |-- pyproject.toml
     |-- docs/
     |   |-- ja/
     |   |-- en/
     |   |-- tutorial/
     |-- src/
     |   |-- cif2x/
     |       |-- __init__.py
     |       |-- main.py
     |       |-- cif2struct.py
     |       |-- struct2qe.py
     |       |-- qe/
     |       |   |-- __init__.py
     |	     |   |-- calc_mode.py
     |	     |   |-- cards.py
     |	     |   |-- content.py
     |	     |   |-- qeutils.py
     |	     |   |-- tools.py
     |       |-- struct2vasp.py
     |       |-- struct2openmx.py
     |       |-- openmx/
     |       |   |-- __init__.py
     |       |   |-- vps_table.py
     |       |-- struct2akaikkr.py
     |       |-- akaikkr/
     |       |   |-- make_input.py
     |       |   |-- read_input.py
     |       |   |-- run_cif2kkr.py
     |       |-- utils.py
     |-- sample/


**Basic usage**

  ``cif2x`` is a tool to generate a set of input files for first-principles calculation software. It takes an input parameter file as a template, and generates parameter items that may vary by materials and calculation conditions from crystallographic data. In the present version, ``cif2x`` supports Quantum ESPRESSO, VASP, OpenMX, and AkaiKKR.

  #. Prepare input parameter file

      First, you need to create an input parameter file in YAML format that describes contents of the input file to be generated for the first-principles calculation software.

  #. Prepare crystal structure files and pseudo-potential files

      The crystal structure data need to be prepared for the target materials. The file format is CIF, POSCAR, xfs, or those supported by pymatgen.

      For Quantum ESPRESSO, the pseudo-potential files and the index file in CSV format need to be placed. Their locations are specified in the input parameter file.

      For VASP, the location of the pseudo-potential files will be specified in a file ``~/.config/.pmgrc.yaml`` or by an environment variable. It may be specified in the input parameter file.

  #. Run command

      Run ``cif2x`` command with the input parameter file and the crystal structure data as arguments. To generate input files for Quantum ESPRESSO, the target option ``-t QE`` should be specified. The option turns to ``-t VASP`` for VASP, ``-t OpenMX`` for OpenMX, and ``-t AkaiKKR`` for AkaiKKR.

      .. code-block:: bash

          $ cif2x -t QE input.yaml material.cif

