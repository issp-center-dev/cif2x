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

.. **Official pages**
.. 
..   - `GitHub repository <https://github.com/issp-center-dev/HTP-tools-dev>`_
.. 
.. **Downloads**
.. 
..   HTP-tools can be downloaded by the following command with git:
.. 
..   .. code-block:: bash
.. 
..     $ git clone -b cif2x git@github.com:issp-center-dev/HTP-tools-dev.git

**Downloads**

  *The source package is not yet publicly available.*
   
**Installation**

  Once the source files are obtained, you can install HTP-tools by running the following command. The required libraries will also be installed automatically at the same time. 

  .. code-block:: bash

     $ cd ./HTP-tools-dev
     $ python3 -m pip install .

  The executable file ``cif2x`` will be installed.
  You may need to add ``--user`` option next to ``install`` keyword above in case you are not allowed to install packages system-wide.


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
     |       |-- utils.py
     |-- sample/


**Basic usage**

  ``cif2x`` is a tool to generate a set of input files for first-principles calculation software. It takes an input parameter file as a template, and generates parameter items that may vary by materials and calculation conditions from crystallographic data. In the present version, ``cif2x`` supports Quantum ESPRESSO, VASP, and OpenMX.

  #. Prepare input parameter file

      First, you need to create an input parameter file in YAML format that describes contents of the input file to be generated for the first-principles calculation software.

  #. Prepare crystal structure files and pseudo-potential files

      The crystal structure data need to be prepared for the target materials. The file format is CIF, POSCAR, xfs, or those supported by pymatgen.

      For Quantum ESPRESSO, the pseudo-potential files and the index file in CSV format need to be placed. Their locations are specified in the input parameter file.

      For VASP, the location of the pseudo-potential files will be specified in a file ``~/.config/.pmgrc.yaml`` or by an environment variable. It may be specified in the input parameter file.

  #. Run command

      Run ``cif2x`` command with the input parameter file and the crystal structure data as arguments. To generate input files for Quantum ESPRESSO, the target option ``-t QE`` should be specified. The option turns to ``-t VASP`` for VASP, and ``-t OpenMX`` for OpenMX.

      .. code-block:: bash

          $ cif2x -t QE input.yaml material.cif

