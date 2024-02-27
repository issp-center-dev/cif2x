Installation and basic usage
================================================================

**Prerequisite**

  A tool to retrieve crystallographic data from databases ``getcif`` included in HTP-tools requires the following programs and libraries:

  - python 3.x
  - pymatgen module
  - ruamel.yaml module
  - mp-api module

**Official pages**

  - GitHub repository `https://github.com/issp-center-dev/getcif <https://github.com/issp-center-dev/getcif>`_

**Downloads**

  getcif can be downloaded by the following command with git:

  .. code-block:: bash

     $ git clone https://github.com/issp-center-dev/getcif.git

**Installation**

  Once the source files are obtained, you can install getcif by running the following command. The required libraries will also be installed automatically at the same time. 

  .. code-block:: bash

     $ cd ./getcif
     $ python3 -m pip install .

  The executable file ``getcif`` will be installed.
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
     |   |-- getcif/
     |       |-- __init__.py
     |       |-- main.py
     |-- sample/


**Basic usage**

  ``getcif`` is a tool to retrieve crystallographic information and other properties of materials from the Materials Project database. Users can search database and obtain information by specifying symmetry, composition, or physical properties of materials.

  #. Getting an API key for Materials Project

     First, you need to register to Materials Project at their `website <https://next-gen.materialsproject.org/>`_ and obtain an API key to access the database.
     The API key (string) is stored in the configuration file using pymatgen tool, or set to an environment variable. It can be stored in a file specified in the input parameter file.

  #. Preparing an input parameter file

     Second, you need to prepare an input parameter file in YAML format that describes the search conditions and the types of data to retrieve. See file format section for further details.

  #. Running the command

     Run ``getcif`` command with the input parameter file as an argument.

     .. code-block:: bash

          $ getcif input.yaml

     Then, it connects to the database and submits queries according to the specified condition. The obtained data are placed in the directory with the material ID for each material. The crystal structure data are stored in the CIF format.
