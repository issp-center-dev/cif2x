.. _sec-getcif-tutorial:

Tutorial
================================================================

In this tutorial, the procedure to use the database query tool ``getcif`` is described for searching and obtaining crystallographic information from databases for the materials science.
It consists of getting an API key, preparing an input parameter file, and running the getcif program.
We will explain the steps along an example of searching and obtaining information for ABO3-type materials provided in the ``docs/tutorial/getcif`` directory.

Getting an API key
----------------------------------------------------------------

In order to access the Materials Project database via API, users need to register to the Materials Project and obtain an API key.
Visit the Materials Project website `https://next-gen.materialsproject.org <https://next-gen.materialsproject.org>`_, create an account and do Login. An API key is automatically generated on registration and shown in the user dashboard. The API key should be kept safe and not shared with others.

The API key is made available to getcif by one of the following ways:

  (a) storing in the pymatgen configuration file by typing in as follows:

        .. code:: bash

          $ pmg config --add PMG_MAPI_KEY <API_KEY>

      or editing the file ``~/.config/.pmgrc`` to include the following:
      
        .. code:: bash

          PMG_MAPI_KEY: <API_KEY>

  (b) setting to an environment variable by:
      
        .. code:: bash

          $ MP_API_KEY="<API_KEY>"
	  $ export MP_API_KEY

  (c) storing the API key to a file located in the directory where getcif is run.
      The default value of the file name is ``materials_project.key``. Otherwise, it is given in the input parameter file. The file name must end with ``.key``.

        .. code:: yaml

          database:
	    api_key_file: materials_project.key

      Comment: it will be recommended to exclude files with ``.key`` as a suffix from version control system. (e.g. for Git, add ``*.key`` in ``.gitignore`` file.)


Prepare an input parameter file
----------------------------------------------------------------

An input parameter file describes search conditions and data items to retrieve from databases.

An example is presented below. It is a text file in YAML format that contains information for accessing the database, search conditions, and types of data to obtain.
See :ref:`file format <sec-getcif-fileformat>` section for the details of specification.

In YAML format, parameters are given in dictionary form as ``keyword: value``, where ``value`` is a scalar such as a number or a string, or a set of values enclosed in ``[ ]`` or listed in itemized form, or a nested dictionary.
For the search conditions and data fields, a list may be given by a space-separated items without brackets as a special notation.

.. literalinclude:: ../../../../tutorial/getcif/input.yaml
   :language: yaml

The input parameter file consists of ``database``, ``option``, ``properties``, and ``fields`` sections.
The ``database`` section describes settings about connecting to databases.
In the example, ``target`` is set to Materials Project, though this term is not considered at present. ``api_key`` can be used to set the API key. The key may also be set in the pymatgen configuration file or in the environment variable. The latter is assumed in the tutorial.

The ``option`` section describes optional settings for the command execution.
``output_dir`` specifies the directory to place the obtained data. The default is the current directory. If ``dry_run`` is set to ``true``, getcif does not connect to the database; instead, it just prints the search conditions and exits. ``dry_run`` may be specified in the command-line option.

The ``properties`` section describes search conditions. They are given in the form of ``keyword: value`` and treated as AND conditions.
In the example, the search condition is specified to find materials with band gap less than or equal to 1.0, stable insulator, having composition formula of ABO3 (where A and B are arbitrary species), that belong to the space group ``Pm-3m`` (perovskite).
The ``band_gap`` takes a pair of values for the lower and upper limits, as well as the description such as ``< 1.0``.
The available terms for specifying search conditions are listed in the Appendix.

The ``fields`` section describes the data items to obtain. It is given as a YAML list, or a space-sparated list.
``structure`` specifies the crystal structure data that will be stored in CIF format.
``band_gap`` specifies the value of band gap, and ``symmetry`` specifies the information on the symmetry. ``material_id`` that refers to the index of material data in the Materials Project, and ``formula_pretty`` that refers to the composition formula are automatically obtained.
The available items are listed in the Appendix, or can be found in the help message of getcif command.
	      
Obtaining data
----------------------------------------------------------------

The program ``getcif`` is executed with the input parameter file (``input.yaml``) as follows.

.. code-block:: bash

  $ getcif input.yaml

Then it connects to the Materials Project database, and obtains the data that match the specified conditions.
The summary including the material IDs, the composition formulas, and other data items is printed to the standard output as follows.

.. literalinclude:: ../../../../tutorial/getcif/output_log.txt
   :language: text

The obtained data are placed in the directory specified by ``output_dir`` with the subdirectories of the material ID for each material.
In this example, seven subdirectories with names from mp-3163 to mp-977455 are created within ``result`` directory, and each subdirectory contains the following files:

  - band_gap
      the value of band gap

  - formula
      the composition formula (that corresponds to the field ``formula_pretty``)

  - structure.cif
      the crystal structure data in CIF format

  - symmetry
      the information about symmetry

If an option ``--dry-run`` is added as a command-line option to ``getcif``,
the program prints the search condition as follows, and exits.
It will be useful for checking the search parameters.
	      
.. literalinclude:: ../../../../tutorial/getcif/output_dryrun.txt
   :language: text
