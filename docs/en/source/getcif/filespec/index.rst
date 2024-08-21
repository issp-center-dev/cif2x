.. _sec-getcif-fileformat:

================================
File format
================================

Input parameter file
================================

An input parameter file describes information to search for crystallographic and other data from Materials Project database by getcif. It should be given in YAML format, and consist of the following sections.

  #. database section: describes information on the database to connect.

  #. option section: describes output directory and other parameters for command execution.

  #. properties section: describes search conditions.

  #. fields section: describes types of data to be retrieved.

     
database
--------------------------------

``target``

  This parameter specifies the database to connected to. At present this parameter is ignored.

``api_key_file`` (default value: ``materials_project.key``)

  This parameter specifies a name of a file that contains the API key to access to the database.
  The suffix of the file name must be ``.key``.
  If the file does not exist or it does not contain a valid value, the API key is obtained from the environment variable ``MP_API_KEY``, or from the parameter ``PMG_MAPI_KEY`` of the pymatgen configuration file in ``~/.config/.pmgrc``.

  The API key file is a text file. A line starting with ``#`` is regarded as a comment. The heading and trailing spaces are ignored. When the file contains more than one line, the API key is taken from the first valid line.
    

option
--------------------------------

This section contains global settings needed for the first-principles calculation software. The available parameters are described in the corresponding sections below.

``output_dir`` (default value: ``""``)

  This parameter specifies the directory name to store the data. The retrieved data are placed in this directory under the subdirectories by the material ID for each material. The default value is the current directory.

``dry_run`` (default value: ``False``)

  When this parameter is set to True, getcif prints the search conditions and exists without connecting to the database. It is useful to check the content of the query.

``symprec`` (default value: ``None``)

  This parameter specifies the tolerance in calculating the symmetry of a crystal structure when the structure data are written to a CIF file. By default, ``None`` is specified, in which case a CIF file is generated without considering symmetry.

  ``symprec`` is a parameter that specifies the tolerance used to determine the symmetry of a crystal structure. When calculating the symmetry of a crystal structure, it is essential to consider the slight displacements of atomic positions and the precision of numerical calculations. ``symprec`` controls the allowable range of these displacements and serves as a threshold for deciding whether a symmetry operation should be applied.

  If ``symprec`` is set to a smaller value (e.g., 0.01), the symmetry determination becomes more stringent, and even minor displacements in the crystal structure may prevent the application of symmetry operations. This can result in the identification of a lower-symmetry space group. Conversely, if ``symprec`` is set to a larger value (e.g., 1.0), the symmetry determination is more lenient, allowing small displacements to be ignored, which may lead to the recognition of a higher-symmetry space group.

  When the ``symmetry`` field is specified in the fields section, the symmetry information determined using the default ``symprec=0.1`` in the Materials Project is obtained and written to a text file (``symmetry``).


properties
--------------------------------

This section defines the search conditions.
The conditions such as the element types, the crystal symmetry, or the values of physical properties are specified in the ``keyword: value`` format. They are treated as AND condition.
The available terms, based on the Materials Project API, conform to the parameters of
the ``materials.summary.search`` method in the mp-api library. The list of terms are summarized in the Appendix, and can be seen by ``getcif --help``.

The format of the parameter values is shown below. It follows the YAML specification with several extension for brief description.

- a number, a string

    describe as-is.

- a boolean value

    describe as ``true`` or ``false``.

- a list of numbers or strings

    describe in the indented style (block style) or in the comma-separated list enclosed by the bracket (flow style) in YAML notation.
    It is also available that it is described as a space-separated list, for example:

    .. code:: yaml

	element: Sr Ti

- a range of numerical value

    described as a list of two numbers such as ``[ min, max ]``, or a pair of two numbers separated by a space as ``min max``. The following formats are also available.

    ``<= max``
      less than or equal to ``max``.

    ``< max``
      less than ``max``. (For a real number, it is equivalent to ``<= max``. For an integer, it is treated as ``<= max-1``.)

    ``>= min``
      more than or equal to ``min``.

    ``> min``
      more than ``min``. (For a real number, it is equivalent to ``>= min``. For an integer, it is treated as ``>= min+1``.)

    ``min ~ max``
      between ``min`` and ``max``.

    N.B.:

      - A space must be placed between the symbol and the number.

      - Due to the YAML syntax that the symbol ``">"`` at the beginning of a term is treated as a special character, ``> min`` and ``>= min`` should be enclosed by quotes as ``"> min"`` and ``">= min"``, respectively.

      - In list notations, ``<= max`` and ``>= min`` are denoted as ``[ None, max ]`` and ``[ min, None ]``, respectively.


- wild card symbols

    The term ``formula`` accepts wild card symbols ``*`` for elements. In this case, the whole value is enclosed by ``" "``. For example,

    .. code:: yaml

	formula: "**O3"

    for :math:`ABO_3`-type materials.


fields
--------------------------------

This section defines the types of data to be retrieved.
A list of types is described in the YAML format, or as a space-sparated strings. In the latter format, it can be given in multiple-line format using the "|" notation of YAML.

The available types of data conform to the ``field`` parameter of the Materials Project API. They are listed in the Appendix, and can be viewd by ``getcif --help``.

The types ``material_id`` and ``formula_pretty`` are retrieved automatically.

The obtained data are placed in the directory specified by ``output_dir`` parameter under the subdirectories of the material_id for each material. Each item is stored as a separate file of the item name. The crystal structure data (``structure``) is stored in a file ``structure.cif`` in CIF format.
