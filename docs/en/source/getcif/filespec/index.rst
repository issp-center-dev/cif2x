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

``api_key`` (default value: None)

  This parameter specifies the API key to access to the database. It is obtained from Materials Project as a registered user.
  If api_key is not specified or the value is empty, the API key is obtained from the environment variable ``MP_API_KEY``, or from the parameter ``PMG_MAPI_KEY`` of the pymatgen configuration file in ``~/.config/.pmgrc``.
    

option
--------------------------------

This section contains global settings needed for the first-principles calculation software. The available parameters are described in the corresponding sections below.

``output_dir`` (default value: ``""``)

  This parameter specifies the directory name to store the data. The retrieved data are placed in this directory under the subdirectories by the material ID for each material. The default value is the current directory.

``dry_run`` (default value: ``False``)

  When this parameter is set to True, getcif prints the search conditions and exists without connecting to the database. It is useful to check the content of the query.

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

    ``< max``
      less than or equal to ``max``.

    ``> min``
      more than or equal to ``min``.

    ``min ~ max``
      between ``min`` and ``max``.

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
