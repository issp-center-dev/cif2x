Command reference
================================================================

getcif
----------------------------------------------------------------

  Retrieve crystallographic and other data from databases.

SYNOPSIS:

  .. code-block:: bash

    getcif [-v][-q] [--dry-run] input_yaml
    getcif -h
    getcif --version

DESCRIPTION:

  This program reads an input parameter file specified by ``input_yaml``, and connects to the database to submit a query and obtain the crystallographic data of materials. 
  It takes the following command line options.

  - ``input_yaml``

    specifies an input parameter file in YAML format.

  - ``-v``

    increases verbosity of the runtime messages. When specified multiple times, the program becomes more verbose.
    
  - ``-q``

    decreases verbosity of the runtime messages. It cancels the effect of ``-v`` option, and when specified multiple times, the program becomes more quiet.

  - ``--dry-run``

    displays search parameters and exits without connecting to the database. It allows to confirm the search conditions. This option supersedes the ``dry_run`` parameter in the input file.

  - ``-h``

    displays help and exits.

  - ``--version``

    displays version information.
