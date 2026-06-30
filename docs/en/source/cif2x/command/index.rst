Command reference
================================================================

cif2x
----------------------------------------------------------------

  Generate input files for first-principles calculation software

SYNOPSIS:

  .. code-block:: bash

    cif2x [-v][-q][--dry-run] -t target input_yaml material.cif
    cif2x [-v][-q][--dry-run] -t target --mp-id ID [--symprec PREC][--api-key-file FILE] input_yaml
    cif2x -h
    cif2x --version

DESCRIPTION:

  This program reads an input parameter file specified by ``input_yaml`` and a crystal data file specified by ``material.cif``, and generates a set of input files for first-principles calculation software. In the current version, the supported software includes Quantum ESPRESSO, VASP, OpenMX, and AkaiKKR.
  It takes the following command line options.

  - ``-v``

    increases verbosity of the runtime messages. When specified multiple times, the program becomes more verbose.
    
  - ``-q``

    decreases verbosity of the runtime messages. It cancels the effect of ``-v`` option, and when specified multiple times, the program becomes more quiet.

  - ``-t`` *target*

    specifies the target first-principles calculation software. The supported software for *target* is listed as follows:

    - ``QE``, ``espresso``, ``quantum_espresso``: generates input files for Quantum ESPRESSO.

    - ``VASP``: generates input files for VASP.

    - ``OpenMX``: generates input files for OpenMX.

    - ``AkaiKKR``: generates input files for AkaiKKR.

    - ``respack``: generates the RESPACK workflow — the Quantum ESPRESSO inputs (``scf``/``nscf``/``bands`` tasks) and the RESPACK control file ``input.in`` (a ``mode: respack`` task). The structure must be the standardized primitive cell (``structure.use_primitive: true``, ``use_ibrav: false``). ``N_wannier`` and the energy windows are user physics supplied in the ``content``/template; initial guesses use SCDM (``N_initial_guess = 0``), targeting the ``respack-wannier-py`` port. The ``nscf`` task must set ``nosym = .true.``/``noinv = .true.`` for ``qe2respack``. Run order: ``scf`` → ``nscf`` → ``qe2respack`` → ``respack``.

  - ``--dry-run``

    prints the generated input files to standard output instead of writing them to disk. This is useful for previewing the result without creating any files or directories.

  - ``input_yaml``

    specifies an input parameter file in YAML format.

  - ``material.cif``

    specifies crystal structure data file, in CIF (Crystallographic Information Framework) format or another format supported by pymatgen. It is optional when ``--mp-id`` is given.

  - ``--mp-id`` *ID*

    fetches the crystal structure for the Materials Project material id *ID* (for example ``mp-149``) directly, instead of reading ``material.cif``. Provide exactly one of ``material.cif`` or ``--mp-id``. The API key is resolved from ``--api-key-file`` (default ``materials_project.key``), then the ``MP_API_KEY`` environment variable, then the ``PMG_MAPI_KEY`` entry in the pymatgen configuration (``~/.config/.pmgrc.yaml``) -- the same way as ``getcif``.

  - ``--symprec`` *PREC*

    symmetry tolerance applied when writing the fetched structure (default ``0.1``, matching the Materials Project). ``--symprec 0`` disables symmetry refinement. Used only with ``--mp-id``.

  - ``--api-key-file`` *FILE*

    file containing the Materials Project API key, one key per non-``#`` line (default ``materials_project.key``). Used only with ``--mp-id``.

  .. note::

     ``--mp-id`` reproduces the ``getcif`` then ``cif2x`` flow by writing the fetched structure to a temporary CIF and reading it back; the Materials Project final structure is used (not a conventionalized cell). Because the structure passes through CIF, site properties such as magnetic moments are not carried over, and disordered (partial-occupancy) structures are rejected by the Quantum ESPRESSO generator -- the same limitations as the two-step flow.

  - ``-h``

    displays help and exits.

  - ``--version``

    displays version information.
