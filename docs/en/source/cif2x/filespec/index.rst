.. _sec-cif2x-fileformat:

======================
File format
======================

Input parameter file
======================

An input parameter file describes information necessary to generate input files for first-principles calculation software by ``cif2x``. It should be given in YAML format, and consist of the following sections.

  1. structure section: describes how to handle crystal structure data.

  2. optional section: describes pseudo-potential files, and symbol definitions for reference feature of YAML.

  3. tasks section: describes contents of input files.

     
structure
---------

  ``use_ibrav`` (default value: ``false``)

    This parameter specifies whether ``ibrav`` parameter is used for Quantum ESPRESSO as the input of the crystal structure. When it is set to ``true``, the lattice is transformed to match the convention of Quantum ESPRESSO, and the lattice parameters ``a``, ``b``, ``c``, ``cosab``, ``cosac``, and ``cosbc`` are written to the input file as needed. The reconstructed structure is checked against the original one with a fixed tolerance of ``0.01``; the ``tolerance`` parameter below does not affect this check.

  ``use_primitive`` (default value: ``false``)

    When it is set to ``true``, the structure is reduced to the primitive cell (by ``get_primitive_structure`` of pymatgen) before the input files are generated. It is required for the ``respack`` target and for the ``bands`` mode of Quantum ESPRESSO, because the high-symmetry k-path is defined in the standardized primitive cell.

  ``tolerance`` (default value: 0.01)

    This parameter specifies the tolerance used in the primitive-cell reduction when ``use_primitive`` is set to ``true``. It has no effect otherwise.

  ``supercell`` (default value: none)

    This parameter specifies the size of supercell, when it is adopted, in the form of  [:math:`n_x`, :math:`n_y`, :math:`n_z`]. For example, to build a supercell doubled along each axis:

    .. code-block:: yaml

      structure:
        supercell: [2, 2, 2]

  ``cell_shape`` (required for molecular structure data)

    When the structure file holds an isolated molecule (for example an XYZ file) instead of a periodic structure, this parameter specifies the box in which the molecule is placed, centered. It takes a scalar ``L`` (a cubic box of side ``L`` angstrom), a list ``[Lx, Ly, Lz]`` (box lengths in angstrom), or ``auto``. With ``auto``, a tentative cubic box of 10 angstrom is used and, for OpenMX, the ``Atoms.UnitVectors`` block is omitted so that OpenMX determines the cell automatically. ``cif2x`` stops with an error when a molecular structure is given without ``cell_shape``.

  ``tol_deg`` (default value: 5.0)

    When the structure data carries magnetic moments (for example an mCIF file), this parameter specifies the angle tolerance in degrees used to decide whether the moments are collinear. If they are not collinear, a noncollinear calculation is enforced for Quantum ESPRESSO even when ``nspin = 2`` is given.


optional
--------
This section contains global settings needed for the first-principles calculation software. The available parameters are described in the corresponding sections below.
    
tasks
-----
This section defines contents of the input files. It is organized as a list of blocks, each corresponding to an input file, to allow for generating a set of input files for an input. The terms described in each block are explained in the following.


  ``mode`` (Quantum ESPRESSO, RESPACK)

    This parameter specifies the type of calculation. For pw.x of Quantum ESPRESSO, ``scf``, ``nscf``, ``relax``, ``vc-relax``, and ``bands`` are supported; the value is written to the ``calculation`` parameter of ``&control``. If an unsupported mode is specified, the settings in ``content`` will be exported as is. For the ``respack`` target, the tasks with ``scf``, ``nscf``, and ``bands`` are generated as Quantum ESPRESSO inputs, and a task with ``respack`` generates the RESPACK input (see :ref:`sec-cif2x-filespec-respack`).

  ``output_file`` (required for Quantum ESPRESSO, OpenMX, AkaiKKR, and RESPACK)

    This parameter specifies the file name of the output. For VASP, the file names are fixed (``INCAR``, ``KPOINTS``, ``POSCAR``, ``POTCAR``) and this parameter is not used.

  ``output_dir``

    This parameter specifies the directory name of the output. The default value is the current directory.

  ``content``

    This parameter describes the content of the output.
    For Quantum ESPRESSO, it contains the namelist data (blocks starting from ``&system``, ``&control``, etc.) in ``namelist`` block, and other card data (such as ``K_POINTS``) as individual blocks. Some card data may take parameters.

  ``template`` (Quantum ESPRESSO, OpenMX, AkaiKKR, RESPACK)

  ``template_dir`` (VASP)

    These parameters specifies the template file and the template directory for the input files, respectively. If they are not given, templates will not be used. The content of the template file is merged with those of ``content``. The entries in the template file will be superseded by those of ``content`` if the entries of the same keys appear both.

  ``precision`` (OpenMX)

    This parameter specifies the set of pseudo-atomic orbitals listed in Tables 1 and 2 of Section 10.6 of the OpenMX manual. It is one of ``quick``, ``standard``, or ``precise``. The default value is ``quick``. Note that it is a parameter of the task itself, not of ``content``.


Specifying parameter set
----------------------------------------
An input parameter may be given a list or range of parameters. In this case, a separate directory is created for every combination of parameters to store the generated input files. A special syntax ``${...}`` is used to specify the parameter set as follows:

- a list: ``${[ A, B, ... ]}``

  a set of parameter values is described as a Python list. Each entry may be a scalar value, or a list of values.

- a range: ``${range(N)}``, ``${range(start, end, step)}``

  a range of parameter is given by the keyword ``range``. The former specifies the values from ``0`` to ``N-1``, and the latter from ``start`` to ``end`` with every ``step``. (If ``step`` is omitted, it is assumed to be ``1``.)
  

Parameters for Quantum ESPRESSO
===============================

The entries of ``optional`` section and ``content`` part of the ``tasks`` section specific to Quantum ESPRESSO are explained below.
In the current version, the ``scf``, ``nscf``, ``relax``, ``vc-relax``, and ``bands`` modes of ``pw.x`` are supported.

optional section
------------------

  ``pp_file``

    This parameter specifies the index file in CSV format that relates the element type and the pseudo-potential file. The first line is a header, and the columns are ``element,pseudopotential,nexclude,orbitals``. The meaning of each column is:

    - ``element``: the element symbol.
    - ``pseudopotential``: the name stem of the pseudo-potential file, i.e. the part of the file name excluding the element symbol and the ``.UPF`` extension. The pseudo-potential file name is assembled as ``{element}.{pseudopotential}.UPF`` (for SOC runs, i.e. noncollinear with spin-orbit coupling, the generated ``ATOMIC_SPECIES`` filename uses the ``rel-`` variant, e.g. ``{element}.rel-{pseudopotential}.UPF``).
    - ``nexclude``: the number of core bands to exclude when counting the number of bands (``nbnd``) for Wannier functions.
    - ``orbitals``: the orbitals to be Wannierized, given as a string of the characters ``s``, ``p``, ``d``, ``f``. It is used to compute ``nbnd`` (see ``_find_nbnd_info`` in ``src/cif2x/struct2qe.py``).

    An example line is given as:

    .. code-block::

      element,pseudopotential,nexclude,orbitals
      Fe,pbe-spn-rrkjus_psl.0.2.1,4,spd

    This line corresponds to the pseudo-potential file ``Fe.pbe-spn-rrkjus_psl.0.2.1.UPF``.

  ``cutoff_file``

    This parameter specifies the index file in CSV format that relates the pseudo-potential file and the cutoff values. The columns are ``pseudofile,ecutwfc,ecutrho``: the name of the pseudo-potential file, the ``ecutwfc`` value, and the ``ecutrho`` value, respectively.

  ``pseudo_dir``

    This parameter specifies the name of the directory that holds pseudo-potential files. It is used when the cutoff values are obtained from the pseudo-potential files.
    It is indenepent from the ``pseudo_dir`` parameter in the input files for Quantum ESPRESSO.
    

content
--------

  **namelist**

  - The lattice specifications in ``&system`` block will be superseded according to ``use_ibrav`` parameter in the ``structure`` section.

    - ``use_ibrav = false``:
      ``ibrav`` is set to ``0``, and the lattice parameters including ``a``, ``b``, ``c``, ``cosab``, ``cosac``, ``cosbc``, ``celldm`` are removed.

    - ``use_ibrav = true``:
      ``ibrab`` is set to the index of Bravais lattices obtained from the crystal structure data. The Structure data will be reconstructed to match the convention of Quantum ESPRESSO.

  - ``nat`` (the number of atoms) and ``ntyp`` (the number of element types) will be superseded by the values obtained from the crystal structure data.

  - The cutoff values ``ecutwfc`` and ``ecutrho`` are obtained from the pseudo-potential files if these parameters are left blank.

  **CELL_PARAMETERS**

  - This block will not be generated if ``use_ibrav`` is set to ``true``. Otherwise, the lattice vectors are exported in units of angstrom.

  - The lattice vectors are always obtained from the crystal structure data; they cannot be overridden from the input parameter file.

  **ATOMIC_SPECIES**

  - This block exports a list of atom species, atomic mass, and the file name of the pseudo-potential data.

  - The information of the atoms are obtained from the crystal structure data. The file names of the pseudo-potential data are referred from the CSV-formatted index file specified by ``pp_file`` parameter. ``pp_file`` is required to generate this block.

  **ATOMIC_POSITIONS**

  - This block exports the atomic species and their fractional coordinates. The coordinate type is given by the ``option`` parameter (default: ``crystal``).

  - When ``ignore_species`` is given to specify an atomic species or a list of species, the values of ``if_pos`` for these species will be set to ``0``. It is used for MD or structure relaxations.

  **K_POINTS**

  - This block exports the information of k points. The type of the output is specified by the ``option`` parameter that takes one of the following:

    - ``gamma``: uses :math:`\Gamma` point.

    - ``crystal``: generates a list of k points in mesh pattern. The mesh width is given by the ``grid`` parameter, or derived from the ``vol_density`` or ``k_resolution`` parameters.

    - ``automatic``: generates a mesh of k points. It is given by the ``grid`` parameter, or derived from the ``vol_density`` or ``k_resolution`` parameters. The shift is obtained from the ``kshifts`` parameter.

    - ``crystal_b``: generates a band path from the high-symmetry k-path of the crystal (``HighSymmKpath`` of pymatgen). The number of points between high-symmetry points is given by ``line_npoints`` (default: 20), and a custom path may be given by ``path`` as a list of label sequences. This option is required for the ``bands`` mode. See :ref:`sec-cif2x-appendix` for details.

  - The mesh width is determined in the following order:

    - the ``grid`` parameter, specified by a list of :math:`n_x, n_y, n_z`, or a scalar value :math:`n`. For the latter, :math:`n_x = n_y = n_z = n` is assumed.
    - derived from the ``vol_density`` parameter.
    - derived from the ``k_resolution`` parameter, whose default value is 0.15.

  - The ``data`` field is honoured only for the band-path options ``crystal_b`` and ``tpiba_b`` (and for an unrecognized ``option``, which is passed through as written). In these cases, the rows given in ``data`` are written as they are instead of being generated. For the other options, and for the other cards above, the data are always generated from the crystal structure.


Parameters for VASP
===============================

The entries of ``optional`` section and ``content`` part of the ``tasks`` section specific to VASP are explained below.

optional
--------

The type and the location of pseudo-potential files are specified.

According to pymatgen, the pseudo-potential files are obtained from 
``PMG_VASP_PSP_DIR``/*functional*/POTCAR.{element}(.gz) or
``PMG_VASP_PSP_DIR``/*functional*/{element}/POTCAR,
where
``PMG_VASP_PSP_DIR`` points to the directory and it is given in the configuration file
``~/.config/.pmgrc.yaml`` or by the environment variable of the same name.
*functional* refers to the type of the pseudo-potential, whose value is predefined as
``POT_GGA_PAW_PBE``, ``POT_LDA_PAW``, etc.


  ``pseudo_functional``

    This parameter specifies the type of the pseudo-potential. The relation to the *functional* value above is defined in the table of pymatgen, for example, by ``PBE`` to ``POT_GGA_PAW_PBE``, or by ``LDA`` to ``POT_LDA_PAW``, or in a similar manner.
    

When the ``pseudo_dir`` parameter is specified, it is used as the directory that holds the pseudo-potential files, ignoring the convention of pymatgen.
    
  ``pseudo_dir``

    This parameter specifies the directory that holds the pseudo-potential files. The paths to the pseudo-potential file turn to ``pseudo_dir``/POTCAR.{element}(.gz), or ``pseudo_dir``/{element}/POTCAR.

tasks
-----

The template files are assumed to be placed in the directory specified by the ``template_dir`` parameter by the names ``INCAR``, ``KPOINTS``, ``POSCAR``, and ``POTCAR``. The missing files will be ignored.



content
-------

  **incar**

  - This block contains parameters described in the INCAR file

  **kpoints**

  - ``type``

    The ``type`` parameter describes how KPOINTS are specified. The default value is ``gamma_automatic``. The following values are allowed, with some types accepting parameters; the parameter names below are the keys to be written under ``kpoints`` (they are case-insensitive). See pymatgen.io.vasp manual for the meaning of each scheme.

    - ``automatic``

      parameter: ``grid`` (a single integer, the subdivision length; default: 1)

    - ``gamma_automatic``

      parameter: ``kpoints`` (a list :math:`[n_x, n_y, n_z]`; default: ``[1, 1, 1]``), ``shift`` (default: ``[0, 0, 0]``)

    - ``monkhorst_automatic``

      parameter: ``kpoints`` (a list :math:`[n_x, n_y, n_z]`; default: ``[2, 2, 2]``), ``shift`` (default: ``[0, 0, 0]``)

    - ``automatic_density``

      parameter: ``grid_density`` (k points per reciprocal atom; default: 0.1), ``force_gamma`` (default: ``false``)

    - ``automatic_gamma_density``

      parameter: ``grid_density`` (default: 0.1)

    - ``automatic_density_by_vol``

      parameter: ``grid_density`` (k points per reciprocal volume; default: 1), ``force_gamma`` (default: ``false``)

    - ``automatic_density_by_lengths``

      parameter: ``length_density`` (a list of three lengths; default: ``[10.0, 10.0, 10.0]``), ``force_gamma`` (default: ``false``)

    - ``automatic_linemode``

      parameter: ``division`` (default: 1), ``path_type`` (corresponding to the ``path_type`` parameter of HighSymmKpath.)

    For example, a :math:`6 \times 6 \times 6` :math:`\Gamma`-centered mesh is written as:

    .. code-block:: yaml

      kpoints:
        type: gamma_automatic
        kpoints: [6, 6, 6]


Parameters for OpenMX
===============================

The entries of ``optional`` section and ``content`` part of the ``tasks`` section specific to OpenMX are explained below.

optional
--------

  ``data_path``

    This parameter specifies the name of directory that holds files for pseudo-atomic orbitals and pseudo-potentials. It corresponds to the ``DATA.PATH`` parameter.

tasks
-----

  ``precision``

    This parameter specifies the set of pseudo-atomic orbitals listed in Tables 1 and 2 of Section 10.6 of the OpenMX manual. It is one of ``quick``, ``standard``, or ``precise``. The default value is ``quick``. It is written at the same level as ``output_file`` and ``template``, not inside ``content``:

    .. code-block:: yaml

      tasks:
        - output_file: openmx.dat
          template: template.in
          precision: precise
          content:
            system.name: NaCl

content
--------

The ``content`` part contains OpenMX keywords, written to the output file as they are (a list value is written as a block ``<Key ... Key>``). The following keywords are generated from the crystal structure data and overwrite the values in the template or ``content``: ``Species.Number``, ``Definition.of.Atomic.Species``, ``Atoms.Number``, ``Atoms.SpeciesAndCoordinates.Unit``, ``Atoms.SpeciesAndCoordinates``, ``Atoms.UnitVectors.Unit``, and ``Atoms.UnitVectors``. ``System.Name`` defaults to the base name of the structure file when it is not given.

Parameters for AkaiKKR
===============================

The entries of ``optional`` section and ``content`` part of the ``tasks`` section specific to AkaiKKR are explained below.

optional
--------

  ``workdir``

    This parameter specifies the directory in which temporal files are stored. If it is not given, ``/tmp`` or the value of the environment variable ``TMPDIR`` is used.


content
--------

The ``content`` part contains the input parameters of AkaiKKR.
A blank is written to the input file for an unspecified parameter, to which the default value defined in AkaiKKR will be assumed.
The parameter values listed below are replaced by the values obtained from the crystal structure data.

- ``brvtyp``, except when it is set to ``aux`` (or a string that contains ``aux``).

- lattice parameters, ``a``, ``c/a``, ``b/a``, ``alpha``, ``beta``, ``gamma``, ``r1``, ``r2``, ``r3``.

- type information, ``ntyp``, ``type``, ``ncmp``, ``rmt``, ``field``, ``mxl``, ``anclr``, ``conc``.

- element information, ``natm``, ``atmicx``, ``atmtyp``.

For ``rmt`` and ``field``, the values specified in the input parameter file will be used only when they are lists having the same number of elements as ``ntyp``.


.. _sec-cif2x-filespec-respack:

Parameters for RESPACK
===============================

The ``respack`` target generates the inputs of a RESPACK workflow: the ``scf`` and ``nscf`` (and optionally ``bands``) runs of Quantum ESPRESSO, followed by the RESPACK control file. The ``structure`` section must have ``use_primitive: true`` and ``use_ibrav: false``, because the high-symmetry k-path written to the RESPACK input is defined in the standardized primitive cell.

tasks
-----

Each task is distinguished by ``mode``:

- ``scf``, ``nscf``, ``bands``: generated exactly as for the Quantum ESPRESSO target, with the parameters described in :ref:`sec-cif2x-fileformat` above (``optional`` with ``pp_file``, ``cutoff_file``, ``pseudo_dir`` and ``content`` with ``namelist`` and the cards). The ``nscf`` task should set ``nosym: true`` and ``noinv: true`` in ``content.namelist.system``, as required by ``qe2respack``; ``cif2x`` warns when they are missing.

- ``respack``: generates the RESPACK control file. ``output_file`` and ``template`` are required. The parameter-set syntax ``${...}`` is not expanded for this task.

template and content (``mode: respack``)
----------------------------------------

The template is a namelist-only file: it may contain only the RESPACK namelists ``&param_chiqw``, ``&param_wannier``, ``&param_interpolation``, ``&param_visualization``, and ``&param_calc_int``. Any other namelist, or any text outside the namelists (such as a k-point block), is rejected, because ``cif2x`` generates the k-path block itself.

``content`` holds the same namelists as mappings; its values are merged over those of the template (the ``content`` value wins for the same key). Namelist and parameter names are case-insensitive.

.. code-block:: yaml

  - mode: respack
    output_file: input.in
    template: respack.in_tmpl
    content:
      param_wannier:
        N_wannier: 3
      param_interpolation:
        dense: [8, 8, 8]

The following conditions are checked, and the run stops with an error when they are not met:

- ``N_wannier`` in ``&param_wannier`` must be a positive integer. It is the number of Wannier functions chosen by the user and is not derived automatically.

- ``Lower_energy_window`` and ``Upper_energy_window`` in ``&param_wannier`` must be numbers with ``Lower_energy_window < Upper_energy_window``.

- ``N_initial_guess`` in ``&param_wannier`` must be ``0`` or omitted; the initial guess is always SCDM, and ``N_initial_guess = 0`` is written.

- ``dense`` in ``&param_interpolation`` must be an integer :math:`n` (used as :math:`n \times n \times n`) or a list of three integers.

- ``reading_sk_format`` in ``&param_interpolation`` must be ``0`` (the default).

The following are generated by ``cif2x`` and written to the output: ``N_sym_points`` in ``&param_interpolation`` (the number of high-symmetry points of each segment of the path), and the block of high-symmetry k-point coordinates that follows the ``&param_interpolation`` namelist, obtained from ``HighSymmKpath`` of pymatgen. The namelists are written in the order ``&param_chiqw``, ``&param_wannier``, ``&param_interpolation``, ``&param_visualization``, ``&param_calc_int``; a namelist is written only when it appears in the template or ``content``.

A complete example is given in ``docs/tutorial/cif2x/respack/`` and in ``sample/cif2x/respack/`` of the repository.
