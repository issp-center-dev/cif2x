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

    This parameter specifies whether ``ibrav`` parameter is used for Quantum ESPRESSO as the input of the crystal structure. When it is set to ``true``, the lattice is transformed to match the convention of Quantum ESPRESSO, and the lattice parameters ``a``, ``b``, ``c``, ``cosab``, ``cosac``, and ``cosbc`` are written to the input file as needed.

  ``tolerance`` (default value: 0.01)

    This parameter specifies the tolerance in the difference between the reconstructed Structure data and the original data when ``use_ibrav`` is set to ``true``.

  ``supercell`` (default value: none)

    This parameter specifies the size of supercell, when it is adopted, in the form of  [:math:`n_x`, :math:`n_y`, :math:`n_z`].
    

optional
--------
This section contains global settings needed for the first-principles calculation software. The available parameters are described in the corresponding sections below.
    
tasks
-----
This section defines contents of the input files. It is organized as a list of blocks, each corresponding to an input file, to allow for generating a set of input files for an input. The terms described in each block are explained in the following.


  ``mode`` (Quantum ESPRESSO)

    This parameter specifies the type of calculation. In the current version, the supported mode includes ``scf`` and ``nscf`` for pw.x of Quantum ESPRESSO. If an unsupported mode is specified, the settings in ``content`` will be exported as is.

  ``output_file`` (Quantum ESPRESSO)

    This parameter specifies the file name of the output.
    
  ``output_dir``

    This parameter specifies the directory name of the output. The default value is the current directory.

  ``content``

    This parameter describes the content of the output.
    For Quantum ESPRESSO, it contains the namelist data (blocks starting from ``&system``, ``&control``, etc.) in ``namelist`` block, and other card data (such as ``K_POINTS``) as individual blocks. Some card data may take parameters.

  ``template`` (Quantum ESPRESSO)

  ``template_dir`` (VASP)

    These parameters specifies the template file and the template directory for the input files, respectively. If they are not given, templates will not be used. The content of the template file is merged with those of ``content``. The entries in the template file will be superseded by those of ``content`` if the entries of the same keys appear both.
    

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
In the current version, ``scf`` mode and ``nscf`` mode of ``pw.x`` are supported.

optional section
------------------

  ``pp_file``

    This parameter specifies the index file in CSV format that relates the element type and the pseudo-potential file. This file contains the following columns: element name, type of pseudo-potential, nexclude, orbitals. An example line is given as:

    .. code-block::

      Fe,pbe-spn-rrkjus_psl.0.2.1,4,spd

    The name of the pseudo-potential file corresponding to the above example reads
    Fe.pbe-spn-rrkjus_psl.0.2.1.UPF .
      
  ``cutoff_file``

    This parameter specifies the index file in CSV format that relates the pseudo-potential file and the cutoff values. This file contains the following columns: name of pseudo-potential file, ``ecutwfc`` value, ``ecutrho`` value.

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

  - The information of the lattice vectors are obtained from the crystal structure data. When the ``data`` field is defined and contains a 3x3 matrix, that value will be used for the set of lattice vectors instead.

  **ATOMIC_SPECIES**

  - This block exports a list of atom species, atomic mass, and the file name of the pseudo-potential data.

  - The information of the atoms are obtained from the crystal structure data. The file names of the pseudo-potential data are referred from the CSV-formatted index file specified by ``pp_list`` parameter.

  - When the ``data`` field is defined and contains the required data, these values will be used instead.

  **ATOMIC_POSITIONS**

  - This block exports the atomic species and their fractional coordinates.

  - When ``ignore_species`` is given to specify an atomic species or a list of species, the values of ``if_pos`` for these species will be set to ``0``. It is used for MD or structure relaxations.

  - When the ``data`` field is defined and contains the required data, these values will be used instead.

  **K_POINTS**

  - This block exports the information of k points. The type of the output is specified by the ``option`` parameter that takes one of the following:

    - ``gamma``: uses :math:`\Gamma` point.

    - ``crystal``: generates a list of k points in mesh pattern. The mesh width is given by the ``grid`` parameter, or derived from the ``vol_density`` or ``k_resolution`` parameters.

    - ``automatic``: generates a mesh of k points. It is given by the ``grid`` parameter, or derived from the ``vol_density`` or ``k_resolution`` parameters. The shift is obtained from the ``kshifts`` parameter.

  - The mesh width is determined in the following order:

    - the ``grid`` parameter, specified by a list of :math:`n_x, n_y, n_z`, or a scalar value :math:`n`. For the latter, :math:`n_x = n_y = n_z = n` is assumed.
    - derived from the ``vol_density`` parameter.
    - derived from the ``k_resolution`` parameter, whose default value is 0.15.

  - When the ``data`` field is defined and contains the required data, these values will be used.


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
    
  ``psuedo_dir``

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

    The ``type`` parameter describes how KPOINTS are specified. The following values are allowed, with some types accepting parameters. See pymatgen.io.vasp manual for further details.

    - ``automatic``

      parameter: ``grid``

    - ``gamma_automatic``

      parameter: ``grid``, ``shift``

    - ``monkhorst_automatic``

      parameter: ``grid``, ``shift``

    - ``automatic_density``

      parameter: ``kppa``, ``force_gamma``

    - ``automatic_gamma_density``

      parameter: ``grid_density``

    - ``automatic_density_by_vol``

      parameter: ``grid_density``, ``force_gamma``

    - ``automatic_density_by_lengths``

      parameter: ``length_density``, ``force_gamma``

    - ``automatic_linemode``

      parameter: ``division``, ``path_type`` (corresponding to the ``path_type`` parameter of HighSymmKpath.)


Parameters for OpenMX
===============================

The entries of ``optional`` section and ``content`` part of the ``tasks`` section specific to OpenMX are explained below.

optional
--------

  ``data_dir``

    This parameter specifies the name of directory that holds files for pseudo-atomic orbitals and pseudo-potentials. It corresponds to the ``DATA.DIR`` parameter.

content
--------

  ``precision``

    This parameter specifies the set of pseudo-atomic orbitals listed in Tables 1 and 2 of Section 10.6 of the OpenMX manual. It is one of ``quick``, ``standard``, or ``precise``. The default value is ``quick``.

