.. _sec-cif2x-tutorial:

Tutorial
================================================================

The procedure to use the input file generator ``cif2x`` for first-principles calculation software consists of preparing an input parameter file, crystal structure data, and pseudo-potential files, and running the program ``cif2x``.
In the current version, the supported software includes Quantum ESPRESSO, VASP, OpenMX, and AkaiKKR.
In this tutorial, we will explain the steps along a sample for Quantum ESPRESSO in ``docs/tutorial/cif2x``.


Prepare an input parameter file
----------------------------------------------------------------

An input parameter file describes the content of input files for the first-principles calculation software.
An example is presented below. It is a text file in YAML format that contains options to crystal structure data, and contents of the input file used as an input for the first-principle calculation. See :ref:`file format <sec-cif2x-fileformat>` section for the details of specification.

In YAML format, parameters are given in dictionary form as ``keyword: value``, where ``value`` is a scalar such as a number or a string, or a set of values enclosed in ``[ ]`` or listed in itemized form, or a nested dictionary.

.. literalinclude:: ../../../../tutorial/cif2x/input.yaml
   :language: yaml


The input parameter file consists of ``structure``, ``optional``, and ``tasks`` sections.
The ``structure`` section specifies options to the crystal structure data.
The ``optional`` section holds global settings concerning the pseudo-potentials.

The ``tasks`` section describes inputs for the first-principles calculations. In case of generating multiple files for a series of calculations, the ``tasks`` section takes a list of parameter sets.
For each set, the calculation type is specified by the ``mode`` parameter: ``scf`` and ``nscf`` are supported as modes, as well as arbitrary modes for generic output.

The content of the output is given in ``content`` section.
The input files of Quantum ESPRESSO are composed of the parts in namelist format of Fortran90 starting from ``&keyword``, and the blocks called cards that start with keywords such as ``K_POINTS`` and end with blank lines. The ``content`` block holds namelist and cards in a form of nested dictionary.
Basically, the specified items are written to the input files as-is, except for several cases. If a keyword is left blank, its value will be obtained form the crystal structure data or other sources.

Besides, templates of the input files can be used. The content of the file given by the ``template`` keyword is considered as input data along with the entries in ``content`` block. When the entries of the same keywords appear both, those of the input parameter files will be used. Therefore, it is possible to use template files and overwrite some entries by the input parameter file as needed.
In the present example, the file (``scf.in_tmpl``) shown below is read as a template, and the entries on cutoff parameters as well as cards of CELL_PARAMETER, ATOMIC_SPECIES, ATOMIC_POSITIONS, K_POINTS are generated from the crystal structure data and pseudo-potential files. It is noted that the values of ``ecutwfc`` and ``ecutrho`` are overwritten by the empty lines.

.. literalinclude:: ../../../../tutorial/cif2x/scf.in_tmpl
   :language: fortran

	      
Generating input files
----------------------------------------------------------------

The program ``cif2x`` is executed with the input parameter file (``input.yaml``) and crystal structure data (``Co3SnS2_nosym.cif``) as follows.

.. code-block:: bash

  $ cif2x -t QE input.yaml Co3SnS2_nosym.cif

The required pseudo-potential files should be placed in the directory ``./pseudo``, and the index file for the pseudo-potential should be prepared as ``./psudo/pp_psl_pbe_rrkjus.csv``.

Run ``cif2x`` and a set of input files for Quantum ESPRESSO will be created. The output file is specified by ``output_file`` parameter of the input parameter file, and stored in the directory given by ``output_dir``. In this example, the input file for SCF calculation is created as ``./scf/scf.in``.


Specifying parameter sets
----------------------------------------------------------------

In some cases, a series of input files should be generated with varying their parameter values. For example, the convergence is examined by modifying the cutoff values or grid resolution of k points. The input parameter can be given a list or a range of values, and the input files for every combination from the choices of parameter values are generated and stored in separate directories. To specify parameter set, a special syntax ``${...}`` is adopted.

.. code-block:: yaml

   content:
     K_POINTS:
       option: automatic
       grid:   ${ [ [4,4,4], [8,8,8], [12,12,12] ] }

When ``K_POINTS`` is given as above, the input files having the ``grid`` value to be ``[4,4,4]``, ``[8,8,8]``, ``[12,12,12]`` will be generated in the sub-directories,  ``4x4x4/``, ``8x8x8/``, ``12x12x12/``, respectively.
