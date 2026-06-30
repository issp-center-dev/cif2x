.. _sec-cif2x-appendix:

================================================================
Extension guide
================================================================

Adding modes of Quantum ESPRESSO
----------------------------------------------------------------

In order to add supports to modes of Quantum ESPRESSO, the mapping between the modes and the transformation classes should be added to ``create_modeproc()`` function in ``src/cif2x/qe/calc_mode.py``.

.. code-block:: python

  def create_modeproc(mode, qe):
      if mode in ["scf", "nscf", "relax", "vc-relax", "bands"]:
          modeproc = QEmode_pw(qe)
      else:
          modeproc = QEmode_generic(qe)
      return modeproc
	

The transformation functionality for each mode is provided as a derived class of ``QEmode_base`` class. This class implements methods ``update_namelist()`` for updating the namelist block, and ``update_cards()`` for generating data of card blocks.
In the current version, two classes are provided: ``QEmode_pw`` class for scf, nscf, relax, vc-relax, and bands calculations of pw.x, and ``QEmode_generic`` class for generating output as-is.

.. code-block:: python

  class QEmode_base:
      def __init__(self, qe):
      def update_namelist(self, content):
      def update_cards(self, content):


For the namelist, the transformation class generates values for blank entries from crystal structure data and other sources. It may also force to set values such as the lattice parameters that are determined from the crystal structure data, or those that must be specified consistently with other parameters. The functions are provided for each mode separately.

For card blocks, a function is provided for each card, and the mapping between the card type and the function is given in the ``card_table`` variable.
The method ``update_cards()`` in the base class picks up and runs the function associated to the card, and updates the content of the card. Of course, a new ``update_cards()`` function may be defined.

.. code-block:: python

    self.card_table = {
        'CELL_PARAMETERS': generate_cell_parameters,
        'ATOMIC_SPECIES': generate_atomic_species,
        'ATOMIC_POSITIONS': generate_atomic_positions,
        'K_POINTS': generate_k_points,
    }

The functions for cards are gathered in ``src/cif2x/qe/cards.py`` with the function names as ``generate_{card name}``. These functions takes parameters for card blocks as argument, and returns a dictionary containing the card name, the options, and the data field.


Supported calculation modes
----------------------------------------------------------------

``QEmode_pw`` generates the structure cards (``CELL_PARAMETERS``,
``ATOMIC_SPECIES``, ``ATOMIC_POSITIONS``, ``K_POINTS``) and sets ``calculation``
from the task ``mode`` for ``scf``, ``nscf``, ``relax``, ``vc-relax``, and
``bands``. The ``&ions`` (for ``relax``/``vc-relax``) and ``&cell`` (for
``vc-relax``) namelists are taken from the template/``content`` as written —
make sure they are present, as ``pw.x`` errors without them.

For ``bands``, set the K_POINTS option to ``crystal_b`` in ``content``; the path
is generated from the crystal symmetry with ``pymatgen``'s high-symmetry k-path,
for example:

.. code-block:: yaml

    tasks:
      - mode: bands
        output_file: bands.in
        content:
          system:
            nbnd: 16        # set explicitly if conduction bands are wanted
          K_POINTS:
            option: crystal_b
            line_npoints: 30   # points generated between high-symmetry points

``line_npoints`` defaults to 20. A high-symmetry path may contain breaks
(disjoint segments); the generator marks them with a 0 in the per-line integer
so ``bands.x`` does not interpolate across the gap. To use a custom path, give a
``path`` list of label sequences (labels must exist in the auto-generated
k-path).

Quantum ESPRESSO computes only the occupied (valence) bands by default, so set
``nbnd`` explicitly in ``content.system`` to include conduction bands in the band
structure.

The high-symmetry coordinates are defined in the **standardized primitive cell**,
so ``bands`` requires the structure to be that cell — use ``use_primitive: true``
(and ``use_ibrav: false``). If the input is a conventional cell, a supercell, or
otherwise non-standard, ``cif2x`` logs a warning ("band path: ... the path may be
incorrect") and the sampled path will not match the written ``CELL_PARAMETERS``.

There is no ``calculation='dos'`` in pw.x. A density-of-states run is the
chain ``scf`` -> ``nscf`` (a dense mesh) -> ``dos`` (a ``dos.x`` input). The
``dos`` task is emitted as-is via ``QEmode_generic``; put the ``&DOS`` namelist
in its template/``content``:

.. code-block:: yaml

    tasks:
      - mode: scf
        output_file: scf.in
        content: { K_POINTS: { option: automatic, grid: [8, 8, 8] } }
      - mode: nscf
        output_file: nscf.in
        content: { K_POINTS: { option: automatic, grid: [16, 16, 16] } }
      - mode: dos
        template: dos.in_tmpl     # contains the &DOS namelist
        output_file: dos.in

Run the generated inputs in order — ``scf``, then ``nscf``, then ``dos.x`` — and
give all three the same ``prefix`` and ``outdir`` so each step can read the
previous step's saved data.


Troubleshooting
----------------------------------------------------------------

Common errors that may occur when running ``cif2x`` and how to resolve them are summarized below.

``cif2x`` validates the target and ``input.yaml`` before generating any files. When the target name or the input is invalid, it logs a single error message and exits with status ``1`` (no Python traceback).

- **unsupported target '...'**

  Reported when the ``-t``/``--target`` value is not one of ``quantum_espresso`` (``qe``, ``espresso``), ``vasp``, ``openmx``, or ``akaikkr``. The message lists the valid choices.

- **task N: 'mode' is required for target 'quantum_espresso'**

  Reported for a Quantum ESPRESSO run when an element of ``tasks`` does not have ``mode``. Specify ``mode`` (such as ``scf`` or ``nscf``) for each task.

- **task N: 'output_file' is required for target '...'**

  Reported when a task that needs an output file name does not have ``output_file``. Specify the output file name by ``output_file`` for each task. (Required for ``quantum_espresso``, ``openmx``, and ``akaikkr``.)

- **task N: unknown key '...' for target '...'**

  Reported when a task contains an unrecognized key, for example a typo such as ``output_fil``. The message lists the keys allowed for the target. The free-form blocks (``content``, ``optional``, ``structure``) are not key-checked.

- **Warnings: pp_file / cutoff_file / pseudo_dir not specified or not found**

  A warning is emitted (and the run continues) when these parameters are not given in the ``optional`` section, or when the specified path does not exist. Check that the paths of the pseudo-potential index file (``pp_file``), the cutoff index file (``cutoff_file``), and the pseudo-potential directory (``pseudo_dir``) are correct. Note that ``cutoff_file`` / ``pseudo_dir`` are warning-and-continue, but ``pp_file`` is effectively required for standard QE generation: ``ATOMIC_SPECIES`` generation (and automatic ``nbnd``) raises a ``RuntimeError`` later in the run if it is not specified.

- **Cutoff values become 0.0**

  When ``ecutwfc`` / ``ecutrho`` are left blank to be obtained automatically, the cutoff value falls back to ``0.0`` if — given a valid ``pp_file`` mapping for the element — the corresponding ``.UPF`` file is not found in ``pseudo_dir`` and the ``cutoff_file`` has no matching entry. Check that the pseudo-potential files and the cutoff index entries are all in place. (Note that a missing element row in ``pp_file`` itself is a hard error (``KeyError``), not a ``0.0`` fallback.)
