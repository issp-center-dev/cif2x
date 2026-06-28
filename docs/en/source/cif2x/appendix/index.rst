.. _sec-cif2x-appendix:

================================================================
Extension guide
================================================================

Adding modes of Quantum ESPRESSO
----------------------------------------------------------------

In order to add supports to modes of Quantum ESPRESSO, the mapping between the modes and the transformation classes should be added to ``create_modeproc()`` function in ``src/cif2x/qe/calc_mode.py``.

.. code-block:: python

  def create_modeproc(mode, qe):
      if mode in ["scf", "nscf"]:
          modeproc = QEmode_pw(qe)
      else:
          modeproc = QEmode_generic(qe)
      return modeproc
	

The transformation functionality for each mode is provided as a derived class of ``QEmode_base`` class. This class implements methods ``update_namelist()`` for updating the namelist block, and ``update_cards()`` for generating data of card blocks.
In the current version, two classes are provided: ``QEmode_pw`` class for scf and nscf calculations of pw.x, and ``QEmode_generic`` class for generating output as-is.

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


Troubleshooting
----------------------------------------------------------------

Common errors that may occur when running ``cif2x`` and how to resolve them are summarized below.

- **mode not specified (RuntimeError)**

  This is raised for a Quantum ESPRESSO run when an element of ``tasks`` does not have ``mode``. Specify ``mode`` (such as ``scf`` or ``nscf``) for each task.

- **output_file not specified (RuntimeError)**

  This is raised for a Quantum ESPRESSO run when an element of ``tasks`` does not have ``output_file``. Specify the output file name by ``output_file`` for each task.

- **Warnings: pp_file / cutoff_file / pseudo_dir not specified or not found**

  A warning is emitted (and the run continues) when these parameters are not given in the ``optional`` section, or when the specified path does not exist. Check that the paths of the pseudo-potential index file (``pp_file``), the cutoff index file (``cutoff_file``), and the pseudo-potential directory (``pseudo_dir``) are correct. Note that ``cutoff_file`` / ``pseudo_dir`` are warning-and-continue, but ``pp_file`` is effectively required for standard QE generation: ``ATOMIC_SPECIES`` generation (and automatic ``nbnd``) raises a ``RuntimeError`` later in the run if it is not specified.

- **Cutoff values become 0.0**

  When ``ecutwfc`` / ``ecutrho`` are left blank to be obtained automatically, the cutoff value falls back to ``0.0`` if — given a valid ``pp_file`` mapping for the element — the corresponding ``.UPF`` file is not found in ``pseudo_dir`` and the ``cutoff_file`` has no matching entry. Check that the pseudo-potential files and the cutoff index entries are all in place. (Note that a missing element row in ``pp_file`` itself is a hard error (``KeyError``), not a ``0.0`` fallback.)
