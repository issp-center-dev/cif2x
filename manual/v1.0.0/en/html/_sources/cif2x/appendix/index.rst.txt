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
