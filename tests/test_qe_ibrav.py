import numpy as np
import pytest

pytest.importorskip("qe_tools")
pytest.importorskip("pymatgen")

from pymatgen.core import Lattice, Structure
from cif2x.struct2qe import Struct2QE


def _set_ibrav(structure):
    obj = Struct2QE.__new__(Struct2QE)
    return obj._set_ibrav_structure(structure)


@pytest.mark.parametrize("params", [
    (4.0, 5.0, 6.0, 70.0, 80.0, 100.0),
    (3.0, 4.0, 5.0, 85.0, 95.0, 110.0),
    (5.0, 5.5, 6.5, 100.0, 110.0, 75.0),
])
def test_ibrav14_triclinic_lattice_is_consistent(params):
    # a generic triclinic cell maps to ibrav=14; the rebuilt lattice must be
    # crystallographically equivalent to the input (internal matches() check),
    # which only holds when v2 = b*(cos gamma, sin gamma, 0).
    latt = Lattice.from_parameters(*params)
    structure = Structure(latt, ["Si"], [[0.0, 0.0, 0.0]])

    new_structure, system = _set_ibrav(structure)
    assert system["ibrav"] == 14
    # the second cell vector's y-component must be b*sin(gamma), not b*cos(gamma)
    b = system["B"]
    gamma = np.radians(latt.gamma)
    v2 = new_structure.lattice.matrix[1]
    assert np.isclose(v2[1], b * np.sin(gamma), atol=1e-6)
    assert not np.isclose(v2[1], b * np.cos(gamma), atol=1e-3)
