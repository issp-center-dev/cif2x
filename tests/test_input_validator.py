import pytest

from cif2x.input_validator import (
    InputValidationError,
    normalize_target,
)


@pytest.mark.parametrize("alias,canonical", [
    ("qe", "quantum_espresso"),
    ("QE", "quantum_espresso"),
    ("espresso", "quantum_espresso"),
    ("quantum_espresso", "quantum_espresso"),
    ("vasp", "vasp"),
    ("VASP", "vasp"),
    ("openmx", "openmx"),
    ("akaikkr", "akaikkr"),
])
def test_normalize_target_aliases(alias, canonical):
    assert normalize_target(alias) == canonical


def test_normalize_target_unknown_lists_choices():
    with pytest.raises(InputValidationError) as ei:
        normalize_target("qe2")
    msg = str(ei.value)
    assert "qe2" in msg
    assert "quantum_espresso" in msg
    assert "vasp" in msg and "openmx" in msg and "akaikkr" in msg
