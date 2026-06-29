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


from cif2x.input_validator import validate_input


def _qe_ok():
    return {"tasks": [{"mode": "scf", "output_file": "scf.in"}]}


def test_validate_ok_qe_does_not_raise():
    validate_input(_qe_ok(), "quantum_espresso")


def test_validate_none_input():
    with pytest.raises(InputValidationError, match="empty"):
        validate_input(None, "vasp")


def test_validate_top_level_not_mapping():
    with pytest.raises(InputValidationError, match="mapping at the top level"):
        validate_input([1, 2], "vasp")


def test_validate_unknown_top_level_key():
    with pytest.raises(InputValidationError, match="unknown top-level key"):
        validate_input({"taskz": [{}]}, "vasp")


def test_validate_tasks_missing():
    with pytest.raises(InputValidationError, match="'tasks' is required"):
        validate_input({"structure": {}}, "vasp")


def test_validate_tasks_not_list():
    with pytest.raises(InputValidationError, match="'tasks' must be a list"):
        validate_input({"tasks": {"mode": "scf"}}, "quantum_espresso")


def test_validate_tasks_empty():
    with pytest.raises(InputValidationError, match="empty"):
        validate_input({"tasks": []}, "vasp")


def test_validate_task_not_mapping():
    with pytest.raises(InputValidationError, match="must be a mapping"):
        validate_input({"tasks": [None]}, "vasp")


def test_validate_qe_missing_mode():
    with pytest.raises(InputValidationError, match="'mode' is required"):
        validate_input({"tasks": [{"output_file": "scf.in"}]}, "quantum_espresso")


def test_validate_qe_missing_output_file():
    with pytest.raises(InputValidationError, match="'output_file' is required"):
        validate_input({"tasks": [{"mode": "scf"}]}, "quantum_espresso")


@pytest.mark.parametrize("val", [None, ""])
def test_validate_output_file_null_or_empty(val):
    with pytest.raises(InputValidationError, match="output_file"):
        validate_input(
            {"tasks": [{"mode": "scf", "output_file": val}]}, "quantum_espresso"
        )


def test_validate_unknown_task_key():
    with pytest.raises(InputValidationError, match="unknown key 'output_fil'"):
        validate_input(
            {"tasks": [{"mode": "scf", "output_file": "x", "output_fil": "y"}]},
            "quantum_espresso",
        )


def test_validate_output_dir_must_be_string():
    with pytest.raises(InputValidationError, match="output_dir"):
        validate_input(
            {"tasks": [{"mode": "scf", "output_file": "x", "output_dir": 3}]},
            "quantum_espresso",
        )


def test_validate_free_form_blocks_pass():
    d = {
        "structure": {"use_ibrav": True, "tol_deg": 5},
        "optional": {"pseudo_dir": "/x", "pp_file": "pp.csv",
                     "pseudo_map": {"Na": "x"}},
        "tasks": [{"mode": "scf", "output_file": "scf.in",
                   "content": {"anything": 1}}],
    }
    validate_input(d, "quantum_espresso")


def test_validate_vasp_requires_nothing():
    validate_input({"tasks": [{"template_dir": "base"}]}, "vasp")


def test_validate_structure_block_must_be_mapping():
    with pytest.raises(InputValidationError, match="'structure' must be a mapping"):
        validate_input({"structure": [1], "tasks": [{"template_dir": "b"}]}, "vasp")


def test_validate_empty_optional_and_structure_blocks_pass():
    # `optional:` / `structure:` with nothing after them parse to None and are
    # accepted as "no entries" (common in samples).
    d = {"structure": None, "optional": None,
         "tasks": [{"output_file": "test.in", "template": "t"}]}
    validate_input(d, "akaikkr")
