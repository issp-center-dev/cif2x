import glob
import os

import pytest
from ruamel.yaml import YAML

from cif2x.input_validator import normalize_target, validate_input

_SAMPLE_GLOB = os.path.join(
    os.path.dirname(__file__), "..", "sample", "cif2x", "*", "*", "input.yaml"
)
_SAMPLES = sorted(glob.glob(_SAMPLE_GLOB))


def test_samples_were_found():
    assert _SAMPLES, f"no sample inputs matched {_SAMPLE_GLOB}"


@pytest.mark.parametrize("path", _SAMPLES)
def test_sample_input_validates(path):
    # sample/cif2x/<target_dir>/<case>/input.yaml -> target_dir is index -3
    target_dir = path.split(os.sep)[-3]
    target = normalize_target(target_dir)
    with open(path) as fp:
        info_dict = YAML(typ="safe").load(fp)
    validate_input(info_dict, target)  # must not raise
