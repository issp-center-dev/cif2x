import json

from cif2x.utils import inflate as utils_inflate
from cif2x.qe.content import Content, inflate as qe_inflate


class DictContent:
    """Minimal stand-in for the Content objects inflate() operates on."""

    def __init__(self, data):
        self.data = data

    def serialize(self):
        return json.dumps(self.data)

    @classmethod
    def deserialize(cls, s):
        return cls(json.loads(s))


def test_utils_inflate_expands_list_placeholder():
    results = utils_inflate(DictContent({"x": "${[30, 40]}"}))
    assert len(results) == 2
    assert sorted(c.data["x"] for _, c in results) == [30, 40]


def test_utils_inflate_expands_range_placeholder():
    results = utils_inflate(DictContent({"x": "${range(3)}"}))
    assert len(results) == 3
    assert sorted(c.data["x"] for _, c in results) == [0, 1, 2]


def test_utils_inflate_without_placeholder_returns_single():
    results = utils_inflate(DictContent({"x": 5}))
    assert len(results) == 1
    key, content = results[0]
    assert key == ""
    assert content.data["x"] == 5


def test_utils_inflate_two_placeholders_cartesian_product():
    results = utils_inflate(DictContent({"a": "${[1, 2]}", "b": "${range(2)}"}))
    combos = sorted((c.data["a"], c.data["b"]) for _, c in results)
    assert combos == [(1, 0), (1, 1), (2, 0), (2, 1)]


def test_qe_inflate_expands_list_placeholder():
    c = Content()
    c.namelist = {"system": {"ecutwfc": "${[30, 40]}"}}
    c.cards = None
    c.textblock = None

    results = qe_inflate(c)
    assert len(results) == 2
    assert sorted(cc.namelist["system"]["ecutwfc"] for _, cc in results) == [30, 40]


def test_qe_inflate_expands_range_placeholder():
    c = Content()
    c.namelist = {"system": {"ecutwfc": "${range(30, 50, 10)}"}}
    c.cards = None
    c.textblock = None

    results = qe_inflate(c)
    assert len(results) == 2
    assert sorted(cc.namelist["system"]["ecutwfc"] for _, cc in results) == [30, 40]
