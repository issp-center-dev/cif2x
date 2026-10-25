import re
import ast
import sys
import math
import itertools
import json
from collections import UserDict

import numpy as np

from cif2x.input_validator import InputValidationError

import logging
logger = logging.getLogger(__name__)


def dryrun_emit(path, text, stream=None):
    """Print a generated input file (path header + content) for --dry-run."""
    stream = stream if stream is not None else sys.stdout
    stream.write("# === {} ===\n".format(path))
    stream.write(text)
    if not text.endswith("\n"):
        stream.write("\n")


def serializer(data):
    return json.dumps(data)

def deserializer(str):
    return json.loads(str)


def _range_values(args):
    """Values of ``${range(...)}``: start up to but excluding end, in steps.

    Computed as start + i*step with the values rounded, so that a float step
    neither includes ``end`` through rounding error (as ``numpy.arange`` can)
    nor produces values such as 0.30000000000000004 in directory names.
    """
    if len(args) == 1:
        start, end, step = 0, args[0], 1
    elif len(args) == 2:
        start, end, step = args[0], args[1], 1
    elif len(args) == 3:
        start, end, step = args
    else:
        raise InputValidationError(
            "parameter set: range() takes 1 to 3 numbers, got {}".format(len(args)))
    if not all(isinstance(a, (int, float)) for a in (start, end, step)):
        raise InputValidationError(
            "parameter set: range() arguments must be numbers, got {!r}".format(args))
    if step == 0:
        raise InputValidationError("parameter set: range() step must not be 0")
    n = math.ceil((end - start) / step - 1.0e-9)
    all_int = all(isinstance(a, int) for a in (start, end, step))
    values = []
    for i in range(max(n, 0)):
        v = start + i * step
        values.append(int(v) if all_int else float(round(v, 10)))
    return values


def inflate(content):
    """Expand ``${...}`` parameter sets in ``content`` into one content per combination.

    ``content`` must provide ``serialize()`` and a classmethod ``deserialize()``
    (JSON round trip). Returns ``[(key, content), ...]`` where ``key`` names the
    sub-directory of the combination (``""`` when nothing is expanded).
    """
    tbl = {}

    def _matcher(m):
        count = len(tbl) + 1
        key = f"${count}"
        matched = m.group(1)
        tbl[key] = _item_parser(matched)
        return key

    def _item_parser(s):
        # 1. range(n) or range(start, end [,step])
        is_range = re.search(r"range\((.*?)\)", s)
        if is_range:
            _a = [_to_number(_t) for _t in is_range.group(1).split(",")]
            return _range_values(_a)

        # 2. [ c1, c2, ... ] (-> list) or  c1, c2, c3 (-> tuple)
        try:
            _a = ast.literal_eval(s)
            if isinstance(_a, (list, tuple)):
                return list(_a)
            return [_a]
        except Exception as e:
            pass

        # 3. [ s1, s2, ... ] where s1, s2, ... are unquoted strings
        is_list = re.search(r"\[(.*?)\]", s)
        if is_list:
            _a = [_to_number(_t) for _t in is_list.group(1).split(",")]
            return _a

        # or parse failed
        raise InputValidationError(
            "parameter set: cannot parse ${{{}}}; use ${{[a, b, ...]}} or "
            "${{range(start, end, step)}}".format(s))

    def _to_number(s):
        s = s.strip()
        try:
            v = int(s)
            return v
        except Exception as e:
            pass
        try:
            v = float(s)
            return v
        except Exception as e:
            pass
        return s

    def _to_string(x):
        if isinstance(x, (int, float, str)):
            return str(x)
        elif isinstance(x, list):
            return "x".join([_to_string(t) for t in x])
        elif isinstance(x, dict):
            return "_".join([_to_string(k) + "-" + _to_string(v) for k,v in x.items()])
        else:
            return str(x)

    # flatten into single string
    ss = content.serialize()

    # find parameter sets that form a whole (JSON string) value: "${...}"
    ss = re.sub(r'"\$\{(.*?)\}"', lambda m: json.dumps(_matcher(m)), ss)

    # a ${...} left over sits inside a longer string, which is not supported
    leftover = re.search(r"\$\{(.*?)\}", ss)
    if leftover:
        raise InputValidationError(
            "parameter set: ${{{}}} must be the whole value of a parameter, "
            "not part of a longer string.".format(leftover.group(1)))

    logger.debug("inflate: tbl={}".format(tbl))

    if len(tbl) == 0:
        infos = [("", content)]
        return infos

    tbl_keys = tbl.keys()
    tbl_vals = tbl.values()

    # expand range lists
    infos = []
    for x in itertools.product(*tbl_vals):
        _map = { k: v for k, v in zip(tbl_keys, x) }

        ss_new = ss
        for from_val, to_val in _map.items():
            ss_new = ss_new.replace(serializer(from_val), serializer(to_val))

        content_new = type(content).deserialize(ss_new)
        key = "_".join([_to_string(_) for _ in x])

        infos += [(key, content_new)]

    return infos


class CaseInsensitiveKey(object):
    def __init__(self, key):
        if isinstance(key, CaseInsensitiveKey):
            self.key = key.key
        else:
            self.key = key
    def __hash__(self):
        return hash(self.key.casefold())
    def __eq__(self, other):
        return self.key.casefold() == other.key.casefold()
    def __str__(self):
        return self.key

class CaseInsensitiveDict(UserDict):
    def __init__(self, mapping=None, /, **kwargs):
        if mapping is not None:
            mapping = { CaseInsensitiveKey(key): value for key, value in mapping.items() }
        else:
            mapping = {}
        if kwargs:
            mapping.update({ CaseInsensitiveKey(key): value for key, value in kwargs.items() })
        super().__init__(mapping)

    def __setitem__(self, key, value):
        super().__setitem__(CaseInsensitiveKey(key), value)

    def __getitem__(self, key):
        value = super().__getitem__(CaseInsensitiveKey(key))
        return value

    def __contains__(self, key):
        return super().__contains__(CaseInsensitiveKey(key))

    def __str__(self):
        s = "{ " + ", ".join([ str(k) + ": " + str(v) for k, v in super().items()]) + " }"
        return s

    def as_dict(self):
        return { str(k): v for k, v in super().items() }
