import re
import itertools
import json
from collections import UserDict

import logging
logger = logging.getLogger(__name__)


def serializer(data):
    return json.dumps(data)

def deserializer(str):
    return json.loads(str)

def inflate(content):
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
            return list(np.arange(*_a))

        # 2. [ c1, c2, ... ] (-> list) or  c1, c2, c3 (-> tuple)
        try:
            _a = ast.literal_eval(s)
            return _a
        except Exception as e:
            pass

        # 3. [ s1, s2, ... ] where s1, s2, ... are unquoted strings
        is_list = re.search(r"\[(.*?)\]", s)
        if is_list:
            _a = [_to_number(_t) for _t in is_list.group(1).split(",")]
            return _a

        # or parse failed
        raise ValueError("unable to parse item: {}".format(s))

    def _to_number(s):
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
        return s.strip()

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

    # find range keywords
    ss = re.sub("\$\{(.*?)\}", _matcher, ss)

    logger.debug("inflate: tbl={}".format(tbl))
    # logger.debug("inflate: ss={}".format(ss))

    if len(tbl) == 0:
        infos = [("", content)]
        return infos

    tbl_keys = tbl.keys()
    tbl_vals = tbl.values()

    # expand range lists
    infos = []
    for x in itertools.product(*tbl_vals):
        _map = { k: v for k, v in zip(tbl_keys, x) }
        # logger.debug("inflate: _map={}".format(_map))

        ss_new = ss
        for from_val, to_val in _map.items():
            ss_new = ss_new.replace(serializer(from_val), serializer(to_val))

        # logger.debug("inflate: ss_new={}".format(ss_new))

        content_new = type(content).deserialize(ss_new)
        key = "_".join([_to_string(_) for _ in x])

        # logger.debug(f"key={key}, content_new={content_new.cards}")

        infos += [(key, content_new)]

    return infos


class CaseInsensitiveKey(object):
    def __init__(self, key):
        # logger.debug("Key: __init__ with type {}, value {}".format(type(key),key))
        if isinstance(key, CaseInsensitiveKey):
            self.key = key.key
        else:
            self.key = key
    def __hash__(self):
        # logger.debug("Key: __hash__ with type {}, value {}".format(type(self.key),self.key))
        return hash(self.key.casefold())
    def __eq__(self, other):
        # logger.debug("Key: __eq__ self={}, other={}({})".format(self.key, other, type(other)))
        return self.key.casefold() == other.key.casefold()
    def __str__(self):
        return self.key

class CaseInsensitiveDict(UserDict):
    def __init__(self, mapping=None, /, **kwargs):
        if mapping is not None:
            # logger.debug("Dict: __init__ with {}".format(mapping))
            mapping = { CaseInsensitiveKey(key): value for key, value in mapping.items() }
        else:
            mapping = {}
        if kwargs:
            # logger.debug("Dict: __init__ with kwargs={}".format(kwargs))
            mapping.update({ CaseInsensitiveKey(key): value for key, value in kwargs.items() })
        super().__init__(mapping)

    def __setitem__(self, key, value):
        # logger.debug("Dict: __setitem__ with key={}, value={}".format(key, value))
        super().__setitem__(CaseInsensitiveKey(key), value)

    def __getitem__(self, key):
        value = super().__getitem__(CaseInsensitiveKey(key))
        # logger.debug("Dict: __getitem__ with key={}, value={}".format(key, value))
        return value

    def __contains__(self, key):
        # logger.debug("Dict: __contains__ {} {}".format(key, type(key)))
        return super().__contains__(CaseInsensitiveKey(key))

    def __str__(self):
        s = "{ " + ", ".join([ str(k) + ": " + str(v) for k, v in super().items()]) + " }"
        return s

    def as_dict(self):
        return { str(k): v for k, v in super().items() }

