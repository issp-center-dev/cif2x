import re
import itertools
import json

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
