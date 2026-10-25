import logging
logger = logging.getLogger(__name__)

def deepupdate(dict1, dict2):
    """
    merge dict2 into dict1; update nested dictionary recursively
    """
    for k,v in dict2.items():
        if isinstance(v, dict) and isinstance(dict1.get(k), dict):
            deepupdate(dict1[k], v)
        else:
            # a mapping replaces a scalar/None (e.g. "system:" left blank in a
            # template-less content) instead of descending into it
            dict1[k] = v

def is_empty_key(tbl, key):
    return key in tbl and tbl[key] is None

