import re
import logging

logger = logging.getLogger(__name__)

def _to_float(s, default=None):
    try:
        t = float(s)
        return t
    except ValueError:
        return default

def _to_int(s, default=None):
    try:
        t = int(s)
        return t
    except ValueError:
        return default

def read_input_file(input_file):
    tokens = tokenize_input_file(input_file)

    params,rest = parse_readin(tokens)
    logger.debug(f"read_input_file: params={params}, rest={rest}")
    if rest:
        fmt,kpath,kdiv = parse_readk(rest)
        params.update({"kpath": kpath, "kdiv": kdiv, "fmt": fmt})

    return params

def tokenize_input_file(input_file):
    with open(input_file, "r") as fp:
        lines = [line.rstrip() for line in fp.readlines() if line[0] not in ["c","C","#"]]
    st = " ".join(lines)
    st = st.strip()

    tokens = re.split(r"(?:\s*,\s*|\s+)", st)
    return tokens

def parse_readin(tokens):
    params = {}
    idx = 0

    params["go"]       = tokens[idx].lower(); idx += 1
    params["potentialfile"]  = tokens[idx]; idx += 1
    params["brvtyp"]   = tokens[idx].lower(); idx += 1

    if "prv" in params["brvtyp"] or "aux" in params["brvtyp"]:
        r11 = _to_float(tokens[idx], 1.0); idx += 1
        r21 = _to_float(tokens[idx], 0.0); idx += 1
        r31 = _to_float(tokens[idx], 0.0); idx += 1

        r12 = _to_float(tokens[idx], 0.0); idx += 1
        r22 = _to_float(tokens[idx], 1.0); idx += 1
        r32 = _to_float(tokens[idx], 0.0); idx += 1

        r13 = _to_float(tokens[idx], 0.0); idx += 1
        r23 = _to_float(tokens[idx], 0.0); idx += 1
        r33 = _to_float(tokens[idx], 1.0); idx += 1

        params["r1"] = [r11, r21, r31]
        params["r2"] = [r12, r22, r32]
        params["r3"] = [r13, r23, r33]

    if "tlt" in params["brvtyp"]:
        angl1 = _to_float(tokens[idx], 0.0); idx += 1
        angl2 = _to_float(tokens[idx], 0.0); idx += 1
        angl3 = _to_float(tokens[idx], 0.0); idx += 1

        params["angl"] = [angl1, angl2, angl3]

    params["a"] = _to_float(tokens[idx]); idx += 1

    if not ("prv" in params["brvtyp"] or "aux" in params["brvtyp"]):
        params["c/a"]    = _to_float(tokens[idx]); idx += 1
        params["b/a"]    = _to_float(tokens[idx]); idx += 1
        params["alpha"]  = _to_float(tokens[idx]); idx += 1
        params["beta"]   = _to_float(tokens[idx]); idx += 1
        params["gamma"]  = _to_float(tokens[idx]); idx += 1

    params["edelt"]    = _to_float(tokens[idx]); idx += 1
    params["ewidth"]   = tokens[idx].lower(); idx += 1

    params["reltyp"]   = tokens[idx].lower(); idx += 1
    params["sdftyp"]   = tokens[idx].lower(); idx += 1
    params["magtyp"]   = tokens[idx].lower(); idx += 1
    params["record"]   = tokens[idx].lower(); idx += 1
    params["outtyp"]   = tokens[idx].lower(); idx += 1
    params["bzqlty"]   = tokens[idx].lower(); idx += 1

    params["maxitr"]   = _to_int(tokens[idx]); idx += 1
    params["pmix"]     = tokens[idx].lower(); idx += 1
    params["ntyp"]     = _to_int(tokens[idx]); idx += 1

    params["type"]     = []
    params["ncmp"]     = []
    params["rmt"]      = []
    params["field"]    = []
    params["mxl"]      = []

    params["anclr"]    = []
    params["conc"]     = []

    for i in range(params["ntyp"]):
        typ = tokens[idx]; idx += 1
        ncmp = _to_int(tokens[idx]); idx += 1
        rmt = _to_float(tokens[idx]); idx += 1
        field = _to_float(tokens[idx]); idx += 1
        mxl = _to_int(tokens[idx]); idx += 1

        anclr = []
        conc = []
        for j in range(ncmp):
            anclr.append(_to_float(tokens[idx])); idx += 1
            conc.append(_to_float(tokens[idx])); idx += 1

        params["type"].append(typ)
        params["ncmp"].append(ncmp)
        params["rmt"].append(rmt)
        params["field"].append(field)
        params["mxl"].append(mxl)

        params["anclr"].append(anclr)
        params["conc"].append(conc)

    params["natm"] = _to_int(tokens[idx]); idx += 1
    params["atmicx"] = []

    for i in range(params["natm"]):
        atmicx1 = tokens[idx].lower(); idx += 1
        atmicx2 = tokens[idx].lower(); idx += 1
        atmicx3 = tokens[idx].lower(); idx += 1
        atmtyp  = tokens[idx]; idx += 1
        params["atmicx"].append([atmicx1, atmicx2, atmicx3, atmtyp])

    if idx >= len(tokens):
        logger.debug(f"end of items {idx}/{len(tokens)}")
    elif tokens[idx].lower() == "end":
        idx += 1
        logger.debug(f"end mark detected")
    else:
        logger.debug("remains:", tokens[idx:])

    rest = tokens[idx:]
    return params, rest

def parse_readk(tokens, fmt=3):
    ntoken = len(tokens)
    idx = 0

    if fmt == 1:
        kval = []
        while idx < ntoken and tokens[idx][0].isdigit():
            v1 = tokens[idx].lower(); idx += 1
            v2 = tokens[idx].lower(); idx += 1
            v3 = tokens[idx].lower(); idx += 1
            kval += [[v1,v2,v3]]
        kdiv = 0

    elif fmt == 2:
        kval = []
        while idx < ntoken and tokens[idx][0].isdigit():
            v1 = tokens[idx].lower(); idx += 1
            v2 = tokens[idx].lower(); idx += 1
            v3 = tokens[idx].lower(); idx += 1
            nv = _to_int(tokens[idx]); idx += 1
            kval += [[v1,v2,v3,nv]]
        kdiv = 0

    elif fmt == 3:
        if idx < ntokens:
            kdiv = _to_int(tokens[idx]); idx += 1
        else:
            kdiv = 0
        kval = []
        while idx < ntoken and tokens[idx][0].isdigit():
            v1 = tokens[idx].lower(); idx += 1
            v2 = tokens[idx].lower(); idx += 1
            v3 = tokens[idx].lower(); idx += 1
            kval += [[v1,v2,v3]]

    else:
        pass

    return fmt, kval, kdiv

