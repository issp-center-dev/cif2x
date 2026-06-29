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
    raw = _read_raw_lines(input_file)
    option_lines, rest_lines = _extract_option_block(raw)

    tokens = _tokenize_lines(rest_lines)

    params,rest = parse_readin(tokens)
    logger.debug(f"read_input_file: params={params}, rest={rest}")
    if rest:
        fmt,kpath,kdiv = parse_readk(rest)
        params.update({"kpath": kpath, "kdiv": kdiv, "fmt": fmt})

    option = parse_option_block(option_lines)
    if option:
        params["option"] = option

    return params

def _read_raw_lines(input_file):
    with open(input_file, "r") as fp:
        return [line.rstrip() for line in fp.readlines() if line[0] not in ["c","C","#"]]

def _extract_option_block(raw):
    # Separate the option block (begin_option ... end_option) from the rest of
    # the input. The markers are case-insensitive. The inner lines are returned
    # separately so they can be parsed back into params["option"] instead of
    # being misread as trailing k-path data. Returns (option_lines, rest_lines)
    # where option_lines is None when no block is present.
    option_lines = None
    rest_lines = []
    depth = 0  # nesting level; an inner option key literally named "option"
               # produces a nested begin_option/end_option pair
    for line in raw:
        token = line.strip().lower()
        if token == "begin_option":
            if depth == 0:
                if option_lines is not None:
                    logger.error("multiple option blocks")
                    raise ValueError("multiple option blocks")
                option_lines = []
            else:
                option_lines.append(line)  # inner sub-block content
            depth += 1
            continue
        if depth > 0:
            if token == "end_option":
                depth -= 1
                if depth > 0:
                    option_lines.append(line)  # inner sub-block content
            else:
                option_lines.append(line)
            continue
        rest_lines.append(line)

    if depth > 0:
        logger.error("unterminated option block (missing end_option)")
        raise ValueError("unterminated option block (missing end_option)")

    return option_lines, rest_lines

def _tokenize_lines(lines):
    st = " ".join(lines)
    st = st.strip()

    tokens = re.split(r"(?:\s*,\s*|\s+)", st)
    return tokens

def tokenize_input_file(input_file):
    raw = _read_raw_lines(input_file)
    _, rest_lines = _extract_option_block(raw)
    return _tokenize_lines(rest_lines)

def parse_option_block(lines):
    # Reconstruct the option block emitted by make_input._make_option_card.
    # `key = value` lines become option[key] = value (a string), and
    # `begin_<name> ... end_<name>` sub-blocks become option[name] = [values],
    # where each value is a whitespace-separated token (a string). The writer
    # stringifies everything, so all parsed values are returned as strings.
    # Markers are matched case-insensitively. Returns None when there is no
    # option block.
    #
    # Contract / known limits:
    # - duplicate keys / sub-blocks: last one wins (matches the writer, which is
    #   driven by a dict).
    # - the writer space-joins list values, so empty/whitespace-only elements are
    #   not representable and are dropped on round-trip.
    # - an unterminated inner `begin_<name>` sub-block raises ValueError, like the
    #   outer block.
    if lines is None:
        return None

    option = {}
    idx = 0
    n = len(lines)
    while idx < n:
        token = lines[idx].strip()
        if not token:
            idx += 1
            continue
        if "=" in token:
            # scalar assignment; checked before the begin_ branch so a scalar
            # key that happens to start with "begin_" is not misread as a
            # list sub-block opener.
            key, value = token.split("=", 1)
            option[key.strip()] = value.strip()
            idx += 1
        elif token.lower().startswith("begin_"):
            # preserve the key spelling as written (round-trip fidelity); the
            # matching end_<name> marker is compared case-insensitively.
            name = token[len("begin_"):].strip()
            end_marker = "end_" + name.lower()
            idx += 1
            values = []
            found_end = False
            while idx < n:
                if lines[idx].strip().lower() == end_marker:
                    found_end = True
                    idx += 1
                    break
                values.extend(lines[idx].split())
                idx += 1
            if not found_end:
                logger.error(f"unterminated option sub-block (missing {end_marker})")
                raise ValueError(f"unterminated option sub-block (missing {end_marker})")
            option[name] = values
        else:
            idx += 1

    return option

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
    # The k-path block is a contiguous run of fixed-width numeric records.
    # A record is consumed only when a full set of non-empty fields is available
    # and the leading field starts with a digit; otherwise parsing stops. We stop
    # (rather than skip-and-resync) on a malformed/empty record because a missing
    # field makes the remaining record boundaries ambiguous to re-align.
    ntoken = len(tokens)
    idx = 0

    if fmt == 1:
        kval = []
        while idx + 2 < ntoken and all(tokens[idx+j] for j in range(3)) and tokens[idx][0].isdigit():
            v1 = tokens[idx].lower(); idx += 1
            v2 = tokens[idx].lower(); idx += 1
            v3 = tokens[idx].lower(); idx += 1
            kval += [[v1,v2,v3]]
        kdiv = 0

    elif fmt == 2:
        kval = []
        while idx + 3 < ntoken and all(tokens[idx+j] for j in range(4)) and tokens[idx][0].isdigit():
            v1 = tokens[idx].lower(); idx += 1
            v2 = tokens[idx].lower(); idx += 1
            v3 = tokens[idx].lower(); idx += 1
            nv = _to_int(tokens[idx]); idx += 1
            kval += [[v1,v2,v3,nv]]
        kdiv = 0

    elif fmt == 3:
        if idx < ntoken:
            kdiv = _to_int(tokens[idx], 0); idx += 1
        else:
            kdiv = 0
        kval = []
        while idx + 2 < ntoken and all(tokens[idx+j] for j in range(3)) and tokens[idx][0].isdigit():
            v1 = tokens[idx].lower(); idx += 1
            v2 = tokens[idx].lower(); idx += 1
            v3 = tokens[idx].lower(); idx += 1
            kval += [[v1,v2,v3]]

    else:
        pass

    return fmt, kval, kdiv

