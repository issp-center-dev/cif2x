import logging
logger = logging.getLogger("make_input")

def make_inputcard(params):

    def _fill_title(s):
        return s.center(60, "-")

    def _fill_value(v, fmt="{}"):
        if v is None or v == "":
            _add_comma = not _fill_value._is_prev_empty
            _fill_value._is_prev_empty = True
            if _add_comma:
                return ",  ,"
            else:
                return ","
        else:
            _fill_value._is_prev_empty = False
            return fmt.format(v)
    _fill_value._is_prev_empty = True

    def _add_title(s):
        return [ "c" + s.center(60, "-") ]

    def _add_separator():
        return [ "c" + ("-" * 60) ]

    def _comment_line(s):
        return "c   " + s

    def _value_line(vv):
        return "    " + "  ".join(vv)

    def _add_values(params, keys, header=True):
        s = []
        if header:
            s.append(_comment_line(" ".join(keys)))
        s.append(_value_line([_fill_value(params.get(k)) for k in keys]))
        return s

    def _make_lattice_vectors(params):
        s = []
        fmt = "{:12.8f}"

        for r in ["r1", "r2", "r3"]:
            vv = params.get(r, [None, None, None])
            s.append(_value_line([_fill_value(v, fmt) for v in vv]))

        if "tlt" in params["brvtyp"]:
            angl = params.get("angl", [None, None, None])
            s.append(_value_line([_fill_value(v, fmt) for v in angl]))

        s.append(_value_line([_fill_value(params.get("a", None), fmt)]))

        return s

    def _make_type_card(params):
        s = []

        s.append(_comment_line(" ".join(["typ", "ncmp", "rmt", "field", "mxl", "[anclr conc]"])))

        ntyp = params.get("ntyp", 0)
        for i in range(ntyp):
            typ   = params["type"][i]
            ncmp  = params["ncmp"][i]
            rmt   = params["rmt"][i]
            field = params["field"][i]
            mxl   = params["mxl"][i]

            s.append(_value_line([_fill_value(v) for v in [typ, ncmp, rmt, field, mxl]]))

            for j in range(ncmp):
                anclr = params["anclr"][i][j]
                conc  = params["conc"][i][j]
                s.append(_value_line([" " * 24] + [_fill_value(v) for v in [anclr, conc]]))

        return s

    def _make_atom_card(params):
        s = []

        s.append(_comment_line(" ".join(["atmicx", "atmtyp"])))

        natm = params.get("natm", 0)
        for i in range(natm):
            atmicx = params["atmicx"][i]
            s.append(_value_line([_fill_value(v) for v in atmicx]))

        return s

    def _make_option_card(params):
        s = []

        if "option" in params:
            options = params["option"]
            s.append("")
            s.append("begin_option")
            for k, v in options.items():
                if isinstance(v, list):
                    v = list(map(str, v))
                    s.append(" begin_{}".format(k))
                    s.append(" " + " ".join(v))
                    s.append(" end_{}".format(k))
                else:
                    v = str(v)
                    s.append(" {} = {}".format(k, v))
            s.append("end_option")
            s.append("")

        return s

    def _make_kpath_card(params):
        s = []

        if "kpath" in params:
            fmt = params.get("fmt", 3)
            if fmt == 3:
                s.append(_value_line([_fill_value(params.get("kdiv", 0))]))
            for kline in params["kpath"]:
                s.append(_value_line([_fill_value(v) for v in kline]))

        return s

    title = params.get("title", "")

    retv = []
    retv += _add_title(title)
    retv += _add_values(params, ["go", "potentialfile"], header=False)
    retv += _add_separator()

    if "aux" in params["brvtyp"] or "prv" in params["brvtyp"]:
        retv += _add_values(params, ["brvtyp"])
        retv += _make_lattice_vectors(params)
    else:
        retv += _add_values(params, ["brvtyp", "a", "c/a", "b/a", "alpha", "beta", "gamma"])

    retv += _add_separator()
    retv += _add_values(params, ["edelt", "ewidth", "reltyp", "sdftyp", "magtyp", "record"])
    retv += _add_separator()
    retv += _add_values(params, ["outtyp", "bzqlty", "maxitr", "pmix"])
    retv += _add_separator()
    retv += _add_values(params, ["ntyp"])
    retv += _add_separator()
    retv += _make_type_card(params)
    retv += _add_separator()
    retv += _add_values(params, ["natm"])
    retv += _add_separator()
    retv += _make_atom_card(params)

    if False:
        retv += _add_separator()
        retv += _make_option_card(params)

    retv += _add_title(title)
    retv += _add_values({"term": "end"}, ["term"], header=False)

    if "spc" in params["go"]:
        retv += _make_kpath_card(params)

    result = "\n".join(retv) + "\n"
    return result
