from cif2x.akaikkr.make_input import make_inputcard


def _params(**extra):
    p = {
        "go": "go", "potentialfile": "pot", "brvtyp": "fcc",
        "ntyp": 1, "type": ["A"], "ncmp": [1], "rmt": [1.0],
        "field": [0.0], "mxl": [2], "anclr": [[1]], "conc": [[1]],
        "natm": 1, "atmicx": [["0", "0", "0", "A"]],
    }
    p.update(extra)
    return p


def test_option_block_emitted_when_present():
    out = make_inputcard(_params(option={"foo": "bar"}))
    assert "begin_option" in out
    assert "foo = bar" in out
    assert "end_option" in out


def test_option_list_value_emitted_as_block():
    out = make_inputcard(_params(option={"orbital_moment": [1, 2, 3]}))
    assert "begin_orbital_moment" in out
    assert "1 2 3" in out
    assert "end_orbital_moment" in out


def test_no_option_block_when_absent():
    out = make_inputcard(_params())
    assert "begin_option" not in out


def test_option_none_does_not_crash_and_emits_no_block():
    # `option:` with no value in YAML parses to None; must not crash
    out = make_inputcard(_params(option=None))
    assert "begin_option" not in out


def test_option_empty_dict_emits_no_block():
    out = make_inputcard(_params(option={}))
    assert "begin_option" not in out


def test_option_block_precedes_terminating_end_record():
    out = make_inputcard(_params(option={"foo": "bar"}))
    # the option block must appear before the closing "end" record
    assert out.index("begin_option") < out.rindex("\n    end\n")


def test_generated_option_block_not_misparsed_by_reader(tmp_path):
    # a generated card reused as a template must not have its option block
    # misread as trailing k-path data (which would inject a spurious kpath)
    from cif2x.akaikkr.read_input import read_input_file

    out = make_inputcard(_params(option={"foo": "bar", "mag": [1, 2, 3]}))
    fn = tmp_path / "x.in"
    fn.write_text(out)
    params = read_input_file(str(fn))
    assert "kpath" not in params


def test_reader_skips_mixed_case_option_markers(tmp_path):
    from cif2x.akaikkr.read_input import read_input_file

    text = (
        "go pot\n"
        "fcc 1.0 1.0 1.0 90 90 90\n"
        "0.01 e 0.0 0.0 0.0 0.0 0.0 0.0 100 0.0 1\n"
        "A 1 1.0 0.0 2 1.0 1.0\n"
        "1\n"
        "0 0 0 A\n"
        "BEGIN_OPTION\n"
        " foo = bar\n"
        "END_OPTION\n"
        "end\n"
    )
    fn = tmp_path / "tmpl.in"
    fn.write_text(text)
    params = read_input_file(str(fn))
    assert "kpath" not in params


def test_reader_raises_on_unterminated_option_block(tmp_path):
    import pytest
    from cif2x.akaikkr.read_input import read_input_file

    out = make_inputcard(_params(option={"foo": "bar"}))
    # drop the end_option line to simulate a malformed template
    broken = out.replace("end_option\n", "")
    fn = tmp_path / "broken.in"
    fn.write_text(broken)
    with pytest.raises(ValueError):
        read_input_file(str(fn))


def test_spc_option_block_does_not_consume_real_kpath(tmp_path):
    # for an spc run, the option block sits before `end` and the real k-path
    # after it; skipping the option block must not eat the trailing k-path
    from cif2x.akaikkr.read_input import read_input_file

    p = _params(go="spc", option={"foo": "bar"},
                kpath=[["0", "0", "0"], ["0.5", "0", "0"]], kdiv=100, fmt=3)
    out = make_inputcard(p)
    fn = tmp_path / "spc.in"
    fn.write_text(out)
    params = read_input_file(str(fn))
    assert params["kdiv"] == 100
    assert params["kpath"] == [["0", "0", "0"], ["0.5", "0", "0"]]


def test_option_block_roundtrip_parsed_back(tmp_path):
    # a generated option block must be parsed back into params["option"] so
    # that the writer->reader round-trip is lossless. Values are returned as
    # strings (the writer stringifies everything).
    from cif2x.akaikkr.read_input import read_input_file

    out = make_inputcard(_params(option={"foo": "bar", "mag": [1, 2, 3]}))
    fn = tmp_path / "rt.in"
    fn.write_text(out)
    params = read_input_file(str(fn))
    assert params["option"] == {"foo": "bar", "mag": ["1", "2", "3"]}


def test_option_block_roundtrip_mixed_case_markers(tmp_path):
    from cif2x.akaikkr.read_input import read_input_file

    out = make_inputcard(_params(option={"foo": "bar", "mag": [1, 2, 3]}))
    out = out.replace("begin_option", "BEGIN_OPTION").replace("end_option", "END_OPTION")
    fn = tmp_path / "rt_mc.in"
    fn.write_text(out)
    params = read_input_file(str(fn))
    assert params["option"] == {"foo": "bar", "mag": ["1", "2", "3"]}


def test_spc_option_and_kpath_roundtrip(tmp_path):
    # an spc input with BOTH an option block (before `end`) and a real numeric
    # k-path (after `end`) must round-trip: option parsed AND kpath/kdiv kept.
    from cif2x.akaikkr.read_input import read_input_file

    p = _params(go="spc", option={"foo": "bar", "mag": [1, 2, 3]},
                kpath=[["0", "0", "0"], ["0.5", "0", "0"]], kdiv=100, fmt=3)
    out = make_inputcard(p)
    fn = tmp_path / "spc_rt.in"
    fn.write_text(out)
    params = read_input_file(str(fn))
    assert params["option"] == {"foo": "bar", "mag": ["1", "2", "3"]}
    assert params["kdiv"] == 100
    assert params["kpath"] == [["0", "0", "0"], ["0.5", "0", "0"]]


def test_option_block_roundtrip_unterminated_still_raises(tmp_path):
    import pytest
    from cif2x.akaikkr.read_input import read_input_file

    out = make_inputcard(_params(option={"foo": "bar"}))
    broken = out.replace("end_option\n", "")
    fn = tmp_path / "broken_rt.in"
    fn.write_text(broken)
    with pytest.raises(ValueError):
        read_input_file(str(fn))


def test_no_option_key_when_block_absent(tmp_path):
    from cif2x.akaikkr.read_input import read_input_file

    out = make_inputcard(_params())
    fn = tmp_path / "noopt.in"
    fn.write_text(out)
    params = read_input_file(str(fn))
    assert "option" not in params


def test_option_subblock_unterminated_raises(tmp_path):
    import pytest
    from cif2x.akaikkr.read_input import read_input_file

    out = make_inputcard(_params(option={"mag": [1, 2, 3]}))
    broken = out.replace(" end_mag\n", "")  # drop the inner end marker
    fn = tmp_path / "bad_inner.in"
    fn.write_text(broken)
    with pytest.raises(ValueError):
        read_input_file(str(fn))


def test_option_duplicate_scalar_key_last_wins():
    # documented policy: a repeated scalar key keeps the last value
    from cif2x.akaikkr.read_input import parse_option_block
    assert parse_option_block(["a = 1", "a = 2"]) == {"a": "2"}


def test_option_empty_list_element_collapses(tmp_path):
    # known limitation: the writer space-joins list values, so empty elements
    # cannot be represented and are dropped on round-trip
    from cif2x.akaikkr.read_input import read_input_file

    out = make_inputcard(_params(option={"mag": ["1", "", "3"]}))
    fn = tmp_path / "empty_el.in"
    fn.write_text(out)
    params = read_input_file(str(fn))
    assert params["option"]["mag"] == ["1", "3"]


def test_empty_option_block_adds_no_option_key(tmp_path):
    import re as _re
    from cif2x.akaikkr.read_input import read_input_file

    out = make_inputcard(_params(option={"foo": "bar"}))
    # collapse the block to an empty begin_option/end_option
    out = _re.sub(r"begin_option\n.*?end_option\n", "begin_option\nend_option\n",
                  out, flags=_re.S)
    fn = tmp_path / "empty_blk.in"
    fn.write_text(out)
    params = read_input_file(str(fn))
    assert "option" not in params


def test_duplicate_top_level_option_block_raises(tmp_path):
    import pytest
    from cif2x.akaikkr.read_input import read_input_file

    out = make_inputcard(_params(option={"foo": "bar"}))
    # inject a second top-level option block
    out = out.replace("end_option\n",
                      "end_option\nbegin_option\n bar = baz\nend_option\n", 1)
    fn = tmp_path / "dup_blk.in"
    fn.write_text(out)
    with pytest.raises(ValueError):
        read_input_file(str(fn))


def test_option_inner_markers_case_insensitive_key_preserved(tmp_path):
    # markers are matched case-insensitively, but the key spelling is preserved
    # verbatim for round-trip fidelity
    from cif2x.akaikkr.read_input import parse_option_block
    parsed = parse_option_block(["BEGIN_MAG", " 1 2 3", "END_MAG"])
    assert parsed == {"MAG": ["1", "2", "3"]}


def test_option_mixed_case_list_key_roundtrips(tmp_path):
    from cif2x.akaikkr.read_input import read_input_file

    out = make_inputcard(_params(option={"Mag": [1, 2, 3]}))
    fn = tmp_path / "mixedkey.in"
    fn.write_text(out)
    params = read_input_file(str(fn))
    assert params["option"] == {"Mag": ["1", "2", "3"]}


def test_option_scalar_key_starting_with_begin_roundtrips(tmp_path):
    # a scalar key that happens to start with "begin_" must be parsed as a
    # key=value assignment, not as a list sub-block opener
    from cif2x.akaikkr.read_input import read_input_file

    out = make_inputcard(_params(option={"begin_mag": "1"}))
    fn = tmp_path / "beginkey.in"
    fn.write_text(out)
    params = read_input_file(str(fn))
    assert params["option"] == {"begin_mag": "1"}


def test_option_inner_key_named_option_roundtrips(tmp_path):
    # an inner list key literally named "option" produces a nested
    # begin_option/end_option pair; the extractor must not close the outer
    # block early
    from cif2x.akaikkr.read_input import read_input_file

    out = make_inputcard(_params(option={"option": [1, 2, 3]}))
    fn = tmp_path / "nested.in"
    fn.write_text(out)
    params = read_input_file(str(fn))
    assert params["option"] == {"option": ["1", "2", "3"]}
    assert "kpath" not in params
