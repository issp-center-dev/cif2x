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
