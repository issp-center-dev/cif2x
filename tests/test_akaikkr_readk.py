from cif2x.akaikkr.read_input import parse_readk


def test_parse_readk_fmt3_reads_kdiv_and_kpoints():
    # fmt 3: leading kdiv, followed by k-point triplets
    fmt, kval, kdiv = parse_readk(["100", "0", "0", "0"], fmt=3)
    assert fmt == 3
    assert kdiv == 100
    assert kval == [["0", "0", "0"]]


def test_parse_readk_fmt3_empty_tokens():
    fmt, kval, kdiv = parse_readk([], fmt=3)
    assert kdiv == 0
    assert kval == []


def test_parse_readk_fmt3_ignores_truncated_triplet():
    # kdiv present, then an incomplete k-point (only 2 of 3 coords) must not crash
    fmt, kval, kdiv = parse_readk(["100", "0", "0"], fmt=3)
    assert kdiv == 100
    assert kval == []


def test_parse_readk_fmt3_empty_token_does_not_crash():
    # double comma in input (e.g. "100,,0 0 0") yields an empty token after split
    fmt, kval, kdiv = parse_readk(["100", "", "0", "0", "0"], fmt=3)
    assert kdiv == 100
    assert kval == []


def test_parse_readk_skips_records_with_empty_interior_field():
    # an empty field anywhere in a record (double comma) marks malformed data; skip it
    assert parse_readk(["0", "", "0"], fmt=1)[1] == []
    assert parse_readk(["0", "0", "", "1"], fmt=2)[1] == []
    assert parse_readk(["100", "0", "", "0"], fmt=3)[1] == []


def test_parse_readk_fmt3_malformed_kdiv_defaults_to_zero():
    fmt, kval, kdiv = parse_readk(["", "0", "0", "0"], fmt=3)
    assert kdiv == 0


def test_parse_readk_reads_multiple_valid_records():
    # several well-formed records in a row must all be returned
    assert parse_readk(["0", "0", "0", "1", "1", "1"], fmt=1)[1] == [
        ["0", "0", "0"], ["1", "1", "1"]]
    assert parse_readk(["0", "0", "0", "2", "1", "1", "1", "3"], fmt=2)[1] == [
        ["0", "0", "0", 2], ["1", "1", "1", 3]]


def test_parse_readk_stops_at_non_numeric_terminator():
    # a trailing terminator token ("end") cleanly ends the k-path block
    fmt, kval, kdiv = parse_readk(["5", "0", "0", "0", "end"], fmt=3)
    assert kdiv == 5
    assert kval == [["0", "0", "0"]]


def test_parse_readk_fmt1_exact_and_truncated():
    fmt, kval, kdiv = parse_readk(["0", "0", "0"], fmt=1)
    assert kval == [["0", "0", "0"]]
    # truncated tail (< 3 coords) must be ignored, not crash
    fmt, kval, kdiv = parse_readk(["0", "0"], fmt=1)
    assert kval == []


def test_parse_readk_fmt2_exact_and_truncated():
    fmt, kval, kdiv = parse_readk(["0", "0", "0", "1"], fmt=2)
    assert kval == [["0", "0", "0", 1]]
    # truncated tail (< 4 tokens) must be ignored, not crash
    fmt, kval, kdiv = parse_readk(["0", "0", "0"], fmt=2)
    assert kval == []
