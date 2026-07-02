import pytest

pytest.importorskip("bs4")

from utils.pp_cutoff import read_pseudo_cutoff


def _write_upf(path, header_attrs="", pp_info=""):
    body = ""
    if pp_info:
        body += f"<PP_INFO>\n{pp_info}\n</PP_INFO>"
    body += f"<PP_HEADER {header_attrs}/>"
    path.write_text(f"<UPF>{body}</UPF>")


def test_zero_header_cutoff_treated_as_missing(tmp_path):
    # ld1.x-generated UPFs carry rho_cutoff="0.0" when unset: the CSV row must
    # get a blank cell (which cif2x reads as "not specified"), not a literal
    # 0.0 that the generator would then reject as "not found".
    upf = tmp_path / "Fe.UPF"
    _write_upf(upf, header_attrs='wfc_cutoff="40.0" rho_cutoff="0.0"')
    assert read_pseudo_cutoff(str(upf)) == [str(upf), 40.0, ""]


def test_missing_wfc_cutoff_returns_none(tmp_path):
    upf = tmp_path / "Fe.UPF"
    _write_upf(upf, header_attrs='rho_cutoff="320.0"')
    assert read_pseudo_cutoff(str(upf)) is None


def test_pp_info_fills_only_missing_fields(tmp_path):
    # A valid header wfc_cutoff must survive a pp_info block whose
    # "Suggested ... wavefunctions" line carries an unset 0.0.
    upf = tmp_path / "Fe.UPF"
    _write_upf(
        upf,
        header_attrs='wfc_cutoff="40.0"',
        pp_info=(" Suggested minimum cutoff for wavefunctions:    0.0 Ry\n"
                 " Suggested minimum cutoff for charge density: 320.0 Ry"),
    )
    assert read_pseudo_cutoff(str(upf)) == [str(upf), 40.0, 320.0]
