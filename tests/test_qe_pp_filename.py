import pytest

pytest.importorskip("qe_tools")
pd = pytest.importorskip("pandas")

from cif2x.struct2qe import Struct2QE


def _qe(is_soc):
    qe = Struct2QE.__new__(Struct2QE)
    qe.pp_list = pd.DataFrame({"pseudopotential": ["pbe-spn-rrkjus_psl.1.0.0"]},
                              index=["Fe"])
    qe.is_soc = is_soc
    return qe


def test_pp_filename_non_soc():
    assert _qe(False)._pp_filename("Fe") == "Fe.pbe-spn-rrkjus_psl.1.0.0.UPF"


def test_pp_filename_soc_has_rel_prefix():
    assert _qe(True)._pp_filename("Fe") == "Fe.rel-pbe-spn-rrkjus_psl.1.0.0.UPF"


def test_cutoff_table_lookup_uses_soc_filename():
    # the cutoff CSV is keyed by the actual (rel-) pseudopotential filename for
    # SOC runs; the lookup must use the same name the generator writes
    qe = _qe(True)
    qe.cutoff_list = pd.DataFrame(
        {"ecutwfc": [40.0], "ecutrho": [320.0]},
        index=["Fe.rel-pbe-spn-rrkjus_psl.1.0.0.UPF"],
    )
    ecutwfc, ecutrho = qe._find_elem_cutoff_from_table("Fe")
    assert ecutwfc == 40.0
    assert ecutrho == 320.0


def test_cutoff_file_lookup_uses_soc_filename(tmp_path):
    pytest.importorskip("bs4")
    qe = _qe(True)
    qe.pseudo_dir = str(tmp_path)
    upf = tmp_path / "Fe.rel-pbe-spn-rrkjus_psl.1.0.0.UPF"
    upf.write_text('<UPF><PP_HEADER wfc_cutoff="40.0" rho_cutoff="320.0"/></UPF>')

    ecutwfc, ecutrho = qe._find_elem_cutoff_from_file("Fe")
    assert ecutwfc == 40.0
    assert ecutrho == 320.0
