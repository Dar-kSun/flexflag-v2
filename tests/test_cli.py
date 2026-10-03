"""`flexflag check`, offline, using committed AlphaFold DB fixtures."""

from pathlib import Path

import pytest
from typer.testing import CliRunner

import flexflag.data.alphafold as alphafold
from flexflag.check import RULE, check, parse_residues
from flexflag.cli import app

FIX = Path(__file__).parent / "fixtures"


@pytest.fixture(autouse=True)
def offline_afdb(monkeypatch):
    monkeypatch.setattr(
        alphafold, "fetch", lambda url, subdir, name, retries=3: FIX / "afdb" / name
    )


def test_parse_residues():
    assert parse_residues("13, 31,35-38,31") == [13, 31, 35, 36, 37, 38]


def test_band_matches_rule_edges():
    rep = check("P00698", list(range(30, 60)))
    edges = RULE["band_edges"]
    expected = sum(rep.pocket_plddt > e for e in edges)
    assert rep.band == expected
    assert 0 < rep.probability < 1


def test_lower_pocket_plddt_means_higher_probability():
    # Lysozyme's first residues (signal peptide) have low pLDDT.
    low = check("P00698", list(range(1, 15)))
    high = check("P00698", list(range(60, 80)))
    assert low.pocket_plddt < high.pocket_plddt
    assert low.probability > high.probability


def test_no_pocket_gives_no_flag():
    out = CliRunner().invoke(app, ["check", "P00698"])
    assert out.exit_code == 0
    assert "No pocket given, so no flag" in out.output


def test_residues_out_of_range_is_an_error():
    out = CliRunner().invoke(app, ["check", "P00698", "--residues", "140-160"])
    assert out.exit_code == 1


def test_pocket_from_local_pdb_file():
    # Adenylate kinase with Ap5A: pocket taken from the committed holo fixture.
    out = CliRunner().invoke(app, ["check", "P69441", "--from-pdb", str(FIX / "1ake_A.cif")])
    assert out.exit_code == 0, out.output
    assert "AP5 in 1ake_A chain A" in out.output
    assert "Pocket-change risk" in out.output


def test_wrong_protein_pdb_is_rejected():
    # Trypsin structure given for adenylate kinase: must not produce a report.
    out = CliRunner().invoke(app, ["check", "P69441", "--from-pdb", str(FIX / "3ptb_A.cif")])
    assert out.exit_code == 1
    assert "no chain in 3ptb_A matches P69441" in out.output


def test_pdb_without_ligand_is_a_clean_error():
    out = CliRunner().invoke(app, ["check", "P69441", "--from-pdb", str(FIX / "4ake_A.cif")])
    assert out.exit_code == 1
    assert "no small-molecule ligand" in out.output


def test_unreadable_residues_are_a_clean_error():
    out = CliRunner().invoke(app, ["check", "P00698", "--residues", "12,abc"])
    assert out.exit_code == 1
    assert "cannot read residue 'abc'" in out.output
