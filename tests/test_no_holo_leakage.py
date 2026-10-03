"""The feature pipeline must never read a PDB entry, the holo structure or the ligand."""

import ast
from pathlib import Path

import gemmi
import pytest

import flexflag.data.alphafold as alphafold
import flexflag.data.apoholo as apoholo
import flexflag.features as features

AFDB_FIX = Path(__file__).parent / "fixtures" / "afdb"
FEATURE_DIR = Path(features.__file__).parent
FORBIDDEN_IMPORTS = {"gemmi", "flexflag.labels", "flexflag.data.apoholo"}


def test_feature_modules_do_not_import_structure_code():
    for path in FEATURE_DIR.glob("*.py"):
        for node in ast.walk(ast.parse(path.read_text())):
            names = []
            if isinstance(node, ast.Import):
                names = [a.name for a in node.names]
            elif isinstance(node, ast.ImportFrom) and node.module:
                names = [node.module] + [f"{node.module}.{a.name}" for a in node.names]
            for name in names:
                assert not any(name == f or name.startswith(f + ".") for f in FORBIDDEN_IMPORTS), (
                    f"{path.name} imports {name}"
                )


def test_features_only_touch_alphafold_db(monkeypatch):
    urls = []

    def fake_fetch(url, subdir, name, retries=3):
        urls.append(url)
        return AFDB_FIX / name

    def trap(*args, **kwargs):
        raise AssertionError("feature pipeline tried to read a PDB structure")

    monkeypatch.setattr(alphafold, "fetch", fake_fetch)
    monkeypatch.setattr(apoholo, "fetch", trap)
    monkeypatch.setattr(apoholo, "structure_path", trap)
    monkeypatch.setattr(gemmi, "read_structure", trap)

    feats = features.compute_features("P00698")

    assert urls and all(u.startswith("https://alphafold.ebi.ac.uk/") for u in urls)
    assert feats["seq_length"] == 147
    assert feats["plddt_mean"] == pytest.approx(93.88, abs=0.05)  # AFDB's own global value
