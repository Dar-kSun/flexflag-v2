from pathlib import Path

import gemmi
import pytest

from flexflag.labels import LabelError, align_indices, site_rmsd

FIX = Path(__file__).parent / "fixtures"
AK_HOLO, AK_APO = FIX / "1ake_A.cif", FIX / "4ake_A.cif"
TRYP_HOLO, TRYP_APO = FIX / "3ptb_A.cif", FIX / "2ptn_A.cif"


def test_known_flexible_protein_scores_high():
    # Adenylate kinase: LID and NMP domains close over Ap5A.
    r = site_rmsd(AK_HOLO, ["A"], AK_APO, ["A"])
    assert r.ligand.startswith("AP5")
    assert r.large_change
    assert r.rmsd == pytest.approx(6.1, abs=0.3)


def test_known_rigid_protein_scores_low():
    # Trypsin: the S1 pocket barely moves when benzamidine binds.
    r = site_rmsd(TRYP_HOLO, ["A"], TRYP_APO, ["A"])
    assert r.ligand.startswith("BEN")  # the calcium ion is skipped, too small
    assert not r.large_change
    assert r.rmsd < 0.5


def test_structure_against_itself_is_zero():
    r = site_rmsd(AK_HOLO, ["A"], AK_HOLO, ["A"])
    assert r.rmsd == pytest.approx(0.0, abs=1e-6)
    assert r.n_matched == r.n_site


def test_renumbered_and_renamed_apo_gives_same_rmsd(tmp_path):
    # Apo and holo entries often number residues differently and use other chain IDs.
    st = gemmi.read_structure(str(AK_APO))
    for ch in st[0]:
        ch.name = "Z"
        for res in ch:
            res.seqid = gemmi.SeqId(res.seqid.num + 100, res.seqid.icode)
    st.setup_entities()
    moved = tmp_path / "renumbered.cif"
    st.make_mmcif_document().write_file(str(moved))

    base = site_rmsd(AK_HOLO, ["A"], AK_APO, ["A"])
    shifted = site_rmsd(AK_HOLO, ["A"], moved, ["Z"])
    assert shifted.apo_chain == "Z"
    assert shifted.rmsd == pytest.approx(base.rmsd, abs=1e-6)


def test_missing_apo_chain_raises():
    with pytest.raises(LabelError):
        site_rmsd(AK_HOLO, ["A"], AK_APO, ["Q"])


def test_apo_structure_has_no_ligand():
    with pytest.raises(LabelError, match="no small-molecule ligand"):
        site_rmsd(AK_APO, ["A"], AK_HOLO, ["A"])


def test_align_indices_skips_gaps_and_mismatches():
    holo = list("MKVLAAGG")
    apo = list("MKVAAGWG")  # L deleted; G->W substitution
    m = align_indices(holo, apo)
    assert m[0] == 0 and m[2] == 2  # M, V
    assert 3 not in m  # L has no partner
    assert m[4] == 3  # A shifts left by one
    assert all(holo[i] == apo[j] for i, j in m.items())


def test_residue_inside_a_peptide_chain_is_not_a_ligand():
    # 5BTR: sirtuin-1 with a Fluor-de-Lys substrate peptide (FDL is a residue of the
    # peptide chain D) and resveratrol (STL). The peptide is out of scope.
    from flexflag.labels import find_ligand, read_model

    ligand, _, chain = find_ligand(read_model(FIX / "5btr_AD.cif"), ["A"])
    assert ligand.name == "STL"
    assert chain == "A"
