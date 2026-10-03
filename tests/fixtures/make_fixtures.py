"""Regenerate the test fixtures: chain A only, first model, no waters.

    python tests/fixtures/make_fixtures.py

Adenylate kinase (E. coli): 4AKE apo (open), 1AKE holo with Ap5A (closed).
Trypsin (bovine): 2PTN apo, 3PTB holo with benzamidine.
Neither pair is in APObind, so the tests check the method against known
biology without touching the training set.
"""

from pathlib import Path

import gemmi

from flexflag.data.apoholo import structure_path

HERE = Path(__file__).parent

for pdb_id in ["1ake", "4ake", "3ptb", "2ptn"]:
    st = gemmi.read_structure(str(structure_path(pdb_id)))
    while len(st) > 1:
        del st[1]
    for ch in list(st[0]):
        if ch.name != "A":
            st[0].remove_chain(ch.name)
    st.remove_waters()
    st.setup_entities()
    doc = st.make_mmcif_document()
    out = HERE / f"{pdb_id}_A.cif"
    doc.write_file(str(out))
    print(out.name, out.stat().st_size // 1024, "KB")
