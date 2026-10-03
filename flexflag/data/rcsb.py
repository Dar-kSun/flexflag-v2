"""RCSB search and data APIs, for building the external PDB-2019+ set."""

import json

import requests

from flexflag.config import CACHE_DIR, HTTP_TIMEOUT_S

SEARCH = "https://search.rcsb.org/rcsbsearch/v2/query"
GRAPHQL = "https://data.rcsb.org/graphql"


def _term(attribute: str, operator: str, value) -> dict:
    return {
        "type": "terminal",
        "service": "text",
        "parameters": {"attribute": attribute, "operator": operator, "value": value},
    }


def search_entries(name: str, terms: list[tuple]) -> list[str]:
    """All entry IDs matching every (attribute, operator, value) term. Cached by name."""
    path = CACHE_DIR / "rcsb" / f"search_{name}.json"
    if path.exists():
        return json.loads(path.read_text())
    query = {
        "query": {"type": "group", "logical_operator": "and", "nodes": [_term(*t) for t in terms]},
        "return_type": "entry",
        "request_options": {"return_all_hits": True},
    }
    resp = requests.post(SEARCH, json=query, timeout=300)
    resp.raise_for_status()
    ids = sorted(r["identifier"].lower() for r in resp.json()["result_set"])
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(ids))
    return ids


ENTRY_QUERY = """query($ids: [String!]!) { entries(entry_ids: $ids) {
  rcsb_id
  rcsb_entry_info { resolution_combined }
  rcsb_accession_info { initial_release_date }
  nonpolymer_entities { nonpolymer_comp { chem_comp { id formula_weight } } }
} }"""


def entry_details(ids: list[str], batch: int = 300) -> dict[str, dict]:
    """Resolution, release date and non-polymer components per entry. Cached per batch."""
    out = {}
    for start in range(0, len(ids), batch):
        chunk = ids[start : start + batch]
        path = CACHE_DIR / "rcsb" / f"entries_{chunk[0]}_{chunk[-1]}_{len(chunk)}.json"
        if path.exists():
            data = json.loads(path.read_text())
        else:
            resp = requests.post(
                GRAPHQL,
                json={"query": ENTRY_QUERY, "variables": {"ids": [i.upper() for i in chunk]}},
                timeout=HTTP_TIMEOUT_S * 2,
            )
            resp.raise_for_status()
            data = resp.json()["data"]["entries"]
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text(json.dumps(data))
        for e in data:
            res = (e["rcsb_entry_info"] or {}).get("resolution_combined") or [None]
            chems = [n["nonpolymer_comp"]["chem_comp"] for n in e["nonpolymer_entities"] or []]
            comps = [(c["id"], c["formula_weight"]) for c in chems]
            out[e["rcsb_id"].lower()] = {
                "resolution": res[0],
                "release": e["rcsb_accession_info"]["initial_release_date"][:10],
                "ligands": comps,
            }
    return out
