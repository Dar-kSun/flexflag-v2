"""E1 of docs/plan-v0.5-overnight.md, step 1: ESM-2 per-residue embeddings, plus the
pocket residues each protein's features are taken over.

    python scripts/11_esm_embed.py [--limit N]

Needs torch and fair-esm (`pip install torch fair-esm`). Uses the GPU if present.
Writes data/cache/esm/<UniProt>.npy (float16, length x 1280) and
results/robustness/pockets.json (AlphaFold residue numbers of each 5 Å pocket).
"""

import argparse
import importlib.util
import json
import os
import time
from concurrent.futures import ThreadPoolExecutor

import numpy as np
import pandas as pd

from flexflag.config import CACHE_DIR, ROOT
from flexflag.evaluate import RESULTS

ESM_DIR = CACHE_DIR / "esm"
POCKETS = RESULTS / "robustness" / "pockets.json"
WINDOW, STRIDE = 1000, 500  # ESM-2 takes at most 1022 residues; longer ones are tiled


def datasets() -> pd.DataFrame:
    a = pd.read_csv(RESULTS / "dataset.csv").assign(dataset="apobind")
    e = pd.read_csv(RESULTS / "external" / "dataset.csv").assign(dataset="external")
    return pd.concat([a, e], ignore_index=True)


def write_pockets(ds: pd.DataFrame) -> None:
    if POCKETS.exists():
        return
    spec = importlib.util.spec_from_file_location("p2", ROOT / "scripts" / "10_p2rank.py")
    p2 = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(p2)

    def one(row):
        try:
            return f"{row.dataset}|{row.uniprot}", sorted(p2.true_pocket(row))
        except Exception:  # counted as missing downstream
            return f"{row.dataset}|{row.uniprot}", None

    with ThreadPoolExecutor(16) as pool:
        pockets = dict(pool.map(one, ds.itertuples()))
    POCKETS.parent.mkdir(parents=True, exist_ok=True)
    POCKETS.write_text(json.dumps(pockets))


def embed_all(ds: pd.DataFrame, limit: int | None) -> None:
    # fair-esm checkpoints are full pickles; torch >= 2.6 refuses them by default.
    os.environ.setdefault("TORCH_FORCE_NO_WEIGHTS_ONLY_LOAD", "1")
    # On a laptop the display shares GPU memory; avoid fragmentation as lengths vary.
    os.environ.setdefault("PYTORCH_CUDA_ALLOC_CONF", "expandable_segments:True")
    import esm
    import torch

    device = "cuda" if torch.cuda.is_available() else "cpu"
    print(f"ESM-2 on {device}" + (f" ({torch.cuda.get_device_name()})" if device == "cuda" else ""))
    model, alphabet = esm.pretrained.esm2_t33_650M_UR50D()
    model = model.eval().to(device)
    if device == "cuda":
        model = model.half()
    convert = alphabet.get_batch_converter()
    cpu_model = None  # float32 copy, made only if a sequence does not fit on the GPU

    def represent(tokens):
        nonlocal cpu_model
        with torch.no_grad():
            if device == "cuda":
                try:
                    return model(tokens.to(device), repr_layers=[33])["representations"][33]
                except torch.cuda.OutOfMemoryError:
                    torch.cuda.empty_cache()
            if cpu_model is None:
                cpu_model = esm.pretrained.esm2_t33_650M_UR50D()[0].eval()
            return cpu_model(tokens, repr_layers=[33])["representations"][33]

    seqs = ds.drop_duplicates("uniprot").set_index("uniprot").sequence
    todo = [
        u
        for u in seqs.sort_values(key=lambda s: s.str.len()).index
        if not (ESM_DIR / f"{u}.npy").exists()
    ]
    if limit:
        todo = todo[:limit]
    ESM_DIR.mkdir(parents=True, exist_ok=True)
    t0 = time.perf_counter()
    for k, u in enumerate(todo, 1):
        seq = seqs[u]
        total = np.zeros((len(seq), 1280), dtype=np.float32)
        count = np.zeros(len(seq), dtype=np.float32)
        starts = (
            list(range(0, max(1, len(seq) - WINDOW + STRIDE), STRIDE)) if len(seq) > WINDOW else [0]
        )
        for s in starts:
            chunk = seq[s : s + WINDOW]
            _, _, tokens = convert([(u, chunk)])
            rep = represent(tokens)
            total[s : s + len(chunk)] += rep[0, 1 : len(chunk) + 1].float().cpu().numpy()
            count[s : s + len(chunk)] += 1
            del rep
        np.save(ESM_DIR / f"{u}.npy", (total / count[:, None]).astype(np.float16))
        if device == "cuda":
            torch.cuda.empty_cache()
        if k % 25 == 0 or k == len(todo):
            print(
                f"{k}/{len(todo)} embedded (last: {len(seq)} aa), {time.perf_counter() - t0:.0f}s",
                flush=True,
            )


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--limit", type=int)
    args = ap.parse_args()
    ds = datasets()
    write_pockets(ds)
    embed_all(ds, args.limit)


if __name__ == "__main__":
    main()
