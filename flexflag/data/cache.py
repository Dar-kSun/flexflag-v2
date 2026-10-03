"""Download-once file cache. Every network fetch in flexflag goes through here."""

import time
from pathlib import Path

import requests

from flexflag.config import CACHE_DIR, HTTP_TIMEOUT_S


class NotFound(Exception):
    """The remote resource does not exist (HTTP 404)."""


def fetch(url: str, subdir: str, name: str, retries: int = 3) -> Path:
    """Return a local path for `url`, downloading it on first use."""
    path = CACHE_DIR / subdir / name
    if path.exists():
        return path
    path.parent.mkdir(parents=True, exist_ok=True)
    for attempt in range(retries):
        try:
            resp = requests.get(url, timeout=HTTP_TIMEOUT_S)
        except requests.RequestException:
            if attempt == retries - 1:
                raise
            time.sleep(2**attempt)
            continue
        if resp.status_code == 404:
            raise NotFound(url)
        if resp.status_code >= 500 and attempt < retries - 1:
            time.sleep(2**attempt)
            continue
        resp.raise_for_status()
        tmp = path.with_suffix(path.suffix + ".part")
        tmp.write_bytes(resp.content)
        tmp.replace(path)
        return path
    raise RuntimeError(f"unreachable: {url}")
