"""Download Hetionet v1.0, the source of every dataset in this project.

Hetionet (Himmelstein et al., eLife 2017; CC0) integrates 29 public resources
into one network keyed by Entrez gene IDs. This project uses:

* ``GiG`` - 147,164 human protein-protein interactions (HI-II-14, the
  Menche et al. 2015 interactome and other binary/literature sources);
* ``DaG`` - disease-gene associations (GWAS Catalog, DISEASES, DisGeNET, DOAF);
* ``DuG`` / ``DdG`` - genes up/down-regulated in the disease (STARGEO);
* ``GpPW`` / ``GpBP`` - pathway (Reactome, WikiPathways, Pathway Interaction
  Database) and GO biological-process membership;
* ``CtD`` / ``CbG`` - approved drugs treating the disease and the proteins
  every DrugBank compound binds.

Both files are pinned to a fixed commit of github.com/hetio/hetionet and
verified with SHA-256.

Usage::

    python src/data_download.py
"""

from __future__ import annotations

import hashlib
import logging
import os
import shutil
import sys
import tempfile
import time
import urllib.request
from dataclasses import dataclass
from pathlib import Path

if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from src.config import load_config

logger = logging.getLogger(__name__)

HETIONET_COMMIT = "8a6cc0c5604d0e2908786e4de38af67b8e46ee4a"


@dataclass(frozen=True)
class RemoteFile:
    filename: str
    url: str
    sha256: str
    description: str


HETIONET_FILES = {
    "nodes": RemoteFile(
        "hetionet-v1.0-nodes.tsv",
        f"https://raw.githubusercontent.com/hetio/hetionet/{HETIONET_COMMIT}/hetnet/tsv/hetionet-v1.0-nodes.tsv",
        "3ad9eb86c031df8bc7d5efefa6edc93562182eb9a6a71e5129d3f18b9d42c50e",
        "47,031 nodes (genes, diseases, compounds, pathways, GO terms, ...)",
    ),
    "edges": RemoteFile(
        "hetionet-v1.0-edges.sif.gz",
        # Large file stored with Git LFS, served from the media host.
        f"https://media.githubusercontent.com/media/hetio/hetionet/{HETIONET_COMMIT}/hetnet/tsv/hetionet-v1.0-edges.sif.gz",
        "611f0411ac666be4e0270e54bbca67407e14720bfc3099311f0f5eec25c5a947",
        "2,250,197 typed edges (source, metaedge, target)",
    ),
}


def sha256sum(path: str | Path) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _fetch(url: str, destination: Path, retries: int = 3, timeout: float = 120) -> None:
    for attempt in range(1, retries + 1):
        try:
            with urllib.request.urlopen(url, timeout=timeout) as response, open(destination, "wb") as out:
                shutil.copyfileobj(response, out)
            return
        except OSError as exc:
            if attempt == retries:
                raise
            logger.warning("Download failed (%s); retry %d/%d", exc, attempt, retries)
            time.sleep(2**attempt)


def download_file(remote: RemoteFile, raw_dir: str | Path, force: bool = False) -> Path:
    """Download ``remote`` into ``raw_dir`` unless a verified copy already exists."""
    target = Path(raw_dir) / remote.filename
    if target.exists() and not force and sha256sum(target) == remote.sha256:
        return target
    target.parent.mkdir(parents=True, exist_ok=True)
    logger.info("Downloading %s", remote.url)
    fd, tmp_name = tempfile.mkstemp(dir=target.parent, suffix=".part")
    os.close(fd)
    tmp = Path(tmp_name)
    try:
        _fetch(remote.url, tmp)
        actual = sha256sum(tmp)
        if actual != remote.sha256:
            raise ValueError(f"Checksum mismatch for {remote.filename}: {actual}")
        tmp.replace(target)
    finally:
        tmp.unlink(missing_ok=True)
    return target


def download_hetionet(raw_dir: str | Path | None = None, force: bool = False) -> dict[str, Path]:
    """Download (or reuse) the Hetionet node and edge files."""
    raw_dir = raw_dir or load_config()["data"]["raw_dir"]
    return {key: download_file(remote, raw_dir, force) for key, remote in HETIONET_FILES.items()}


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(name)s: %(message)s")
    for name, path in download_hetionet(force="--force" in sys.argv).items():
        print(f"{name:6s} {path}  ({path.stat().st_size / 1e6:.1f} MB)")
