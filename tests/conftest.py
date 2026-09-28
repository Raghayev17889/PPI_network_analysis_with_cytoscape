"""Shared fixtures: a miniature Hetionet with a toy disease."""

from __future__ import annotations

import matplotlib

matplotlib.use("Agg")

import pandas as pd
import pytest

NODES = [
    ("Disease::DOID:0001", "toy disease", "Disease"),
    ("Compound::DB0001", "Drugamab", "Compound"),
    ("Compound::DB0002", "Otherol", "Compound"),
    ("Pathway::PW1", "Toy signalling", "Pathway"),
    ("Pathway::PW2", "Toy metabolism", "Pathway"),
    ("Biological Process::GO:1", "toy process", "Biological Process"),
] + [(f"Gene::{i}", f"G{i}", "Gene") for i in range(1, 11)]

EDGES = [
    ("Gene::1", "GiG", "Gene::2"),
    ("Gene::2", "GiG", "Gene::3"),
    ("Gene::3", "GiG", "Gene::1"),
    ("Gene::3", "GiG", "Gene::4"),
    ("Gene::4", "GiG", "Gene::5"),
    ("Gene::6", "GiG", "Gene::7"),
    ("Disease::DOID:0001", "DaG", "Gene::1"),
    ("Disease::DOID:0001", "DaG", "Gene::2"),
    ("Disease::DOID:0001", "DaG", "Gene::6"),
    ("Disease::DOID:0001", "DuG", "Gene::4"),
    ("Disease::DOID:0001", "DdG", "Gene::9"),
    ("Gene::1", "GpPW", "Pathway::PW1"),
    ("Gene::2", "GpPW", "Pathway::PW1"),
    ("Gene::7", "GpPW", "Pathway::PW2"),
    ("Gene::3", "GpBP", "Biological Process::GO:1"),
    ("Compound::DB0001", "CtD", "Disease::DOID:0001"),
    ("Compound::DB0001", "CbG", "Gene::2"),
    ("Compound::DB0002", "CbG", "Gene::8"),
    ("Gene::1", "Gr>G", "Gene::2"),  # ignored metaedge
]


@pytest.fixture
def hetionet_files(tmp_path):
    nodes = tmp_path / "nodes.tsv"
    pd.DataFrame(NODES, columns=["id", "name", "kind"]).to_csv(nodes, sep="\t", index=False)
    edges = tmp_path / "edges.sif.gz"
    pd.DataFrame(EDGES, columns=["source", "metaedge", "target"]).to_csv(edges, sep="\t", index=False)
    return nodes, edges
