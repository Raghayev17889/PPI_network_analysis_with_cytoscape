from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
from scipy.stats import hypergeom

from src.enrichment import benjamini_hochberg, clean_term, cluster_labels, enrich, enrich_clusters

BACKGROUND = {str(i) for i in range(1000)}
SETS = {
    "Target pathway [PW1]": {str(i) for i in range(10)},
    "Unrelated [PW2]": {str(i) for i in range(500, 520)},
    "Tiny [PW3]": {"1", "2"},
    "Partly outside [PW4]": {str(i) for i in range(995, 1010)},
}


def test_benjamini_hochberg():
    p = [0.01, 0.04, 0.03, 0.5]
    np.testing.assert_allclose(benjamini_hochberg(p), [0.04, 0.04 * 4 / 3, 0.04 * 4 / 3, 0.5])
    assert benjamini_hochberg([]).size == 0
    assert (benjamini_hochberg([0.9, 0.95]) <= 1).all()


def test_enrich_matches_scipy():
    query = [str(i) for i in range(8)] + ["600", "601"]
    table = enrich(query, SETS, BACKGROUND, min_size=5, symbols={"0": "GENE0"})
    assert table["term"].tolist() == ["Target pathway [PW1]"]  # PW2 has no overlap, PW3 too small
    row = table.iloc[0]
    assert row["overlap"] == 8 and row["set_size"] == 10 and row["query_size"] == 10
    assert row["p_value"] == pytest.approx(hypergeom.sf(7, 1000, 10, 10))
    assert row["fold_enrichment"] == pytest.approx((8 / 10) / (10 / 1000))
    assert row["genes"].startswith("1, 2") and "GENE0" in row["genes"]
    # FDR accounts for every tested set (PW1, PW2 and PW4 have >= 5 background members).
    assert row["fdr"] == pytest.approx(min(row["p_value"] * 3, 1.0))


def test_enrich_empty():
    assert enrich(["999999"], SETS, BACKGROUND).empty


def test_enrich_clusters_and_labels():
    labels = {**{str(i): 1 for i in range(8)}, **{str(i): 2 for i in range(700, 710)}, "5": 1, "x": 0}
    result = enrich_clusters(labels, {"pathway": SETS}, BACKGROUND, fdr=0.05)
    assert set(result["cluster"]) == {1}
    assert result.columns[:2].tolist() == ["cluster", "source"]
    assert cluster_labels(result) == {1: "Target pathway"}
    assert clean_term("Translation [PC7_8331]") == "Translation"
    long = pd.DataFrame({"cluster": [3], "source": ["pathway"], "term": ["A" * 60 + " [X]"], "p_value": [0.01]})
    assert len(cluster_labels(long, max_length=20)[3]) == 20
