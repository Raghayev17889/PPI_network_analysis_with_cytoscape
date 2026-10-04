"""Over-representation analysis of pathways and GO biological processes.

For a query gene list (a cluster, or the whole disease network) each gene set
is tested with the one-sided hypergeometric test

    P(X >= k),  X ~ Hypergeometric(N, K, n)

where ``N`` is the background (proteins of the interactome that have at least
one annotation), ``K`` the gene-set size within the background, ``n`` the
query size within the background and ``k`` the overlap. P-values are
corrected with the Benjamini-Hochberg false discovery rate.
"""

from __future__ import annotations

import logging

import numpy as np
import pandas as pd
from scipy.stats import hypergeom

logger = logging.getLogger(__name__)


def benjamini_hochberg(p_values: np.ndarray | list[float]) -> np.ndarray:
    """Benjamini-Hochberg adjusted p-values (q-values), monotone and capped at 1."""
    p = np.asarray(p_values, dtype=float)
    n = p.size
    if n == 0:
        return p
    order = np.argsort(p)
    ranked = p[order] * n / np.arange(1, n + 1)
    adjusted = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.minimum(adjusted, 1.0)
    return out


def enrich(
    query: list[str] | set[str],
    gene_sets: dict[str, set[str]],
    background: set[str],
    min_size: int = 5,
    max_size: int = 500,
    symbols: dict[str, str] | None = None,
) -> pd.DataFrame:
    """Hypergeometric enrichment of ``query`` in every gene set.

    Gene sets are restricted to ``background`` and kept if their size lies in
    ``[min_size, max_size]``. Returns a table sorted by p-value with columns
    ``term``, ``overlap``, ``set_size``, ``query_size``, ``background_size``,
    ``fold_enrichment``, ``p_value``, ``fdr`` and ``genes``.
    """
    background = set(background)
    query_in = set(query) & background
    n_background, n_query = len(background), len(query_in)
    rows = []
    for term, members in gene_sets.items():
        members_in = members & background
        size = len(members_in)
        if size < min_size or size > max_size:
            continue
        hits = query_in & members_in
        if not hits:
            continue
        rows.append(
            {
                "term": term,
                "overlap": len(hits),
                "set_size": size,
                "query_size": n_query,
                "background_size": n_background,
                "fold_enrichment": (len(hits) / n_query) / (size / n_background) if n_query else np.nan,
                "p_value": float(hypergeom.sf(len(hits) - 1, n_background, size, n_query)),
                "genes": ", ".join(sorted((symbols or {}).get(g, g) for g in hits)),
            }
        )
    columns = ["term", "overlap", "set_size", "query_size", "background_size", "fold_enrichment", "p_value",
               "fdr", "genes"]
    if not rows:
        return pd.DataFrame(columns=columns)
    table = pd.DataFrame(rows)
    # FDR over every tested set, including those without overlap (they have p = 1).
    tested = sum(min_size <= len(m & background) <= max_size for m in gene_sets.values())
    p_all = np.concatenate([table["p_value"].to_numpy(), np.ones(tested - len(table))])
    table["fdr"] = benjamini_hochberg(p_all)[: len(table)]
    return table.sort_values(["p_value", "term"]).reset_index(drop=True)[columns]


def enrich_clusters(
    labels: dict[str, int],
    collections: dict[str, dict[str, set[str]]],
    background: set[str],
    symbols: dict[str, str] | None = None,
    fdr: float = 0.05,
    min_size: int = 5,
    max_size: int = 500,
) -> pd.DataFrame:
    """Significant terms (``fdr`` threshold) for every cluster and every collection."""
    tables = []
    for cluster in sorted(c for c in set(labels.values()) if c > 0):
        genes = [g for g, c in labels.items() if c == cluster]
        for source, gene_sets in collections.items():
            table = enrich(genes, gene_sets, background, min_size, max_size, symbols)
            table = table[table["fdr"] <= fdr]
            if len(table):
                tables.append(table.assign(cluster=cluster, source=source))
    if not tables:
        return pd.DataFrame(columns=["cluster", "source", "term", "overlap", "fdr"])
    result = pd.concat(tables, ignore_index=True)
    first = ["cluster", "source"]
    return result[first + [c for c in result.columns if c not in first]]


def clean_term(term: str) -> str:
    """Drop the ``[identifier]`` suffix: ``'Apoptosis [PC7_1]'`` -> ``'Apoptosis'``."""
    return term.rsplit(" [", 1)[0]


def cluster_labels(enrichment: pd.DataFrame, source: str = "pathway", max_length: int = 40) -> dict[int, str]:
    """A short functional label per cluster: its most significant term in ``source``."""
    labels = {}
    subset = enrichment[enrichment["source"] == source].sort_values(["cluster", "p_value"])
    for cluster, group in subset.groupby("cluster"):
        term = clean_term(group.iloc[0]["term"])
        labels[int(cluster)] = term if len(term) <= max_length else term[: max_length - 1].rstrip() + "…"
    return labels
