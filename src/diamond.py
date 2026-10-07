"""DIAMOnD: disease module detection by connectivity significance.

Implementation of Ghiassian, Menche & Barabási, *PLoS Comput Biol* 11,
e1004120 (2015). Known disease genes are an incomplete sample of the true
disease module. DIAMOnD grows the module one protein at a time, always adding
the protein whose number of links to the current module is least likely to
arise by chance:

    p(k, k_s) = P(X >= k_s),   X ~ Hypergeometric(N, s, k)

where ``N`` is the number of proteins in the interactome, ``s`` the current
module size, ``k`` the candidate's degree and ``k_s`` its links into the
module. With ``alpha > 1`` links to the original seeds count ``alpha`` times.
"""

from __future__ import annotations

import logging

import networkx as nx
import numpy as np
import pandas as pd
from scipy.special import gammaln

logger = logging.getLogger(__name__)


def hypergeom_log_sf(x_min: np.ndarray, total: int, successes: int, draws: np.ndarray) -> np.ndarray:
    """``log P(X >= x_min)`` for ``X ~ Hypergeometric(total, successes, draws)``, vectorised.

    Sums the probability mass function over ``x_min..min(successes, draws)``
    in log space using log-gamma functions. This is exact, avoids
    underflow for the tiny p-values DIAMOnD produces, and is far faster than
    calling ``scipy.stats.hypergeom.sf`` thousands of times per iteration.
    """
    x_min = np.asarray(x_min, dtype=np.int64)
    draws = np.asarray(draws, dtype=np.int64)
    upper = np.minimum(successes, draws)
    lengths = np.maximum(upper - x_min + 1, 0)
    result = np.full(x_min.shape, -np.inf)
    valid = lengths > 0
    if not valid.any():
        return result

    lengths_v = lengths[valid]
    seg = np.repeat(np.arange(lengths_v.size), lengths_v)
    starts = np.cumsum(lengths_v) - lengths_v
    x = x_min[valid][seg] + (np.arange(lengths_v.sum()) - starts[seg])
    n = draws[valid][seg]

    def log_choose(a, b):
        return gammaln(a + 1) - gammaln(b + 1) - gammaln(a - b + 1)

    log_pmf = log_choose(successes, x) + log_choose(total - successes, n - x) - log_choose(total, n)
    peak = np.maximum.reduceat(log_pmf, starts)
    sums = np.add.reduceat(np.exp(log_pmf - peak[seg]), starts)
    result[valid] = peak + np.log(sums)
    return np.minimum(result, 0.0)


def diamond(graph: nx.Graph, seeds: list[str], n_genes: int = 200, alpha: int = 1) -> pd.DataFrame:
    """Rank ``n_genes`` proteins to add to the seed module.

    Parameters
    ----------
    graph:
        Interactome.
    seeds:
        Known disease genes (those absent from ``graph`` are ignored).
    n_genes:
        Number of DIAMOnD iterations, i.e. proteins to add.
    alpha:
        Weight of links to seed genes (``1`` = original unweighted DIAMOnD).

    Returns
    -------
    DataFrame with one row per added protein, in order: ``rank``, ``gene``,
    ``degree``, ``links_to_module``, ``p_value``.
    """
    if alpha < 1:
        raise ValueError("alpha must be >= 1")
    nodes = list(graph.nodes())
    index = {n: i for i, n in enumerate(nodes)}
    seed_idx = sorted({index[s] for s in seeds if s in index})
    if not seed_idx:
        raise ValueError("None of the seed genes is in the interactome")

    adjacency = nx.to_scipy_sparse_array(graph, nodelist=nodes, weight=None, format="csr")
    n_total = len(nodes)
    degree = np.asarray(adjacency.sum(axis=1)).ravel().astype(np.int64)

    in_module = np.zeros(n_total, dtype=bool)
    in_module[seed_idx] = True
    seed_mask = in_module.copy()
    # links_to_module[i] = number of module members adjacent to i
    links = np.asarray(adjacency[:, seed_idx].sum(axis=1)).ravel().astype(np.int64)
    seed_links = links.copy()
    n_seeds = len(seed_idx)
    module_size = n_seeds

    rows = []
    for rank in range(1, n_genes + 1):
        candidates = np.flatnonzero((links > 0) & ~in_module)
        if candidates.size == 0:
            logger.warning("DIAMOnD stopped after %d genes: no connected candidates left", rank - 1)
            break
        extra = (alpha - 1) * seed_links[candidates]
        k_s = links[candidates] + extra
        k = degree[candidates] + extra
        s = module_size + (alpha - 1) * n_seeds
        log_p = hypergeom_log_sf(k_s, n_total + (alpha - 1) * n_seeds, s, k)
        # Lowest p-value; ties broken by more links, then by node order.
        order = np.lexsort((candidates, -k_s, log_p))
        best = int(candidates[order[0]])
        rows.append((rank, nodes[best], int(degree[best]), int(links[best]), float(np.exp(log_p[order[0]]))))

        in_module[best] = True
        module_size += 1
        neighbours = adjacency.indices[adjacency.indptr[best] : adjacency.indptr[best + 1]]
        links[neighbours] += 1

    result = pd.DataFrame(rows, columns=["rank", "gene", "degree", "links_to_module", "p_value"])
    logger.info("DIAMOnD added %d genes to %d seeds", len(result), n_seeds)
    return result


def diamond_cross_validation(
    graph: nx.Graph,
    seeds: list[str],
    n_genes: int = 200,
    folds: int = 5,
    seed: int = 42,
) -> pd.DataFrame:
    """How well does DIAMOnD recover held-out disease genes?

    Seeds are split into ``folds`` parts. For each fold DIAMOnD runs on the
    remaining seeds and we count how many held-out seeds appear among the
    first ``n_genes`` additions, compared with the number expected if the
    additions were random proteins (outside the training seeds).

    Returns one row per fold with ``recovered``, ``expected``,
    ``fold_enrichment`` and the DIAMOnD ranks at which held-out genes were
    found (``recovered_ranks``, for recovery curves).
    """
    seeds = sorted({s for s in seeds if s in graph})
    rng = np.random.default_rng(seed)
    order = rng.permutation(len(seeds))
    rows = []
    for fold in range(folds):
        held_out = {seeds[i] for i in order[fold::folds]}
        training = [s for s in seeds if s not in held_out]
        added = diamond(graph, training, n_genes=n_genes)
        hits = added[added["gene"].isin(held_out)]
        recovered = len(hits)
        background = graph.number_of_nodes() - len(training)
        expected = len(added) * len(held_out) / background
        rows.append(
            {
                "fold": fold + 1,
                "training_seeds": len(training),
                "held_out_seeds": len(held_out),
                "background": background,
                "recovered": recovered,
                "expected": expected,
                "fold_enrichment": recovered / expected if expected else float("nan"),
                "recovered_ranks": hits["rank"].tolist(),
            }
        )
    return pd.DataFrame(rows)
