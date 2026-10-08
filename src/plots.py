"""Static figures (matplotlib) that complement the Cytoscape network views."""

from __future__ import annotations

from contextlib import contextmanager
from pathlib import Path
from typing import Iterator

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from src.cytoscape_viz import CLUSTER_COLORS, INK, cluster_color

PALETTE = CLUSTER_COLORS
STYLE = {
    "figure.facecolor": INK["surface"],
    "axes.facecolor": INK["surface"],
    "savefig.facecolor": INK["surface"],
    "axes.edgecolor": "#c3c2b7",
    "axes.labelcolor": INK["secondary"],
    "axes.titlecolor": INK["primary"],
    "axes.titlesize": 12,
    "axes.titleweight": "bold",
    "axes.titlelocation": "left",
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.grid": True,
    "grid.color": "#e1e0d9",
    "grid.linewidth": 0.6,
    "xtick.color": INK["muted"],
    "ytick.color": INK["muted"],
    "xtick.labelcolor": INK["secondary"],
    "ytick.labelcolor": INK["secondary"],
    "legend.frameon": False,
    "legend.fontsize": 9,
    "lines.linewidth": 1.8,
}
CLASS_COLORS = {
    "hub-bottleneck": PALETTE[7],
    "non-hub bottleneck": PALETTE[3],
    "hub non-bottleneck": PALETTE[0],
    "peripheral": "#c3c2b7",
}


@contextmanager
def plot_style() -> Iterator[None]:
    with plt.rc_context(STYLE):
        yield


def save(fig: plt.Figure, path: str | Path, dpi: int = 160) -> Path:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    return path


def plot_module_significance(random_lcc: np.ndarray, observed: int, z_score: float, n_genes: int) -> plt.Figure:
    """Histogram of random LCC sizes with the observed seed LCC."""
    with plot_style():
        fig, ax = plt.subplots(figsize=(6.4, 4.0))
        bins = np.arange(random_lcc.min() - 1, max(random_lcc.max(), observed) + 3, 3)
        ax.hist(random_lcc, bins=bins, color=PALETTE[0], alpha=0.85, label="degree-matched random gene sets")
        ax.axvline(observed, color=PALETTE[7], linewidth=2.2, label=f"Alzheimer's seed genes ({observed})")
        ax.set(xlabel=f"Largest connected component (of {n_genes} genes)", ylabel="Random sets")
        ax.set_title(f"The disease genes form a connected module (z = {z_score:.1f})")
        ax.legend(loc="upper left")
        ax.grid(axis="x", visible=False)
    return fig


def plot_diamond_recovery(cv: pd.DataFrame, n_genes: int) -> plt.Figure:
    """Cumulative recovery of held-out seed genes along the DIAMOnD ranking, per fold."""
    with plot_style():
        fig, ax = plt.subplots(figsize=(6.4, 4.2))
        steps = np.arange(1, n_genes + 1)
        for i, row in enumerate(cv.itertuples()):
            ranks = np.asarray(row.recovered_ranks, dtype=int)
            recovered = np.array([(ranks <= s).sum() for s in steps]) / row.held_out_seeds
            ax.plot(steps, recovered, color=PALETTE[0], alpha=0.35 + 0.1 * i, linewidth=1.4,
                    label="DIAMOnD (5 folds)" if i == 0 else None)
        expected = steps / cv["background"].mean()  # random picks recover held-out genes at this rate
        ax.plot(steps, expected, color=INK["muted"], linestyle=(0, (5, 2)), label="random expectation")
        ax.set(xlabel="DIAMOnD iteration", ylabel="Held-out disease genes recovered")
        ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda v, _: f"{v:.0%}"))
        total = cv["recovered"].sum() / cv["expected"].sum()
        ax.set_title(f"DIAMOnD recovers held-out Alzheimer's genes ({total:.1f}× random)")
        ax.legend(loc="upper left")
    return fig


def plot_centrality_scatter(table: pd.DataFrame, n_labels: int = 12) -> plt.Figure:
    """Degree vs betweenness coloured by hub/bottleneck class, key proteins labelled."""
    with plot_style():
        fig, ax = plt.subplots(figsize=(7.0, 5.0))
        for cls, colour in CLASS_COLORS.items():
            subset = table[table["centrality_class"] == cls]
            ax.scatter(subset["degree"], subset["betweenness"] + 1e-5, s=26 if cls != "peripheral" else 14,
                       color=colour, label=f"{cls} ({len(subset)})", edgecolors=INK["surface"], linewidths=0.5,
                       zorder=3 if cls != "peripheral" else 2)
        top = table.sort_values("betweenness", ascending=False).head(n_labels)
        for row in top.itertuples():
            ax.annotate(row.symbol, (row.degree, row.betweenness + 1e-5), xytext=(4, 3), textcoords="offset points",
                        fontsize=8.5, color=INK["primary"])
        ax.axvline(table.attrs.get("hub_degree_threshold", np.nan), color=INK["muted"], linewidth=1, linestyle=":")
        ax.axhline(table.attrs.get("bottleneck_betweenness_threshold", np.nan), color=INK["muted"], linewidth=1,
                   linestyle=":")
        ax.set(xscale="log", yscale="log", xlabel="Degree (interactions in the AD network)",
               ylabel="Betweenness centrality")
        ax.set_title("Hubs and bottlenecks of the Alzheimer's network")
        ax.legend(loc="lower right")
    return fig


def plot_robustness(curves: pd.DataFrame) -> plt.Figure:
    """Largest-component size as proteins are removed by different strategies."""
    styles = {"betweenness": (PALETTE[7], "-"), "degree": (PALETTE[0], "-"), "random": (INK["muted"], (0, (5, 2)))}
    names = {"betweenness": "highest betweenness first", "degree": "highest degree first", "random": "random order"}
    with plot_style():
        fig, ax = plt.subplots(figsize=(6.4, 4.2))
        for strategy, group in curves.groupby("strategy", sort=False):
            colour, line = styles[strategy]
            ax.plot(group["fraction_removed"], group["lcc_fraction"], color=colour, linestyle=line,
                    label=names[strategy])
        ax.set(xlabel="Proteins removed", ylabel="Largest connected component", ylim=(0, 1.02))
        ax.xaxis.set_major_formatter(plt.FuncFormatter(lambda v, _: f"{v:.0%}"))
        ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda v, _: f"{v:.0%}"))
        ax.set_title("Removing bottlenecks fragments the network")
        ax.legend(loc="lower left")
    return fig


def plot_cluster_enrichment(enrichment: pd.DataFrame, cluster_names: dict[int, str], per_cluster: int = 3,
                            source: str = "pathway") -> plt.Figure:
    """Top terms per cluster as horizontal bars of -log10(FDR), coloured by cluster."""
    from src.enrichment import clean_term

    subset = enrichment[enrichment["source"] == source].sort_values(["cluster", "p_value"])
    top = subset.groupby("cluster").head(per_cluster)
    with plot_style():
        fig, ax = plt.subplots(figsize=(8.2, 0.34 * len(top) + 1.4))
        labels, values, colours, y = [], [], [], []
        for i, row in enumerate(top.itertuples()):
            term = clean_term(row.term)
            labels.append((term[:55] + "…") if len(term) > 56 else term)
            values.append(-np.log10(row.fdr))
            colours.append(cluster_color(int(row.cluster)))
            y.append(i)
        ax.barh(y, values, color=colours, height=0.72)
        ax.set_yticks(y, labels, fontsize=8.5)
        ax.invert_yaxis()
        ax.grid(axis="y", visible=False)
        ax.set(xlabel="-log10 FDR")
        handles = [plt.Rectangle((0, 0), 1, 1, color=cluster_color(c)) for c in sorted(top["cluster"].unique())]
        ax.legend(handles, [f"{c}. {cluster_names.get(int(c), '')}" for c in sorted(top["cluster"].unique())],
                  loc="lower right", fontsize=8)
        ax.set_title("Top enriched pathways per functional cluster")
    return fig
