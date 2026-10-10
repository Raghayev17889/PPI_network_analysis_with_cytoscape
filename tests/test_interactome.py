from __future__ import annotations

import pandas as pd

from src.interactome import build_interactome, induced_subgraph, interactome_summary, map_genes, symbols_to_entrez


def _interactions(pairs):
    return pd.DataFrame(pairs, columns=["gene_a", "gene_b"])


def test_build_interactome_cleans_and_keeps_lcc():
    graph = build_interactome(
        _interactions([("1", "2"), ("2", "1"), ("2", "3"), ("3", "3"), ("8", "9")]),
        symbols={"1": "A", "2": "B"},
    )
    assert set(graph) == {"1", "2", "3"}
    assert graph.number_of_edges() == 2
    assert graph.nodes["1"]["symbol"] == "A"
    assert graph.nodes["3"]["symbol"] == "3"  # falls back to the ID
    full = build_interactome(_interactions([("1", "2"), ("8", "9")]), largest_component=False)
    assert full.number_of_nodes() == 4


def test_map_genes_and_summary():
    graph = build_interactome(_interactions([("1", "2"), ("2", "3")]))
    mapping = map_genes(graph, ["3", "1", "99", 1])
    assert mapping.mapped == ["1", "3"] and mapping.missing == ["99"]
    assert mapping.coverage == 2 / 3
    summary = interactome_summary(graph)
    assert summary["proteins"] == 3 and summary["interactions"] == 2 and summary["max_degree"] == 2


def test_exclude_hub_before_component():
    # Removing the hub "0" splits the star; the largest remaining piece is 1-2.
    graph = build_interactome(_interactions([("0", "1"), ("0", "3"), ("0", "4"), ("1", "2")]), exclude=["0"])
    assert set(graph) == {"1", "2"}
    assert symbols_to_entrez(["UBC", "NOPE"], {"7316": "UBC"}) == ["7316"]


def test_induced_subgraph_keeps_parent_order_and_attributes():
    import networkx as nx

    graph = nx.Graph()
    graph.add_node("z", symbol="Z")
    graph.add_edges_from([("z", "a", {"w": 1}), ("a", "m"), ("m", "q")])
    sub = induced_subgraph(graph, {"m", "a", "z"})
    assert list(sub) == ["z", "a", "m"]
    assert sub.nodes["z"]["symbol"] == "Z" and sub.edges["z", "a"]["w"] == 1
    assert not sub.has_node("q")


def test_build_interactome_is_identical_across_processes():
    import subprocess
    import sys

    script = (
        "import networkx as nx, pandas as pd; from src.interactome import build_interactome;"
        "g = nx.gnm_random_graph(300, 900, seed=1);"
        "df = pd.DataFrame([(str(u), str(v)) for u, v in g.edges()], columns=['gene_a', 'gene_b']);"
        "h = build_interactome(df); print(list(h)[:50], list(h.edges())[:50])"
    )
    outputs = {
        subprocess.run([sys.executable, "-c", script], capture_output=True, text=True, check=True,
                       env={"PYTHONHASHSEED": str(seed), "PYTHONPATH": "."}).stdout
        for seed in (1, 2, 3)
    }
    assert len(outputs) == 1
