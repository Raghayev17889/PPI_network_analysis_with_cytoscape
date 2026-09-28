from __future__ import annotations

import hashlib

import pytest

from src.data_download import RemoteFile, download_file, sha256sum
from src.hetionet import (
    extract_disease_data,
    load_disease_data,
    load_edges,
    load_nodes,
    read_gmt,
    save_disease_data,
    strip_kind,
    write_gmt,
)


@pytest.fixture
def toy(hetionet_files):
    nodes, edges = hetionet_files
    return extract_disease_data(load_nodes(nodes), load_edges(edges), "DOID:0001")


def test_strip_kind():
    assert strip_kind("Gene::5663") == "5663"
    assert strip_kind("Disease::DOID:10652") == "DOID:10652"
    assert strip_kind("plain") == "plain"


def test_load_edges_filters_metaedges(hetionet_files):
    _, edges = hetionet_files
    assert "Gr>G" not in set(load_edges(edges)["metaedge"])
    assert "Gr>G" in set(load_edges(edges, metaedges=None)["metaedge"])


def test_extract_disease_data(toy):
    assert toy.disease_name == "toy disease"
    assert toy.seed_genes == ["1", "2", "6"]
    assert toy.genes_with_evidence("upregulated") == ["4"]
    assert toy.genes_with_evidence("downregulated") == ["9"]
    assert len(toy.interactions) == 6
    assert toy.symbols["3"] == "G3"
    assert toy.pathways == {"Toy signalling [PW1]": {"1", "2"}, "Toy metabolism [PW2]": {"7"}}
    assert toy.go_processes == {"toy process [GO:1]": {"3"}}
    assert toy.drugs["name"].tolist() == ["Drugamab"]
    assert toy.drug_targets[["entrez", "drug"]].values.tolist() == [["2", "Drugamab"]]
    assert toy.druggable == {"2", "8"}
    assert toy.summary()["associated_genes"] == 3


def test_unknown_disease(hetionet_files):
    nodes, edges = hetionet_files
    with pytest.raises(KeyError):
        extract_disease_data(load_nodes(nodes), load_edges(edges), "DOID:9999")


def test_roundtrip(tmp_path, toy):
    directory = save_disease_data(toy, tmp_path / "processed")
    loaded = load_disease_data(directory)
    assert loaded.seed_genes == toy.seed_genes
    assert loaded.pathways == toy.pathways
    assert loaded.druggable == toy.druggable
    assert loaded.summary() == toy.summary()


def test_gmt_roundtrip(tmp_path):
    sets = {"A pathway [PW9]": {"10", "2", "33"}, "B [X]": {"1"}}
    path = write_gmt(sets, tmp_path / "sets.gmt")
    assert path.read_text().splitlines()[0] == "A pathway [PW9]\tPW9\t2\t10\t33"
    assert read_gmt(path) == sets


def test_download_verifies_checksum(tmp_path):
    source = tmp_path / "src.txt"
    source.write_bytes(b"hello")
    good = RemoteFile("f.txt", source.as_uri(), hashlib.sha256(b"hello").hexdigest(), "test")
    path = download_file(good, tmp_path / "raw")
    assert sha256sum(path) == good.sha256
    bad = RemoteFile("g.txt", source.as_uri(), "0" * 64, "test")
    with pytest.raises(ValueError):
        download_file(bad, tmp_path / "raw")
    assert sorted(p.name for p in (tmp_path / "raw").iterdir()) == ["f.txt"]
