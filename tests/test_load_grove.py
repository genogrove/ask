# SPDX-License-Identifier: GPL-3.0-or-later
"""load_grove builds a .gg once, then deserializes it. Runs only where pygenogrove
is installed (CI skips); uses a file:// resource so it needs no network."""

from __future__ import annotations

import hashlib

import pytest

pg = pytest.importorskip("pygenogrove")

from ask import resources
from ask.resources import Resource

GFF3 = (
    "##gff-version 3\n"
    "chr1\tH\tgene\t1000\t2000\t.\t+\t.\tID=g1;gene_name=AAA\n"
    "chr1\tH\ttranscript\t1000\t2000\t.\t+\t.\tID=t1;Parent=g1\n"
    "chr1\tH\texon\t1000\t1500\t.\t+\t.\tID=e1;Parent=t1\n"
).encode()


def _register(monkeypatch, tmp_path) -> str:
    src = tmp_path / "mini.gff3"
    src.write_bytes(GFF3)
    sha = hashlib.sha256(GFF3).hexdigest()
    monkeypatch.setattr(resources, "_CACHE", tmp_path / "cache")
    monkeypatch.setitem(resources.RESOURCES, "_mini", Resource("_mini", src.as_uri(), sha))
    return sha


def _gg_path(tmp_path, sha):
    return tmp_path / "cache" / "groves" / f"{sha}.{resources._GROVE_SCHEMA}.gg"


def test_load_grove_builds_caches_and_deserializes(monkeypatch, tmp_path) -> None:
    sha = _register(monkeypatch, tmp_path)

    g = resources.load_grove("_mini")        # first call builds + serializes
    assert g.size() == 3
    assert _gg_path(tmp_path, sha).exists()

    g2 = resources.load_grove("_mini")       # second call deserializes the .gg
    assert g2.size() == 3
    # edges survive the roundtrip: gene -> transcript (contains)
    gene = next(k for k in g2.intersect(pg.GenomicCoordinate("*", 1500, 1500), "chr1")
                if k.data["type"] == "gene")
    assert [n.data["type"] for n in g2.get_neighbors(gene)] == ["transcript"]


def test_load_grove_rebuilds_unreadable_cache(monkeypatch, tmp_path) -> None:
    sha = _register(monkeypatch, tmp_path)
    resources.load_grove("_mini")            # build the cache
    _gg_path(tmp_path, sha).write_bytes(b"not a grove")  # corrupt it

    g = resources.load_grove("_mini")        # self-heals: rebuild rather than crash
    assert g.size() == 3
