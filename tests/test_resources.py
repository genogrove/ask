# SPDX-License-Identifier: GPL-3.0-or-later
"""Resolve = download + checksum-verify + cache. The mismatch branch is the
trust boundary, so it gets the most attention. Uses ``file://`` URLs so the test
needs no network."""

from __future__ import annotations

import hashlib

import pytest

from ask import resources
from ask.resources import Resource


def _catalog_entry(monkeypatch, tmp_path, name, payload, sha256):
    """Point _CACHE at a tmp dir and register a file:// resource for ``name``."""
    src = tmp_path / f"{name}.gff3.gz"
    src.write_bytes(payload)
    monkeypatch.setattr(resources, "_CACHE", tmp_path / "cache")
    monkeypatch.setitem(
        resources.RESOURCES, name, Resource(name, src.as_uri(), sha256)
    )


def test_resolve_verifies_and_caches(monkeypatch, tmp_path) -> None:
    payload = b"##gff-version 3\nchr1\t.\tgene\t1\t100\t.\t+\t.\tID=g1\n"
    sha = hashlib.sha256(payload).hexdigest()
    _catalog_entry(monkeypatch, tmp_path, "good", payload, sha)

    path = resources.resolve("good")
    assert path.read_bytes() == payload
    assert sha in str(path)  # content-addressed
    # second call is a cache hit: same path, no re-download
    assert resources.resolve("good") == path


def test_resolve_rejects_checksum_mismatch(monkeypatch, tmp_path) -> None:
    _catalog_entry(monkeypatch, tmp_path, "bad", b"tampered", "00" * 32)

    with pytest.raises(ValueError, match="checksum mismatch"):
        resources.resolve("bad")
    # nothing committed to the cache on mismatch
    assert not list((tmp_path / "cache").rglob("*.gff3.gz"))


def test_resolve_unknown_name_raises() -> None:
    with pytest.raises(KeyError):
        resources.resolve("not.curated")
