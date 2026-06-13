# Changelog

All notable changes to this project are documented here. The format is based on
[Keep a Changelog](https://keepachangelog.com/en/1.0.0/), and this project adheres
to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added
- Initial project skeleton: `uv` + hatchling packaging, GPL-3.0-or-later license,
  module layout (`cli`, `llm`, `sandbox`, `registry`, `prompts/system.md`), CLI
  smoke tests, and a CI workflow stub.
- **`pygenogrove` API surface in the codegen system prompt** (`prompts/system.md`):
  the bound `Grove` / `BedGrove` / `GffGrove` surface (insert / intersect / flanking),
  the directed-edge graph overlay, file readers, and serialization — with the
  coordinate-convention rules (closed `Interval` vs half-open BED vs 1-based GFF) and
  a worked 2-hop connected-interval example. Documented against `pygenogrove` 0.2.0.
- **Build pinning for Level 2 reproducibility** (`registry.py`): `BuildPin`,
  `verify_pygenogrove_build()` (raises on version drift, returns the C++ engine
  version), and `build_manifest()` for run provenance. `pygenogrove` is pinned to the
  immutable commit `1a9c975` (tag `v0.2.0`) in `pyproject.toml` and mirrored here.
- **Pin-drift guard test** (`tests/test_registry_pins.py`): asserts the `pyproject.toml`
  `==<version>` pin and the `[tool.uv.sources]` `rev` both match `registry.PYGENOGROVE`,
  and that the pin is a full immutable commit SHA — so the Level 2 "all three agree"
  check fails CI on drift instead of relying on manual review. Parses `pyproject.toml`
  with regexes (no `tomllib`) so it runs on the py3.9 floor
  ([#2](https://github.com/genogrove/ask/pull/2)).
- **Out-of-process sandbox** (`sandbox.py`): runs untrusted model-generated Python in an
  isolated subprocess with parent/OS-enforced hard guarantees (stripped env, `setrlimit`
  CPU/memory/no-write/fd caps, whole-session wall-clock kill, byte-capped output) plus
  in-child defense-in-depth (import allowlist with the network/exec primitives scrubbed,
  read-only `open` restricted to registry data roots). The hard boundary is the parent/OS
  layer; an OS-level backend (seccomp/namespaces) is the documented next step for
  adversarial robustness. Covered by 22 isolation tests
  ([#3](https://github.com/genogrove/ask/pull/3)).
