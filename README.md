# genogrove ask

**Natural-language interface for [genogrove](https://github.com/genogrove/genogrove) — ask plain-English questions over connected genomic intervals, powered by [pygenogrove](https://github.com/genogrove/pygenogrove).**

> ⚠️ **Status: pre-alpha skeleton.** The project structure, packaging, and module
> boundaries are in place; the end-to-end question→answer loop is not implemented
> yet. See [Roadmap](#roadmap).

---

## What this is

genogrove stores genomic annotations as a *connected* structure: intervals indexed
in per-chromosome B+ trees, with a directed graph overlay linking related keys
(exon→transcript, breakpoint→mate, enhancer→gene). Relational questions that would
otherwise require a brittle `intersect | awk | sort | join` pipeline become a single
traversal of that structure.

`genogrove ask` puts a natural-language front end on top of it. You ask a question in
plain English; the tool translates it into Python that drives `pygenogrove`, runs that
Python in a sandbox, and prints the result.

```console
$ genogrove ask "Which transcripts share an exon with the variant at chr7:55,191,822?"
```

## How it works (architecture)

The design is **path B**: the LLM generates Python that calls into the genogrove C++
engine through the `pygenogrove` bindings. There is no fixed query algebra and no
plan interpreter — the bindings *are* the interface the model targets.

```
  question (plain English)
        │
        ▼
  ┌───────────────┐   generated Python (uses pygenogrove)
  │  llm.py       │ ─────────────────────────────────────┐
  │  Anthropic    │                                       │
  │  claude-opus  │                                       ▼
  └───────────────┘                              ┌──────────────────┐
        ▲                                         │  sandbox.py      │
        │  curated tool/schema + resource context │  restricted exec │
  ┌───────────────┐                               └──────────────────┘
  │  registry.py  │  pinned datasets + builds              │
  │  (Level 2     │ ◀─────────────────────────────────────┘
  │  reproducib.) │                                    results → stdout
  └───────────────┘
```

| Module | Responsibility |
|---|---|
| `cli.py` | Thin CLI wrapper: parse the question, orchestrate llm → sandbox, print results. |
| `llm.py` | Anthropic client + the code-generation prompt. Default model: `claude-opus-4-8`. |
| `sandbox.py` | Run the generated Python under restrictions (no network, allowlisted imports, resource caps). **Security-critical.** |
| `registry.py` | Curated registry of datasets (pinned URLs + checksums) and pinned `pygenogrove`/`genogrove` builds. |
| `prompts/system.md` | The system prompt that teaches the model the `pygenogrove` surface and the rules for generated code. |

**Dependency direction is one-way:** `ask → pygenogrove`. `pygenogrove` stays a lean,
stable bindings layer with no LLM dependency, so `pip install pygenogrove` never drags
in an LLM SDK.

## Reproducibility

The project commits to **Level 2** reproducibility: a *curated resource registry* with
pinned dataset versions (URL + checksum) and pinned library builds. Given the same
question and the same registry snapshot, a run is reproducible. Open-web resource
discovery (Level 3) is explicitly out of scope.

## Installation

This project uses [`uv`](https://docs.astral.sh/uv/).

```console
$ git clone https://github.com/genogrove/ask
$ cd ask
$ uv sync                  # creates the venv and resolves deps (incl. pygenogrove from git)
$ export ANTHROPIC_API_KEY=sk-ant-...
$ uv run genogrove-ask --help
```

`pygenogrove` is resolved from its GitHub repository (see `[tool.uv.sources]` in
`pyproject.toml`); it is not yet on PyPI.

## The `genogrove ask` surface

In the genogrove paper and docs the command is written `genogrove ask <question>`.
That is a thin alias over the `genogrove-ask` console script this package installs —
the application layer lives here, separate from the core C++ CLI, because it has a
different release cadence and audience.

## Roadmap

- [ ] `llm.py` — Anthropic codegen loop (adaptive thinking, structured tool surface)
- [x] `sandbox.py` — restricted execution of generated Python (out-of-process
      isolation; OS-level hardening backend tracked as a follow-up)
- [ ] `registry.py` — curated dataset + pinned-build registry
- [x] `prompts/system.md` — flesh out the `pygenogrove` API surface and codegen rules
      (done for the pinned v0.4.0 build; the dataset-gated "Available resources" block
      lands with the registry item above)
- [ ] End-to-end hero query (≥ 2-hop connected-interval question) for the paper demo

## License

GPL-3.0-or-later. See [LICENSE](LICENSE).
