# genogrove ask

**Natural-language interface for [genogrove](https://github.com/genogrove/genogrove) — ask plain-English questions over connected genomic intervals, powered by [pygenogrove](https://github.com/genogrove/pygenogrove).**

> ⚠️ **Status: pre-alpha.** Packaging, the sandbox, the dataset registry, and the
> GFF→Grove data layer are implemented and tested (see [Try it](#try-it-the-gff--grove-model));
> the end-to-end question→answer loop (`llm.py` + CLI wiring) is not yet. See [Roadmap](#roadmap).

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
  │ resources.py  │  pinned datasets + builds              │
  │  (Level 2     │ ◀─────────────────────────────────────┘
  │  reproducib.) │                                    results → stdout
  └───────────────┘
```

| Module | Responsibility |
|---|---|
| `cli.py` | Thin CLI wrapper: parse the question, orchestrate llm → sandbox, print results. |
| `llm.py` | Anthropic client + the code-generation prompt. Default model: `claude-opus-4-8`. |
| `sandbox.py` | Run the generated Python under restrictions (no network, allowlisted imports, resource caps). **Security-critical.** |
| `resources.py` | Curated registry of datasets (pinned URLs + checksums, resolved + sha256-verified into a local cache) and pinned `pygenogrove`/`genogrove` builds. |
| `gff.py` | Loads GFF/GENCODE into a universal `Grove` — gene/transcript/exon keys, `contains`/`first_exon`/`next` edges, CDS folded onto exons. |
| `prompts/system.md` | The system prompt that teaches the model the `pygenogrove` surface, the GENCODE grove model, and the rules for generated code. |

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

## Try it: the GFF → Grove model

The annotation-loading layer (`ask.gff.load_gff`) is implemented and tested. Running it
needs a local `pygenogrove` build (it's a C++/htslib extension, not yet on PyPI).

**macOS** (Homebrew):

```console
$ brew install htslib cmake
$ git clone https://github.com/genogrove/pygenogrove      # sibling checkout
$ python3 -m venv .venv && source .venv/bin/activate
$ CMAKE_PREFIX_PATH=/opt/homebrew CMAKE_ARGS="-DCMAKE_PREFIX_PATH=/opt/homebrew/opt/htslib" \
    pip install ./pygenogrove pytest
```

(On Linux, `pygenogrove` builds htslib from source — see its `.github/scripts`.)

Run the loader tests against the real bindings (they `importorskip` when pygenogrove is absent):

```console
$ PYTHONPATH=src pytest tests/test_gff.py -q
```

Load a GENCODE locus and query it (`PYTHONPATH=src python try.py`):

```python
import pygenogrove as pg
from ask import gff, resources

path = resources.resolve("gencode.human")   # downloads + sha256-verifies v50 (~70 MB) once, then caches
g = gff.load_gff(path, region=("chr7", 55_000_000, 55_300_000))   # streams the gzip (~20s), builds the locus

# Which gene overlaps chr7:55,191,822 ?  (1-based -> 0-based closed = 55,191,821)
q = pg.GenomicCoordinate("*", 55_191_821, 55_191_821)
for k in g.intersect(q, "chr7"):
    if k.data["type"] == "gene":
        print(k.data["name"], k.data["id"], k.data["biotype"])
# EGFR      ENSG00000146648.23  protein_coding
# EGFR-AS1  ENSG00000224057.3   lncRNA
```

Then traverse: `get_neighbors(gene)` gives transcripts (`contains`), and `first_exon` →
`next` walks a transcript's splice chain, each exon carrying its `cds` range. The full
schema is in [`prompts/system.md`](src/ask/prompts/system.md) under "The GENCODE Grove model".

> The region load streams the whole annotation gzip (no tabix index), ~20s. A cached
> serialized `.gg` per resource (build once, `deserialize` in ms) is the planned fix.

## The `genogrove ask` surface

In the genogrove paper and docs the command is written `genogrove ask <question>`.
That is a thin alias over the `genogrove-ask` console script this package installs —
the application layer lives here, separate from the core C++ CLI, because it has a
different release cadence and audience.

## Roadmap

- [ ] `llm.py` — Anthropic codegen loop (adaptive thinking, structured tool surface)
- [x] `sandbox.py` — restricted execution of generated Python (out-of-process
      isolation; OS-level hardening backend tracked as a follow-up)
- [x] `resources.py` — curated dataset catalog + pinned-build registry (`resolve` +
      sha256 cache; GENCODE v50 pinned; `load_grove` builds + caches a serialized `.gg`
      per resource, deserialized fast on reuse). _Remaining: runtime "Available
      resources" prompt injection._
- [x] `gff.py` — load GFF/GENCODE into the universal `Grove` (hierarchy + splice-chain
      edges, CDS folded onto exons)
- [x] `prompts/system.md` — `pygenogrove` API surface + the GENCODE grove model and
      codegen rules (pinned v0.6.2 build)
- [ ] End-to-end hero query (≥ 2-hop connected-interval question) for the paper demo

## License

GPL-3.0-or-later. See [LICENSE](LICENSE).
