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

This project uses [`uv`](https://docs.astral.sh/uv/). `pygenogrove` is resolved from
its GitHub repository (see `[tool.uv.sources]`) and **built from source** — it's a
C++/htslib extension, not yet on PyPI — so you need a compiler, CMake, and htslib first.

```console
# prerequisites (macOS / Homebrew). On Linux: your package manager, or see
# pygenogrove's .github/scripts/install-htslib-linux.sh
$ brew install uv htslib cmake

$ git clone https://github.com/genogrove/ask
$ cd ask
# `env VAR=… cmd` works in bash, zsh AND fish; it points CMake at htslib for the build.
$ env CMAKE_PREFIX_PATH=/opt/homebrew CMAKE_ARGS="-DCMAKE_PREFIX_PATH=/opt/homebrew/opt/htslib" uv sync
$ uv run python -c "import pygenogrove as pg; print(pg.__version__)"   # -> 0.6.2
```

The natural-language loop additionally needs `ANTHROPIC_API_KEY` set — but that part
isn't built yet, and everything in **Try it** below runs without it.

## Try it: the GFF → Grove model

The data layer (`ask.gff` + `ask.resources`) is implemented and tested — **no Claude /
API key needed.** After the `uv sync` above, run the loader tests against the real
bindings (they `importorskip`, so they actually run here rather than skip):

```console
$ uv run --extra dev pytest tests/test_gff.py tests/test_load_grove.py -q
```

(`--extra dev` pulls in `pytest`, which lives in the optional `dev` dependencies.)

Load a GENCODE locus and query it. Save this as `query.py`, then `uv run python query.py`:

```python
import pygenogrove as pg
from ask import gff, resources

# Quick single-locus load (~20s; streams the gzip once, no whole-genome build):
path = resources.resolve("gencode.human")   # downloads + sha256-verifies v50 (~70 MB) once, then caches
g = gff.load_gff(path, region=("chr7", 55_000_000, 55_300_000))

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

For repeated use, `resources.load_grove("gencode.human")` builds the **whole-genome** grove
once (a few minutes, a few GB RAM), caches it as a serialized `.gg`, and `deserialize`s in
well under a second on every later call. `load_gff(region=…)` above is the lighter path for
a one-off locus.

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
