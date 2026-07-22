# SPDX-License-Identifier: GPL-3.0-or-later
"""Command-line entry point for genogrove ask.

A thin wrapper: parse the question, then orchestrate the three stages — generate
Python (:mod:`ask.llm`), execute it under restrictions (:mod:`ask.sandbox`), and
print the result. The host resolves each dataset to a serialized ``.gg`` and
injects its path as a variable; the generated code only deserializes and queries.
"""

from __future__ import annotations

import argparse
import json
import re
import sys
import time
from pathlib import Path

from ask import __version__, llm, resources, sandbox

# Default Anthropic model for code generation. Opus is the most capable tier and
# the connected-interval reasoning here is the paper's headline contribution, so
# we do not downgrade by default.
DEFAULT_MODEL = "claude-opus-4-8"

# The base annotation grove. The regulatory (enhancer→gene) layer is augmented onto it
# per cohort, on demand — see the rE2G helpers in ``ask.resources``.
_BASE = "gencode.human"

# When a question needs enhancers but names no tissue, augment with this cohort and say so.
DEFAULT_COHORT = "EFO:0005726"  # LNCaP clone FGC (prostate cancer) — the flagship cohort

# A question wants the regulatory layer if it mentions it. Cheap intent check (no extra LLM
# call): gene/transcript/exon questions stay on plain GENCODE; only these pull a cohort.
_ENHANCER_HINT = re.compile(r"enhanc|regulat|cis-?reg|\bccre\b|\bcre\b|regulatory element", re.I)


def _wants_enhancers(question: str) -> bool:
    """True if the question asks about the regulatory layer (enhancers / regulation)."""
    return bool(_ENHANCER_HINT.search(question or ""))


def build_parser() -> argparse.ArgumentParser:
    """Construct the argument parser for the ``genogrove-ask`` command."""
    parser = argparse.ArgumentParser(
        prog="genogrove-ask",
        description="Ask plain-English questions over connected genomic intervals.",
    )
    parser.add_argument("question", nargs="?", help="The natural-language question to answer.")
    parser.add_argument(
        "--model",
        default=DEFAULT_MODEL,
        help=f"Anthropic model to use for code generation (default: {DEFAULT_MODEL}).",
    )
    parser.add_argument(
        "--format",
        choices=("text", "bed", "tsv", "json"),
        default="text",
        help="Output format (default: text — an aligned table showing every field, led by "
             "the summary line; bed/tsv/json are machine formats). Scalar answers ignore this.",
    )
    parser.add_argument(
        "--show-code",
        action="store_true",
        help="Print the generated Python before running it.",
    )
    parser.add_argument(
        "--init",
        action="store_true",
        help="Download the dataset grove(s) now (the pinned ~90 MB .gg) and exit, so the "
             "first real query is instant.",
    )
    parser.add_argument(
        "-i", "--interactive",
        action="store_true",
        help="Interactive session: keep the grove(s) open across questions (the ~200 ms "
             "open is paid once, then queries are sub-ms). One question per line.",
    )
    parser.add_argument(
        "--cohort",
        action="append",
        metavar="NAME",
        help="Load the enhancer→gene layer for this ENCODE-rE2G cohort (name or ontology id, "
             "e.g. 'LNCaP' or 'EFO:0005726'; repeatable for several). Enables enhancer "
             "queries. Omit and enhancers still load for an enhancer question, defaulting to "
             f"the flagship cohort. See --list-cohorts.",
    )
    parser.add_argument(
        "--list-cohorts",
        action="store_true",
        help="List the available ENCODE-rE2G cohorts (name, ontology id, type, replicates) "
             "and exit.",
    )
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    return parser


def _resolve_cohorts(specs):
    """Map ``--cohort`` specs to ``{cohort name: [accessions]}`` via ``re2g_cohorts``.

    Each spec matches a cohort by exact ontology id or case-insensitive substring of its
    name; the most-replicated match wins (the catalog is sorted by replicate count). Raises
    ``SystemExit`` with a pointer to ``--list-cohorts`` on no match.
    """
    catalog = resources.re2g_cohorts()
    chosen = {}
    for spec in specs:
        s = spec.strip().lower()
        match = next((c for c in catalog if c["ontology_id"].lower() == s), None) or \
            next((c for c in catalog if s in c["name"].lower()), None)
        if match is None:
            raise SystemExit(f"genogrove-ask: no cohort matches {spec!r} — see --list-cohorts")
        chosen[match["name"]] = match["accessions"]
    return chosen


def _list_cohorts() -> None:
    """Print the available rE2G cohorts (most-replicated first) to stdout."""
    print(f"{'ontology id':16}  {'reps':>4}  {'type':16}  name")
    for c in resources.re2g_cohorts():
        print(f"{c['ontology_id']:16}  {c['n_replicates']:>4}  {c['type']:16}  {c['name']}")


def _grove_context(cohorts):
    """Resolve the grove to (resources_block, code_preamble, data_paths, note).

    ``cohorts`` is ``{name: [accessions]}`` (empty → plain GENCODE). A non-empty set augments
    GENCODE with those cohorts' rE2G enhancer→gene edges (built/cached on first use) and adds
    the regulatory layer to the prompt. One handle, ``GENCODE_HUMAN``, is injected either way;
    the query opens it lazily with ``pg.GroveView.open``. ``note`` (for the user) names the
    cohorts, or is empty for a plain-GENCODE run.
    """
    var = "GENCODE_HUMAN"
    if cohorts:
        gg = str(resources.ensure_augmented_grove(_BASE, cohorts))
        labels = "; ".join(cohorts)
        block = (
            f"- `{var}` (str): a **combined grove** — GENCODE gene/transcript/exon structure "
            f"PLUS the ENCODE-rE2G enhancer→gene layer for cohort(s): **{labels}**. Open it "
            f"lazily with `g = pg.GroveView.open({var})`.\n"
            f"  - **Enhancer nodes** carry `{{'type':'enhancer','class':...}}` and are spatially "
            f"indexed, so `intersect` at a variant/region returns them alongside genes.\n"
            f"  - **Edges** (filter by `rel` with `get_neighbors_if`): `regulates` "
            f"(enhancer→gene) and its reverse `regulated_by` (gene→enhancer, for 'enhancers of "
            f"a gene'). Each carries `byCohort`: `{{'<cohort>': {{'score': <rE2G confidence>, "
            f"'n': <replicates supporting it>}}}}`. Rank/threshold on `score`; use `n` for "
            f"confidence; a link with `class=='promoter'` and ~0 distance is a self-promoter.\n"
            f"  - GENCODE structure is unchanged (see \"The GENCODE Grove model\" above): a "
            f"gene node reached via `regulated_by` still has its `first_exon`/`next` chain.\n"
            f"  - **Return the evidence, not bare intervals.** Each enhancer result record must "
            f"carry the connection: `type` (`\"enhancer\"`), `class`, the connected gene's "
            f"`name`, and from the relevant `byCohort` entry the `score` (put it in a `score` "
            f"field) and `n` and cohort label. Give it a descriptive `name` too, e.g. "
            f"`f\"enh:{{cls}}→{{gene}}\"`. A `.`-only interval with no score/target is not an "
            f"answer. Sort by `n` then `score` so the confident links lead."
        )
        note = f"Enhancers loaded for cohort(s): {labels}."
    else:
        gg = str(resources.ensure_all_grove(_BASE))
        block = (
            f"- `{var}` (str): path to the GENCODE grove "
            f"({resources.RESOURCES[_BASE].description}) — open it lazily with "
            f"`g = pg.GroveView.open({var})`. A **located** query (a variant at chr7:55191822) "
            f"reads just that locus; a **genome-wide / gene-name** query works from the same "
            f"handle. Query-only: `intersect`, `flanking`, `get_neighbors`, `get_edges`, "
            f'`get_neighbors_if`. See "The GENCODE Grove model" above for the node/edge structure.'
        )
        note = ""
    return block, f"{var} = {json.dumps(gg)}\n", [gg], note


def _render(text: str, fmt: str) -> str:
    """Render the generated code's stdout. Non-JSON lines (the agent's ``label: value``
    summary) **lead**, then the JSONL feature records become the chosen ``fmt`` table."""
    records, passthrough = [], []
    for line in text.splitlines():
        s = line.strip()
        if not s:
            continue
        try:
            obj = json.loads(s)
        except ValueError:
            obj = None
        if isinstance(obj, dict):
            records.append(obj)
        else:
            passthrough.append(line)  # a summary / scalar line — shown before the table
    out = list(passthrough)
    if records:
        out.append(_format_records(records, fmt))
    return "\n".join(p for p in out if p) + "\n"


def _record_columns(records: list[dict]) -> list[str]:
    """Column order across (possibly heterogeneous) records: coordinates, then identity,
    then edge evidence, then anything else — union of keys, stable order."""
    order = ["chrom", "start", "end", "strand", "name", "type", "class",
             "score", "n", "cohort", "target", "id", "biotype"]
    cols = [k for k in order if any(k in r for r in records)]
    for r in records:  # append any keys the agent used that aren't in the preferred order
        for k in r:
            if k not in cols:
                cols.append(k)
    return cols


def _format_records(records: list[dict], fmt: str) -> str:
    if fmt == "json":
        return "\n".join(json.dumps(r) for r in records)
    cols = _record_columns(records)
    if fmt in ("text", "tsv"):
        if fmt == "tsv":
            rows = ["\t".join(cols)]
            rows += ["\t".join(str(r.get(c, "")) for c in cols) for r in records]
            return "\n".join(rows)
        # text: an aligned, human-readable table (every field visible), padded per column.
        width = {c: max(len(c), max((len(str(r.get(c, ""))) for r in records), default=0)) for c in cols}
        fmt_row = lambda vals: "  ".join(str(v).ljust(width[c]) for c, v in zip(cols, vals))
        rows = [fmt_row(cols)]
        rows += [fmt_row([r.get(c, "") for c in cols]) for r in records]
        return "\n".join(rows)
    # BED: 0-based closed -> half-open (end + 1); host owns the conversion, once.
    rows = ["#chrom\tstart\tend\tname\tscore\tstrand"]
    for r in records:
        if "start" not in r or "end" not in r:
            rows.append(json.dumps(r))  # not an interval record; emit verbatim
            continue
        rows.append("\t".join(str(v) for v in (
            r.get("chrom", "."), r["start"], int(r["end"]) + 1,
            r.get("name") or r.get("id") or ".", r.get("score", "."),
            r.get("strand", "."),
        )))
    return "\n".join(rows)


def _pygenogrove_site_dir() -> str:
    """The site-packages dir holding ``pygenogrove``, for the sandbox's sys.path."""
    import pygenogrove

    f = Path(pygenogrove.__file__).resolve()
    return str(f.parent.parent if f.name == "__init__.py" else f.parent)


def _cohorts_for(args):
    """Decide which cohorts to load: explicit ``--cohort``, else the default cohort when a
    one-shot question asks about enhancers, else none (plain GENCODE). Returns
    ``{cohort name: [accessions]}`` (empty for a plain-GENCODE run)."""
    if args.cohort:
        return _resolve_cohorts(args.cohort)
    if not args.interactive and _wants_enhancers(args.question):
        return _resolve_cohorts([DEFAULT_COHORT])  # enhancer question, no tissue named
    return {}


def _prepare(cohorts) -> None:
    """Print a first-run notice if the grove(s) this run needs aren't cached yet."""
    if not resources._all_grove_gg(_BASE).exists():
        print(f"Fetching {_BASE} grove (first run only: a pinned ~90 MB .gg)…", file=sys.stderr)
    if cohorts and not resources.augmented_grove_path(_BASE, cohorts).exists():
        print(f"Building the enhancer grove for {', '.join(cohorts)} "
              "(first run: ENCODE download + augment)…", file=sys.stderr)


def _answer(question, *, system_prompt, preamble, args, execute):
    """Translate one question to code, run it via ``execute(script)``, and render.

    ``execute`` is a ``script -> SandboxResult`` callable (``sandbox.run`` for one-shot,
    ``Worker.submit`` for interactive). Returns ``(rendered_stdout, error_msg, gen_s, exec_s)``
    — exactly one of stdout/error is non-empty; the two times split code-gen from execution.
    """
    t0 = time.perf_counter()
    code = llm.generate_query(question, system_prompt, model=args.model)
    gen_s = time.perf_counter() - t0
    if args.show_code:
        print("# --- generated code ---", file=sys.stderr)
        print(code, file=sys.stderr)
    # JSONL is the output contract, so guarantee `json` is importable even if the
    # generated code forgets the import (it's already in the allowlist).
    t1 = time.perf_counter()
    result = execute("import json\n" + preamble + code)
    exec_s = time.perf_counter() - t1
    if result.returncode != 0 or result.timed_out:
        return "", (result.stderr.strip() or "(the generated code failed with no output)"), gen_s, exec_s
    rendered = _render(result.stdout, args.format)
    if not rendered.strip():
        return "", "(the generated code produced no output)", gen_s, exec_s
    return rendered, "", gen_s, exec_s


def _interactive(args, *, system_prompt, preamble, data_paths, site_dir) -> int:
    """Warm-worker REPL: open the grove(s) once, then answer questions until EOF/'exit'."""
    worker = sandbox.Worker(data_paths=data_paths, extra_syspath=[site_dir])
    print("genogrove-ask interactive — one question per line; Ctrl-D or 'exit' to quit.",
          file=sys.stderr)
    try:
        while True:
            try:
                question = input("ask> ").strip()
            except (EOFError, KeyboardInterrupt):
                print(file=sys.stderr)
                break
            if not question:
                continue
            if question in ("exit", "quit"):
                break
            try:
                out, err, gen_s, exec_s = _answer(question, system_prompt=system_prompt,
                                                  preamble=preamble, args=args, execute=worker.submit)
            except Exception as exc:  # e.g. an LLM error — keep the session alive
                print(f"genogrove-ask: {exc}", file=sys.stderr)
                continue
            if err:
                print(err, file=sys.stderr)
            else:
                sys.stdout.write(out)
                sys.stdout.flush()
            print(f"({gen_s + exec_s:.2f}s  llm {gen_s:.2f}s · grove {exec_s:.3f}s)", file=sys.stderr)
    finally:
        worker.close()
    return 0


def main(argv: list[str] | None = None) -> int:
    """Parse arguments and run the end-to-end loop. Returns a process exit code."""
    parser = build_parser()
    args = parser.parse_args(argv)

    if args.list_cohorts:
        _list_cohorts()
        return 0

    if args.init:  # prime the base grove ahead of first use, then exit
        try:
            _prepare({})
            resources.ensure_all_grove(_BASE)
        except Exception as exc:
            print(f"genogrove-ask: {exc}", file=sys.stderr)
            return 1
        print("Ready.", file=sys.stderr)
        return 0

    if not args.question and not args.interactive:
        parser.print_help()
        return 0

    try:
        # Pick the cohort(s) this run needs, then make the grove local (base .gg download +,
        # if a cohort is needed, the augmented build). The query opens it via GroveView.
        cohorts = _cohorts_for(args)
        _prepare(cohorts)
        resources_block, preamble, data_paths, note = _grove_context(cohorts)
        site_dir = _pygenogrove_site_dir()
        system_prompt = llm.build_system_prompt(resources_block)
    except SystemExit:
        raise  # a clean --cohort resolution error already carries its message
    except Exception as exc:  # surface a clean message, not a traceback
        print(f"genogrove-ask: {exc}", file=sys.stderr)
        return 1

    if note:  # tell the user which cohort's enhancers are in play (esp. the default)
        default_used = not args.cohort
        print(note + (" (default — pass --cohort to choose another; --list-cohorts to see them)"
                      if default_used else ""), file=sys.stderr)

    if args.interactive:  # warm worker: grove open paid once for the whole session
        return _interactive(args, system_prompt=system_prompt, preamble=preamble,
                            data_paths=data_paths, site_dir=site_dir)

    try:  # one-shot: a fresh sandbox per invocation
        out, err, _gen_s, _exec_s = _answer(
            args.question, system_prompt=system_prompt, preamble=preamble, args=args,
            execute=lambda s: sandbox.run(s, data_paths=data_paths, extra_syspath=[site_dir]),
        )
    except Exception as exc:
        print(f"genogrove-ask: {exc}", file=sys.stderr)
        return 1
    if err:
        print(err, file=sys.stderr)
        return 1
    sys.stdout.write(out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
