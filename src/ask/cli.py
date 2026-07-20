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
from pathlib import Path

from ask import __version__, llm, resources, sandbox

# Default Anthropic model for code generation. Opus is the most capable tier and
# the connected-interval reasoning here is the paper's headline contribution, so
# we do not downgrade by default.
DEFAULT_MODEL = "claude-opus-4-8"

# Datasets exposed to a query. One curated entry for now; the loop generalizes to
# the whole catalog once more resources are added.
_DATASETS = ("gencode.human",)


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
        choices=("bed", "tsv", "json"),
        default="bed",
        help="Output format for results (default: bed). Scalar answers ignore this.",
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
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    return parser


def _var_name(resource_name: str) -> str:
    """A Python identifier for a dataset's injected path variable."""
    return re.sub(r"\W", "_", resource_name).upper()


def _dataset_context(names):
    """Resolve datasets to (resources_block, code_preamble, data_paths).

    Inject one handle per dataset: the local path to its prebuilt grove ``.gg``. Every
    query — located or genome-wide — opens it with ``pg.GroveView.open`` and pages in only
    the blocks it touches, so there's no tabix GFF, no region-build, and no whole-grove
    load. ``ensure_all_grove`` has already made the ``.gg`` local (see ``main``).
    """
    block, preamble, data_paths = [], [], []
    for name in names:
        gg = str(resources._all_grove_gg(name))  # local .gg (downloaded/built in main)
        var = _var_name(name)
        desc = resources.RESOURCES[name].description
        block.append(
            f'- `{var}` (str): path to the `{name}` grove ({desc}) — open it lazily with '
            f'`g = pg.GroveView.open({var})`. GroveView pages in only the blocks a query '
            f"touches, so a **located** query (e.g. a variant at chr7:55191822) reads just that "
            f"locus, and a **genome-wide / gene-name** query works from the same handle — no "
            f"region to pick, no whole-grove load. Query-only: `intersect`, `flanking`, "
            f"`get_neighbors`, `get_edges`, `get_neighbors_if` (no `insert`/`serialize`). See "
            f'"The GENCODE Grove model" above for the node/edge structure.'
        )
        preamble.append(f"{var} = {json.dumps(gg)}")
        data_paths.append(gg)
    return "\n".join(block), "\n".join(preamble) + "\n", data_paths


def _render(text: str, fmt: str) -> str:
    """Render the generated code's stdout. JSONL feature records become ``fmt``;
    non-JSON lines (a scalar ``label: value``) pass through unchanged."""
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
            passthrough.append(line)
    out = [_format_records(records, fmt)] if records else []
    out.extend(passthrough)
    return "\n".join(p for p in out if p) + "\n"


def _format_records(records: list[dict], fmt: str) -> str:
    if fmt == "json":
        return "\n".join(json.dumps(r) for r in records)
    if fmt == "tsv":
        cols = list(records[0])
        rows = ["\t".join(cols)]
        rows += ["\t".join(str(r.get(c, "")) for c in cols) for r in records]
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


def _ensure_groves(names) -> None:
    """Make each dataset's grove ``.gg`` local (download the pinned one, or build once),
    printing a one-line notice on a real first-run fetch."""
    pending = [n for n in names if not resources._all_grove_gg(n).exists()]
    if pending:
        print(f"Fetching {', '.join(pending)} grove (first run only: a pinned ~90 MB .gg)…",
              file=sys.stderr)
    for name in names:
        resources.ensure_all_grove(name)


def _answer(question, *, system_prompt, preamble, args, execute):
    """Translate one question to code, run it via ``execute(script)``, and render.

    ``execute`` is a ``script -> SandboxResult`` callable (``sandbox.run`` for one-shot,
    ``Worker.submit`` for interactive). Returns ``(rendered_stdout, error_msg)`` — exactly
    one is non-empty.
    """
    code = llm.generate_query(question, system_prompt, model=args.model)
    if args.show_code:
        print("# --- generated code ---", file=sys.stderr)
        print(code, file=sys.stderr)
    # JSONL is the output contract, so guarantee `json` is importable even if the
    # generated code forgets the import (it's already in the allowlist).
    result = execute("import json\n" + preamble + code)
    if result.returncode != 0 or result.timed_out:
        return "", (result.stderr.strip() or "(the generated code failed with no output)")
    rendered = _render(result.stdout, args.format)
    if not rendered.strip():
        return "", "(the generated code produced no output)"
    return rendered, ""


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
                out, err = _answer(question, system_prompt=system_prompt, preamble=preamble,
                                   args=args, execute=worker.submit)
            except Exception as exc:  # e.g. an LLM error — keep the session alive
                print(f"genogrove-ask: {exc}", file=sys.stderr)
                continue
            if err:
                print(err, file=sys.stderr)
            else:
                sys.stdout.write(out)
                sys.stdout.flush()
    finally:
        worker.close()
    return 0


def main(argv: list[str] | None = None) -> int:
    """Parse arguments and run the end-to-end loop. Returns a process exit code."""
    parser = build_parser()
    args = parser.parse_args(argv)

    if args.init:  # prime the grove(s) ahead of first use, then exit
        try:
            _ensure_groves(_DATASETS)
        except Exception as exc:
            print(f"genogrove-ask: {exc}", file=sys.stderr)
            return 1
        print("Ready.", file=sys.stderr)
        return 0

    if not args.question and not args.interactive:
        parser.print_help()
        return 0

    try:
        # Every query reads the grove via GroveView, so make it local up front (pinned
        # .gg download on first run; cached after).
        _ensure_groves(_DATASETS)
        resources_block, preamble, data_paths = _dataset_context(_DATASETS)
        site_dir = _pygenogrove_site_dir()
        system_prompt = llm.build_system_prompt(resources_block)
    except Exception as exc:  # surface a clean message, not a traceback
        print(f"genogrove-ask: {exc}", file=sys.stderr)
        return 1

    if args.interactive:  # warm worker: grove open paid once for the whole session
        return _interactive(args, system_prompt=system_prompt, preamble=preamble,
                            data_paths=data_paths, site_dir=site_dir)

    try:  # one-shot: a fresh sandbox per invocation
        out, err = _answer(
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
