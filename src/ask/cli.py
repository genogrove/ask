# SPDX-License-Identifier: GPL-3.0-or-later
"""Command-line entry point for genogrove ask.

A thin wrapper: parse the question, then orchestrate the three stages — generate
Python (:mod:`ask.llm`), execute it under restrictions (:mod:`ask.sandbox`), and
print the result. The host resolves each dataset to a serialized ``.gg`` and
injects its path as a variable; the generated code only deserializes and queries.
"""

from __future__ import annotations

import argparse
import inspect
import json
import re
import sys
from pathlib import Path

from ask import __version__, gff, llm, resources, sandbox

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
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    return parser


def _var_name(resource_name: str) -> str:
    """A Python identifier for a dataset's injected path variable."""
    return re.sub(r"\W", "_", resource_name).upper()


def _dataset_context(names):
    """Resolve datasets to (resources_block, code_preamble, data_paths).

    Inject the bgzip+tabix GFF path (for fast region reads via the injected
    ``build_grove`` helper) and the whole-genome grove path (for genome-wide
    queries, built lazily). Whitelist every path for the sandbox.
    """
    block, preamble, data_paths = [], [], []
    for name in names:
        gff_path = str(resources.indexed_path(name))  # bgzip+tabix; built once
        all_gg = str(resources._all_grove_gg(name))    # path only — built lazily if referenced
        var = _var_name(name)
        desc = resources.RESOURCES[name].description
        block.append(
            f'- `{var}` (str): path to a tabix-indexed `{name}` GFF ({desc}) Build a '
            f'region-restricted Grove with `g = build_grove({var}, "chr7:55000000-55300000")` — '
            f"the region is a tabix string (**1-based inclusive**), and only features overlapping "
            f"it are loaded (fast). Pick a region covering what the query needs: a point for "
            f"\"what overlaps here\", a gene's span for its full exon/CDS structure.\n"
            f"- `{var}_ALL` (str): path to the **whole-genome** `{name}` Grove — "
            f"`g = pg.Grove.deserialize({var}_ALL)`. Use it ONLY for genome-wide queries or "
            f"gene-name lookups with no known locus (it reads everything; slower). Prefer "
            f"`build_grove` with a region whenever the query names a locus.\n"
            f'Both groves have the structure in "The GENCODE Grove model" above.'
        )
        preamble.append(f"{var} = {json.dumps(gff_path)}")
        preamble.append(f"{var}_ALL = {json.dumps(all_gg)}")
        data_paths += [gff_path, gff_path + ".tbi", all_gg]
    helper = inspect.getsource(gff.build_grove)  # pure-pygenogrove region loader, runs in the sandbox
    return "\n".join(block), helper + "\n\n" + "\n".join(preamble) + "\n", data_paths


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


def main(argv: list[str] | None = None) -> int:
    """Parse arguments and run the end-to-end loop. Returns a process exit code."""
    parser = build_parser()
    args = parser.parse_args(argv)

    if not args.question:
        parser.print_help()
        return 0

    try:
        unindexed = [n for n in _DATASETS if not resources.is_indexed(n)]
        if unindexed:
            print(f"Fetching {', '.join(unindexed)} (first run only: downloads the "
                  "indexed dataset, a few hundred MB)…", file=sys.stderr)
        resources_block, preamble, data_paths = _dataset_context(_DATASETS)
        site_dir = _pygenogrove_site_dir()

        system_prompt = llm.build_system_prompt(resources_block)
        code = llm.generate_query(args.question, system_prompt, model=args.model)
        if args.show_code:
            print("# --- generated code ---", file=sys.stderr)
            print(code, file=sys.stderr)

        # Genome-wide queries reference <VAR>_ALL — build that whole-genome grove
        # lazily, only now that we know the query needs it. Located queries skip it.
        for name in _DATASETS:
            if f"{_var_name(name)}_ALL" in code:
                print(f"Building whole-genome grove for {name} (genome-wide query)…", file=sys.stderr)
                resources.ensure_all_grove(name)

        # JSONL is the output contract, so guarantee `json` is importable even if
        # the generated code forgets the import (it's already in the allowlist).
        script = "import json\n" + preamble + code
        result = sandbox.run(script, data_paths=data_paths, extra_syspath=[site_dir])
    except Exception as exc:  # surface a clean message, not a traceback
        print(f"genogrove-ask: {exc}", file=sys.stderr)
        return 1

    if result.returncode != 0 or result.timed_out:
        print(result.stderr.strip() or "(the generated code failed with no output)", file=sys.stderr)
        return 1
    rendered = _render(result.stdout, args.format)
    if not rendered.strip():
        print("(the generated code produced no output)", file=sys.stderr)
        return 1
    sys.stdout.write(rendered)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
