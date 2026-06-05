# SPDX-License-Identifier: GPL-3.0-or-later
"""Command-line entry point for genogrove ask.

This is a deliberately thin wrapper: it parses the question and options, then
orchestrates the three stages — generate Python (:mod:`genogrove_ask.llm`),
execute it under restrictions (:mod:`genogrove_ask.sandbox`), and print the
result. The orchestration itself is not implemented yet (see Roadmap in README).
"""

from __future__ import annotations

import argparse
import sys

from genogrove_ask import __version__

# Default Anthropic model for code generation. Opus is the most capable tier and
# the connected-interval reasoning here is the paper's headline contribution, so
# we do not downgrade by default.
DEFAULT_MODEL = "claude-opus-4-8"


def build_parser() -> argparse.ArgumentParser:
    """Construct the argument parser for the ``genogrove-ask`` command."""
    parser = argparse.ArgumentParser(
        prog="genogrove-ask",
        description="Ask plain-English questions over connected genomic intervals.",
    )
    parser.add_argument(
        "question",
        nargs="?",
        help="The natural-language question to answer.",
    )
    parser.add_argument(
        "--model",
        default=DEFAULT_MODEL,
        help=f"Anthropic model to use for code generation (default: {DEFAULT_MODEL}).",
    )
    parser.add_argument(
        "--show-code",
        action="store_true",
        help="Print the generated Python before running it.",
    )
    parser.add_argument(
        "--version",
        action="version",
        version=f"%(prog)s {__version__}",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    """Parse arguments and dispatch. Returns a process exit code."""
    parser = build_parser()
    args = parser.parse_args(argv)

    if not args.question:
        parser.print_help()
        return 0

    # TODO(roadmap): wire the end-to-end loop:
    #   1. code = llm.generate_query(args.question, model=args.model)
    #   2. if args.show_code: print(code)
    #   3. result = sandbox.run(code)
    #   4. print(result)
    print(
        "genogrove ask is not implemented yet — this is a pre-alpha skeleton.\n"
        f"Would answer: {args.question!r} (model={args.model}).",
        file=sys.stderr,
    )
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
