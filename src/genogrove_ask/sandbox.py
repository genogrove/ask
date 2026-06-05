# SPDX-License-Identifier: GPL-3.0-or-later
"""Sandboxed execution of LLM-generated Python.

**Security-critical.** The code passed here is produced by a language model and
must be treated as untrusted. The intended approach (to be implemented):

* Run in a separate process (``subprocess``), never in-process ``exec``.
* No network access; an allowlist of imports (``pygenogrove`` and a small set of
  stdlib/data modules); wall-clock, memory, and output-size limits.
* Read access only to registry-resolved, pinned data paths
  (see :mod:`genogrove_ask.registry`).
* Capture stdout/stderr and a structured result; surface failures to the caller
  rather than letting them escape.

Until this is implemented, :func:`run` raises rather than executing anything —
do not relax that without the restrictions above in place.
"""

from __future__ import annotations

from dataclasses import dataclass


@dataclass
class SandboxResult:
    """Outcome of running generated code in the sandbox."""

    stdout: str
    stderr: str
    returncode: int


def run(code: str, *, timeout_s: float = 30.0) -> SandboxResult:
    """Execute ``code`` under restrictions and return its captured output.

    Raises until the sandbox is implemented — generated code must never run
    without the isolation described in this module's docstring.
    """
    raise NotImplementedError("sandbox.run is not implemented yet")
