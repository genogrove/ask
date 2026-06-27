# SPDX-License-Identifier: GPL-3.0-or-later
"""Code generation via the Anthropic API.

Given a natural-language question, ask Claude to emit Python that drives
``pygenogrove`` to compute the answer. The model is given the ``pygenogrove`` API
surface and the curated resource context (see :mod:`genogrove_ask.resources`) via
the system prompt in ``prompts/system.md``.

Implementation notes for when this is built out:

* Use the official ``anthropic`` SDK (``from anthropic import Anthropic``), not
  raw HTTP.
* Default to ``claude-opus-4-8`` with adaptive thinking
  (``thinking={"type": "adaptive"}``); do not set ``temperature`` / ``budget_tokens``
  (both rejected on this model tier).
* Constrain the output to runnable Python — structured outputs
  (``output_config={"format": ...}``) or a fenced-code convention, then strip.
* The generated code is untrusted: it must only ever run through
  :mod:`genogrove_ask.sandbox`, never ``exec``'d directly here.
"""

from __future__ import annotations

DEFAULT_MODEL = "claude-opus-4-8"


def generate_query(question: str, *, model: str = DEFAULT_MODEL) -> str:
    """Translate a natural-language ``question`` into ``pygenogrove`` Python.

    Returns the generated Python source as a string. The caller is responsible
    for executing it via the sandbox.
    """
    raise NotImplementedError("llm.generate_query is not implemented yet")
