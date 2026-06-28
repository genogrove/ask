# SPDX-License-Identifier: GPL-3.0-or-later
"""Code generation via the Anthropic API.

Given a natural-language question, ask Claude to emit Python that drives
``pygenogrove`` to compute the answer. The model is given the ``pygenogrove`` API
surface, the grove model, and the curated resource context (see
:mod:`ask.resources`) via the system prompt in ``prompts/system.md``.

The generated code is untrusted: it must only ever run through :mod:`ask.sandbox`,
never ``exec``'d here.
"""

from __future__ import annotations

import re
from pathlib import Path

DEFAULT_MODEL = "claude-opus-4-8"

_SYSTEM_MD = Path(__file__).with_name("prompts") / "system.md"
# The system prompt's runtime-injected datasets block (the TODO placeholder).
_RESOURCES_HEADING = "## Available resources"


def build_system_prompt(resources_block: str, output_format: str = "bed") -> str:
    """The codegen system prompt: ``system.md`` with the resources block injected.

    ``resources_block`` replaces the placeholder under "Available resources" — it
    names the variables holding each dataset path and what they are.
    ``output_format`` is the requested result format (default ``bed``); it's stated
    so the model honors the user's ``--format`` choice (see the Rules).
    """
    text = _SYSTEM_MD.read_text(encoding="utf-8")
    head, _, _tail = text.partition(_RESOURCES_HEADING)
    return (
        f"{head}{_RESOURCES_HEADING}\n\n{resources_block.strip()}\n"
        f"\n## Output format\n\nEmit the results as `{output_format}` (see the Rules).\n"
    )


def generate_query(question: str, system_prompt: str, *, model: str = DEFAULT_MODEL) -> str:
    """Translate ``question`` into ``pygenogrove`` Python via Claude.

    Returns the generated Python source. The caller runs it through the sandbox;
    nothing is executed here. Raises ``RuntimeError`` if the model declines.
    """
    import anthropic  # lazy: keeps the module importable without the SDK/key

    client = anthropic.Anthropic()  # reads ANTHROPIC_API_KEY from the environment
    response = client.messages.create(
        model=model,
        max_tokens=16000,  # room for adaptive thinking + a small program; non-streaming-safe
        thinking={"type": "adaptive"},
        system=system_prompt,
        messages=[{"role": "user", "content": question}],
    )
    if response.stop_reason == "refusal":
        raise RuntimeError("the model declined to answer this question")
    text = "".join(b.text for b in response.content if b.type == "text")
    return _strip_code_fence(text)


def _strip_code_fence(text: str) -> str:
    """Return the Python inside a ```...``` block, or the text as-is if unfenced.

    system.md asks for a bare program, but models sometimes fence it anyway.
    """
    m = re.search(r"```(?:python)?\n(.*?)```", text, re.DOTALL)
    return (m.group(1) if m else text).strip() + "\n"
