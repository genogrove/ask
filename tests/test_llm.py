# SPDX-License-Identifier: GPL-3.0-or-later
"""Pure codegen helpers — no API key or pygenogrove needed."""

from ask import llm
from ask.cli import _var_name


def test_build_system_prompt_injects_resources_block():
    prompt = llm.build_system_prompt("- `GENCODE_HUMAN` (str): a grove path")
    assert "## Available resources" in prompt
    assert "GENCODE_HUMAN" in prompt
    assert "## The GENCODE Grove model" in prompt  # earlier sections preserved
    assert "TODO: injected at runtime" not in prompt  # placeholder dropped
    assert prompt.index("GENCODE Grove model") < prompt.index("Available resources")


def test_strip_code_fence():
    assert llm._strip_code_fence("```python\nprint(1)\n```") == "print(1)\n"
    assert llm._strip_code_fence("```\nx = 1\n```") == "x = 1\n"
    assert llm._strip_code_fence("print(2)") == "print(2)\n"  # unfenced passes through


def test_var_name():
    assert _var_name("gencode.human") == "GENCODE_HUMAN"
