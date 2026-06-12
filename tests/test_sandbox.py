# SPDX-License-Identifier: GPL-3.0-or-later
"""Isolation tests for :mod:`genogrove_ask.sandbox`.

These assert the guarantees in the module docstring actually hold: allowed code
runs, dangerous imports are blocked, network has no reachable path, the
known-internal posix escape is closed, and the parent-enforced caps (wall-clock,
output size, no-write, restricted reads, stripped env) behave as specified.
"""

from __future__ import annotations

import os
from pathlib import Path

import pytest

from genogrove_ask import sandbox
from genogrove_ask.sandbox import ALLOWED_IMPORTS, SandboxResult, run

pytestmark = pytest.mark.skipif(
    os.name != "posix", reason="sandbox uses POSIX process-group / rlimit isolation"
)


# --- happy path ------------------------------------------------------------- #


def test_prints_stdout() -> None:
    res = run("print('hello world')")
    assert isinstance(res, SandboxResult)
    assert res.stdout.strip() == "hello world"
    assert res.returncode == 0
    assert not res.timed_out and not res.truncated


def test_allowed_compute_import_works() -> None:
    res = run("import math\nprint(math.factorial(5))")
    assert res.returncode == 0
    assert res.stdout.strip() == "120"


def test_user_exception_is_nonzero_with_traceback() -> None:
    res = run("raise ValueError('boom')")
    assert res.returncode != 0
    assert "ValueError" in res.stderr and "boom" in res.stderr


# --- import allowlist / network -------------------------------------------- #


@pytest.mark.parametrize(
    "module",
    ["os", "posix", "socket", "subprocess", "ctypes", "shutil", "urllib", "importlib"],
)
def test_dangerous_imports_blocked(module: str) -> None:
    res = run(f"import {module}\nprint('REACHED')")
    assert res.returncode != 0, f"{module} should not import"
    assert "REACHED" not in res.stdout
    assert "blocked" in res.stderr.lower()


def test_network_has_no_reachable_path() -> None:
    # socket is the base of all networking; if it cannot be imported, nothing
    # downstream (urllib/http/ssl) can open a connection.
    res = run("import socket\ns = socket.socket()")
    assert res.returncode != 0
    assert "blocked" in res.stderr.lower()


def test_exotic_posix_reference_is_closed() -> None:
    # The import machinery preloads `posix` and holds an internal reference to
    # it; verify no module in sys.modules still exposes a posix-like `_os`.
    code = (
        "import sys\n"
        "hit = None\n"
        "for name, mod in list(sys.modules.items()):\n"
        "    o = getattr(mod, '_os', None)\n"
        "    if o is not None and getattr(o, 'system', None) is not None:\n"
        "        hit = name\n"
        "print('FOUND:' + hit if hit else 'NONE')\n"
    )
    res = run(code)
    assert res.returncode == 0
    assert res.stdout.strip() == "NONE", res.stdout


def test_allowlist_excludes_dangerous_modules() -> None:
    assert "pygenogrove" in ALLOWED_IMPORTS
    for bad in ("os", "socket", "subprocess", "ctypes", "sys", "importlib"):
        assert bad not in ALLOWED_IMPORTS


# --- filesystem ------------------------------------------------------------- #


def test_file_write_denied() -> None:
    res = run("open('/tmp/ggask_should_not_exist.txt', 'w')")
    assert res.returncode != 0
    assert "read-only" in res.stderr.lower()


def test_read_outside_roots_denied() -> None:
    res = run("open('/etc/hosts').read()")  # no data_paths -> no readable roots
    assert res.returncode != 0
    assert "restricted" in res.stderr.lower()


def test_read_inside_registry_root_allowed(tmp_path: Path) -> None:
    data = (tmp_path / "peaks.txt").resolve()
    data.write_text("chr1\t100\t200\n", encoding="utf-8")
    res = run(
        f"print(open({str(data)!r}).read().strip())",
        data_paths={"peaks": data},
    )
    assert res.returncode == 0, res.stderr
    assert res.stdout.strip() == "chr1\t100\t200"


def test_read_sibling_outside_root_denied(tmp_path: Path) -> None:
    # A path that escapes the allowed file via ".." must be refused.
    allowed = (tmp_path / "allowed.txt").resolve()
    allowed.write_text("ok\n", encoding="utf-8")
    secret = (tmp_path / "secret.txt").resolve()
    secret.write_text("nope\n", encoding="utf-8")
    res = run(
        f"print(open({str(allowed.parent / 'x' / '..' / 'secret.txt')!r}).read())",
        data_paths={"allowed": allowed},
    )
    assert res.returncode != 0
    assert "restricted" in res.stderr.lower()


# --- caps ------------------------------------------------------------------- #


def test_wall_clock_timeout_kills_runaway() -> None:
    res = run("while True:\n    pass", timeout_s=1)
    assert res.timed_out
    assert res.returncode != 0
    assert "wall-clock" in res.stderr.lower()


def test_output_is_capped() -> None:
    res = run("print('x' * 1_000_000)", output_cap=1000)
    assert res.truncated
    assert len(res.stdout.encode("utf-8")) <= 1000


# --- environment ------------------------------------------------------------ #


def test_child_env_strips_secrets(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setenv("ANTHROPIC_API_KEY", "sk-ant-secret")
    monkeypatch.setenv("AWS_SECRET_ACCESS_KEY", "also-secret")
    env = sandbox._child_env()
    assert "ANTHROPIC_API_KEY" not in env
    assert "AWS_SECRET_ACCESS_KEY" not in env
    assert "sk-ant-secret" not in env.values()
