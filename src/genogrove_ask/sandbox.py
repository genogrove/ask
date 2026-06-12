# SPDX-License-Identifier: GPL-3.0-or-later
"""Sandboxed execution of LLM-generated Python.

**Security-critical.** The code passed to :func:`run` is produced by a language
model and is treated as untrusted. It is executed out-of-process under several
layers of restriction.

Threat model
------------
The realistic adversary here is *LLM-generated code* — code the model wrote to
answer a genomics question, possibly buggy or accidentally dangerous (a runaway
loop, an unintended network call, a stray file write), not a human deliberately
crafting an escape. The guarantees below are sized to that model.

What is enforced by the parent / the OS (a child cannot lift these)
-------------------------------------------------------------------
* **Out-of-process.** The code runs in a separate interpreter (``subprocess``),
  never via in-process ``exec``/``eval``.
* **Stripped environment.** The child gets a minimal env — no secrets such as
  ``ANTHROPIC_API_KEY`` are reachable even if the in-child guards are defeated.
* **Resource caps** (POSIX ``setrlimit`` in a pre-exec hook): CPU seconds,
  address space (memory), ``RLIMIT_FSIZE = 0`` (no file *writes* of any size),
  and an open-file-descriptor cap. These can only be *lowered* by the child.
* **Wall-clock kill.** The child runs in its own session; on timeout the whole
  process group is killed.
* **Output-size cap.** stdout/stderr are read with a byte cap so a flood cannot
  exhaust the parent's memory.

In-child guards (strong against generated code; defense-in-depth)
-----------------------------------------------------------------
A bootstrap prelude prepended to the code:

* installs an **import allowlist** (``pygenogrove`` + a small compute-only set),
* **scrubs** the dangerous primitives that the interpreter preloads (``posix``,
  ``marshal``) and **nulls** the internal ``_frozen_importlib_external._os``
  reference, which together close both the straightforward (``import os`` /
  ``socket`` / ``subprocess``) and the known-internal (``posix.system`` via the
  import machinery) escape paths — so **network has no reachable path**, and
* replaces ``open`` with a **read-only** variant restricted to registry-resolved
  data roots.

Residual risk
-------------
The in-child guards run in the same interpreter as the untrusted code, so a
*determined* adversary with arbitrary object-graph gymnastics is out of scope for
them; the parent/OS layer (env, rlimits, session-kill, output cap) is the hard
boundary. True isolation against an adversary needs an OS-level backend
(seccomp-bpf / network + mount namespaces / a container / an unprivileged user in
a jail); that is the documented next step and the architecture here keeps that
backend pluggable (it would wrap the same ``subprocess`` invocation).
"""

from __future__ import annotations

import json
import os
import subprocess
import sys
import tempfile
import threading
from collections.abc import Iterable, Mapping
from dataclasses import dataclass
from pathlib import Path

try:  # POSIX-only; resource caps are skipped (with a weaker guarantee) elsewhere.
    import resource
except ImportError:  # pragma: no cover - non-POSIX (e.g. Windows)
    resource = None  # type: ignore[assignment]


# Compute-only stdlib the generated code may use. Deliberately excludes anything
# that does I/O, networking, subprocess, or dynamic import (os, io helpers,
# socket, subprocess, ctypes, gzip/csv/pickle, importlib, ...). pygenogrove does
# its own file reading, so file/codec modules are intentionally absent. Verified
# (see tests) not to transitively pull a network/exec primitive back in.
_COMPUTE_MODULES = (
    "math",
    "itertools",
    "collections",
    "functools",
    "operator",
    "heapq",
    "bisect",
    "re",
    "json",
)

#: Top-level modules the generated code is permitted to import.
ALLOWED_IMPORTS = frozenset(("pygenogrove",) + _COMPUTE_MODULES)

#: Default caps.
DEFAULT_TIMEOUT_S = 30.0
DEFAULT_OUTPUT_CAP = 256 * 1024  # bytes, per stream
_DEFAULT_MEMORY_BYTES = 2 * 1024 * 1024 * 1024  # 2 GiB address space
_MAX_OPEN_FDS = 64


@dataclass
class SandboxResult:
    """Outcome of running generated code in the sandbox."""

    stdout: str
    stderr: str
    returncode: int
    timed_out: bool = False
    #: True if either stream hit :data:`DEFAULT_OUTPUT_CAP` and was truncated.
    truncated: bool = False


# --------------------------------------------------------------------------- #
# The in-child bootstrap. Prepended to the untrusted code; runs first and
# installs the import allowlist / read-only open before the code executes.
# --------------------------------------------------------------------------- #

_BOOTSTRAP = '''\
import sys as _sys, builtins as _builtins

_ALLOW = set(_ALLOW_JSON)
_ROOTS = tuple(_ROOTS_JSON)

# 1) Pre-warm the allowlist so no later import needs the file-import machinery
#    (which we are about to disable). A module that is not installed is skipped;
#    importing it from user code then fails cleanly via the normal import error.
for _m in _ALLOW:
    try:
        __import__(_m)
    except Exception:
        pass

# 2) Scrub the dangerous primitives (preloaded, or pulled in while warming a
#    trusted module), then null the internal os reference the import machinery
#    holds. Popping a name forces any later `import` of it back through the
#    allowlist guard below, which refuses it. Together these remove every
#    straightforward and known-internal path to posix/os/socket.
for _m in (
    "posix", "marshal", "os", "socket", "_socket",
    "subprocess", "_posixsubprocess", "ctypes", "_ctypes",
):
    _sys.modules.pop(_m, None)
try:
    import _frozen_importlib_external as _ext
    _ext._os = None
except Exception:
    pass

# 3) Import allowlist: anything already cached (interpreter internals + the
#    pre-warmed allowlist) is fine; any *new* top-level import must be allowed.
class _Guard:
    def find_spec(self, name, path=None, target=None):
        if name in _sys.modules or name.split(".")[0] in _ALLOW:
            return None
        raise ImportError("import %r is blocked in the sandbox" % name)

_sys.meta_path.insert(0, _Guard())

# 4) Read-only open, restricted to registry-resolved data roots. Writes are
#    refused here (and RLIMIT_FSIZE=0 enforces it at the OS level regardless).
_real_open = _builtins.open

def _norm(p):
    # Canonicalize an absolute path with pure string ops — no os/posixpath (which
    # would re-import os). Collapses "." and ".." so "/root/../etc" cannot escape;
    # a relative path has no leading "/" and so matches no (absolute) root.
    p = str(p)
    parts = []
    for seg in p.split("/"):
        if seg in ("", "."):
            continue
        if seg == "..":
            if parts:
                parts.pop()
            continue
        parts.append(seg)
    return ("/" if p.startswith("/") else "") + "/".join(parts)

def _guarded_open(file, mode="r", *args, **kwargs):
    if any(flag in mode for flag in ("w", "a", "x", "+")):
        raise PermissionError("the sandbox is read-only")
    target = _norm(file)
    if not any(target == r or target.startswith(r + "/") for r in _ROOTS):
        raise PermissionError(
            "reads are restricted to registry data paths: %r" % (file,)
        )
    return _real_open(file, mode, *args, **kwargs)

_builtins.open = _guarded_open

# ----- end bootstrap; untrusted code follows -----
'''


def _build_script(code: str, roots: list[str]) -> str:
    """Prepend the bootstrap (with the allowlist and data roots baked in)."""
    header = (
        f"_ALLOW_JSON = {json.dumps(sorted(ALLOWED_IMPORTS))}\n"
        f"_ROOTS_JSON = {json.dumps(roots)}\n"
    )
    return header + _BOOTSTRAP + "\n" + code


def _child_env() -> dict[str, str]:
    """A minimal environment for the child — no inherited secrets.

    Only the few variables needed for a Python interpreter to run with sane
    text handling are passed; notably ``ANTHROPIC_API_KEY`` and everything else
    in the parent environment are dropped.
    """
    env = {
        "PATH": "/usr/bin:/bin",
        "LC_ALL": "C.UTF-8",
        "LANG": "C.UTF-8",
        # Hardening flags for the child interpreter.
        "PYTHONDONTWRITEBYTECODE": "1",
        "PYTHONNOUSERSITE": "1",
    }
    return env


def _apply_limits(timeout_s: float) -> None:  # pragma: no cover - runs in child
    """Pre-exec hook (POSIX): drop into a new session and cap resources."""
    os.setsid()  # own process group, so the parent can kill the whole tree
    if resource is None:
        return
    cpu = int(timeout_s) + 1
    _setrlimit(resource.RLIMIT_CPU, cpu)
    _setrlimit(resource.RLIMIT_AS, _DEFAULT_MEMORY_BYTES)
    _setrlimit(resource.RLIMIT_FSIZE, 0)  # no file writes of any size
    _setrlimit(resource.RLIMIT_NOFILE, _MAX_OPEN_FDS)


def _setrlimit(which: int, value: int) -> None:  # pragma: no cover - runs in child
    try:
        soft, hard = resource.getrlimit(which)
        new_hard = value if hard == resource.RLIM_INFINITY else min(value, hard)
        resource.setrlimit(which, (min(value, new_hard), new_hard))
    except (ValueError, OSError):
        pass  # best effort; some limits are not settable on every platform


def _drain(stream, cap: int, sink: list[bytes], flags: dict[str, bool]) -> None:
    """Read ``stream`` into ``sink`` up to ``cap`` bytes, then keep draining."""
    total = 0
    while True:
        chunk = stream.read(65536)
        if not chunk:
            break
        if total < cap:
            room = cap - total
            sink.append(chunk[:room])
            total += min(room, len(chunk))
            if total >= cap:
                flags["truncated"] = True
        # past the cap we keep reading (and discarding) so the child does not
        # block on a full pipe.


def _normalize_roots(data_paths: Mapping[str, object] | Iterable[object] | None) -> list[str]:
    if data_paths is None:
        return []
    values = data_paths.values() if isinstance(data_paths, Mapping) else data_paths
    return [str(Path(p).resolve()) for p in values]


def run(
    code: str,
    *,
    data_paths: Mapping[str, object] | Iterable[object] | None = None,
    timeout_s: float = DEFAULT_TIMEOUT_S,
    output_cap: int = DEFAULT_OUTPUT_CAP,
) -> SandboxResult:
    """Execute untrusted ``code`` under isolation and return its captured output.

    ``data_paths`` are the registry-resolved paths the code is allowed to read
    (a mapping ``{name: path}`` or a plain iterable of paths); all other reads
    and every write are refused. See the module docstring for the full set of
    guarantees and the residual risk.
    """
    roots = _normalize_roots(data_paths)
    script = _build_script(code, roots)

    with tempfile.TemporaryDirectory(prefix="ggask-sandbox-") as tmp:
        script_path = Path(tmp) / "query.py"
        script_path.write_text(script, encoding="utf-8")

        proc = subprocess.Popen(
            [sys.executable, "-I", "-S", str(script_path)],
            stdin=subprocess.DEVNULL,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            cwd=tmp,
            env=_child_env(),
            preexec_fn=(lambda: _apply_limits(timeout_s)) if os.name == "posix" else None,
        )

        out_chunks: list[bytes] = []
        err_chunks: list[bytes] = []
        flags = {"truncated": False}
        readers = [
            threading.Thread(target=_drain, args=(proc.stdout, output_cap, out_chunks, flags)),
            threading.Thread(target=_drain, args=(proc.stderr, output_cap, err_chunks, flags)),
        ]
        for t in readers:
            t.start()

        timed_out = False
        try:
            proc.wait(timeout=timeout_s)
        except subprocess.TimeoutExpired:
            timed_out = True
            _kill(proc)
            proc.wait()
        finally:
            for t in readers:
                t.join()

    stdout = b"".join(out_chunks).decode("utf-8", "replace")
    stderr = b"".join(err_chunks).decode("utf-8", "replace")
    if timed_out:
        note = f"\n[sandbox] killed after exceeding the {timeout_s:g}s wall-clock limit"
        stderr = (stderr + note) if stderr else note.lstrip()
    return SandboxResult(
        stdout=stdout,
        stderr=stderr,
        returncode=proc.returncode if proc.returncode is not None else -1,
        timed_out=timed_out,
        truncated=flags["truncated"],
    )


def _kill(proc: subprocess.Popen) -> None:
    """Kill the child's whole process group (best effort)."""
    try:
        if os.name == "posix":
            os.killpg(os.getpgid(proc.pid), 9)
        else:  # pragma: no cover - non-POSIX
            proc.kill()
    except (ProcessLookupError, PermissionError, OSError):
        try:
            proc.kill()
        except OSError:
            pass
