# SPDX-License-Identifier: GPL-3.0-or-later
"""Enable ``python -m ask`` as an alias for the console script."""

from ask.cli import main

if __name__ == "__main__":
    raise SystemExit(main())
