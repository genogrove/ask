# SPDX-License-Identifier: GPL-3.0-or-later
"""Enable ``python -m genogrove_ask`` as an alias for the console script."""

from genogrove_ask.cli import main

if __name__ == "__main__":
    raise SystemExit(main())
