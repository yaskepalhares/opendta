"""Entry point of the packaged application (PyInstaller).

Without arguments it opens the interface; with a file (.do or .dta), it opens the
interface and runs or loads the file. The command-line options
(--console, -b) still work for anyone who calls the executable inside
the bundle from the Terminal.
"""

import os
import sys
from pathlib import Path


def _persistent_mplconfig() -> None:
    """PyInstaller's matplotlib run-time hook points MPLCONFIGDIR to a new temporary folder on
    every launch, so matplotlib scans all system fonts again each time (seconds on macOS). The
    application is a fixed folder, so a per-user cache stays valid between launches."""
    if sys.platform == "darwin":
        base = Path.home() / "Library" / "Caches" / "OpenDTA"
    elif sys.platform == "win32":
        base = Path(os.environ.get("LOCALAPPDATA") or Path.home() / "AppData" / "Local") / "OpenDTA"
    else:
        base = Path(os.environ.get("XDG_CACHE_HOME") or Path.home() / ".cache") / "opendta"
    target = base / "matplotlib"
    try:
        target.mkdir(parents=True, exist_ok=True)
    except OSError:
        return          # keep the temporary folder
    os.environ["MPLCONFIGDIR"] = str(target)


if getattr(sys, "frozen", False):
    _persistent_mplconfig()

from opendta.cli import main  # noqa: E402 - after MPLCONFIGDIR

if __name__ == "__main__":
    # older Finder versions pass -psn_0_12345 when opening the app; it is not our argument
    sys.argv = [a for a in sys.argv if not a.startswith("-psn_")]
    sys.exit(main())
