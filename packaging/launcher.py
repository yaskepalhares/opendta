"""Entry point of the packaged application (PyInstaller).

Without arguments it opens the interface; with a file (.do or .dta), it opens the
interface and runs or loads the file. The command-line options
(--console, -b) still work for anyone who calls the executable inside
the bundle from the Terminal.
"""

import sys

from opendta.cli import main

if __name__ == "__main__":
    # older Finder versions pass -psn_0_12345 when opening the app; it is not our argument
    sys.argv = [a for a in sys.argv if not a.startswith("-psn_")]
    sys.exit(main())
