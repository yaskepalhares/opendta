"""Mata function library ([M-5]): the table every Mata call looks up.

The functions live in modules by topic and enter the table when the module is imported
(libbase.lib); this module imports them all, so `from .library import LIBRARY` gives the
complete table.
"""

from __future__ import annotations

from . import (  # noqa: F401 (they fill LIBRARY)
               advanced,
               lib_elementwise,
               lib_linalg,
               lib_matrix,
               lib_misc,
               lib_strings,
               moptimize,
               stata_api,
)
from .libbase import LIBRARY

__all__ = ["LIBRARY"]
