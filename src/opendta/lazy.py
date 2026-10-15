"""Modules imported on first use.

scipy.stats and matplotlib take most of OpenDTA's start-up time (more than a second
in the packaged application). Modules that need them only when a command runs use
lazy("scipy.stats") instead of a module-level import, so the window opens before
they are loaded."""

from __future__ import annotations

import importlib
from types import ModuleType


class _LazyModule:
    __slots__ = ("_name", "_mod")

    def __init__(self, name: str):
        self._name = name
        self._mod: ModuleType | None = None

    def __getattr__(self, attr: str):
        mod = self._mod
        if mod is None:
            mod = self._mod = importlib.import_module(self._name)
        return getattr(mod, attr)

    def __repr__(self) -> str:
        return f"<lazy module {self._name!r}>"


def lazy(name: str):
    """A stand-in for `import name` that imports the module on the first attribute access."""
    return _LazyModule(name)
