"""Comandos do OpenDTA. Importar este pacote registra todos eles."""

from . import data, delimited, display, filecmd, files, infile, progcmd, inspect, label, macro, program, scalar, system  # noqa: F401
from .registry import REGISTRY, lookup  # noqa: F401
