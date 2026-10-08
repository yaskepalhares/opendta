"""Comandos do OpenDTA. Importar este pacote registra todos eles."""

from . import adocmd, data, delimited, display, filecmd, files, infile, preserve, progcmd, inspect, label, macro, program, scalar, syntaxcmd, system  # noqa: F401
from .registry import REGISTRY, lookup  # noqa: F401
