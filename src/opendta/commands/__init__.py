"""Comandos do OpenDTA. Importar este pacote registra todos eles."""

from . import adocmd, data, delimited, display, egen, filecmd, files, helpcmd, summarize, tabulate, infile, inspect, logcmd, matrix, preserve, progcmd, label, macro, program, scalar, syntaxcmd, system  # noqa: F401
from .registry import REGISTRY, lookup  # noqa: F401
