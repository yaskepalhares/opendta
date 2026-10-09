"""Comandos do OpenDTA. Importar este pacote registra todos eles."""

from . import adocmd, collapse, combine, data, datatools, delimited, display, egen, filecmd, files, helpcmd, matacmd, reshape, summarize, tabulate, infile, inspect, logcmd, matrix, preserve, progcmd, label, macro, program, scalar, syntaxcmd, system  # noqa: F401
from .registry import REGISTRY, lookup  # noqa: F401
