"""Comandos do OpenDTA. Importar este pacote registra todos eles."""

from . import adocmd, collapse, combine, data, datatools, delimited, display, egen, estcmds, filecmd, files, graphcmd, helpcmd, matacmd, reshape, summarize, tabulate, infile, inspect, logcmd, matrix, preserve, progcmd, label, macro, program, scalar, syntaxcmd, system, uicmd, nptests, immediate  # noqa: F401
from .registry import REGISTRY, lookup  # noqa: F401
