"""Comandos do OpenDTA. Importar este pacote registra todos eles."""

from . import adocmd, data, delimited, display, filecmd, files, helpcmd, infile, inspect, logcmd, matrix, preserve, progcmd, label, macro, program, scalar, syntaxcmd, system  # noqa: F401
from .registry import REGISTRY, lookup  # noqa: F401
