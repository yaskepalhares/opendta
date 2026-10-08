"""Comandos do OpenDTA. Importar este pacote registra todos eles."""

from . import data, delimited, display, files, infile, inspect, label, macro, program, scalar, system  # noqa: F401
from .registry import REGISTRY, lookup  # noqa: F401
