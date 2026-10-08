"""Comandos do OpenDTA. Importar este pacote registra todos eles."""

from . import data, display, inspect, label, macro, program, scalar, system  # noqa: F401
from .registry import REGISTRY, lookup  # noqa: F401
