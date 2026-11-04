"""Time-series estimation and test commands of phase 18, loaded on first use (start-up time)."""

from __future__ import annotations

from .registry import lazy_commands

lazy_commands("tsunit", "dfgls", "pperron")
lazy_commands("prais", "prais")


lazy_commands("varpost", "varlmar", "varnorm", "varstable", "varwle", "varbasic", "fcast")
lazy_commands("vec", "vec", "vecrank")
lazy_commands("svar", "svar")
lazy_commands("archmodel", "arch")
lazy_commands("arfima", "arfima")
lazy_commands("ucm", "ucm")
lazy_commands("sspace", "sspace", "dfactor")
lazy_commands("mgarch", "mgarch")


