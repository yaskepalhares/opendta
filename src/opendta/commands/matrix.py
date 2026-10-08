"""matrix ([P] matrix). Completado na fase 2e."""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from ..core.errors import StataError


@dataclass
class Matrix:
    data: np.ndarray
    rownames: list[str] = field(default_factory=list)
    colnames: list[str] = field(default_factory=list)

    @property
    def rows(self) -> int:
        return int(self.data.shape[0])

    @property
    def cols(self) -> int:
        return int(self.data.shape[1])


def take_matrix(s, text):
    raise StataError(198, "matrix support arrives in phase 2e")


def ereturn_matrix(s, sub, rest):
    raise StataError(198, "matrix support arrives in phase 2e")
