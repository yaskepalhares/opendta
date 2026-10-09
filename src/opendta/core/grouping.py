"""Grupos definidos por variáveis-chave (by(), collapse, contract, merge...).

`group_ids` numera os grupos 0..k-1 na ordem crescente das chaves, como o
sort os deixaria: números antes de missing (. < .a < ... < .z) e strings em
ordem de bytes. Observações fora da máscara recebem -1.
"""

from __future__ import annotations

import numpy as np

from .dataset import Dataset


def key_columns(ds: Dataset, names: list[str]) -> list[np.ndarray]:
    """Colunas float64 que ordenam como as variáveis (strings viram postos)."""
    cols = []
    for name in names:
        v = ds.get(name)
        if v.is_string:
            _, inv = np.unique(v.raw.astype(str), return_inverse=True)
            cols.append(inv.astype(np.float64).ravel())
        else:
            cols.append(np.asarray(v.data, dtype=np.float64))
    return cols


def group_ids(ds: Dataset, names: list[str], mask: np.ndarray | None = None
              ) -> tuple[np.ndarray, int]:
    """(id do grupo por observação, nº de grupos)."""
    n = ds.nobs
    mask = np.ones(n, dtype=bool) if mask is None else np.asarray(mask, dtype=bool)
    ids = np.full(n, -1, dtype=np.int64)
    if not mask.any():
        return ids, 0
    if not names:
        ids[mask] = 0
        return ids, 1
    mat = np.column_stack([c[mask] for c in key_columns(ds, names)])
    _, inv = np.unique(mat, axis=0, return_inverse=True)
    inv = inv.ravel()
    ids[mask] = inv
    return ids, int(inv.max()) + 1


def first_index(ids: np.ndarray, k: int) -> np.ndarray:
    """Índice da primeira observação de cada grupo (na ordem dos dados)."""
    first = np.full(k, -1, dtype=np.int64)
    sel = np.flatnonzero(ids >= 0)
    # percorre de trás para a frente: a última escrita é a primeira obs.
    first[ids[sel[::-1]]] = sel[::-1]
    return first


def joint_keys(a: list[np.ndarray], b: list[np.ndarray]) -> tuple[np.ndarray, np.ndarray, int]:
    """Ids comuns para chaves de dois conjuntos de dados (merge, joinby).
    Cada lista traz uma coluna por variável-chave, já comparáveis entre si
    (números como float64; strings como objetos str)."""
    na = len(a[0]) if a else 0
    cols = []
    for ca, cb in zip(a, b):
        if ca.dtype == object or cb.dtype == object:
            both = np.concatenate([ca.astype(str), cb.astype(str)])
            _, inv = np.unique(both, return_inverse=True)
            cols.append(inv.astype(np.float64).ravel())
        else:
            cols.append(np.concatenate([ca, cb]).astype(np.float64))
    if not cols:
        return np.zeros(na, dtype=np.int64), np.zeros(len(b[0]) if b else 0, dtype=np.int64), 1
    mat = np.column_stack(cols)
    if len(mat) == 0:
        return np.zeros(0, dtype=np.int64), np.zeros(0, dtype=np.int64), 0
    _, inv = np.unique(mat, axis=0, return_inverse=True)
    inv = inv.ravel().astype(np.int64)
    return inv[:na], inv[na:], int(inv.max()) + 1
