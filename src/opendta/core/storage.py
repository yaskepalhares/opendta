"""Armazenamento nativo das variáveis numéricas.

Cada variável fica na memória no próprio tipo, com 1 byte por valor em byte,
2 em int, 4 em long e float, 8 em double, como no Stata. Os missing usam a
mesma codificação do formato .dta, o que torna ler e gravar arquivos uma
cópia direta:

    byte   int8     .  = 101,         .a..z = 102..127
    int    int16    .  = 32741,       .a..z = 32742..32767
    long   int32    .  = 2147483621,  .a..z = 2147483622..2147483647
    float  float32  .  = bits 0x7f000000, .a..z = + k * 0x800
    double float64  os mesmos códigos usados no resto do OpenDTA (missing.py)

Os cálculos continuam em float64: `decode` produz o vetor float64 com os
códigos de missing do OpenDTA e `encode` faz o caminho inverso.
"""

from __future__ import annotations

import numpy as np

from . import missing as M

DTYPE = {"byte": np.dtype(np.int8), "int": np.dtype(np.int16), "long": np.dtype(np.int32),
         "float": np.dtype(np.float32), "double": np.dtype(np.float64)}
INT_MISS = {"byte": 101, "int": 32741, "long": 2147483621}
INT_RANGES = {"byte": (-127, 100), "int": (-32767, 32740), "long": (-2147483647, 2147483620)}
FLOAT_MISS_BITS = 0x7F000000
FLOAT_MISS_STEP = 0x00000800
FLOAT_MAX = 1.7014117331926443e38   # maior float não missing (bits 0x7effffff)

_DBL_BASE = np.uint64(0x7FE0000000000000)
_DBL_SHIFT = np.uint64(40)                 # passo de 0x0000010000000000 entre .a, .b, ...
_CODES = np.array([M.SYSMISS] + [M.EXTENDED[c] for c in "abcdefghijklmnopqrstuvwxyz"])


def _missing_index(v: np.ndarray, miss: np.ndarray) -> np.ndarray:
    """0 para '.', 1..26 para .a..z (valores desconhecidos viram '.')."""
    bits = v[miss].view(np.uint64)
    k = ((bits - _DBL_BASE) >> _DBL_SHIFT).astype(np.int64)
    known = (k >= 0) & (k <= 26) & (bits == _DBL_BASE + (k.astype(np.uint64) << _DBL_SHIFT))
    return np.where(known, k, 0)


def encode(values, vtype: str) -> np.ndarray:
    """float64 (códigos do OpenDTA) -> vetor no tipo nativo.

    Valores fora da faixa do tipo viram '.', e não inteiros são truncados em
    direção a zero, como em fit_numeric; normalmente os valores já chegam
    ajustados."""
    v = np.asarray(values, dtype=np.float64)
    if vtype == "double":
        return np.array(v, dtype=np.float64, copy=True)
    miss = (v >= M.SYSMISS) | (v != v)
    if vtype == "float":
        over = ~miss & (np.abs(v) > FLOAT_MAX)
        out = np.where(miss | over, 0.0, v).astype(np.float32)
        allmiss = miss | over
        if allmiss.any():
            idx = np.zeros(len(v), dtype=np.int64)
            if miss.any():
                idx[miss] = _missing_index(np.where(v != v, M.SYSMISS, v), miss)
            bits = out.view(np.uint32)
            bits[allmiss] = (FLOAT_MISS_BITS + idx[allmiss] * FLOAT_MISS_STEP).astype(np.uint32)
        return out
    lo, hi = INT_RANGES[vtype]
    over = ~miss & ((v < lo) | (v > hi))
    allmiss = miss | over
    out = np.where(allmiss, 0.0, v).astype(DTYPE[vtype])
    if allmiss.any():
        idx = np.zeros(len(v), dtype=np.int64)
        if miss.any():
            idx[miss] = _missing_index(np.where(v != v, M.SYSMISS, v), miss)
        out[allmiss] = (INT_MISS[vtype] + idx[allmiss]).astype(DTYPE[vtype])
    return out


def decode(raw: np.ndarray, vtype: str) -> np.ndarray:
    """Vetor no tipo nativo -> float64 com os códigos de missing do OpenDTA."""
    if vtype == "double":
        return np.array(raw, dtype=np.float64, copy=True)
    v = raw.astype(np.float64)
    if vtype == "float":
        bits = raw.view(np.uint32)
        # só os positivos acima do maior float válido; com sinal (bit 31) é número
        miss = (bits >= FLOAT_MISS_BITS) & (bits < 0x80000000)
        if miss.any():
            k = ((bits[miss] - FLOAT_MISS_BITS) // FLOAT_MISS_STEP).astype(np.int64).clip(0, 26)
            v[miss] = _CODES[k]
        return v
    base = INT_MISS[vtype]
    miss = raw >= base
    if miss.any():
        v[miss] = _CODES[(raw[miss].astype(np.int64) - base).clip(0, 26)]
    return v


def decode_one(x, vtype: str) -> float:
    if vtype == "double":
        return float(x)
    return float(decode(np.array([x], dtype=DTYPE[vtype]), vtype)[0])


def encode_one(x: float, vtype: str):
    return encode(np.array([x], dtype=np.float64), vtype)[0]


def missing_raw(vtype: str):
    """Valor nativo de '.' no tipo."""
    return encode_one(M.SYSMISS, vtype)
