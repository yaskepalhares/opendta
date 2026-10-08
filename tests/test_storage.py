"""Armazenamento nativo (core/storage.py) e a classe Variable."""

import numpy as np
import pytest

from opendta.core import missing as M
from opendta.core import storage as S
from opendta.core.dataset import Variable

CODES = [M.SYSMISS] + [M.EXTENDED[c] for c in "abcdefghijklmnopqrstuvwxyz"]


@pytest.mark.parametrize("vtype,vals", [
    ("byte", [-127, 0, 100, -1]),
    ("int", [-32767, 0, 32740, -5]),
    ("long", [-2147483647, 0, 2147483620, -7]),
    ("float", [-1.7014117331926443e38, -0.1791305, 0.0, 1.5, 1.7014117331926443e38]),
    ("double", [-8.98846567431158e307, -0.5, 0.0, 1e-300]),
])
def test_round_trip_values_and_missing(vtype, vals):
    v = np.array(vals + CODES, dtype=np.float64)
    if vtype == "float":
        v[:len(vals)] = v[:len(vals)].astype(np.float32)
    raw = S.encode(v, vtype)
    assert raw.dtype == S.DTYPE[vtype]
    back = S.decode(raw, vtype)
    assert np.array_equal(back.view(np.int64), v.view(np.int64))


def test_disk_codes():
    assert S.encode(np.array([M.SYSMISS, M.EXTENDED["z"]]), "byte").tolist() == [101, 127]
    assert S.encode(np.array([M.EXTENDED["a"]]), "long").tolist() == [2147483622]
    bits = S.encode(np.array([M.EXTENDED["b"]]), "float").view(np.uint32)[0]
    assert bits == 0x7F000000 + 2 * 0x800


def test_out_of_range_becomes_missing():
    raw = S.encode(np.array([101.0, -128.0, 2.7, -2.7, np.nan]), "byte")
    assert S.decode(raw, "byte").tolist()[2:4] == [2.0, -2.0]
    assert all(M.is_missing(x) for x in S.decode(raw, "byte")[[0, 1, 4]])


def test_variable_api():
    v = Variable("x", "byte", np.array([1.0, 2.0, M.SYSMISS]))
    assert v.raw.dtype == np.int8 and len(v) == 3
    assert v.value(1) == 2.0 and M.is_missing(v.value(2))
    v.set_value(0, 50)
    assert v.data.tolist()[0] == 50
    v.vtype = "double"                   # promoção mantém os valores
    assert v.raw.dtype == np.float64 and v.value(0) == 50
    with pytest.raises(ValueError):
        v.data[0] = 3                     # data é só leitura
    s = Variable("s", "str3", ["a", "bb", ""])
    s.set_value(2, "c")
    assert s.data.tolist() == ["a", "bb", "c"]
