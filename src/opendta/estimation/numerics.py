"""set numerics stata|precise (OpenDTA).

stata (default): computations follow Stata, including its rounding
(normal equations in regress, D.x stored as float, ml stopping criteria).

precise: the same estimates, computed more exactly:
  * regress solved by QR decomposition (the error grows with κ(X), not κ(X)²);
  * D.x and S.x computed and kept in double;
  * ml and optimize() take one more Newton step after converging;
  * optimize()'s numerical derivatives with Richardson extrapolation;
  * generate, egen, input... create double variables when the default type is float.
compare.py and the format tests use stata.
"""

_STATE = {"precise": False}


def set_precise(on: bool) -> None:
    _STATE["precise"] = bool(on)


def precise() -> bool:
    return _STATE["precise"]
