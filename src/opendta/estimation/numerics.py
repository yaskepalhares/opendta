"""set numerics stata|precise (OpenDTA).

stata (padrão): as contas seguem o Stata, inclusive nos arredondamentos
(equações normais no regress, D.x guardado em float, critérios de parada do ml).

precise: as mesmas estimativas, calculadas com mais exatidão:
  * regress resolvido por decomposição QR (o erro cresce com κ(X), não κ(X)²);
  * D.x e S.x calculados e mantidos em double;
  * ml e optimize() dão um passo de Newton a mais depois de convergir;
  * derivadas numéricas do optimize() com extrapolação de Richardson;
  * generate, egen, input... criam variáveis double quando o tipo padrão é float.
O compare.py e os testes de formato usam stata.
"""

_STATE = {"precise": False}


def set_precise(on: bool) -> None:
    _STATE["precise"] = bool(on)


def precise() -> bool:
    return _STATE["precise"]
