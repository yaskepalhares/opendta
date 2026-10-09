* Fase 4: optimize() do Mata (registro de iterações e resultados).
mata
void q(todo, p, v, g, H)
{
    v = -(p[1] - 1)^2 - (p[2] + 2)^2 - p[1]*p[2]
}
S = optimize_init()
optimize_init_evaluator(S, &q())
optimize_init_params(S, (0, 0))
p = optimize(S)
p
optimize_result_value(S)
optimize_result_converged(S)
optimize_result_iterations(S)
optimize_result_gradient(S)
optimize_result_V(S)
void lg(todo, b, y, X, v, g, H)
{
    real colvector xb, pr
    xb = X * b'
    pr = invlogit(xb)
    v = sum(y :* ln(pr) + (1 :- y) :* ln(1 :- pr))
    if (todo >= 1) g = ((y - pr)' * X)
}
y = (0 \ 0 \ 1 \ 1 \ 0 \ 1 \ 1 \ 1)
X = ((1 \ 2 \ 3 \ 4 \ 5 \ 6 \ 7 \ 8), J(8, 1, 1))
T = optimize_init()
optimize_init_evaluator(T, &lg())
optimize_init_evaluatortype(T, "d1")
optimize_init_argument(T, 1, y)
optimize_init_argument(T, 2, X)
optimize_init_params(T, (0, 0))
b = optimize(T)
b
sqrt(diagonal(optimize_result_V(T)))'
end
