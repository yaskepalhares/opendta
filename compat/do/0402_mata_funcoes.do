* Fase 4: funções do usuário, referências, structs e ponteiros no Mata.
mata
real scalar dobro(real scalar x)
{
    return(2 * x)
}
dobro(21)
void soma1(real scalar x)
{
    x = x + 1
}
y = 5
soma1(y)
y
real scalar fat(real scalar n)
{
    if (n <= 1) return(1)
    return(n * fat(n - 1))
}
fat(6)
function opc(a, | b)
{
    if (args() == 1) return(a)
    return(a + b)
}
opc(1)
opc(1, 2)
string scalar saud(string scalar nome)
{
    return("Olá, " + nome)
}
saud("Ana")
struct ponto {
    real scalar x, y
}
real scalar norma()
{
    struct ponto scalar p
    p.x = 3
    p.y = 4
    return(sqrt(p.x^2 + p.y^2))
}
norma()
a = 1
p = &a
*p = 7
a
mata drop dobro()
end
