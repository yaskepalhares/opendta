* Mata: classes (members, methods, static members, constructors, inheritance,
* virtual methods, super and this)
mata:
class point {
    real scalar x, y
    static real scalar count
    void set()
    real scalar norm()
    virtual string scalar describe()
    void new()
}
void point::new() {
    count = (count == . ? 1 : count + 1)
}
void point::set(real scalar a, real scalar b) {
    x = a
    this.y = b
}
real scalar point::norm() {
    return(sqrt(x^2 + y^2))
}
string scalar point::describe() {
    return("point " + strofreal(norm()))
}
class point3 extends point {
    real scalar z
    virtual string scalar describe()
    real scalar norm()
}
real scalar point3::norm() {
    return(sqrt(super.norm()^2 + z^2))
}
string scalar point3::describe() {
    return("point3 " + strofreal(norm()))
}
real scalar test() {
    class point scalar p
    class point3 scalar q
    p.set(3, 4)
    q.set(1, 2)
    q.z = 2
    printf("%s | %s | %g\n", p.describe(), q.describe(), p.count)
    return(p.norm() + q.norm())
}
test()
p = point()
p.set(6, 8)
p.norm()
p.x
end
