from sympy import Poly, interpolate, factor, divisors, primitive
from sympy.abc import x
from itertools import product as cart_prod

def factorize_kronecker(poly: Poly):
    res = []
    con, poly = primitive(poly)

    if con > 1: res.append(con)
    while True:
        temp_res = kronecker(poly)
        if temp_res is None:
            res.append(poly)
            return res
        else:
            res.append(temp_res[0])
            poly = poly.div(temp_res[0])[0]

def kronecker(poly: Poly):
    poly_factors_in_point = lambda poly, p: divisors(poly.subs(x, p))
    n = poly.degree()

    for i in range(n // 2 + 1):
        if poly.subs(x, i) == 0:
            return Poly(x - i, x), 1

    U = poly_factors_in_point(poly, 0)
    for i in range(1, n // 2 + 1):
        M = poly_factors_in_point(poly, i)
        U = cart_prod(U, M)

        for seq in U:
            points = [(x, y) for x, y in zip(range(i+1), seq)]
            poly_g = Poly(interpolate(points, x), x)

            if poly_g.degree() == i and poly.div(poly_g)[1] == 0:
                return poly_g, i

#poly = Poly(x ** 5  + x ** 4 - 21 * x ** 3 - 45 * x ** 2, x)
poly = Poly(2 * x ** 4 + 14 * x ** 3 + 12 * x ** 2 - 56 * x - 80, x)

print("Polynomial: ", poly)
print(poly.factor_list(), '\n')
print(factorize_kronecker(poly))
print()