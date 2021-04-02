from sympy import lcm, LM, LT, Symbol, expand, reduced, monic, groebner
from sympy.abc import x, y ,z
from itertools import product
from collections import deque

def s_polynomial(f, g):
    return expand(lcm(LM(f), LM(g)) * (1 / LT(f) * f - 1 / LT(g) * g))

def buchberger(F, need_reduce=True):
    G, pairs = list(F), deque()

    for i, f1 in enumerate(F):
        for f2 in F[i + 1:]:
            pairs.append((f1, f2))

    while pairs:
        f1, f2 = pairs.popleft()
        s = s_polynomial(f1, f2)
        h = reduced(s, G)[1]

        if h != 0:
            pairs.extend(product(G, [h]))
            G.append(h)

    if need_reduce:
        N = []
        flag = True
        while flag:
            flag = False
            for i, g in enumerate(G):
                temp = reduced(g, G[:i] + G[i + 1:])

                if temp[1] != 0:
                    N.append(monic(temp[1]))
                else:
                    G.remove(g)
                    flag = True
                    N.clear()
                    break
        return N
    return G

F = [x * y - x * z + y ** 2, y * z - x ** 2 + x ** 2 * y, x - x * y + y]

print("Initial system:", F, '\n')
G = buchberger(F, True)

print("Groebner basis is:", G)
print(groebner(F))