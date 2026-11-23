"""Exact checks for the bioctonions C(x)O (Section 4).

* gamma_O is neither an automorphism nor an anti-automorphism;
* dim Sym(sigma_0) = 2, dim Sym(sigma_1) = 8;
* centre and left nucleus of C(x)O are C*1 (real dimension 2);
* the (1,3) entry of the linearised Jordan identity at x=d[11], y=a[12], z=b[23], w=1[33]
  equals -(1/8)(d,a,b)  (Lemma 4.1), for random rational d in Sym(sigma), a, b;
* the Jordan identity holds in H3(C(x)O, sigma_0) and fails in H3(C(x)O, sigma_1).
"""
import random
from fractions import Fraction as F
import sympy as sp
from bioct_lib import *

random.seed(7)
def rnd_o(): return [F(random.randint(-3, 3)) for _ in range(8)]
def rnd_b(): return (rnd_o(), rnd_o())
def basis16():
    out = []
    for part in range(2):
        for k in range(8):
            a, b = ozero(), ozero()
            (a if part == 0 else b)[k] = F(1)
            out.append((a, b))
    return out
def vec(x): return x[0] + x[1]
B16 = basis16()

# gamma_O
Ie1 = elt(im={1: 1}); I = elt(im={0: 1}); e1 = elt(re={1: 1})
lhs = gamma(bmul(I, e1)); r1 = bmul(gamma(e1), gamma(I)); r2 = bmul(gamma(I), gamma(e1))
print("gamma_O(I e1) =", show(lhs), "; gamma(e1)gamma(I) =", show(r1), "; gamma(I)gamma(e1) =", show(r2))
assert lhs != r1 and lhs != r2

# Sym dimensions
for name, sig in [("sigma_0", sigma0), ("sigma_1", sigma1)]:
    M = sp.Matrix([[vec(sig(b))[i] - vec(b)[i] for b in B16] for i in range(16)])
    print(f"dim Sym({name}) =", 16 - M.rank())
    # involution check on random pairs
    for _ in range(20):
        x, y = rnd_b(), rnd_b()
        assert sig(bmul(x, y)) == bmul(sig(y), sig(x)) and sig(sig(x)) == x
print("sigma_0 and sigma_1 are involutions (20 random pairs each)")

# centre and left nucleus
def linear_rows(fun):
    # rows of the linear map c -> fun(c) (stacked), c in R^16
    cols = [vec_all(fun(b)) for b in B16]
    return sp.Matrix(cols).T
def vec_all(lst):
    out = []
    for t in lst: out += vec(t)
    return out
comm = linear_rows(lambda c: [badd(bmul(c, x), bscale(F(-1), bmul(x, c))) for x in B16])
assoc = lambda u, v, w: badd(bmul(bmul(u, v), w), bscale(F(-1), bmul(u, bmul(v, w))))
nuc = linear_rows(lambda c: [assoc(c, x, y) for x in B16 for y in B16])
print("real dim of centre:", 16 - comm.rank(), "; real dim of left nucleus:", 16 - nuc.rank())

# Lemma 4.1 identity
J = jprod
def lin(x, y, z, w):
    tot = None
    for (p, q, r) in [(x, y, z), (y, z, x), (z, x, y)]:
        t = msub(J(p, J(J(q, r), w)), J(J(q, r), J(p, w)))
        tot = t if tot is None else [[badd(tot[i][j], t[i][j]) for j in range(3)] for i in range(3)]
    return tot
one = elt(re={0: 1})
for name, sig in [("sigma_0", sigma0), ("sigma_1", sigma1)]:
    for _ in range(5):
        d0 = rnd_b(); d = bscale(F(1, 2), badd(d0, sig(d0)))   # random element of Sym(sigma)
        a, b = rnd_b(), rnd_b()
        D = lin(herm({(0, 0): d}, sig), herm({(0, 1): a}, sig), herm({(1, 2): b}, sig), herm({(2, 2): one}, sig))
        assert D[0][2] == bscale(F(-1, 8), assoc(d, a, b))
print("Lemma 4.1: (1,3) entry = -(1/8)(d,a,b) for 5 random (d,a,b) under each involution")
d, a, b = Ie1, elt(re={2: 1}), elt(re={3: 1})
print("(I e1, e2, e3) =", show(assoc(d, a, b)))

# Jordan identity
def rnd_herm(sig):
    ent = {}
    for i in range(3):
        d0 = rnd_b(); ent[(i, i)] = bscale(F(1, 2), badd(d0, sig(d0)))
    for (i, j) in [(0, 1), (0, 2), (1, 2)]:
        ent[(i, j)] = rnd_b()
    return herm(ent, sig)
for name, sig in [("sigma_0", sigma0), ("sigma_1", sigma1)]:
    zero = 0
    for _ in range(3):
        X, Y = rnd_herm(sig), rnd_herm(sig)
        zero += iszero(jordan_defect(X, Y))
    print(f"Jordan identity on 3 random pairs in H3(C(x)O, {name}): holds in {zero} of 3")
