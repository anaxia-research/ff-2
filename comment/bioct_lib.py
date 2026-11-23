"""Exact arithmetic in C(x)O and 3x3 Hermitian matrices over it."""
from fractions import Fraction as F
import itertools

# octonion table e_i e_{i+1} = e_{i+3} (indices 1..7 mod 7)
OT = {}
for i in range(1, 8):
    OT[(0, i)] = (1, i); OT[(i, 0)] = (1, i); OT[(i, i)] = (-1, 0)
OT[(0, 0)] = (1, 0)
for i in range(7):
    a, b, c = i % 7 + 1, (i + 1) % 7 + 1, (i + 3) % 7 + 1
    for x, y, z in [(a, b, c), (b, c, a), (c, a, b)]:
        OT[(x, y)] = (1, z); OT[(y, x)] = (-1, z)

def omul(x, y):
    z = [F(0)] * 8
    for i in range(8):
        if x[i] == 0: continue
        for j in range(8):
            if y[j] == 0: continue
            s, k = OT[(i, j)]; z[k] += s * x[i] * y[j]
    return z
def oadd(x, y): return [a + b for a, b in zip(x, y)]
def oneg(x): return [-a for a in x]
def oconj(x): return [x[0]] + [-a for a in x[1:]]
def ozero(): return [F(0)] * 8
def obasis(i): v = ozero(); v[i] = F(1); return v

# C(x)O element: pair (a, b) meaning a + I b
def bmul(x, y):
    a, b = x; c, d = y
    return (oadd(omul(a, c), oneg(omul(b, d))), oadd(omul(a, d), omul(b, c)))
def badd(x, y): return (oadd(x[0], y[0]), oadd(x[1], y[1]))
def bscale(t, x): return ([t * u for u in x[0]], [t * u for u in x[1]])
def bzero(): return (ozero(), ozero())
def sigma0(x): return (oconj(x[0]), oconj(x[1]))           # C-linear
def sigma1(x): return (oconj(x[0]), oneg(oconj(x[1])))     # C-antilinear
def gamma(x):  # the map of CMCA2025 eq. (2.5): conj z0 - z_k e_k
    a, b = x
    return ([a[0]] + [-u for u in a[1:]], [-b[0]] + [-u for u in b[1:]])

def mmul(X, Y):
    Z = [[bzero() for _ in range(3)] for _ in range(3)]
    for i in range(3):
        for j in range(3):
            for k in range(3):
                Z[i][j] = badd(Z[i][j], bmul(X[i][k], Y[k][j]))
    return Z
def jprod(X, Y):
    A, B = mmul(X, Y), mmul(Y, X)
    return [[bscale(F(1, 2), badd(A[i][j], B[i][j])) for j in range(3)] for i in range(3)]
def msub(X, Y): return [[badd(X[i][j], bscale(F(-1), Y[i][j])) for j in range(3)] for i in range(3)]
def iszero(X): return all(all(u == 0 for u in X[i][j][0] + X[i][j][1]) for i in range(3) for j in range(3))
def herm(entries, sig):
    """entries: dict {(i,j): element} for i<=j; fills (j,i) with sig."""
    X = [[bzero() for _ in range(3)] for _ in range(3)]
    for (i, j), v in entries.items():
        X[i][j] = v
        if i != j: X[j][i] = sig(v)
    return X
def jordan_defect(X, Y):
    X2 = jprod(X, X)
    return msub(jprod(jprod(X2, Y), X), jprod(X2, jprod(Y, X)))
def elt(re=None, im=None):
    a, b = ozero(), ozero()
    for k, c in (re or {}).items(): a[k] = F(c)
    for k, c in (im or {}).items(): b[k] = F(c)
    return (a, b)
def show(x):
    names = ['1'] + ['e%d' % i for i in range(1, 8)]
    t = []
    for part, pre in [(x[0], ''), (x[1], 'I')]:
        for k, c in enumerate(part):
            if c != 0: t.append(f"{c}{'*' if pre or k else ''}{pre}{names[k] if k else ('' if pre else '')}")
    return ' + '.join(t) if t else '0'
