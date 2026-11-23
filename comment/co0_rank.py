"""Rank of Co_0 on the 196560 minimal Leech vectors.

Builds explicit elements of Co_0 (permutations from M24, sign changes on Golay
codewords, and Conway's element xi_T), checks that each preserves the set of
minimal vectors, and computes:
  * the orbit of one minimal vector u (transitivity),
  * the orbits of the stabilizer of u, using random Schreier generators.
Since a subgroup has at least as many orbitals as the whole group, and the
inner product <u,v> already separates 7 classes, finding exactly 7 orbits of
the stabilizer proves that Co_0 has rank 7 on Lambda_4 (and rank 4 on pairs +-v).
"""
import numpy as np
from leech_lib import *

C = golay(); V = leech_min(C); idx = index_of(V)
codeset = {c.tobytes() for c in C}
N = len(V)

def perm_matrix(pi):          # coordinate permutation: new[pi[i]] = old[i]
    M = np.zeros((24, 24)); 
    for i in range(24): M[pi[i], i] = 1
    return M

def preserves_code(pi):
    return all(c[np.argsort(pi)].tobytes() in codeset for c in C)

def inv23(x): return pow(x, p - 2, p)
def mob(f):
    return [f(x) for x in range(24)]
alpha = mob(lambda x: INF if x == INF else (x + 1) % p)
beta  = mob(lambda x: INF if x == INF else (2 * x) % p)
gamma = mob(lambda x: 0 if x == INF else (INF if x == 0 else (-inv23(x)) % p))
def delta_fn(x, swap):
    if x in (0, INF): return x
    inQ = x in Q
    if inQ != swap: return (pow(x, 3, p) * inv23(9)) % p
    return (9 * pow(x, 3, p)) % p
cands = {'alpha': alpha, 'beta': beta, 'gamma': gamma,
         'delta': mob(lambda x: delta_fn(x, False)), 'delta2': mob(lambda x: delta_fn(x, True))}
gens = []
for name, pi in cands.items():
    ok = sorted(pi) == list(range(24)) and preserves_code(np.array(pi))
    print(f"permutation {name}: bijective and preserves the Golay code: {ok}")
    if ok:
        gens.append((name, perm_matrix(pi)))

# sign change on an octad
octad = C[C.sum(1) == 8][0]
gens.append(('eps_octad', np.diag(1 - 2 * octad.astype(float))))

# Conway's xi_T: sextet of a tetrad T
T = list(np.nonzero(octad)[0][:4])
tetrads = [T]
rest = set(range(24)) - set(T)
octads = [set(np.nonzero(c)[0]) for c in C[C.sum(1) == 8]]
while rest:
    a = min(rest)
    # the tetrad of a: the octad through T and a minus T
    o = [O for O in octads if set(T) <= O and a in O][0]
    tet = sorted(o - set(T)); tetrads.append(tet); rest -= set(tet)
assert len(tetrads) == 6
H = 0.5 * np.ones((4, 4)) - np.eye(4)
found = False
for signset in [[], T]:
    for Hs in [H, -H]:
        M = np.zeros((24, 24))
        for tet in tetrads:
            M[np.ix_(tet, tet)] = Hs
        D = np.ones(24); D[signset] = -1
        M = np.diag(D) @ M
        pi = as_perm(M, V, idx)
        if pi is not None:
            print("xi_T found (sign set size %d, block %s): preserves Lambda_4" % (len(signset), '+H' if Hs is H else '-H'))
            gens.append(('xi', M)); found = True; break
    if found: break
assert found

perms = []
for name, M in gens:
    assert np.allclose(M @ M.T, np.eye(24)), name
    pi = as_perm(M, V, idx); assert pi is not None, name
    perms.append(pi)
print("generators:", [g[0] for g in gens])

# orbit of u = V[0] with Schreier tree
u = 0
parent = -np.ones(N, dtype=np.int64); plabel = -np.ones(N, dtype=np.int64)
parent[u] = u; frontier = [u]; seen = 1
while frontier:
    nf = []
    for g, pi in enumerate(perms):
        img = pi[np.array(frontier)]
        new = img[parent[img] < 0]
        new, first = np.unique(new, return_index=True)
        src = np.array(frontier)[parent[pi[np.array(frontier)]] < 0][first]
        parent[new] = src; plabel[new] = g
        nf.extend(new.tolist())
    frontier = nf
orbit = int((parent >= 0).sum())
print("orbit of one minimal vector under the group generated:", orbit)
assert orbit == N

def compose(a, b):  # a after b
    return a[b]
def inverse(a):
    r = np.empty_like(a); r[a] = np.arange(len(a)); return r
ident = np.arange(N)
def transversal(j):  # element mapping u to j
    word = []
    while j != u:
        word.append(plabel[j]); j = parent[j]
    t = ident
    for g in reversed(word):
        t = compose(perms[g], t)
    return t

def orbits(gs):
    lab = np.arange(N)
    while True:
        old = lab.copy()
        for s in gs:
            lab = np.minimum(lab, lab[s])
            np.minimum.at(lab, s, lab.copy())
        lab = lab[lab]
        if np.array_equal(lab, old):
            break
    _, lab = np.unique(lab, return_inverse=True)
    return lab.max() + 1, lab

rng = np.random.default_rng(1)
def random_element(L=60):
    h = ident
    for g in rng.integers(0, len(perms), L):
        h = compose(perms[g], h)
    return h

ip = (V @ V[u]) // 8
classes = {int(k): int((ip == k).sum()) for k in sorted(set(ip.tolist()))}
print("inner-product classes w.r.t. u:", classes)

stab = []
for it in range(12):
    h = random_element()
    s = compose(inverse(transversal(h[u])), h)
    assert s[u] == u
    stab.append(s)
    k, lab = orbits(stab)
    print(f"{len(stab)} random stabilizer element(s): {k} orbits")
    if k == 7:
        break
assert k == 7
# each orbit is exactly one inner-product class
for c in range(k):
    vals = set(ip[lab == c].tolist()); assert len(vals) == 1
print("the 7 orbits of the stabilizer of u are exactly the 7 inner-product classes")
print("=> Co_0 has rank 7 on Lambda_4; on pairs {+-v} the rank is 4 (|<u,v>| = 4, 2, 1, 0)")
