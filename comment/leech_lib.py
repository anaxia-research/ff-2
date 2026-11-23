"""Golay code and Leech minimal vectors (coordinates scaled by sqrt 8, norm 32)."""
import itertools
import numpy as np

p = 23
Q = sorted({(x * x) % p for x in range(1, p)})
INF = 23

def golay():
    rows = []
    for s in range(p):
        v = [0] * 24
        for q in Q:
            v[(q + s) % p] = 1
        v[INF] = sum(v[:23]) % 2
        rows.append(v)
    rows.append([1] * 24)
    # GF(2) basis
    basis, ints = [], []
    for r in rows:
        x = int(''.join(map(str, r)), 2)
        y = x
        for b in ints:
            y = min(y, y ^ b)
        if y:
            ints.append(y); basis.append(r)
            ints.sort(reverse=True)
    assert len(basis) == 12
    B = np.array(basis, dtype=np.int64)
    C = np.array([np.array(c) @ B % 2 for c in itertools.product([0, 1], repeat=12)], dtype=np.int8)
    return C

def leech_min(C):
    wt = C.sum(1)
    vecs = []
    for o in C[wt == 8]:
        idx = np.nonzero(o)[0]
        for signs in itertools.product([1, -1], repeat=8):
            if signs.count(-1) % 2 == 0:
                v = np.zeros(24, dtype=np.int8); v[idx] = 2 * np.array(signs); vecs.append(v)
    for i, j in itertools.combinations(range(24), 2):
        for si, sj in itertools.product([4, -4], repeat=2):
            v = np.zeros(24, dtype=np.int8); v[i] = si; v[j] = sj; vecs.append(v)
    for c in C:
        for j in range(24):
            v = np.ones(24, dtype=np.int8); v[j] = -3
            vecs.append(v * (1 - 2 * c.astype(np.int8)))
    V = np.array(vecs, dtype=np.int64)
    assert len(V) == 196560 and np.all((V * V).sum(1) == 32)
    return V

def index_of(V):
    return {r.tobytes(): i for i, r in enumerate(V.astype(np.int8))}

def as_perm(M, V, idx):
    """Permutation of the rows of V induced by x -> M x (M rational 24x24); None if V not preserved."""
    W = V @ M.T
    Wi = np.rint(W)
    if not np.allclose(W, Wi):
        return None
    Wi = Wi.astype(np.int8)
    out = np.empty(len(V), dtype=np.int64)
    for k, r in enumerate(Wi):
        j = idx.get(r.tobytes())
        if j is None:
            return None
        out[k] = j
    return out
