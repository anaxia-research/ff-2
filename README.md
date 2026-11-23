# Checks for the 27+3, GUT, Dixon–Rosenfeld and Monstrous M-theory refutations

Programs that recheck the finite arithmetic in two papers by Deep Bhattacharjee. Every proof in the papers is analytic; the programs only confirm numbers.

- `null-sym-appendix/`: the further obstructions appendix of "Residual-Dimension Obstructions for Null-Projected Supersymmetric Yang–Mills Theories" (f4(−20) signature, spinor branchings, symmetric squares, root counts).
- `comment/`: "Obstructions in Dixon–Rosenfeld and Monstrous Models" (version 2 of "Derivations, Jordan Structures and Leech Vectors in Dixon–Rosenfeld and Monstrous Models"): derivation algebras of C⊗H, C⊗O, C⊗H⊗O; Jordan closure over C⊗H and C⊗O; Leech minimal vectors, the design property, and the rank of Co_0 on them.

| Program in `comment/` | What it checks |
|---|---|
| `dixon_checks.py` (with `dixon_checks_lib.py`) | the map γ, involutions of C⊗H, derivation dimensions 6, 28, 34, the spans of the D_{a,b} |
| `dixon_der_exact.py` | dim der(H⊗O) = 17 and the centre of C⊗H⊗O, by exact null spaces |
| `bioct_checks.py` (with `bioct_lib.py`) | bioctonions in exact arithmetic: γ_O, Sym dimensions, centre and nucleus, the associator identity, the Jordan identity for σ_0 and its failure for σ_1 |
| `leech_co0.py` | Golay code, the 196560 minimal vectors, inner-product counts, the 11-design property, the subset sums |
| `co0_rank.py` (with `leech_lib.py`) | explicit elements of Co_0; transitivity on the minimal vectors; the seven orbits of the stabilizer |
| `arithmetic_checks.py` | traces, root bounds and dimension identities |

Run each script from its folder with `python3 <script>.py`; it prints the output stored in the matching `.log`. Requires Python 3 with numpy, sympy and mpmath. `dixon_der_exact.py` and `leech_co0.py` take under a minute; the rest finish in seconds.

Author: Deep Bhattacharjee. Licence: see `LICENSE`.
