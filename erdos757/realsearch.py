"""Exhaustive search over realizable (4,5)-sets via ordered midpoint assignments.

Points are labelled 0..n-1 in increasing order. Every AP has a unique middle (weak Sidon), and
the two extreme points are never middles (Ma-Tang Lemma 2.5). So a (4,5)-set is described by
choosing, for each interior point i, either a pair (a,b) with a < i < b and x_i = (x_a+x_b)/2,
or 'none'. We DFS over i = 1..n-2, pruning with:
  * a lower bound on the number of APs needed to reach the target tau
    (Henning-Yeo: 17 tau <= 5n + 3m, valid since AP hypergraphs are linear and F7-free),
  * linearity (two APs share <= 1 point),
  * realizability: an increasing solution must exist (LP),
  * genericity: at a random interior solution, no forced coincidence breaking (4,5)
    (forced coincidences persist in every completion, so this pruning is sound).
At a leaf we compute the actual AP hypergraph of the generic realization and its tau.
Returns the maximum tau over all (4,5)-sets of size n (target mode: whether tau >= T exists).
"""
import sys, time, itertools, random
import numpy as np
from fractions import Fraction
from scipy.optimize import linprog
import sympy as sp
from ortools.sat.python import cp_model

def tau_cp(n, edges):
    m = cp_model.CpModel(); x = [m.NewBoolVar("") for _ in range(n)]
    for e in edges: m.AddBoolOr([x[v] for v in e])
    m.Minimize(sum(x)); s = cp_model.CpSolver(); s.parameters.num_workers = 1; s.Solve(m)
    return int(round(s.ObjectiveValue()))

_T4 = {}
def ok45(x):
    n = len(x)
    if n < 4: return True
    if n not in _T4: _T4[n] = np.array(list(itertools.combinations(range(n), 4))).T
    I, J, K, L = _T4[n]; a, b, c, d = x[I], x[J], x[K], x[L]
    D = np.abs(np.stack([a-b, a-c, a-d, b-c, b-d, c-d], 1)); D.sort(1)
    return bool(((np.diff(D, axis=1) != 0).sum(1) >= 4).all())

class Search:
    def __init__(self, n, target, rng_seed=1):
        self.n = n; self.T = target; self.rng = random.Random(rng_seed)
        self.nodes = 0; self.found = None
        # minimal number of APs for tau >= T, from 17 tau <= 5n + 3m
        self.m_min = max(0, -(-(17*target - 5*n)//3))
        self.max_none = (n - 2) - self.m_min

    def realize(self, eqs):
        """eqs: list of (i,a,b). Return generic integer increasing solution or None."""
        n = self.n
        A = sp.zeros(len(eqs), n)
        for r, (i, a, b) in enumerate(eqs):
            A[r, i] = 2; A[r, a] = -1; A[r, b] = -1
        N = A.nullspace() if eqs else [sp.eye(n)[:, j] for j in range(n)]
        N = sp.Matrix.hstack(*N); d = N.shape[1]
        Nf = np.array(N.tolist(), dtype=float)
        # LP: maximize t s.t. x_{j+1}-x_j >= t, x = N c, |c| <= 1
        G = Nf[1:] - Nf[:-1]
        cobj = np.zeros(d+1); cobj[-1] = -1
        A_ub = np.hstack([-G, np.ones((n-1, 1))]); b_ub = np.zeros(n-1)
        res = linprog(cobj, A_ub=A_ub, b_ub=b_ub, bounds=[(-1, 1)]*d + [(None, 1)], method="highs")
        if res.status != 0 or -res.fun <= 1e-9: return None
        c0 = res.x[:d]; gap = -res.fun
        # exact generic point: c = round(c0*S) + small random integers
        S = 10**9
        for _ in range(3):
            c = [int(round(v*S)) + self.rng.randint(-int(gap*S/(10*d*max(1,np.abs(Nf).max()))), int(gap*S/(10*d*max(1,np.abs(Nf).max())))) for v in c0]
            xs = N * sp.Matrix(c)
            den = sp.ilcm(*[sp.fraction(v)[1] for v in xs]) if xs else 1
            xi = [int(v*den) for v in xs]
            if all(xi[j] < xi[j+1] for j in range(n-1)):
                return np.array(xi, dtype=object)
        return None

    def run(self):
        self.t0 = time.time()
        self.dfs(1, [], set(), 0)
        return self.found

    def dfs(self, i, eqs, used_pairs, nones):
        if self.found: return
        self.nodes += 1
        n = self.n
        if i == n - 1:
            x = self.realize(eqs)
            if x is None or not ok45(x): return
            S = {int(v): j for j, v in enumerate(x)}; E = []
            for p in range(n):
                for q in range(p+1, n):
                    s = int(x[p]) + int(x[q])
                    if s % 2 == 0 and s//2 in S: E.append((p, S[s//2], q))
            t = tau_cp(n, E)
            if t >= self.T: self.found = (t, [int(v) for v in x], E)
            return
        # option: i is not a midpoint
        if nones < self.max_none:
            self.dfs(i+1, eqs, used_pairs, nones+1)
        for a in range(0, i):
            for b in range(i+1, n):
                if (a, b) in used_pairs: continue
                e = {a, i, b}
                if any(len(e & {p, q, r}) >= 2 for (q, p, r) in eqs): continue
                neq = eqs + [(i, a, b)]
                x = self.realize(neq)
                if x is None or not ok45(x): continue
                self.dfs(i+1, neq, used_pairs | {(a, b)}, nones)
                if self.found: return

if __name__ == "__main__":
    n, T = int(sys.argv[1]), int(sys.argv[2])
    s = Search(n, T); t0 = time.time(); r = s.run()
    print(f"n={n} target tau>={T} (m_min={s.m_min}, max_none={s.max_none}): {'FOUND tau=%d set=%s' % (r[0], r[1]) if r else 'NONE EXISTS'}  nodes={s.nodes}  {time.time()-t0:.1f}s", flush=True)
