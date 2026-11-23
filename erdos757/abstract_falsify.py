"""Stage 5 falsification: search the ABSTRACT class for tau > 3n/7.
Class: 3-uniform, linear, F7-free, hereditarily sparse (|E(S)| <= |S|-2 for all S spanning an edge).
Sparsity test (exact): for each edge e, the edge multiset E + 2 extra copies of e must have a
system of distinct representatives (each edge -> one of its vertices, injective)  [Hall].
Objective: tau (CP-SAT) + small LP-relaxation tiebreak. Simulated annealing, many restarts.
"""
import sys, random, math, time, itertools, json
import numpy as np
from scipy.optimize import linprog
from ortools.sat.python import cp_model
sys.path.insert(0, '.')
from c5b_n_search import f7_copy

def sdr(edges):
    match = {}  # vertex -> edge index
    def aug(i, seen):
        for v in edges[i]:
            if v in seen: continue
            seen.add(v)
            if v not in match or aug(match[v], seen):
                match[v] = i; return True
        return False
    return all(aug(i, set()) for i in range(len(edges)))

def sparse(E):
    for e in E:
        if not sdr(E + [e, e]): return False
    return True

def linear_with(E, t):
    return all(len(set(t) & set(f)) <= 1 for f in E)

def tau(n, E):
    m = cp_model.CpModel(); x = [m.NewBoolVar("") for _ in range(n)]
    for e in E: m.AddBoolOr([x[v] for v in e])
    m.Minimize(sum(x)); s = cp_model.CpSolver(); s.parameters.num_workers = 1; s.Solve(m)
    return int(round(s.ObjectiveValue()))

def lp(n, E):
    if not E: return 0.0
    A = np.zeros((len(E), n))
    for k, e in enumerate(E):
        for v in e: A[k, v] = -1
    return linprog(np.ones(n), A_ub=A, b_ub=-np.ones(len(E)), bounds=(0, 1), method="highs").fun

def valid_add(E, t):
    return linear_with(E, t) and sparse(E + [t]) and f7_copy(E + [t]) is None

def anneal(n, iters, seed):
    rng = random.Random(seed)
    T = list(itertools.combinations(range(n), 3))
    E = []
    for _ in range(4000):
        if len(E) >= n - 2: break
        t = rng.choice(T)
        if t not in E and valid_add(E, t): E.append(t)
    def score(E): t = tau(n, E); return t, t + 0.3 * lp(n, E)
    ct, cs = score(E); best = (ct, list(E)); Temp = 0.4
    for it in range(iters):
        F = list(E)
        if F and (len(F) >= n - 2 or rng.random() < 0.5):
            F.pop(rng.randrange(len(F)))
        for _ in range(40):
            t = rng.choice(T)
            if t not in F and valid_add(F, t): F.append(t); break
        if F == E: continue
        nt, ns = score(F)
        if ns >= cs or rng.random() < math.exp((ns - cs) / Temp):
            E, ct, cs = F, nt, ns
            if ct > best[0]: best = (ct, list(E))
        Temp = max(0.03, Temp * 0.999)
    return best

if __name__ == "__main__":
    n0, n1, iters, restarts, budget = map(int, sys.argv[1:6])
    t0 = time.time()
    for n in range(n0, n1 + 1):
        if time.time() - t0 > budget: break
        bt, bE = -1, None
        for s in range(restarts):
            t, E = anneal(n, iters, s)
            if t > bt: bt, bE = t, E
            if 7 * bt > 3 * n: break
        flag = "**ABSTRACT COUNTEREXAMPLE**" if 7 * bt > 3 * n else ""
        print(f"n={n} best tau={bt} 3n/7={3*n/7:.2f} m={len(bE)} {flag} ({time.time()-t0:.0f}s)", flush=True)
        if flag: json.dump({"n": n, "tau": bt, "edges": bE}, open(f"abstract_cex_n{n}.json", "w"))
