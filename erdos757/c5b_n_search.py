#!/usr/bin/env python3
"""
Erdős #757 / Tao constant C_5b: search for a counterexample to the transversal conjecture

  Conjecture. Every linear, F7-free 3-uniform hypergraph H on n vertices with
  |E(S)| <= |S| - 2 for every vertex set S (|S| >= 2) has tau(H) <= 3n/7.

(Ma-Tang 2026 show the 3-AP hypergraph of every (4,5)-set has these three properties,
so the conjecture would give C_5b = 4/7.)  By Henning-Yeo (17 tau <= 5n + 3m for linear
F7-free 3-graphs) the conjecture holds automatically for n <= 15; n = 16 is the first
open case.

The script asks CP-SAT for a linear 3-graph with tau >= floor(3n/7) + 1 and m <= n - 2.
Sparsity on subsets and F7-freeness are added lazily: each time the solver returns a
candidate, the script finds a violated subset / an F7 copy, adds a cut, and re-solves.
  * INFEASIBLE  -> proof that no counterexample exists for this n (given the solver is correct)
  * COUNTEREXAMPLE -> printed and saved; it satisfies every condition (re-checked by brute force)

usage:  python c5b_n_search.py 16 --workers 8 --hours 24
requires: pip install ortools
"""
import argparse, itertools, json, sys, time
from ortools.sat.python import cp_model


def sparse_violation(n, E):
    Es = [frozenset(e) for e in E]
    for s in range(3, n + 1):
        for S in itertools.combinations(range(n), s):
            Ss = set(S)
            if sum(1 for e in Es if e <= Ss) > s - 2:
                return S
    return None


def f7_copy(E):
    """Return 5 edges forming F7 = {c a1 a2},{c b1 b2},{c d1 d2},{a? b? d?},{a? b? d?} (the AG(2,3)-minus-2-points
    configuration), or None."""
    Es = [tuple(e) for e in E]
    S = {frozenset(e) for e in Es}
    by_v = {}
    for e in Es:
        for v in e:
            by_v.setdefault(v, []).append(e)
    for c, inc in by_v.items():
        for e1, e2, e3 in itertools.combinations(inc, 3):
            P = [tuple(x for x in e if x != c) for e in (e1, e2, e3)]
            if len({*P[0], *P[1], *P[2]}) != 6:
                continue
            for i, j, k in itertools.product(range(2), repeat=3):
                f = frozenset((P[0][i], P[1][j], P[2][k]))
                g = frozenset((P[0][1 - i], P[1][1 - j], P[2][1 - k]))
                if f in S and g in S:
                    return [e1, e2, e3, tuple(f), tuple(g)]
    return None


def tau_brute(n, E):
    for k in range(n + 1):
        for T in itertools.combinations(range(n), k):
            Ts = set(T)
            if all(Ts & set(e) for e in E):
                return k


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("n", type=int)
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--hours", type=float, default=24.0)
    ap.add_argument("--no-f7", action="store_true", help="drop the F7-free condition (calibration only)")
    a = ap.parse_args()
    n = a.n
    k = (3 * n) // 7            # need tau >= k+1  (> 3n/7)
    deadline = time.time() + 3600 * a.hours

    m = cp_model.CpModel()
    T = list(itertools.combinations(range(n), 3))
    x = {t: m.NewBoolVar(str(t)) for t in T}
    for t1, t2 in itertools.combinations(T, 2):          # linear
        if len(set(t1) & set(t2)) == 2:
            m.AddBoolOr([x[t1].Not(), x[t2].Not()])
    m.Add(sum(x.values()) <= n - 2)                       # global sparsity
    deg = [sum(x[t] for t in T if v in t) for v in range(n)]
    for v in range(n - 1):
        m.Add(deg[v] >= deg[v + 1])                       # symmetry breaking (relabel by degree)
    m.Add(deg[n - 1] >= 1)                                # isolated vertex => smaller counterexample
    for K in itertools.combinations(range(n), k):         # tau >= k+1
        rest = [v for v in range(n) if v not in K]
        m.AddBoolOr([x[t] for t in itertools.combinations(rest, 3)])
    # cheap up-front sparsity cuts on all 4..7-vertex sets (the rest are added lazily)
    for s in range(4, min(n, 7) + 1):
        for S in itertools.combinations(range(n), s):
            m.Add(sum(x[t] for t in itertools.combinations(S, 3)) <= s - 2)

    it, cuts = 0, 0
    log = open(f"c5b_n{n}_log.txt", "a")
    def say(s):
        print(s, flush=True); log.write(s + "\n"); log.flush()
    say(f"n={n}: searching for linear{'' if a.no_f7 else ' F7-free'} sparse 3-graph with tau >= {k+1} (> 3n/7 = {3*n/7:.3f})")
    while True:
        left = deadline - time.time()
        if left <= 0:
            say("STOPPED: time limit reached (inconclusive)"); return
        s = cp_model.CpSolver()
        s.parameters.num_workers = a.workers
        s.parameters.max_time_in_seconds = left
        r = s.Solve(m)
        it += 1
        if r == cp_model.INFEASIBLE:
            say(f"RESULT n={n}: INFEASIBLE -> no counterexample exists (iterations={it}, cuts={cuts})"); return
        if r not in (cp_model.OPTIMAL, cp_model.FEASIBLE):
            say(f"STOPPED: solver status {s.StatusName(r)} (inconclusive)"); return
        E = [t for t in T if s.Value(x[t])]
        S = sparse_violation(n, E)
        if S is not None:
            m.Add(sum(x[t] for t in itertools.combinations(S, 3)) <= len(S) - 2); cuts += 1
            if cuts % 50 == 0: say(f"  {cuts} cuts so far")
            continue
        if not a.no_f7:
            F = f7_copy(E)
            if F is not None:
                m.Add(sum(x[t] for t in F) <= 4); cuts += 1
                if cuts % 50 == 0: say(f"  {cuts} cuts so far")
                continue
        t = tau_brute(n, E)
        assert t >= k + 1, "solver/encoding mismatch"
        say(f"RESULT n={n}: COUNTEREXAMPLE FOUND  tau={t}  m={len(E)}  edges={E}")
        json.dump({"n": n, "tau": t, "edges": E}, open(f"c5b_counterexample_n{n}.json", "w"))
        return


if __name__ == "__main__":
    main()
