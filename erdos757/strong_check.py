# Test the strengthened form  7*tau <= 2n + m + 2  on all linear, hereditarily sparse 3-graphs with n vertices
import sys, itertools, time
from ortools.sat.python import cp_model
from c5b_n_search import f7_copy
from abstract_falsify import sparse
def find(n, t, f7free, tl=280):
    # tau >= t and m <= 7t - 2n - 3  (violates 7tau <= 2n+m+2)
    mmax = 7*t - 2*n - 3
    if mmax < 1: return "trivially none"
    m = cp_model.CpModel(); T = list(itertools.combinations(range(n),3)); x = {e: m.NewBoolVar("") for e in T}
    for a,b in itertools.combinations(T,2):
        if len(set(a)&set(b))==2: m.AddBoolOr([x[a].Not(),x[b].Not()])
    m.Add(sum(x.values()) <= min(mmax, n-2))
    for K in itertools.combinations(range(n), t-1):
        rest=[v for v in range(n) if v not in K]; m.AddBoolOr([x[e] for e in itertools.combinations(rest,3)])
    for s in range(4, n+1):
        for S in itertools.combinations(range(n), s): m.Add(sum(x[e] for e in itertools.combinations(S,3)) <= s-2)
    t0=time.time()
    while time.time()-t0 < tl:
        sv=cp_model.CpSolver(); sv.parameters.max_time_in_seconds=tl; sv.parameters.num_workers=2
        r=sv.Solve(m)
        if r==cp_model.INFEASIBLE: return "NONE"
        if r not in (cp_model.OPTIMAL,cp_model.FEASIBLE): return "TIMEOUT"
        E=[e for e in T if sv.Value(x[e])]
        F=f7_copy(E) if f7free else None
        if F: m.Add(sum(x[e] for e in F)<=4); continue
        return ("FOUND", E)
    return "TIMEOUT"
for n in range(int(sys.argv[2]) if len(sys.argv)>2 else 3, int(sys.argv[1])+1):
    for t in range(1, n):
        if 7*t - 2*n - 3 < 1: continue
        if 7*t - 2*n - 3 > n - 2 and False: pass
        r = find(n, t, f7free=False)
        if r not in ("NONE","trivially none"): print(f"n={n} tau>={t}: {r}", flush=True)
    print(f"n={n} checked", flush=True)
