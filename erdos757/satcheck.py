# Is there a linear 3-graph on n vertices with |E(S)|<=|S|-2 for every S (|S|>=3), and tau > k ?
import sys, itertools, time
from ortools.sat.python import cp_model
def check(n, k, extra_f7=False, timeout=600):
    m=cp_model.CpModel()
    T=list(itertools.combinations(range(n),3)); e={t:m.NewBoolVar(str(t)) for t in T}
    # linearity
    for t1,t2 in itertools.combinations(T,2):
        if len(set(t1)&set(t2))==2: m.AddBoolOr([e[t1].Not(),e[t2].Not()])
    # hereditary sparsity
    for s in range(3,n+1):
        for S in itertools.combinations(range(n),s):
            Sset=set(S)
            m.Add(sum(e[t] for t in itertools.combinations(S,3)) <= s-2)
    # tau > k: every k-set misses some edge
    for K in itertools.combinations(range(n),k):
        rest=[v for v in range(n) if v not in K]
        m.AddBoolOr([e[t] for t in itertools.combinations(rest,3)])
    # symmetry breaking: vertex 0 in an edge... (light) 
    s=cp_model.CpSolver(); s.parameters.max_time_in_seconds=timeout; s.parameters.num_workers=8
    r=s.Solve(m)
    if r in (cp_model.OPTIMAL,cp_model.FEASIBLE):
        return "FOUND",[t for t in T if s.Value(e[t])]
    return s.StatusName(r),None
if __name__=="__main__":
    for n in range(int(sys.argv[1]),int(sys.argv[2])+1):
        k=(3*n)//7   # tau>k  <=> tau > 3n/7 (since tau integer, tau>=k+1 > 3n/7 iff k+1>3n/7)
        if not (k+1>3*n/7): k+=1
        t=time.time(); st,E=check(n,k)
        print(f"n={n}: exists sparse linear 3-graph with tau>={k+1} (>3n/7={3*n/7:.2f})? {st} {E if E else ''} ({time.time()-t:.0f}s)",flush=True)
