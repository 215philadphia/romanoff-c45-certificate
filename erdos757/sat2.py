import sys, itertools, time
from ortools.sat.python import cp_model
def sparse_violation(n, E):
    # find S with |E(S)| > |S|-2 : check unions of connected edge subsets via exhaustive over vertex subsets spanned (n<=16: brute force over subsets)
    Es=[set(e) for e in E]
    best=None
    for s in range(3,n+1):
        for S in itertools.combinations(range(n),s):
            Ss=set(S); c=sum(1 for e in Es if e<=Ss)
            if c>s-2: return S
    return None
def run(n, k, tlimit):
    m=cp_model.CpModel()
    T=list(itertools.combinations(range(n),3)); e={t:m.NewBoolVar("") for t in T}
    for t1,t2 in itertools.combinations(T,2):
        if len(set(t1)&set(t2))==2: m.AddBoolOr([e[t1].Not(),e[t2].Not()])
    m.Add(sum(e.values())<=n-2)
    deg=[sum(e[t] for t in T if v in t) for v in range(n)]
    for v in range(n-1): m.Add(deg[v]>=deg[v+1])          # symmetry breaking
    m.Add(deg[n-1]>=1)
    for K in itertools.combinations(range(n),k):
        rest=[v for v in range(n) if v not in K]
        m.AddBoolOr([e[t] for t in itertools.combinations(rest,3)])
    t0=time.time(); it=0
    while time.time()-t0<tlimit:
        it+=1
        s=cp_model.CpSolver(); s.parameters.max_time_in_seconds=tlimit-(time.time()-t0); s.parameters.num_workers=2
        r=s.Solve(m)
        if r==cp_model.INFEASIBLE: return "NO COUNTEREXAMPLE (proved)", it
        if r not in (cp_model.OPTIMAL,cp_model.FEASIBLE): return "TIMEOUT", it
        E=[t for t in T if s.Value(e[t])]
        S=sparse_violation(n,E)
        if S is None: return ("COUNTEREXAMPLE", E), it
        m.Add(sum(e[t] for t in itertools.combinations(S,3))<=len(S)-2)
    return "TIMEOUT", it
if __name__=="__main__":
    n=int(sys.argv[1]); k=(3*n)//7
    t=time.time(); r,it=run(n,k,float(sys.argv[2]))
    print(f"n={n}: sparse linear 3-graph with tau>={k+1} (>3n/7={3*n/7:.2f})? {r} | CEGAR iterations={it} ({time.time()-t:.0f}s)",flush=True)
