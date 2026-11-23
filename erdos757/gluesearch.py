"""Broader gluing search (Stage 5 falsification, abstract class).
Blocks: every known extremal (4,5)-set hypergraph (SK n=6..21 certificates; all tight m=n-2).
Constructions for each ordered block pair (A,B), including A=B:
  D2: disjoint union + 2 extra edges  (N=nA+nB, budget m<=N-2 -> exactly 2 extra)
  S1: identify one vertex of A with one of B + 1 extra edge  (N=nA+nB-1, budget 1)
  S2: identify two vertices (no extra edges; budget 0)
Extra edges must keep the hypergraph linear, F7-free and hereditarily sparse.
A construction is only explored if its best conceivable tau could exceed 3N/7.
tau(base + new edges) is decided exactly from the family of ALL optimal transversals of the base:
  1 edge raises tau  iff no optimal transversal meets it;
  2 edges raise tau  iff no optimal transversal meets both (they can raise it by at most 2;
  we compute the final tau with CP-SAT for every hit).
"""
import itertools, glob, re, json, time, sys
import numpy as np
from ortools.sat.python import cp_model
from tools import aps
from abstract_falsify import tau, sparse, linear_with
from c5b_n_search import f7_copy

def H_of(A):
    A = sorted(A); idx = {a: i for i, a in enumerate(A)}
    return [tuple(sorted(idx[v] for v in e)) for e in aps(A)]

def all_opt(n, E, t, cap=200000):
    m = cp_model.CpModel(); x = [m.NewBoolVar("") for _ in range(n)]
    for e in E: m.AddBoolOr([x[v] for v in e])
    m.Add(sum(x) == t)
    class CB(cp_model.CpSolverSolutionCallback):
        def __init__(s): super().__init__(); s.sols = []
        def on_solution_callback(s):
            s.sols.append([s.Value(v) for v in x])
            if len(s.sols) >= cap: s.StopSearch()
    s = cp_model.CpSolver(); s.parameters.enumerate_all_solutions = True; s.parameters.num_workers = 1
    cb = CB(); s.Solve(m, cb)
    return np.array(cb.sols, dtype=bool)

blocks = []
for f in sorted(glob.glob('/home/claude/ext5/sk/c5b-erdos757-small-n/certificates/n*.txt'), key=lambda f: int(re.findall(r'n(\d+)', f)[0])):
    A = [int(v) for v in open(f).read().split()]; E = H_of(A); n = len(A)
    blocks.append((n, E, tau(n, E)))
print("blocks:", [(n, t) for n, E, t in blocks], flush=True)

def try_extra(N, G, budget, label):
    t0 = tau(N, G)
    if budget == 0:
        return [(t0, G)] if 7 * t0 > 3 * N else []
    if 7 * (t0 + budget) <= 3 * N: return []
    O = all_opt(N, G, t0)                       # rows: optimal transversals
    if budget == 1:
        never = [v for v in range(N) if not O[:, v].any()]
        hits = []
        for e in itertools.combinations(never, 3):
            if not linear_with(G, e): continue
            GG = G + [e]
            if sparse(GG) and f7_copy(GG) is None:
                t = tau(N, GG)
                if 7 * t > 3 * N: hits.append((t, GG))
        return hits
    cands = []
    for e in itertools.combinations(range(N), 3):
        if not linear_with(G, e): continue
        hit = O[:, list(e)].any(1)
        cands.append((e, hit))
    hits = []
    if budget >= 1:
        for e, h in cands:
            if not h.any():
                GG = G + [e]
                if sparse(GG) and f7_copy(GG) is None:
                    t = tau(N, GG)
                    if 7 * t > 3 * N: hits.append((t, GG))
    if budget >= 2 and not hits:
        H = np.array([h for e, h in cands]); E2 = [e for e, h in cands]
        for i in range(len(E2)):
            both = (H[i] & H[i+1:]).any(1)
            for j in np.nonzero(~both)[0] + i + 1:
                e1, e2 = E2[i], E2[j]
                if len(set(e1) & set(e2)) >= 2: continue
                GG = G + [e1, e2]
                if not (sparse(GG) and f7_copy(GG) is None): continue
                t = tau(N, GG)
                if 7 * t > 3 * N: hits.append((t, GG))
    return hits

def d2_factored(nA, EA, tA, nB, EB, tB, G):
    N = nA + nB
    OA = all_opt(nA, EA, tA); OB = all_opt(nB, EB, tB)
    def hA(P): return OA[:, list(P)].any(1) if P else np.zeros(len(OA), bool)
    def hB(Q): return OB[:, [q - nA for q in Q]].any(1) if Q else np.zeros(len(OB), bool)
    cands = []
    for e in itertools.combinations(range(N), 3):
        P = [v for v in e if v < nA]; Q = [v for v in e if v >= nA]
        if not P or not Q: continue               # edges inside one block break that block's tightness
        if not linear_with(G, e): continue
        cands.append((e, hA(P), hB(Q)))
    # a pair can only survive if each member has an uncoverable side (lemma: else the partner's other side covers)
    cands = [c for c in cands if not c[1].any() or not c[2].any()]
    hits = []
    for (e1, a1, b1), (e2, a2, b2) in itertools.combinations(cands, 2):
        if (a1 & a2).any() or (a1.any() and b2.any()) or (b1.any() and a2.any()) or (b1 & b2).any(): continue
        GG = G + [e1, e2]
        if not (sparse(GG) and f7_copy(GG) is None): continue
        t = tau(N, GG)
        if 7 * t > 3 * N: hits.append((t, GG))
    return hits

t_start = time.time(); report = []
if __name__ == '__main__':
    START = tuple(map(int, sys.argv[1:3])) if len(sys.argv) > 2 else (0, 0)
    for (iA, (nA, EA, tA)), (iB, (nB, EB, tB)) in itertools.combinations_with_replacement(list(enumerate(blocks)), 2):
        if (nA, nB) < START: continue
        # D2
        N = nA + nB
        if 7 * (tA + tB + 2) > 3 * N:
            G = EA + [tuple(v + nA for v in e) for e in EB]
            h = d2_factored(nA, EA, tA, nB, EB, tB, G)
            report.append(("D2", nA, nB, len(h)));
            if h: print("HIT D2", nA, nB, h[0], flush=True)
        # S1 / S2
        for share in (1, 2):
            N = nA + nB - share; budget = (N - 2) - (len(EA) + len(EB))
            if budget < 0 or 7 * (tA + tB + budget) <= 3 * N: continue
            cnt = 0; nh = 0
            for U in itertools.combinations(range(nA), share):
                for W in itertools.permutations(range(nB), share):
                    mp = {}; nxt = nA
                    for v in range(nB):
                        if v in W: mp[v] = U[W.index(v)]
                        else: mp[v] = nxt; nxt += 1
                    G = list(EA)
                    ok = True
                    for e in EB:
                        f = tuple(sorted(mp[v] for v in e))
                        if not linear_with(G, f): ok = False; break
                        G.append(f)
                    if not ok or not sparse(G) or f7_copy(G) is not None: continue
                    cnt += 1
                    h = try_extra(N, G, budget, f"S{share}")
                    if h: nh += 1; print(f"HIT S{share}", nA, nB, U, W, h[0][0], flush=True)
            report.append((f"S{share}", nA, nB, cnt, nh))
        print(f"done pair ({nA},{nB})  {time.time()-t_start:.0f}s", flush=True)
    json.dump(report, open("gluesearch_report.json", "w"))
    print("TOTAL HITS:", sum(r[-1] for r in report))
    