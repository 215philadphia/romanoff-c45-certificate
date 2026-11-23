import itertools, json, time, sys
import numpy as np
from gluesearch import blocks, all_opt, d2_factored, try_extra
from abstract_falsify import tau, sparse, linear_with
from c5b_n_search import f7_copy
from ortools.sat.python import cp_model
def in_opt(n, E, t, v):
    m = cp_model.CpModel(); x = [m.NewBoolVar("") for _ in range(n)]
    for e in E: m.AddBoolOr([x[u] for u in e])
    m.Add(sum(x) == t); m.Add(x[v] == 1); s = cp_model.CpSolver(); s.parameters.num_workers = 1
    return s.Solve(m) in (cp_model.OPTIMAL, cp_model.FEASIBLE)
fam = [all_opt(n, E, t) for n, E, t in blocks]
flex = [bool(O.any(0).all()) for O in fam]
print("flexible:", [(b[0], f) for b, f in zip(blocks, flex)], flush=True)
t0 = time.time(); hits = []; stats = {"D2_skipped_by_lemma": 0, "D2_run": 0, "S1_bases": 0, "S1_edges_tested": 0, "S2_skipped": 0, "S2_run": 0}
START = tuple(map(int, sys.argv[1:3])) if len(sys.argv) > 2 else (0, 0)
for (iA, (nA, EA, tA)), (iB, (nB, EB, tB)) in itertools.combinations_with_replacement(list(enumerate(blocks)), 2):
    if (nA, nB) < START: continue
    OA, OB = fam[iA], fam[iB]
    # D2
    N = nA + nB
    if 7 * (tA + tB + 2) > 3 * N:
        if flex[iA] and flex[iB]: stats["D2_skipped_by_lemma"] += 1
        else:
            stats["D2_run"] += 1
            G = EA + [tuple(v + nA for v in e) for e in EB]
            h = d2_factored(nA, EA, tA, nB, EB, tB, G); hits += [("D2", nA, nB, x) for x in h]
    # S1 (budget 1)
    N = nA + nB - 1; budget = (N - 2) - (len(EA) + len(EB))
    if budget >= 0 and 7 * (tA + tB + budget) > 3 * N:
        for u in range(nA):
            for w in range(nB):
                mp = {v: (u if v == w else (nA + v if v < w else nA + v - 1)) for v in range(nB)}
                G = list(EA); ok = True
                for e in EB:
                    f = tuple(sorted(mp[v] for v in e))
                    if not linear_with(G, f): ok = False; break
                    G.append(f)
                if not ok: continue
                rowsA = OA[OA[:, u]]; rowsB = OB[OB[:, w]]
                if len(rowsA) and len(rowsB):          # tau(G0)=tA+tB-1, optimal = TA u TB both through shared vertex
                    tG = tA + tB - 1
                    covA = rowsA.any(0); covB = rowsB.any(0)
                    never = [v for v in range(nA) if not covA[v]] + [mp[v] for v in range(nB) if v != w and not covB[v]]
                else:
                    tG = tau(N, G); never = [v for v in range(N) if not in_opt(N, G, tG, v)]
                stats["S1_bases"] += 1
                if 7 * (tG + budget) <= 3 * N: continue
                if not (sparse(G) and f7_copy(G) is None): continue
                for e in itertools.combinations(never, 3):
                    if not linear_with(G, e): continue
                    GG = G + [e]; stats["S1_edges_tested"] += 1
                    if sparse(GG) and f7_copy(GG) is None:
                        t = tau(N, GG)
                        if 7 * t > 3 * N: hits.append(("S1", nA, nB, (u, w, e, t))); print("HIT S1", nA, nB, u, w, e, t, flush=True)
    # S2 (budget 0)
    N = nA + nB - 2; budget = (N - 2) - (len(EA) + len(EB))
    if budget >= 0 and 7 * (tA + tB + budget) > 3 * N:
        if flex[iA] and flex[iB] and 7 * (tA + tB - 1 + budget) <= 3 * N: stats["S2_skipped"] += 1
        else:
            stats["S2_run"] += 1
            covA = OA.any(0); covB = OB.any(0)
            for U in itertools.combinations(range(nA), 2):
                for W in itertools.permutations(range(nB), 2):
                    # lemma: if some shared pair (u,w) has u in a min transversal of A and w in one of B,
                    # gluing saves >=1, so tau <= tA+tB-1 (+0 budget); skip when that cannot beat 3N/7
                    if any(covA[u] and covB[w] for u, w in zip(U, W)) and 7 * (tA + tB - 1 + budget) <= 3 * N: continue
                    mp = {}; nxt = nA
                    for v in range(nB):
                        if v in W: mp[v] = U[W.index(v)]
                        else: mp[v] = nxt; nxt += 1
                    G = list(EA); ok = True
                    for e in EB:
                        f = tuple(sorted(mp[v] for v in e))
                        if not linear_with(G, f): ok = False; break
                        G.append(f)
                    if not ok: continue
                    t = tau(N, G)
                    if 7 * t > 3 * N and sparse(G) and f7_copy(G) is None:
                        hits.append(("S2", nA, nB, (U, W, t))); print("HIT S2", nA, nB, U, W, t, flush=True)
    print(f"pair ({nA},{nB}) done {time.time()-t0:.0f}s", flush=True)
print("STATS", stats); print("TOTAL HITS", len(hits))
json.dump({"stats": stats, "hits": [str(h) for h in hits]}, open("glue_fast_report.json", "w"))
