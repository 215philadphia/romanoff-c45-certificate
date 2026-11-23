"""Executed checks for the new lemmas.
(1) Lemma A: connected linear 3-graphs with max degree <= 2 have sigma = 2n+m+2-7tau >= 2.
    Checked on random instances, and the dual identity tau = m - nu(G) checked directly.
(2) Chvatal-Hanson for Delta=3: e(G) <= 3nu + floor(nu/2), checked on random subcubic graphs.
"""
import random, itertools, networkx as nx, sys
sys.path.insert(0, '.')
from s_falsify import tau, connected

def rand_deg2(n, rng):
    deg = [0]*n; E = []
    for _ in range(2000):
        t = tuple(sorted(rng.sample(range(n), 3)))
        if t in E or any(deg[v] >= 2 for v in t): continue
        if any(len(set(t) & set(f)) > 1 for f in E): continue
        E.append(t)
        for v in t: deg[v] += 1
    return E

def nu_dual(E):
    G = nx.Graph(); G.add_nodes_from(range(len(E)))
    for i, j in itertools.combinations(range(len(E)), 2):
        if set(E[i]) & set(E[j]): G.add_edge(i, j)
    return G, len(nx.max_weight_matching(G, maxcardinality=True))

rng = random.Random(1); checked = 0; minsig = 99
for trial in range(4000):
    n = rng.randint(3, 16); E = rand_deg2(n, rng)
    used = sorted({v for e in E for v in e}); idx = {v: i for i, v in enumerate(used)}
    E = [tuple(idx[v] for v in e) for e in E]; n = len(used)
    if not E or not connected(n, E): continue
    t = tau(n, E); G, nu = nu_dual(E)
    assert t == len(E) - nu, (E, t, nu)
    sig = 2*n + len(E) + 2 - 7*t; minsig = min(minsig, sig); assert sig >= 2, E
    checked += 1
print(f"Lemma A: {checked} connected Delta<=2 instances, tau = m - nu(dual) on all, min sigma = {minsig}")

bad = 0; cnt = 0
for trial in range(3000):
    k = rng.randint(2, 14); G = nx.gnm_random_graph(k, rng.randint(1, 2*k), seed=rng.randrange(10**9))
    G.remove_edges_from([e for e in list(G.edges()) if G.degree(e[0]) > 3 or G.degree(e[1]) > 3])
    if max((d for _, d in G.degree()), default=0) > 3 or G.number_of_edges() == 0: continue
    nu = len(nx.max_weight_matching(G, maxcardinality=True)); cnt += 1
    if G.number_of_edges() > 3*nu + nu//2: bad += 1
print(f"Chvatal-Hanson (Delta<=3): {cnt} graphs, violations = {bad}")
