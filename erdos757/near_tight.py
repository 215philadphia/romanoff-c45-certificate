"""Catalog connected, linear, hereditarily sparse 3-graphs with sigma = 2n+m+2-7tau <= 1
(the components that block the vertex-deletion reduction). Heuristic (annealing), deduped up to isomorphism.
Usage: python near_tight.py nmin nmax seconds_per_n
"""
import sys, random, itertools, time, json, networkx as nx
sys.path.insert(0, '.')
from abstract_falsify import sparse, linear_with
from s_falsify import tau, connected

def inc(n, E):
    G = nx.Graph()
    for v in range(n): G.add_node(('v', v), b=0)
    for i, e in enumerate(E):
        G.add_node(('e', i), b=1)
        for v in e: G.add_edge(('v', v), ('e', i))
    return G

def search(n, secs, seed, found):
    rng = random.Random(seed); T = list(itertools.combinations(range(n), 3)); t0 = time.time()
    nm = nx.algorithms.isomorphism.categorical_node_match('b', 0)
    while time.time() - t0 < secs:
        E = []
        for _ in range(3000):
            t = rng.choice(T)
            if t not in E and linear_with(E, t) and sparse(E + [t]): E.append(t)
            if len(E) >= n - 2: break
        cur = 7*tau(n, E) - len(E) if connected(n, E) else -99
        for it in range(300):
            F = list(E); r = rng.random()
            if r < 0.5 and F: F.pop(rng.randrange(len(F)))
            if r > 0.3:
                for _ in range(50):
                    t = rng.choice(T)
                    if t not in F and linear_with(F, t) and sparse(F + [t]): F.append(t); break
            if not connected(n, F): continue
            s = 7*tau(n, F) - len(F)
            if s >= cur or rng.random() < 0.05:
                E, cur = F, s
                sig = 2*n + 2 - s
                if sig <= 1:
                    G = inc(n, F)
                    if not any(nx.is_isomorphic(G, H, node_match=nm) for (_, _, H) in found.get(n, [])):
                        found.setdefault(n, []).append((sig, sorted(F), G))
                        print(f"n={n} m={len(F)} tau={(s+len(F))//7} sigma={sig} E={sorted(F)}", flush=True)

if __name__ == "__main__":
    found = {}
    for n in range(int(sys.argv[1]), int(sys.argv[2]) + 1):
        search(n, float(sys.argv[3]), n, found)
        print(f"n={n} done: {len(found.get(n, []))} near-tight classes found", flush=True)
    json.dump({n: [(s, E) for (s, E, _) in v] for n, v in found.items()}, open(f"near_tight_out_{sys.argv[1]}.json", "w"))
