"""Heuristic search for a counterexample to Conjecture S:
connected, linear, hereditarily sparse 3-graph with 7*tau > 2n + m + 2.
No F7-freeness imposed (S is claimed without it). Simulated annealing on edge sets.
Score = 7*tau - m (n fixed); target > 2n + 2.  Usage: python s_falsify.py n seconds seed
"""
import sys, random, itertools, time, json
sys.path.insert(0, '.')
from abstract_falsify import sparse, linear_with

def tau(n, E):
    best = [n]
    def rec(chosen, k):
        if k >= best[0]: return
        for e in E:
            if not (chosen & ((1 << e[0]) | (1 << e[1]) | (1 << e[2]))):
                for v in e: rec(chosen | (1 << v), k + 1)
                return
        best[0] = k
    rec(0, 0); return best[0]

def connected(n, E):
    if not E: return False
    adj = {v: set() for v in range(n)}
    for e in E:
        for v in e: adj[v] |= set(e)
    seen = {E[0][0]}; st = [E[0][0]]
    while st:
        u = st.pop()
        for w in adj[u]:
            if w not in seen: seen.add(w); st.append(w)
    return len(seen) == n

def run(n, secs, seed):
    rng = random.Random(seed); T = list(itertools.combinations(range(n), 3))
    best = -99; bestE = None; t0 = time.time()
    while time.time() - t0 < secs:
        E = []
        for _ in range(3000):
            t = rng.choice(T)
            if t not in E and linear_with(E, t) and sparse(E + [t]): E.append(t)
            if len(E) >= n - 2: break
        cur = 7 * tau(n, E) - len(E) if connected(n, E) else -99
        for it in range(400):
            F = list(E); r = rng.random()
            if r < 0.5 and F: F.pop(rng.randrange(len(F)))
            if r > 0.3:
                for _ in range(50):
                    t = rng.choice(T)
                    if t not in F and linear_with(F, t) and sparse(F + [t]): F.append(t); break
            if not connected(n, F): continue
            s = 7 * tau(n, F) - len(F)
            if s >= cur or rng.random() < 0.05:
                E, cur = F, s
                if s > best:
                    best, bestE = s, list(F)
                    if best > 2 * n + 2:
                        print("COUNTEREXAMPLE", n, bestE, flush=True)
                        json.dump(bestE, open(f"S_cex_{n}_{seed}.json", "w")); return
    print(f"n={n} best 7tau-m={best} target>{2*n+2} slack={2*n+2-best}", flush=True)

if __name__ == "__main__":
    run(int(sys.argv[1]), float(sys.argv[2]), int(sys.argv[3]))
