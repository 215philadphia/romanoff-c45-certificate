"""Independent exact verifier/optimizer for Romanoff-constant upper bounds.

Bound (obstruction argument): for a finite odd prime set P,
  M = prod q, T = lcm ord_q(2), nu(a) = #{r mod T: gcd(a-2^r, M)=1},
  dbar(A) <= sum_a min(1/(2M), nu(a)/(T phi(M) log 2)).
Structure used: by CRT, a <-> (a mod q)_q. For each q, a mod q is either one of the
ord_q(2) powers of 2 (removes the class r = e mod ord_q) or one of the q - ord_q other
residues (removes nothing). Primes are split into a 'seed' whose orders divide T0, and
clusters of primes sharing an order l with gcd(l, T0) = 1 and the l's pairwise coprime.
Then nu = nu_seed * prod_l nu_l exactly (CRT on r).
All arithmetic exact (integers / Fractions); log 2 bounded below by a truncated series.
"""
import sys, json
from fractions import Fraction
from math import gcd, lcm, prod
from collections import defaultdict
import numpy as np


def ordp(q):
    x, d = 2 % q, 1
    while x != 1:
        x = x * 2 % q; d += 1
    return d


def seed_histogram(primes):
    """Exact histogram {nu: #a mod prod(primes)} for r mod T0=lcm orders.
    State: vector c[x] = number of surviving r mod T0 with r = x mod L, where L is the
    lcm of orders of primes not yet processed (projection). States are merged up to
    cyclic rotation of x (choices are translation-equivariant)."""
    ords = {q: ordp(q) for q in primes}
    T0 = lcm(*ords.values())
    # process primes with largest orders first (different order from Griego's)
    order = sorted(primes, key=lambda q: (-ords[q], q))
    Ls = [lcm(*[ords[q] for q in order[i:]]) if order[i:] else 1 for i in range(len(order) + 1)]
    def canon(v):
        L = len(v)
        if L == 1: return (int(v[0]),)
        best = None
        for s in range(L):
            w = tuple(int(t) for t in np.roll(v, -s))
            if best is None or w < best: best = w
        return best
    states = {canon(np.full(Ls[0], T0 // Ls[0], dtype=np.int64)): 1}
    for i, q in enumerate(order):
        d, L, nL = ords[q], Ls[i], Ls[i + 1]
        nxt = defaultdict(int)
        for key, w in states.items():
            v = np.array(key, dtype=np.int64)
            # non-power residues: no removal, just project
            out = v.reshape(L // nL, nL).sum(axis=0) if nL < L else v
            nxt[canon(out)] += w * (q - d)
            for s in range(d):
                u = v.copy(); u[s::d] = 0
                out = u.reshape(L // nL, nL).sum(axis=0) if nL < L else u
                nxt[canon(out)] += w
        states = nxt
        print(f"  seed step q={q} ord={d} L {L}->{nL} states={len(states)}", file=sys.stderr)
    hist = defaultdict(int)
    for key, w in states.items():
        assert len(key) == 1
        hist[key[0]] += w
    assert sum(hist.values()) == prod(primes)
    return dict(hist), T0


def cluster_hist(l, qs):
    """primes qs all of order l: nu_l = l - #distinct chosen exponents. Exact histogram."""
    # DP over primes: state = number of distinct exponents chosen so far (by symmetry)
    st = {0: 1}
    for q in qs:
        n = defaultdict(int)
        for k, w in st.items():
            n[k] += w * (q - l)          # non-power residue
            n[k] += w * k                # power residue hitting an already-chosen exponent
            if k < l: n[k + 1] += w * (l - k)  # new exponent
        st = n
    return {l - k: w for k, w in st.items()}


def log2_lower(terms=60):
    return 2 * sum(Fraction(1, (2 * j + 1) * 3 ** (2 * j + 1)) for j in range(terms))


def bound(primes, seed_cache={}):
    ords = {q: ordp(q) for q in primes}
    # seed = primes whose order shares a factor with some other non-prime-cluster structure:
    # define clusters = primes with prime order l such that l is coprime to all other orders
    by = defaultdict(list)
    for q in primes: by[ords[q]].append(q)
    cl = {}
    for o, qs in by.items():
        others = [ords[p] for p in primes if ords[p] != o]
        if o > 2 and all(o % p for p in range(2, int(o ** .5) + 1)) and all(x % o for x in others):
            cl[o] = sorted(qs)
    seed = sorted(q for q in primes if ords[q] not in cl)
    key = tuple(seed)
    if key not in seed_cache:
        import subprocess
        out = subprocess.run(["./seed"] + [str(q) for q in seed], capture_output=True, text=True, check=True).stdout.split("\n")
        T0 = int(out[0].split()[1])
        h = {int(a): int(b) for a, b in (l.split() for l in out[1:] if l.strip())}
        assert sum(h.values()) == prod(seed)
        seed_cache[key] = (h, T0)
    sh, T0 = seed_cache[key]
    T = T0
    mult = {1: 1}
    for l, qs in sorted(cl.items()):
        ch = cluster_hist(l, qs)
        new = defaultdict(int)
        for a, wa in mult.items():
            for b, wb in ch.items():
                new[a * b] += wa * wb
        mult = new; T *= l
    M = prod(primes)
    phi = prod(q - 1 for q in primes)
    assert T == lcm(*ords.values()) and sum(mult.values()) * sum(sh.values()) == M
    import bisect
    ms = sorted(mult); W = [mult[m] for m in ms]
    P0 = [0]; P1 = [0]
    for m, w in zip(ms, W):
        P0.append(P0[-1] + w); P1.append(P1[-1] + w * m)
    L2 = log2_lower()
    cap = Fraction(1, 2 * M)
    D = T * phi * L2          # nu/D is the per-class prime-count bound
    s = Fraction(0)
    for nu0, w0 in sh.items():
        if nu0 == 0: continue
        # nu0*m/D <= cap  <=>  m <= cap*D/nu0
        thr = cap * D / nu0
        k = bisect.bisect_right(ms, thr)   # ms[:k] use the nu bound (exact compare via Fraction)
        while k > 0 and Fraction(ms[k-1]) > thr: k -= 1
        while k < len(ms) and Fraction(ms[k]) <= thr: k += 1
        s += w0 * (Fraction(nu0) * P1[k] / D + cap * (P0[-1] - P0[k]))
    hist = mult
    return s, dict(M=M, T=T, phi=phi, entries=len(hist), seed_entries=len(sh), seed=seed, clusters=cl)


if __name__ == "__main__":
    P = json.loads(sys.argv[1])
    s, info = bound(P)
    print("primes", len(P), "| T =", info["T"], "| hist entries", info["entries"])
    print("M =", info["M"])
    print("upper bound sum =", f"{float(s):.15f}", " exact-decimal:", str(s.numerator * 10**30 // s.denominator))
