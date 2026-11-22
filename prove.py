"""Exact rigorous upper bound with upward bucketing.
For prime set = SEED (orders | T0) + clusters (prime orders l, coprime to T0 and to each other):
  dbar(A) <= sum_classes min(1/(2M), nu/(T phi L2)),  L2 <= log 2 (rational lower bound),
  nu = nu0 * prod_l nu_l.  We replace prod_l (nu_l/l) by an upper bucket value K/G
  (rounded UP after every cluster multiply), which only increases each term because
  min(cap, c*nu) is nondecreasing in nu. All weights exact integers; final sum exact rational.
usage: python3 prove.py clusters.json G
"""
import sys, json, subprocess, math
from math import lcm, prod
from fractions import Fraction
import numpy as np

def ordp(q):
    x, d = 2 % q, 1
    while x != 1:
        x = x * 2 % q; d += 1
    return d
def isprime(n):
    if n < 2: return False
    i = 2
    while i * i <= n:
        if n % i == 0: return False
        i += 1
    return True

SEED = [3,5,7,11,13,17,19,29,31,37,41,61,73,257]
C = {int(k): sorted(v) for k, v in json.load(open(sys.argv[1])).items()}
G = int(sys.argv[2])
# ---- hypothesis checks ----
seed_ords = [ordp(q) for q in SEED]; T0 = lcm(*seed_ords)
allp = SEED + sum(C.values(), [])
assert len(set(allp)) == len(allp) and all(isprime(q) for q in allp)
for l, qs in C.items():
    assert isprime(l) and math.gcd(l, T0) == 1, l
    assert all(ordp(q) == l for q in qs), l
assert len(set(C)) == len(C)
# ---- seed histogram (independent C++ generator) ----
import os
if not os.path.exists("seedY.txt"):
    open("seedY.txt","w").write(subprocess.run(["./seed"] + list(map(str, SEED)), capture_output=True, text=True, check=True).stdout)
out = open("seedY.txt").read().split("\n")
assert int(out[0].split()[1]) == T0
seed = [(int(a), int(b)) for a, b in (l.split() for l in out[1:] if l.strip())]
assert sum(w for _, w in seed) == prod(SEED)
# ---- cluster multipliers with upward bucketing ----
def cluster(l, qs):
    st = {0: 1}
    for q in qs:
        n = {}
        for k, w in st.items():
            n[k] = n.get(k, 0) + w * (q - l) + w * k
            if k < l: n[k + 1] = n.get(k + 1, 0) + w * (l - k)
        st = n
    return [(l - k, w) for k, w in st.items()]
keys = np.array([G], dtype=np.int64); W = np.array([1], dtype=object)
for l, qs in sorted(C.items()):
    ch = cluster(l, qs)
    assert sum(w for _, w in ch) == prod(qs)
    nk = np.concatenate([(keys * v + (l - 1)) // l for v, _ in ch])        # exact ceil(K*v/l)
    nw = np.concatenate([W * w for _, w in ch])
    o = np.argsort(nk, kind="stable"); nk = nk[o]; nw = nw[o]
    starts = np.flatnonzero(np.r_[True, nk[1:] != nk[:-1]])
    keys = nk[starts]; W = np.add.reduceat(nw, starts)
    print(f"  cluster l={l} qs={qs}: buckets={len(keys)}", file=sys.stderr, flush=True)
Mcl = prod(sum(C.values(), []))
assert sum(W) == Mcl
M = prod(allp); phi = prod(q - 1 for q in allp)
L2 = 2 * sum(Fraction(1, (2 * j + 1) * 3 ** (2 * j + 1)) for j in range(60))
assert L2 < Fraction(693147180559945309418, 10**21)  # partial sums of a positive series < log 2
cap = Fraction(1, 2 * M)
D = G * T0 * phi * L2                    # per-class term = nu0*K / D
P0 = [0]; P1 = [0]
for k, w in zip(keys.tolist(), W.tolist()):
    P0.append(P0[-1] + w); P1.append(P1[-1] + w * k)
import bisect
kl = keys.tolist()
s = Fraction(0)
for nu0, w0 in seed:
    if nu0 == 0: continue
    thr = cap * D / nu0                  # K <= thr  -> uncapped term
    i = bisect.bisect_right(kl, math.floor(thr))
    s += w0 * (Fraction(nu0) * P1[i] / D + cap * (P0[-1] - P0[i]))
print("primes:", len(allp), "| clusters:", len(C), "| G =", G, "| buckets:", len(keys))
print("rigorous upper bound on dbar(A):", f"{float(s):.15f}")
print("exact (first 40 digits):", str(s.numerator * 10**40 // s.denominator))
print("below Yoo record 0.490180063290061:", s < Fraction(490180063290061, 10**15))
