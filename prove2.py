"""Rigorous bound via mult.cpp (upward-rounded keys and fixed-point weights) + exact rational assembly.
usage: python3 prove2.py clusters.json G"""
import sys, json, math, subprocess
from fractions import Fraction
from math import lcm, prod

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

import os
SEEDP = [int(x) for x in os.environ.get("SEEDP","3,5,7,11,13,17,19,29,31,37,41,61,73,257").split(",")]
C = {int(k): sorted(v) for k, v in json.load(open(sys.argv[1])).items()}
G = int(sys.argv[2])
T0 = lcm(*[ordp(q) for q in SEEDP]); assert T0 == 5040
allp = SEEDP + sum(C.values(), [])
assert len(set(allp)) == len(allp) and all(isprime(q) for q in allp)
for l, qs in C.items():
    assert isprime(l) and math.gcd(l, T0) == 1 and all(ordp(q) == l for q in qs), l
lines = open(os.environ.get("SEEDFILE","seedY.txt")).read().split("\n")
assert int(lines[0].split()[1]) == T0
seed = [(int(a), int(b)) for a, b in (x.split() for x in lines[1:] if x.strip())]
M0 = prod(SEEDP); assert sum(w for _, w in seed) == M0
ONE = 1 << 62
def cluster(l, qs):
    st = {0: 1}
    for q in qs:
        n = {}
        for k, w in st.items():
            n[k] = n.get(k, 0) + w * (q - l) + w * k
            if k < l: n[k + 1] = n.get(k + 1, 0) + w * (l - k)
        st = n
    Q = prod(qs); assert sum(st.values()) == Q
    return [(l - k, -(-w * ONE // Q)) for k, w in sorted(st.items())]   # ceil
phiM = prod(Fraction(q - 1, q) for q in allp)
L2 = 2 * sum(Fraction(1, (2 * j + 1) * 3 ** (2 * j + 1)) for j in range(60))
c = 1 / (T0 * phiM * L2)
thr = {nu0: math.floor(Fraction(G) / (2 * nu0 * c)) for nu0, _ in seed if nu0 > 0}
inp = [f"{G} {len(C)}"]
for l, qs in sorted(C.items()):
    ch = cluster(l, qs)
    inp.append(f"{l} {len(ch)} " + " ".join(f"{v} {p}" for v, p in ch))
ts = sorted(set(thr.values()))
inp.append(f"{len(ts)} " + " ".join(map(str, ts)))
r = subprocess.run(["./mult"], input="\n".join(inp), capture_output=True, text=True, check=True)
out = r.stdout.split("\n")
n = int(out[0].split()[1]); tot = int(out[1].split()[1])
S = {int(a): (int(b), int(d)) for a, b, d in (x.split() for x in out[2:] if x.strip())}
s = Fraction(0)
for nu0, w0 in seed:
    if nu0 == 0: continue
    S0, S1 = S[thr[nu0]]
    s += Fraction(w0, M0) * (nu0 * c * Fraction(S1, G * ONE) + Fraction(tot - S0, 2 * ONE))
print(f"primes: {len(allp)} | clusters: {len(C)} | G = {G} | entries: {n} | weight overshoot: {float(Fraction(tot, ONE) - 1):.3e}")
print("rigorous upper bound:", f"{float(s):.15f}")
print("exact (first 40 digits):", str(s.numerator * 10**40 // s.denominator))
