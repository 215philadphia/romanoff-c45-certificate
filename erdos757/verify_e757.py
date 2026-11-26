#!/usr/bin/env python3
"""
Verification script for: J. Clark, "A short proof that f(18)=11, f(20)=12 and f(22)=13 in
Erdos Problem #757", Zenodo, version 3 of doi:10.5281/zenodo.24547831.
These values were first established by other means in hturcotte1/Erdos-757 (GitHub, 2026-10-17)
and fortunto2/erdos-757-sidon45 (doi:10.5281/zenodo.23494596); this script checks the short proof.

f(n) = min over (4,5)-sets A with |A| = n of h(A), the size of the largest Sidon subset of A.

In a (4,5)-set every 4 elements have 5 distinct differences, so a Sidon violation a+d=b+c with
4 distinct elements is impossible; the only obstructions are 3-term APs (a+c = 2b). Hence
h(A) = n - tau(H_A), where H_A is the 3-uniform hypergraph of 3-APs in A.

Part 1 (upper bounds, exact computation):
  - Ma-Tang's 14-point set A_base (arXiv:2602.23282, Sec. 5.1) is a (4,5)-set with h = 8.
  - an explicit 6-point (4,5)-set with h = 4 and 8-point (4,5)-set with h = 5.
  - Ma-Tang Prop. 3.7 (subadditivity f(m+n) <= f(m)+f(n)) then gives
    f(18) <= f(14)+f(4) <= 11, f(20) <= f(14)+f(6) <= 12 and f(22) <= f(14)+f(8) <= 13.
  - As an extra check, the gluing of Ma-Tang Lemma 3.6 is carried out explicitly and the
    resulting 20- and 22-point sets are verified directly.

Part 2 (lower bounds, arithmetic of the bound):
  Relies on (not re-proved here):
   (a) H_A is linear and F7-free (as used by Ma-Tang to apply Henning-Yeo);
   (b) every connected component C of H_A has m_C <= n_C - 2 edges;
   (c) Henning-Yeo 2021, Theorem 6: 17 tau <= 4n + 4m + h for such H, h = # single-edge components,
       and 17 tau < 4n + 4m if H has no exceptional vertex (none of its components is E3, H10, H10').
  Per component this gives tau_C <= floor((8 n_C - 9)/17), except tau_C <= 1 for E3 and <= 4 for H10, H10'.
  The script maximises the total over all ways to split n vertices into components and
  compares with the old bound floor((8n-6)/17) and with floor(3n/7).
"""
from itertools import combinations
from fractions import Fraction
import sys

def is45(A):
    return all(len({abs(x - y) for x, y in combinations(B, 2)}) >= 5
               for B in combinations(A, 4))

def is_sidon(S):
    sums = set()
    for i, a in enumerate(S):
        for b in S[i:]:
            s = a + b
            if s in sums:
                return False
            sums.add(s)
    return True

def h_exact(A):
    """Largest Sidon subset, by exhaustive search from the top size down."""
    A = sorted(A)
    for k in range(len(A), 0, -1):
        for S in combinations(A, k):
            if is_sidon(S):
                return k, S
    return 0, ()

def three_aps(A):
    S = set(A)
    return [(a, (a + c) // 2, c) for a, c in combinations(sorted(A), 2)
            if (a + c) % 2 == 0 and (a + c) // 2 in S]

def tau_exact(A):
    idx = {a: i for i, a in enumerate(sorted(A))}
    edges = [frozenset(idx[v] for v in e) for e in three_aps(A)]
    best = [len(A)]
    def rec(chosen, k):
        if k >= best[0]:
            return
        for e in edges:
            if not (e & chosen):
                for v in e:
                    rec(chosen | {v}, k + 1)
                return
        best[0] = k
    rec(frozenset(), 0)
    return best[0]

def glue(A, B):
    """Ma-Tang Lemma 3.6: C = A u (qB + t), with q not in R and t - D_A > max(D_A, q D_B).
    Uses integer q (scaled) so C stays in Z."""
    A = [a - min(A) for a in A]; B = [b - min(B) for b in B]
    DA, DB = max(A), max(B)
    R = {Fraction(a - a2, b - b2) for a in A for a2 in A for b in B for b2 in B if b != b2}
    q = 1
    while Fraction(q) in R:
        q += 1
    t = max(DA, q * DB) + DA + 1
    return sorted(A + [q * b + t for b in B])

def main():
    ok = True
    def check(cond, msg):
        nonlocal ok
        print(("PASS " if cond else "FAIL ") + msg)
        ok &= bool(cond)

    print("== Part 1: explicit sets ==")
    A14 = [0, 136, 200, 243, 246, 249, 272, 286, 298, 323, 400, 528, 596, 1056]  # Ma-Tang Sec 5.1
    A8 = [5, 11, 16, 21, 24, 37, 41, 58]
    # 6-point set: found by search below if not given
    A6 = None
    for S in combinations(range(0, 40), 6):
        if S[0] != 0:
            break
        if is45(S) and h_exact(S)[0] == 4:
            A6 = list(S); break
    check(A6 is not None, f"found a 6-point (4,5)-set with h=4: {A6}")
    for name, A, target in [("A6", A6, 4), ("A8", A8, 5), ("A14 (Ma-Tang)", A14, 8)]:
        k, S = h_exact(A)
        check(is45(A), f"{name} is a (4,5)-set (n={len(A)})")
        check(k == target, f"{name}: h = {k} (expected {target}); a largest Sidon subset: {list(S)}")
        check(len(A) - tau_exact(A) == k, f"{name}: h = n - tau(3-AP hypergraph) = {len(A)} - {tau_exact(A)}")

    print("\n== Part 1b: explicit glued sets (Ma-Tang Lemma 3.6) ==")
    for name, X, Y, target in [("C20 = A14 + A6", A14, A6, 12), ("C22 = A14 + A8", A14, A8, 13)]:
        C = glue(X, Y)
        check(is45(C), f"{name} is a (4,5)-set (n={len(C)})")
        hC = len(C) - tau_exact(C)
        check(hC <= target, f"{name}: h = {hC} <= {target}")
        print(f"     {name} = {C}")

    print("\n== Part 1c: 18-point set (v3) ==")
    A4 = [1, 2, 3, 6]                    # Ma-Tang: a (4,5)-set that is not Sidon
    check(is45(A4) and h_exact(A4)[0] == 3, "A4 = {1,2,3,6} is a (4,5)-set with h = 3")
    C18 = glue(A14, A4)
    check(is45(C18), f"C18 = A14 + A4 is a (4,5)-set (n={len(C18)})")
    check(len(C18) - tau_exact(C18) == 11, f"C18: h = {len(C18) - tau_exact(C18)} = 11")
    print(f"     C18 = {C18}")

    print("\n== Part 2: lower-bound arithmetic (Henning-Yeo Thm 6, with its strict clause) ==")
    # Every component C of H(A) is linear, F7-free and has m_C <= n_C - 2 (its edges are the
    # 3-APs of the (4,5)-set V(C); Ma-Tang Lemmas 2.4, 2.5, 5.2). Thm 6 applied to C:
    #   E3: tau = 1.  H10, H10' (n = 10): 17 tau <= 4n + 4m = 68, tau <= 4.
    #   any other component has no exceptional vertex: 17 tau < 4n + 4m <= 8n - 8.
    def comp_bound(n, strict):
        if n == 1: return 0
        if n == 2: return None           # no 3-uniform component on 2 vertices
        if n == 3: return 1              # E3
        if not strict: return (8 * n - 8) // 17
        b = (8 * n - 9) // 17
        return max(b, 4) if n == 10 else b
    def dp(strict, N):
        best = [0] + [-1] * N
        for n in range(1, N + 1):
            for c in range(1, n + 1):
                b = comp_bound(c, strict)
                if b is None or best[n - c] < 0: continue
                best[n] = max(best[n], best[n - c] + b)
        return best
    N = 60
    weak, best = dp(False, N), dp(True, N)
    print("  n  tau<= (Thm6)  tau<= (Thm6 strict)  Ma-Tang floor((8n-6)/17)  => f(n) >=")
    for n in range(4, N + 1):
        old = (8 * n - 6) // 17
        if n in (18, 20, 22) or best[n] < old:
            print(f" {n:3d}      {weak[n]:3d}            {best[n]:3d}                  {old:3d}              {n - best[n]}")
    check(best[18] == 7, "f(18) >= 11")
    check(best[20] == 8, "f(20) >= 12")
    check(best[22] == 9, "f(22) >= 13")

    print("\nRESULT:", "ALL CHECKS PASSED" if ok else "SOME CHECKS FAILED")
    sys.exit(0 if ok else 1)

if __name__ == "__main__":
    main()
