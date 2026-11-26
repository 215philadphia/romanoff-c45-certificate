# Proof progress on Conjecture S (session of 2026-11-22)

> **REFUTED (2026-11-26).** Conjecture S (7τ ≤ 2n + m + 2) and the target τ ≤ 3n/7 are false: the 65-point (4,5)-set C003 has m = 63 and τ = 28 (7τ = 196 > 195), and it gives c* ≤ 37/65 < 4/7, so C_5b = 4/7 is false too. Kept for the record; see RESULTS.md.

Notation: H is 3-uniform, linear, and hereditarily sparse (|E(S)| <= |S| - 2 for every S spanning an edge).
For connected H, write the slack as sigma(H) = 2n + m + 2 - 7 tau(H). Conjecture S says sigma(H) >= 0.
Labels: PROVEN (written proof below), EXECUTED (computer check), OBSERVED (heuristic search, not exhaustive).

## Lemma A (PROVEN). If H is connected and linear with max degree <= 2, then sigma(H) >= 2.
Sparsity is not needed.

*Proof.* Let G be the dual graph. G has one node for each edge of H, and one edge for each degree-2 vertex of H, joining the two H-edges that contain it.
- G is simple, because H is linear.
- G has max degree <= 3, because an H-edge has 3 vertices.
- Write n1 and n2 for the numbers of degree-1 and degree-2 vertices of H. Then n = n1 + n2, 3m = n1 + 2 n2 and |E(G)| = n2, so n = 3m - n2.

Claim: tau(H) = m - nu(G).
- Upper bound. Take a maximum matching M of G. Its |M| vertices cover 2|M| edges of H. Add one vertex from each remaining H-edge. The total is m - |M|.
- Lower bound. A minimal transversal is a cover of the nodes of G by G-edges and single nodes, where a degree-1 vertex covers one node. The G-edges of a minimal cover form vertex-disjoint stars. A star with s edges covers s + 1 nodes, and a single vertex covers 1 node. So m - t equals the number of stars, which is at most nu(G).

Therefore sigma = 2(3m - n2) + m + 2 - 7(m - nu) = 7 nu + 2 - 2|E(G)|.

By the Chvátal–Hanson theorem (1976), a graph with max degree 3 and matching number nu has at most 3 nu + floor(nu/2) edges. So 2|E(G)| <= 7 nu, and sigma >= 2. ∎

EXECUTED (`lemma_checks.py`):
- The identity tau = m - nu(G) and sigma >= 2 hold on 3,947 random connected instances.
- Chvátal–Hanson holds on 2,774 random graphs with max degree <= 3.

## Lemma D (PROVEN). Vertex-deletion inequality for a minimal counterexample.
Let H be a connected counterexample to S with n + m minimal, so 7 tau(H) >= 2n + m + 3. Let y be a vertex of degree d. In H - y (delete y and its d edges), let C_1, ..., C_r be the components that have edges, and let I be the set of isolated vertices. Then

    d + 2|I| + sum_i sigma(C_i) <= 2 + 2r.

*Proof.* We have tau(H) <= 1 + sum_i tau(C_i). Each C_i is connected, linear and hereditarily sparse, with a smaller n + m than H, so sigma(C_i) >= 0. Also sum_i n_i = n - 1 - |I| and sum_i m_i = m - d. This gives

    7 tau(H) <= 7 + 2(n - 1 - |I|) + (m - d) + 2r - sum_i sigma(C_i).

Comparing with 7 tau(H) >= 2n + m + 3 gives the inequality. ∎

## Corollary (PROVEN). Structure of a minimal counterexample H.
1. H has max degree >= 3 (Lemma A).
2. Every vertex y of degree >= 3 satisfies sum_i (2 - sigma(C_i)) >= d(y) + 2|I_y| - 2 >= 1. So H - y has a component with sigma <= 1, and by Lemma A that component also has max degree >= 3.
3. If H - y is connected (r = 1), then d(y) + sigma(C_1) + 2|I_y| <= 4. So d(y) is 3 or 4, and C_1 is near-tight.

Together with the earlier bound n >= 11 (exhaustive check for n <= 10), this replaces the naive "delete a vertex of degree >= 5" reduction. Every high-degree vertex of a minimal counterexample must separate off near-tight pieces.

## Computational evidence from this session
- **No counterexample to S for n = 11..18 (OBSERVED, `s_falsify.py`, 15 min annealing per n, no F7-freeness imposed).** The best slack found was 3, 1, 1, 0, 1, 2, 2, 3 for n = 11..18. The slack-0 case at n = 14 matches H14.
- **Near-tight pieces are common (OBSERVED, `near_tight.py`, `near_tight_catalog.json`).** Connected instances with sigma = 1, up to isomorphism, found by heuristic search:

| n | sigma = 1 classes found |
|---|---|
| 10 | 6 |
| 12 | 134 |
| 13 | 2 |
| 15 | 6 |
| 8, 9, 11, 14 | none found |

  The only sigma = 0 instances seen are F7 and H14.

**Consequence for strategy.** A Henning–Yeo-style deficiency with a *small* exceptional family will not work at the sigma <= 1 level, because sigma = 1 is too common. The tight family (sigma = 0) may still be small. The next step should therefore be:
- **(L5') Prove the corollary's configuration is reducible.** Use flexibility: when a near-tight component C of H - y has a minimum transversal that meets the edges at y, cover those edges from inside C instead of taking y. This saves the +7 for y.
- **(L6') Characterize sigma = 0.** Show that the tight connected members are exactly the F7-chains (n = 7k, m = 7k - 2, tau = 3k).

## Still open
- L5' and L6' above.
- Exhaustive confirmation at n = 11 (only heuristic so far).
