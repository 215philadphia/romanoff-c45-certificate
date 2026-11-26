# Proof sketch: C_5b = 4/7 (Erdős #757). Status: incomplete, with one explicit GAP

> **REFUTED (2026-11-26).** Conjecture S (7τ ≤ 2n + m + 2) and the target τ ≤ 3n/7 are false: the 65-point (4,5)-set C003 has m = 63 and τ = 28 (7τ = 196 > 195), and it gives c* ≤ 37/65 < 4/7, so C_5b = 4/7 is false too. Kept for the record; see RESULTS.md.

Every step is labeled PROVEN, EXECUTED (computer-verified), or GAP. This is single-model work and has not been independently reviewed.

## 0. Reduction (Ma–Tang 2026, cited)
For a (4,5)-set A, the 3-AP hypergraph H(A) has three properties:
- It is linear.
- It is F7-free.
- It satisfies |E(S)| <= |S| - 2 for every vertex set S spanning an edge. This is Lemma 2.5 applied to the subset S, which is itself a (4,5)-set.

Also h(A) = n - tau(H(A)), and C_5b = inf f(n)/n. So it suffices to prove:

**Conjecture C.** If H is linear and hereditarily sparse (|E(S)| <= |S| - 2), then tau(H) <= 3n/7.

## 1. A strengthened, induction-ready form
**Conjecture S.** If H is linear, hereditarily sparse and connected, then 7 tau(H) <= 2n + m + 2.

- S implies C (PROVEN). Apply S to each component (n_i, m_i). Each component has m_i <= n_i - 2, so 7 tau = sum 7 tau_i <= sum (2 n_i + m_i + 2) <= sum 3 n_i = 3n.
- S is false without connectivity (PROVEN). Two disjoint copies of the 14-point hypergraph give 84 > 82. So the +2 must be counted per component.
- S needs sparsity (PROVEN). AG(2,3) minus a point (n=8, m=8, tau=4) violates S, but it has 8 > 6 edges.
- S does not seem to need F7-freeness (EXECUTED). F7 (n=7, m=5, tau=3) attains S with equality.
- Equality cases known (EXECUTED): F7 (7,5,3), the 14-point Ma–Tang block (14,12,6), and the 21-point SK block (21,19,9). Every other known extremal block satisfies S strictly.
- S verified exhaustively for all n <= 10 (EXECUTED, CP-SAT, strong_check.py). This check used the stronger global form with +2 even for disconnected H.

S has the same shape as the Chvátal–McDiarmid bound 4 tau <= n + m and the Henning–Yeo bound 17 tau <= 5n + 3m. This makes a minimal-counterexample proof plausible.

## 2. Minimal-counterexample skeleton for S
Let H be a counterexample to S with n + m minimal.

| Step | Statement | Status |
|---|---|---|
| L1 | n >= 11 | EXECUTED (SAT, n <= 10) |
| L2 | H has no isolated vertex (H connected, n >= 2) | PROVEN (definition) |
| L3 | Every proper connected sub-hypergraph satisfies S | PROVEN (minimality plus heredity of the class) |
| L4 | Gluing lemmas: H is not two flexible blocks joined by two edges, nor two blocks sharing a vertex that lies in a minimum transversal on both sides | PROVEN for blocks satisfying S (see RESULTS.md) |
| L5 | Low-degree reduction: H has no vertex of degree <= 2 / no "removable" configuration | **GAP** |
| L6 | Discharging: if no reducible configuration exists, then sparsity (m <= n - 2, hereditary) contradicts 7 tau > 2n + m + 2 | **GAP** |

**Why L5 is the core difficulty (PROVEN).** The naive reduction deletes a vertex v of degree d. It loses 1 from tau but lowers the right-hand side 2n + m + 2c by 2 + d (when no component splits off). It closes only if d >= 5. But hereditary sparsity forces average degree < 3. So plain vertex deletion cannot work. One needs to delete larger configurations (an edge together with its neighbourhood), or use a deficiency term as in Henning–Yeo 2018/2021.

**Suggested route.** Copy the Henning–Yeo deficiency method. Prove 7 tau <= 2n + m + 2 - def(H), where def counts copies of the equality blocks F7, H14 and H21. Then reduce configurations around vertices of degree 1–2, using sparsity to bound how often they occur. The Henning–Yeo 2021 paper (requested from the authors) contains exactly this machinery for 17 tau <= 5n + 3m.

## 3. What would refute the approach
- A connected, linear, hereditarily sparse hypergraph with 7 tau > 2n + m + 2. The first open size is n = 11. Such an example would refute S but not necessarily C.
- A hereditarily sparse hypergraph with tau > 3n/7. This would refute C, and would also need a realizability check before it says anything about C_5b.
