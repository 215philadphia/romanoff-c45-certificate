# Erdős #757 / Tao constant C_5b: working notes (AI-assisted, single-model)

## Status update, 2026-11-26 (supersedes the notes below where they conflict)

**Prior work found after these notes were written.** A prior-art check on 2026-11-26 found results that settle
much of what is listed as open below:

- **Exact values.** f(16)=10, f(18)=11, f(20)=12, f(22)=13 (and f(23)=14, f(27)=16) were established by
  computer-assisted certificates in the public repository
  [hturcotte1/Erdos-757](https://github.com/hturcotte1/Erdos-757) (results dated 2026-10-17 onward). The
  corresponding lower bounds were obtained independently in
  [fortunto2/erdos-757-sidon45](https://github.com/fortunto2/erdos-757-sidon45),
  Zenodo doi:[10.5281/zenodo.23494596](https://doi.org/10.5281/zenodo.23494596) (October 2026).
  So "n = 16 is the decisive size" below is closed: no (4,5)-set with at most 24 points beats 4/7.
- **Upper bound.** A self-similar family posted by an anonymous reader on Erdős Problem a Day (October 2026),
  verified in hturcotte1/Erdos-757, gives c* <= 5/9 (reported, not yet refereed). Its first member below 4/7
  is a 65-point set, C003.

**Both conjectures in this branch are refuted by C003.** Independent check (2026-11-26):
C003 is a (4,5)-set with n = 65, m = 63 three-term APs, tau = 28, so h = 37.
- (4,5)-property: exhaustive check of all equal-difference pairs, exact integer arithmetic, 0 violations (`check_c003.py`).
- h >= 37: the stored 37-element witness is a subset of C003 and is Sidon (all pairwise sums distinct), checked directly (`tau_highs.py`).
- h <= 37, i.e. tau >= 28: proved optimal by two independent exact solvers, OR-Tools CP-SAT (`check_c003.py`)
  and HiGHS MILP via scipy (`tau_highs.py`), both tau = 28 with matching lower bound. The hypergraph is a single
  connected component, too large for the plain branch-and-bound used elsewhere here; no hand-checkable certificate.
- Consequence, independent of the 5/9 family: by Ma–Tang Lemma 3.6 (gluing) f(65k) <= 37k, and f is subadditive,
  so by Fekete c* = inf f(n)/n <= 37/65 ~ 0.5692 < 4/7 ~ 0.5714. (Construction due to the anonymous reader; only the check is ours.)
- "tau <= 3n/7 on the class (linear, F7-free, |E(S)| <= |S|-2)" fails: 7*28 = 196 > 195 = 3*65.
  Its consequence C_5b = 4/7 is therefore not available by this route (and c* <= 37/65 < 4/7 is verified below; 5/9 is reported).
- Conjecture S, 7*tau <= 2n + m + 2 (PROOF_SKETCH.md), fails: 196 > 2*65 + 63 + 2 = 195.

**My published note.** Zenodo concept doi:[10.5281/zenodo.24547831](https://doi.org/10.5281/zenodo.24547831),
current version v3, doi:[10.5281/zenodo.24571312](https://doi.org/10.5281/zenodo.24571312). Versions 1 and 2
wrongly said no earlier statement of f(20), f(22) was found; v3 corrects this, credits the work above, and
gives a short, computer-free proof of f(18)=11, f(20)=12, f(22)=13 from Henning–Yeo (Graphs Combin. 2021)
Theorem 6 including its strict clause, with Ma–Tang Lemmas 2.4, 2.5, 5.2 and 3.6. It claims no new values.

**Still open (as far as checked):** the constant itself, 9/17 <= c* <= 5/9 (reported). The repository
hturcotte1/Erdos-757 proposes 11*tau <= 3n + 2m on the class above (proved there by computation only for
n <= 27); a general proof would give c* >= 6/11.

---

## Original notes (2026-11-22)

Known range: 9/17 <= C_5b <= 4/7. Prior art: Ma–Tang arXiv:2602.23282, the SK small-n certificates linked from Tao's 5b page, and the erdosproblemaday.com #757 report.

## Conjecture under test (not proved, not found in the literature) — REFUTED, see status update above
Every linear, F7-free 3-uniform hypergraph with |E(S)| <= |S|-2 for every vertex set S has tau <= 3n/7.
If true, it implies C_5b = 4/7, via Ma–Tang's lemmas.

## Executed results
- The conjecture holds automatically for n <= 15 (Henning–Yeo's 17tau <= 5n+3m with m <= n-2). n = 16 is the first open size.
- Transversal flexibility of the known extremal blocks (SK certificates, n = 6..21): every vertex lies in some minimum transversal, except for the blocks n = 7, 9, 11. See glue_run_part1.txt.
- Gluing search over all 91 pairs of these blocks: disjoint union + 2 edges, one shared vertex + 1 edge, two shared vertices. Linearity, F7-freeness and hereditary sparsity were enforced throughout. **0 hypergraphs with tau > 3n/7.**
  Lemmas used to skip cases:
  1. Two flexible blocks plus two cross edges can never raise tau.
  2. Sharing a vertex pair (u,w), where u is in some minimum transversal of A and w is in some minimum transversal of B, saves at least 1.
- Exhaustive realizable search (realsearch.py) agrees with the known exact values for n = 8, 9. It is too slow for n >= 10.
- The abstract SAT search (c5b_n_search.py) is too slow at n = 11 on 2 cores.

## Open (as of 2026-11-22; see status update above)
n = 16 is the decisive size. The sharpness examples of Henning–Yeo 2021 (paywalled) have been requested from the authors.
