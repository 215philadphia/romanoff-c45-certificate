# Erdős #757 / Tao constant C_5b: working notes (2026-11-22, AI-assisted, single-model)

Known range: 9/17 <= C_5b <= 4/7. Prior art: Ma–Tang arXiv:2602.23282, the SK small-n certificates linked from Tao's 5b page, and the erdosproblemaday.com #757 report.

## Conjecture under test (not proved, not found in the literature)
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

## Open
n = 16 is the decisive size. The sharpness examples of Henning–Yeo 2021 (paywalled) have been requested from the authors.
