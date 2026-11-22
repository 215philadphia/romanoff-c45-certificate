# Romanoff's constant: upper-density bound below 0.490179107442766

Let A = {n odd : n = p + 2^k, p prime, k ≥ 0}. This repository certifies

    upper density of A  <  0.490179107442766

using a finite obstruction set of 113 primes (`primes_113.txt`). The previous record was 0.490180063290061 (Yoo 2026); this improves it by about 9.6 × 10⁻⁷.

## Inequality used
This is the same obstruction inequality used by Habsieger–Roblot, Chen–Dai–Li, Griego and Yoo. With M = ∏q, T = lcm ord_q(2), and ν(a) = #{r mod T : gcd(a − 2^r, M) = 1},

    dbar(A) ≤ Σ_{a mod M} min( 1/(2M), ν(a) / (T φ(M) log 2) ).

## Structure of the computation
- **Seed.** The primes {3, 5, 7, 11, 13, 17, 19, 29, 31, 37, 41, 61, 73, 257} have orders dividing T0 = 5040. Their exact ν-histogram is computed by `seed.cpp`, a dynamic program over r mod L that projects as primes are eliminated and merges states up to cyclic rotation (Booth canonical form). It has 2,104 entries for the 13-prime Griego seed; the totals are checked against M0.
- **Clusters.** Each remaining prime has a prime order l with gcd(l, 5040) = 1, grouped by l (90 clusters with distinct l). By the Chinese remainder theorem on r, ν = ν_seed · ∏ ν_l, where ν_l = l − (number of distinct exponents hit). The cluster histograms are exact integers.
- **Rigorous rounding (`mult.cpp`).**
  - The product ∏(ν_l / l) is stored as an integer key K on a grid of step 1/G. After each cluster it is rounded **up**: K → ⌈K·v/l⌉, in exact 128-bit arithmetic.
  - Probabilities are stored as 62-bit fixed point and also rounded **up**: W → ⌈W·p/2^62⌉.
  - Every term min(c, x) is nondecreasing in x and nonnegative, so both roundings can only raise the computed sum. The result is therefore a rigorous upper bound.
- **Assembly (`prove2.py`).** The final sum is assembled in exact rational arithmetic. log 2 is replaced by the lower bound 2·Σ_{j<60} 1/((2j+1)·3^(2j+1)).

## Checks
- `romanoff.py`, an exact unrounded prover, reproduces Griego's 0.4902494078111542556868594292… to 30 digits.
- `romanoff.py` also reproduces Yoo's bound as 0.4901800632900604571…; `clusters_yoo38.json` is included as a calibration input.
- On the 53-prime set (`clusters_53.json`, G = 5·10⁶), the exact-weight prover `prove.py` gives 0.490179502482281 and the fixed-point prover gives 0.490179502485397, an agreement of 3 × 10⁻¹². The gap comes from the upward weight rounding.
- `brute.py` enumerates every residue class on a small instance (M = 10,101,945, T = 3036) and confirms the seed-times-clusters factorization: the histograms are identical.

## Reproduce
    ./verify.sh        # g++ and python3 with numpy; ~6 min, ~3 GB RAM

The expected output is in `EXPECTED_OUTPUT.txt`. The SHA-256 of the generated `seedY.txt` is
6f5bca5058744cad8d2c913ae2569a728dc9598785f64057821b3d8f1c299310.

## Provenance
The prime set was found with a float-valued greedy search over primes whose multiplicative order of 2 is prime. The prime set, the code and this text were prepared with AI assistance (Claude, Anthropic). The verification was run by the AI in its own environment; it has not yet been reproduced independently.
