import json, sys, numpy as np
from itertools import combinations
from scipy.optimize import milp, LinearConstraint, Bounds
A = json.load(open(sys.argv[1]))['set']; n = len(A); S = set(A); idx = {a: i for i, a in enumerate(A)}
E = [(idx[a], idx[(a+c)//2], idx[c]) for a, c in combinations(sorted(A), 2) if (a+c) % 2 == 0 and (a+c)//2 in S]
M = np.zeros((len(E), n))
for r, e in enumerate(E): M[r, list(e)] = 1
res = milp(c=np.ones(n), constraints=LinearConstraint(M, lb=1, ub=np.inf), integrality=np.ones(n), bounds=Bounds(0, 1),
           options={"time_limit": 500, "mip_rel_gap": 0})
x = np.round(res.x).astype(int)
assert all(x[list(e)].sum() >= 1 for e in E), "returned set is not a transversal"
print('HiGHS status:', res.status, res.message, '| tau =', int(x.sum()), '| dual bound:', getattr(res, 'mip_dual_bound', None))
# also check the claimed Sidon witness of size 37 directly
W = json.load(open(sys.argv[1]))['witness_sidon_subset']
sums = [a+b for i, a in enumerate(W) for b in W[i:]]
print('witness size', len(W), 'subset of A:', set(W) <= S, 'Sidon:', len(sums) == len(set(sums)))
