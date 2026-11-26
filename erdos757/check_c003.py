import json, sys, time
from itertools import combinations
from ortools.sat.python import cp_model
A = json.load(open(sys.argv[1]))['set']; n = len(A)
t = time.time()
# (4,5)-set check: every 4 points give >=5 distinct differences <=> no x-y = z-w equalities among
# a 4-set beyond what 5 distinct allows; check directly but efficiently: count pairs with equal differences
from collections import defaultdict
diffs = defaultdict(list)
for i, j in combinations(range(n), 2):
    diffs[abs(A[i]-A[j])].append((i, j))
bad = 0
for d, prs in diffs.items():
    if len(prs) < 2: continue
    for (a, b), (c, e) in combinations(prs, 2):
        S = {a, b, c, e}
        if len(S) == 3: continue          # 3-AP: one repeated difference inside 3 points is allowed in a 4-set
        # four distinct points with two equal differences: need <5 distinct among 6 -> check
        pts = sorted(A[k] for k in S)
        if len({abs(x-y) for x, y in combinations(pts, 2)}) < 5: bad += 1
# 4-sets containing a 3-AP plus another point: covered by enumerating every equal-difference pair, incl. shared point
for d, prs in diffs.items():
    for (a, b), (c, e) in combinations(prs, 2):
        S = {a, b, c, e}
        if len(S) == 3:
            for k in range(n):
                if k in S: continue
                pts = sorted(A[x] for x in S | {k})
                if len({abs(x-y) for x, y in combinations(pts, 2)}) < 5: bad += 1
print('(4,5)-set violations:', bad, f'({time.time()-t:.1f}s)')
S = set(A); idx = {a: i for i, a in enumerate(A)}
E = [(idx[a], idx[(a+c)//2], idx[c]) for a, c in combinations(sorted(A), 2) if (a+c) % 2 == 0 and (a+c)//2 in S]
print('n =', n, 'm =', len(E))
m = cp_model.CpModel(); x = [m.NewBoolVar('') for _ in range(n)]
for e in E: m.AddBoolOr([x[v] for v in e])
m.Minimize(sum(x)); s = cp_model.CpSolver(); s.parameters.max_time_in_seconds = 300; s.parameters.num_workers = 4
r = s.Solve(m)
print('status', s.StatusName(r), 'tau =', s.ObjectiveValue(), 'lower bound', s.BestObjectiveBound(), '=> h =', n - int(s.ObjectiveValue()))
tau = int(s.ObjectiveValue())
print('3n/7 =', 3*n/7, '| tau > 3n/7:', 7*tau > 3*n, '| h/n =', (n-tau)/n, '< 4/7:', 7*(n-tau) < 4*n, '| Conj S (7tau <= 2n+m+2):', 7*tau, '<=', 2*n+len(E)+2, 7*tau <= 2*n+len(E)+2)
