import itertools, sympy as sp, random
from fractions import Fraction
from tools import aps, is45
def realize(n, edges, trials=1):
    """edges: list of 3-tuples of vertex ids. Try every middle assignment; solve x_a+x_c=2x_b with x0=0,x1=1 fixed
    plus generic free params; return list of valid (4,5) integer realizations whose AP hypergraph == edges."""
    out=[]
    target={frozenset(e) for e in edges}
    for mids in itertools.product(range(3), repeat=len(edges)):
        xs=sp.symbols(f'x0:{n}')
        eqs=[]
        for e,mi in zip(edges,mids):
            b=e[mi]; a,c=[v for v in e if v!=b]
            eqs.append(sp.Eq(xs[a]+xs[c],2*xs[b]))
        sol=sp.solve(eqs, xs, dict=True)
        if not sol: continue
        sol=sol[0]
        free=[x for x in xs if x not in sol]
        for _ in range(trials):
            sub={f:random.randint(-10**6,10**6) for f in free}
            vals=[sp.nsimplify(sp.sympify(sol.get(x,x)).subs(sub)) for x in xs]
            den=sp.ilcm(*[sp.fraction(v)[1] for v in vals])
            iv=[int(v*den) for v in vals]
            if len(set(iv))<n: continue
            E={frozenset(iv.index(t) for t in ap) for ap in aps(iv)}
            if E==target and is45(iv): out.append((mids,iv))
    return out
