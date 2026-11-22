# Brute-force check of nu-histogram factorization on a small instance.
import numpy as np, subprocess
from math import lcm, prod
from collections import Counter
def ordp(q):
    x,d=2%q,1
    while x!=1: x=x*2%q; d+=1
    return d
seedP=[3,5,7]; cl={11:[23,89],23:[47]}
P=seedP+[23,89,47]; M=prod(P); T=lcm(*[ordp(q) for q in P])
r=np.arange(T)
pw={q:np.array([pow(2,int(x),q) for x in r]) for q in P}
hist=Counter()
for a0 in range(0,M,20000):
    a=np.arange(a0,min(M,a0+20000))
    ok=np.ones((len(a),T),bool)
    for q in P: ok &= (a[:,None]%q)!=pw[q][None,:]
    hist.update(ok.sum(1).tolist())
# formula
out=subprocess.run(["./seed"]+list(map(str,seedP)),capture_output=True,text=True).stdout.split("\n")
sh={int(x):int(y) for x,y in (l.split() for l in out[1:] if l.strip())}
def cluster(l,qs):
    st={0:1}
    for q in qs:
        n={}
        for k,w in st.items():
            n[k]=n.get(k,0)+w*(q-l)+w*k
            if k<l: n[k+1]=n.get(k+1,0)+w*(l-k)
        st=n
    return {l-k:w for k,w in st.items()}
F=dict(sh)
for l,qs in cl.items():
    c=cluster(l,qs); G=Counter()
    for a,wa in F.items():
        for b,wb in c.items(): G[a*b]+=wa*wb
    F=G
print("M=",M,"T=",T,"brute entries",len(hist),"formula entries",len({k:v for k,v in F.items() if v}))
print("histograms identical:", dict(hist)=={k:v for k,v in F.items() if v})
