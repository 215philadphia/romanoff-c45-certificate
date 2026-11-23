from itertools import combinations
def aps(A):
    S=set(A); A=sorted(A); E=[]
    for i,a in enumerate(A):
        for c in A[i+1:]:
            if (a+c)%2==0 and (a+c)//2 in S: E.append((a,(a+c)//2,c))
    return E
def is45(A):
    for B in combinations(A,4):
        if len({abs(x-y) for x,y in combinations(B,2)})<5: return False
    return True
def tau(n, edges):
    # edges: list of frozensets of vertex indices; branch on unhit edge
    best=[n]
    def rec(chosen, k):
        if k>=best[0]: return
        for e in edges:
            if not (e & chosen):
                for v in e: rec(chosen|{v}, k+1)
                return
        best[0]=k
    rec(frozenset(),0); return best[0]
