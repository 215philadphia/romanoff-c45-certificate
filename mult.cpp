// Rigorous multiplier distribution with upward rounding (keys and weights).
// stdin:  G  nclusters  then per cluster: l  nv  (v_i  p_i)*   [p_i = ceil(prob * 2^62), as u64]
//         nthr  then thresholds t_j (int64, ascending not required)
// stdout: "n <entries>" "tot <sum W>" then per threshold: "<t> <S0> <S1>" where
//   S0 = sum_{K<=t} W,  S1 = sum_{K<=t} W*K   (exact integer sums of upward-rounded weights)
// Keys: K -> ceil(K*v/l)  (exact). Weights: W -> ceil(W*p / 2^62) (exact u128 math).
#include <bits/stdc++.h>
using namespace std;
typedef unsigned long long u64; typedef unsigned __int128 u128; typedef long long i64;
static string s128(u128 x){if(!x)return "0";string s;while(x){s+=char('0'+int(x%10));x/=10;}reverse(s.begin(),s.end());return s;}
int main(){
  i64 G; int nc; cin>>G>>nc;
  vector<i64> K{G}; vector<u64> W{(u64)1<<62};
  for(int c=0;c<nc;c++){
    i64 l; int nv; cin>>l>>nv; vector<i64> v(nv); vector<u64> p(nv);
    for(int i=0;i<nv;i++) cin>>v[i]>>p[i];
    // each branch is sorted (monotone map); k-way merge with dedupe
    vector<vector<i64>> bk(nv); vector<vector<u64>> bw(nv);
    for(int i=0;i<nv;i++){
      bk[i].resize(K.size()); bw[i].resize(K.size());
      for(size_t j=0;j<K.size();j++){
        bk[i][j]=(i64)(((u128)K[j]*(u128)v[i]+(u128)(l-1))/(u128)l);
        u128 prod=(u128)W[j]*(u128)p[i];
        bw[i][j]=(u64)((prod+(((u128)1<<62)-1))>>62);
      }
    }
    vector<i64> nk; vector<u64> nw; nk.reserve(K.size()*nv); nw.reserve(K.size()*nv);
    vector<size_t> pos(nv,0);
    while(true){
      int best=-1;
      for(int i=0;i<nv;i++) if(pos[i]<bk[i].size() && (best<0||bk[i][pos[i]]<bk[best][pos[best]])) best=i;
      if(best<0) break;
      i64 key=bk[best][pos[best]]; u64 w=bw[best][pos[best]]; pos[best]++;
      if(!nk.empty()&&nk.back()==key) nw.back()+=w; else {nk.push_back(key); nw.push_back(w);}
    }
    K.swap(nk); W.swap(nw);
    fprintf(stderr,"cluster %d (l=%lld): entries=%zu\n",c,l,K.size());
  }
  int nt; cin>>nt; vector<i64> t(nt); for(auto&x:t)cin>>x;
  vector<int> ord(nt); iota(ord.begin(),ord.end(),0); sort(ord.begin(),ord.end(),[&](int a,int b){return t[a]<t[b];});
  vector<u128> S0(nt),S1(nt); u128 a0=0,a1=0; size_t j=0;
  for(int oi:ord){ while(j<K.size()&&K[j]<=t[oi]){a0+=W[j];a1+=(u128)W[j]*(u128)K[j];j++;} S0[oi]=a0;S1[oi]=a1; }
  u128 tot=0; for(u64 w:W) tot+=w;
  printf("n %zu\ntot %s\n",K.size(),s128(tot).c_str());
  for(int i=0;i<nt;i++) printf("%lld %s %s\n",t[i],s128(S0[i]).c_str(),s128(S1[i]).c_str());
}
