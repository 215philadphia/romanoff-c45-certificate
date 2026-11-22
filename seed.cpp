// Independent seed-histogram generator (own implementation).
// usage: ./seed q1 q2 ...   -> prints "T0 <t>" then lines "nu count" (count as decimal u128)
// State: c[x] = #surviving r mod T0 with r = x mod L (L = lcm of orders still to process),
// merged up to cyclic rotation (Booth least rotation). Greedy order: next prime is the one
// whose removal makes the remaining lcm smallest.
#include <bits/stdc++.h>
using namespace std;
typedef unsigned __int128 u128;
static int ordp(int q){int x=2%q,d=1;while(x!=1){x=x*2%q;d++;}return d;}
static long long L_of(const vector<int>&o){long long L=1;for(int d:o)L=lcm(L,(long long)d);return L;}
static vector<uint16_t> booth(const vector<uint16_t>&s){
  int n=s.size(); if(n<=1) return s;
  vector<int> f(2*n,-1); int k=0;
  for(int j=1;j<2*n;j++){
    uint16_t sj=s[j%n]; int i=f[j-k-1];
    while(i!=-1 && sj!=s[(k+i+1)%n]){ if(sj<s[(k+i+1)%n]) k=j-i-1; i=f[i]; }
    if(sj!=s[(k+i+1)%n]){ if(sj<s[k%n]) k=j; f[j-k]=-1; } else f[j-k]=i+1;
  }
  vector<uint16_t> r(n); for(int i=0;i<n;i++) r[i]=s[(k+i)%n]; return r;
}
struct H{size_t operator()(const vector<uint16_t>&v)const{uint64_t h=1469598103934665603ull;for(auto x:v){h^=x;h*=1099511628211ull;}return h;}};
static string str(u128 x){if(!x)return "0";string s;while(x){s+=char('0'+int(x%10));x/=10;}reverse(s.begin(),s.end());return s;}
int main(int argc,char**argv){
  vector<int> P; for(int i=1;i<argc;i++) P.push_back(atoi(argv[i]));
  vector<int> O; for(int q:P) O.push_back(ordp(q));
  long long T0=L_of(O);
  // greedy order
  vector<int> idx(P.size()); iota(idx.begin(),idx.end(),0);
  vector<int> order;
  while(!idx.empty()){
    int best=-1; long long bestL=LLONG_MAX;
    for(size_t t=0;t<idx.size();t++){
      vector<int> rest; for(size_t u=0;u<idx.size();u++) if(u!=t) rest.push_back(O[idx[u]]);
      long long L=L_of(rest);
      if(L<bestL||(L==bestL&&O[idx[t]]<O[idx[best]])){bestL=L;best=t;}
    }
    order.push_back(idx[best]); idx.erase(idx.begin()+best);
  }
  vector<long long> Ls(order.size()+1);
  for(size_t i=0;i<=order.size();i++){vector<int> r;for(size_t j=i;j<order.size();j++)r.push_back(O[order[j]]);Ls[i]=L_of(r);}
  unordered_map<vector<uint16_t>,u128,H> st;
  st[booth(vector<uint16_t>(Ls[0],(uint16_t)(T0/Ls[0])))]=1;
  for(size_t i=0;i<order.size();i++){
    int q=P[order[i]], d=O[order[i]]; long long L=Ls[i], nL=Ls[i+1];
    unordered_map<vector<uint16_t>,u128,H> nx; nx.reserve(st.size()*4);
    for(auto&kv:st){
      const auto&v=kv.first; u128 w=kv.second;
      vector<uint16_t> out(nL,0);
      for(long long x=0;x<L;x++) out[x%nL]+=v[x];
      nx[booth(out)]+=w*(u128)(q-d);
      for(int s=0;s<d;s++){
        vector<uint16_t> o2(nL,0);
        for(long long x=0;x<L;x++) if(x%d!=s) o2[x%nL]+=v[x];
        nx[booth(o2)]+=w;
      }
    }
    st.swap(nx);
    fprintf(stderr,"step q=%d ord=%d L %lld->%lld states=%zu\n",q,d,L,nL,st.size());
  }
  map<int,u128> hist; u128 tot=0;
  for(auto&kv:st){ hist[kv.first[0]]+=kv.second; tot+=kv.second; }
  u128 M=1; for(int q:P) M*=q;
  if(tot!=M){fprintf(stderr,"TOTAL MISMATCH\n");return 1;}
  printf("T0 %lld\n",T0);
  for(auto&kv:hist) printf("%d %s\n",kv.first,str(kv.second).c_str());
}
