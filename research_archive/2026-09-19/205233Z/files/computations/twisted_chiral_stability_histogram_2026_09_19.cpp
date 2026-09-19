// Exact projective energy and strict one-spin-local-maximum histograms.
// Full even-order sign matrices only: every local field is odd, hence nonzero.
#include <algorithm>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <vector>
int main() {
  int n; if (!(std::cin >> n) || n < 2 || n > 20 || n%2) return 2;
  int a[20][20]{};
  for(int i=0;i<n;i++) for(int j=0;j<n;j++) if(!(std::cin>>a[i][j])) return 2;
  for(int i=0;i<n;i++) for(int j=0;j<n;j++)
    if(a[i][j]!=a[j][i] || (i==j ? a[i][j]!=0 : std::abs(a[i][j])!=1)) return 2;
  const int bound=n*(n-1)/2;
  std::vector<uint64_t> positive(bound+1), stable(bound+1);
  int x[20], field[20]{}; std::fill(x,x+n,1);
  int energy=0;
  for(int i=0;i<n;i++) for(int j=0;j<n;j++) field[i]+=a[i][j];
  for(int i=0;i<n;i++) for(int j=i+1;j<n;j++) energy+=a[i][j];
  const uint64_t states=uint64_t(1)<<(n-1);
  uint64_t zeros=0; int minimum=bound, maximum=-bound;
  for(uint64_t k=0;k<states;k++) {
    if(k) {
      int v=1+__builtin_ctzll(k), old=x[v];
      energy-=2*old*field[v];
      for(int j=0;j<n;j++) if(j!=v) field[j]-=2*old*a[j][v];
      x[v]=-old;
    }
    minimum=std::min(minimum,energy); maximum=std::max(maximum,energy);
    if(energy==0) ++zeros;
    if(energy>0) {
      ++positive[energy]; bool local=true;
      // Coordinate zero is fixed only for enumeration, NOT for stability.
      for(int i=0;i<n;i++) if(x[i]*field[i]<=0) {local=false;break;}
      if(local) ++stable[energy];
    }
  }
  std::cout<<"{\"order\":"<<n<<",\"projective_states\":"<<states
    <<",\"minimum\":"<<minimum<<",\"maximum\":"<<maximum
    <<",\"cap\":"<<std::max(-minimum,maximum)<<",\"zero_count\":"<<zeros;
  auto emit=[&](const char* name,const std::vector<uint64_t>& hist) {
    std::cout<<",\""<<name<<"\":{"; bool first=true;
    for(int e=1;e<=bound;e++) if(hist[e]) {
      if(!first) std::cout<<","; first=false; std::cout<<"\""<<e<<"\":"<<hist[e];
    } std::cout<<"}";
  };
  emit("positive_energy_histogram",positive);emit("positive_stable_histogram",stable);
  std::cout<<"}\n";
}
