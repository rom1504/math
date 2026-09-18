#include <algorithm>
#include <bit>
#include <cstdint>
#include <iostream>
#include <vector>
int main(){constexpr int n=8,m=21,E=28;std::vector<uint64_t> sm(128);std::vector<int> base(128);
for(uint64_t xm=0;xm<128;++xm){uint64_t s=0;int e=0,q=0;for(int j=1;j<n;++j)q+=(xm>>(j-1)&1)?-1:1;for(int i=1;i<n;++i)for(int j=i+1;j<n;++j,++e){int p=((xm>>(i-1)&1)?-1:1)*((xm>>(j-1)&1)?-1:1);q+=p;if(p<0)s|=uint64_t(1)<<e;}sm[xm]=s;base[xm]=q;}
std::array<unsigned long long,4> best;best.fill(~0ull);std::array<uint64_t,4> bm{};std::array<int,4>bq{};
for(uint64_t a=0;a<(uint64_t(1)<<m);++a){int aw=std::popcount(a),Q=0;std::array<unsigned long long,4> sums{};for(int x=0;x<128;++x){long long q=base[x]-2*aw+4*std::popcount(a&sm[x]);Q=std::max(Q,int(std::abs(q)));long long z=q*q,p=z;for(int k=0;k<4;++k){sums[k]+=p;p*=z;}}
for(int k=0;k<4;++k)if(sums[k]<best[k]){best[k]=sums[k];bm[k]=a;bq[k]=Q;}}
for(int k=0;k<4;++k)std::cout<<"moment "<<2*(k+1)<<" minavg="<<best[k]/128<<" mask="<<bm[k]<<" Q="<<bq[k]<<'\n';}
