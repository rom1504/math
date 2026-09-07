#include <array>
#include <bit>
#include <cstdint>
#include <iostream>
#include <random>
#include <unordered_map>
#include <vector>

struct Key { std::array<int,6> t; bool operator==(Key const&o)const{return t==o.t;} };
struct Hash { size_t operator()(Key const&k)const { size_t h=0; for(int v:k.t) h=h*1000003u+std::hash<int>{}(v); return h;} };

int main(int argc,char**argv){
  constexpr int n=8, m=21;
  uint64_t limit=argc>1?std::stoull(argv[1]):(uint64_t(1)<<m);
  bool random=argc>3;
  std::vector<uint64_t> sm(1<<(n-1)); std::vector<int> base(sm.size());
  for(uint64_t xm=0;xm<sm.size();++xm){uint64_t s=0;int e=0,q=0;for(int j=1;j<n;++j)q+=(xm>>(j-1)&1)?-1:1;
    for(int i=1;i<n;++i)for(int j=i+1;j<n;++j,++e){int p=((xm>>(i-1)&1)?-1:1)*((xm>>(j-1)&1)?-1:1);q+=p;if(p<0)s|=uint64_t(1)<<e;}sm[xm]=s;base[xm]=q;}
  std::unordered_map<Key,std::pair<int,uint64_t>,Hash> seen; std::mt19937_64 gen(1);
  for(uint64_t it=0;it<limit;++it){uint64_t mask=random?(gen()&((uint64_t(1)<<m)-1)):it;int aw=std::popcount(mask),Q=0;
    for(size_t x=0;x<sm.size();++x){int q=base[x]-2*aw+4*std::popcount(mask&sm[x]);Q=std::max(Q,std::abs(q));}
    int A[n][n]{};for(int j=1;j<n;++j)A[0][j]=A[j][0]=1;int e=0;for(int i=1;i<n;++i)for(int j=i+1;j<n;++j,++e)A[i][j]=A[j][i]=(mask>>e&1)?-1:1;
    int P[n][n],T[n][n];for(int i=0;i<n;++i)for(int j=0;j<n;++j)P[i][j]=A[i][j];Key key{};
    for(int p=2;p<=8;++p){for(int i=0;i<n;++i)for(int j=0;j<n;++j){int z=0;for(int k=0;k<n;++k)z+=P[i][k]*A[k][j];T[i][j]=z;}for(int i=0;i<n;++i)for(int j=0;j<n;++j)P[i][j]=T[i][j];if(p>=3){int tr=0;for(int i=0;i<n;++i)tr+=P[i][i];key.t[p-3]=tr;}}
    auto [pos,ins]=seen.emplace(key,std::pair{Q,mask});if(!ins&&pos->second.first!=Q && (argc < 3 || pos->second.first==10 || Q==10)){std::cout<<"FOUND Q1="<<pos->second.first<<" mask1="<<pos->second.second<<" Q2="<<Q<<" mask2="<<mask<<" traces=";for(int v:key.t)std::cout<<v<<',';std::cout<<'\n';return 0;}
  }
  std::cout<<"none; unique spectra="<<seen.size()<<'\n';
}
