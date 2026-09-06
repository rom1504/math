// Exact integer scan of all switching-gauged signings, n <= 8 by default.
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <map>
#include <vector>
int main(int argc,char**argv){
 int top=argc>1?std::stoi(argv[1]):8;
 for(int n=3;n<=top;++n){
  int ns=1<<(n-1); std::vector<std::vector<int>> x(ns,std::vector<int>(n,1));
  for(int u=0;u<ns;++u)for(int i=1;i<n;++i)x[u][i]=(u>>(i-1)&1)?-1:1;
  std::vector<std::pair<int,int>> ed;
  for(int i=1;i<n;++i)for(int j=i+1;j<n;++j)ed.push_back({i,j});
  int ne=ed.size(); if(ne>=63)return 2;
  std::vector<std::vector<int>> col(ne,std::vector<int>(ns));
  for(int e=0;e<ne;++e)for(int u=0;u<ns;++u)col[e][u]=x[u][ed[e].first]*x[u][ed[e].second];
  std::vector<std::vector<int>> dot(ns,std::vector<int>(ns));
  for(int u=0;u<ns;++u)for(int v=0;v<ns;++v){int s=0;for(int i=0;i<n;++i)s+=x[u][i]*x[v][i];dot[u][v]=std::abs(s);}
  std::vector<int> q(ns),sg(ne,1);
  for(int u=0;u<ns;++u)for(int i=0;i<n;++i)for(int j=i+1;j<n;++j)q[u]+=x[u][i]*x[u][j];
  int best=n*n;std::map<int,uint64_t> hist,witness;uint64_t count=0;uint64_t first=0;
  std::map<std::vector<int>,std::map<int,uint64_t>> spectra;
  for(uint64_t k=0;k<(1ULL<<ne);++k){
   if(k){int e=__builtin_ctzll(k);sg[e]=-sg[e];for(int u=0;u<ns;++u)q[u]+=2*sg[e]*col[e][u];}
   int cap=0;for(int v:q)cap=std::max(cap,std::abs(v));
   if(cap>best)continue;
   if(cap<best){best=cap;hist.clear();witness.clear();spectra.clear();count=0;first=k^(k>>1);}
   ++count;int ext=n*n;
   for(int a=0;a<ns;++a){int v=0;for(int u=0;u<ns;++u){v=std::max(v,std::abs(q[u])+dot[a][u]);if(v>=ext)break;}ext=std::min(ext,v);}
   if(!hist.count(ext-best))witness[ext-best]=k^(k>>1);
   ++hist[ext-best];
   std::vector<int> profile(best+1);for(int v:q)++profile[std::abs(v)];
   if(!spectra[profile].count(ext-best))spectra[profile][ext-best]=k^(k>>1);
  }
  std::cout<<"n="<<n<<" M="<<best<<" count="<<count<<" first_gray="<<first<<" E_hist=";
  for(auto [e,c]:hist)std::cout<<e<<":"<<c<<",";
  std::cout<<" witness_gray=";for(auto [e,g]:witness)std::cout<<e<<":"<<g<<",";
  std::cout<<" sharp_target="<<1.5*best/n<<"\n"<<std::flush;
  int mixed=0;for(const auto& [p,w]:spectra)if(w.size()>1)++mixed;
  std::cout<<"spectra="<<spectra.size()<<" mixed_E_spectra="<<mixed;
  for(const auto& [p,w]:spectra)if(w.size()>1){std::cout<<" example_profile=";for(int v:p)std::cout<<v<<",";std::cout<<" E_gray=";for(auto[e,g]:w)std::cout<<e<<":"<<g<<",";break;}
  std::cout<<"\n"<<std::flush;
 }
}
