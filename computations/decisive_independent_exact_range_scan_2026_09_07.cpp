// Exhaustive switching-gauged scan of BOTH original cap and half-width.
// Integer energies and comparisons only. Output is JSON lines.
#include <algorithm>
#include <cstdint>
#include <iostream>
#include <map>
#include <vector>
int main(int argc,char**argv){
 const int top=argc>1?std::stoi(argv[1]):8;
 if(top>9)return 2;
 for(int n=3;n<=top;++n){
  const int ns=1<<(n-1);
  std::vector<std::vector<int16_t>> x(ns,std::vector<int16_t>(n,1));
  for(int u=0;u<ns;++u)for(int i=1;i<n;++i)x[u][i]=(u>>(i-1)&1)?-1:1;
  std::vector<std::pair<int,int>> ed;
  for(int i=1;i<n;++i)for(int j=i+1;j<n;++j)ed.emplace_back(i,j);
  const int ne=ed.size(); if(ne>=63)return 2;
  std::vector<std::vector<int16_t>> col(ne,std::vector<int16_t>(ns));
  for(int e=0;e<ne;++e)for(int u=0;u<ns;++u)col[e][u]=x[u][ed[e].first]*x[u][ed[e].second];
  std::vector<int16_t> q(ns);std::vector<int>sg(ne,1);
  for(int u=0;u<ns;++u)for(int i=0;i<n;++i)for(int j=i+1;j<n;++j)q[u]+=x[u][i]*x[u][j];
  int bestcap=n*n,bestwidth=2*n*n,mincapamongwidth=n*n,maxcapamongwidth=0;
  uint64_t countcap=0,countwidth=0,firstwidth=0,firstcap=0;
  std::map<int,uint64_t> width_cap_hist;
  const int childbits=(n-2)*(n-3)/2;
  std::vector<unsigned char> childstatus(1ULL<<childbits,0);
  std::vector<uint64_t> childwitness(1ULL<<childbits,0);
  std::vector<uint64_t> noncap_width_parents;
  for(uint64_t k=0;k<(1ULL<<ne);++k){
   if(k){int e=__builtin_ctzll(k);sg[e]=-sg[e];for(int u=0;u<ns;++u)q[u]+=2*sg[e]*col[e][u];}
   int lo=q[0],hi=q[0];for(int v:q){lo=std::min(lo,v);hi=std::max(hi,v);}
   const int cap=std::max(hi,-lo),width=hi-lo;
   if(cap<bestcap){bestcap=cap;countcap=0;firstcap=k^(k>>1);}if(cap==bestcap)++countcap;
   if(width<bestwidth){bestwidth=width;countwidth=0;firstwidth=k^(k>>1);mincapamongwidth=n*n;maxcapamongwidth=0;width_cap_hist.clear();std::fill(childstatus.begin(),childstatus.end(),0);noncap_width_parents.clear();}
   if(width==bestwidth){++countwidth;mincapamongwidth=std::min(mincapamongwidth,cap);maxcapamongwidth=std::max(maxcapamongwidth,cap);++width_cap_hist[cap];
    // Delete vertex zero and gauge vertex one to all positive.
    int cm=0,bit=0,e=n-2;
    for(int i=2;i<n;++i)for(int j=i+1;j<n;++j,++e,++bit)
     if(sg[e]*sg[i-2]*sg[j-2]<0)cm|=1<<bit;
    childstatus[cm]|=1;childwitness[cm]=k^(k>>1);
    // Known exact cap minima through nine; an independent full cap scan is
    // still performed above and printed, so this target is checked at EOF.
    const int knownM[]={0,0,1,3,4,4,5,9,10,12};
    if(cap==knownM[n])childstatus[cm]|=2;
    else noncap_width_parents.push_back(k^(k>>1));
   }
  }
  std::cout<<"{\"n\":"<<n<<",\"M\":"<<bestcap<<",\"twice_W\":"<<bestwidth<<",\"M_count\":"<<countcap<<",\"W_count\":"<<countwidth<<",\"first_M_gray\":"<<firstcap<<",\"first_W_gray\":"<<firstwidth<<",\"min_cap_of_W_minimizers\":"<<mincapamongwidth<<",\"max_cap_of_W_minimizers\":"<<maxcapamongwidth<<",\"W_cap_hist\":{";
  bool first=true;for(auto[c,k]:width_cap_hist){if(!first)std::cout<<",";first=false;std::cout<<"\""<<c<<"\":"<<k;}
  uint64_t children=0,badchildren=0,firstbad=0,badparent=0;
  for(uint64_t c=0;c<childstatus.size();++c)if(childstatus[c]){++children;if(childstatus[c]==1){if(!badchildren){firstbad=c;badparent=childwitness[c];}++badchildren;}}
  std::map<int,uint64_t> repair_hist;
  uint64_t no_repair_witness=0;
  for(auto code:noncap_width_parents){
   int A[9][9]{};for(int i=0;i<n;++i)for(int j=i+1;j<n;++j)A[i][j]=A[j][i]=1;
   for(int e=0;e<ne;++e)if(code>>e&1)A[ed[e].first][ed[e].second]=A[ed[e].second][ed[e].first]=-1;
   int repairs=0;
   for(int del=0;del<n;++del){
    std::vector<int> vv;for(int i=0;i<n;++i)if(i!=del)vv.push_back(i);
    int cm=0,bit=0;
    for(int i=1;i<n-1;++i)for(int j=i+1;j<n-1;++j,++bit)
     if(A[vv[i]][vv[j]]*A[vv[0]][vv[i]]*A[vv[0]][vv[j]]<0)cm|=1<<bit;
    repairs+=(childstatus[cm]&2)!=0;
   }
   ++repair_hist[repairs];if(!repairs)no_repair_witness=code;
  }
  std::cout<<"},\"width_child_classes\":"<<children<<",\"width_children_without_cap_opt_extension\":"<<badchildren<<",\"first_bad_child\":"<<firstbad<<",\"bad_parent_gray\":"<<badparent<<",\"noncap_parent_repairable_vertex_hist\":{";
  first=true;for(auto[c,k]:repair_hist){if(!first)std::cout<<",";first=false;std::cout<<"\""<<c<<"\":"<<k;}
  std::cout<<"},\"no_single_row_repair_witness\":"<<no_repair_witness<<"}\n"<<std::flush;
 }
}
