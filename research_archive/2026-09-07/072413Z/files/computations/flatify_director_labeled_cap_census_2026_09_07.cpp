// Exhaustive labeled, first-row-positive cap census. Integer arithmetic only.
// Build: g++ -O3 -std=c++17 computations/flatify_director_labeled_cap_census_2026_09_07.cpp -o /home/math/quadra/tmp/flatify_director_labeled_cap_census_2026_09_07
// Arguments: order (3..8). Output: JSON with every minimizing Gray-code mask.
#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <map>
#include <utility>
#include <vector>

int main(int argc,char** argv) {
  if(argc!=2) return 2;
  const int n=std::atoi(argv[1]);
  if(n<3 || n>8) return 2;
  const int states=1<<(n-1);
  std::vector<std::pair<int,int>> edges;
  for(int i=1;i<n;++i) for(int j=i+1;j<n;++j) edges.emplace_back(i,j);
  const unsigned total=1U<<edges.size();
  std::vector<std::vector<int>> products(edges.size(),std::vector<int>(states));
  std::vector<int> energy(states,0);
  for(int s=0;s<states;++s) {
    std::vector<int> x(n,1);
    for(int i=1;i<n;++i) x[i]=1-2*((s>>(i-1))&1);
    for(int i=0;i<n;++i) for(int j=i+1;j<n;++j) energy[s]+=x[i]*x[j];
    for(unsigned e=0;e<edges.size();++e)
      products[e][s]=x[edges[e].first]*x[edges[e].second];
  }
  int best=n*n;
  std::vector<unsigned> minimizers;
  std::map<int,unsigned> histogram;
  for(unsigned step=0;step<total;++step) {
    const unsigned gray=step^(step>>1);
    if(step) {
      const unsigned e=__builtin_ctz(step);
      const int delta=(gray&(1U<<e))?-2:2;
      for(int s=0;s<states;++s) energy[s]+=delta*products[e][s];
    }
    int cap=0;
    for(int value:energy) cap=std::max(cap,std::abs(value));
    ++histogram[cap];
    if(cap<best) {best=cap;minimizers.clear();}
    if(cap==best) minimizers.push_back(gray);
  }
  std::sort(minimizers.begin(),minimizers.end());
  std::cout<<"{\"n\":"<<n<<",\"signings\":"<<total<<",\"M\":"<<best<<",\"histogram\":{";
  bool first=true;
  for(const auto& kv:histogram) {
    if(!first) std::cout<<",";
    first=false;
    std::cout<<"\""<<kv.first<<"\":"<<kv.second;
  }
  std::cout<<"},\"minimizers\":[";
  for(unsigned i=0;i<minimizers.size();++i) {if(i) std::cout<<",";std::cout<<minimizers[i];}
  std::cout<<"]}\n";
}
