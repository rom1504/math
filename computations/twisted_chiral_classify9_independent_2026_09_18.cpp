#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <unordered_map>
#include <unordered_set>
#include <vector>

// Independently exhaust all 2^28 root-gauged order-nine signings. Independent
// half-edge signed-sum tables evaluate q_A(x) directly, with no solver and no parent machinery.
// Input: initial representative count followed by their 9x9 matrices.
constexpr int n=9,m=28,total_spins=256,target=12;
using Matrix=std::vector<std::vector<int>>;
std::unordered_map<uint32_t,int> orbit;
std::vector<Matrix> representatives;
std::vector<uint64_t> counts,orbit_sizes;
const uint32_t all_bits=(uint32_t(1)<<m)-1;

void add_class(const Matrix& a) {
  int cls=representatives.size();representatives.push_back(a);counts.push_back(0);
  std::vector<int> p(n);std::iota(p.begin(),p.end(),0);
  std::unordered_set<uint32_t> this_orbit;
  do {
    uint32_t code=0;int bit=0;
    for(int i=1;i<n;i++)for(int j=i+1;j<n;j++,bit++)
      if(a[p[0]][p[i]]*a[p[0]][p[j]]*a[p[i]][p[j]]<0)code|=uint32_t(1)<<bit;
    this_orbit.insert(code);this_orbit.insert(code^all_bits);
  }while(std::next_permutation(p.begin(),p.end()));
  for(auto code:this_orbit)if(!orbit.emplace(code,cls).second){std::cerr<<"overlapping input classes\n";std::exit(3);}
  orbit_sizes.push_back(this_orbit.size());
  std::cout<<"{\"kind\":\"class\",\"class\":"<<cls<<",\"root_orbit_size\":"<<this_orbit.size()<<",\"matrix\":[";
  for(int i=0;i<n;i++){std::cout<<(i?",":"")<<"[";for(int j=0;j<n;j++)std::cout<<(j?",":"")<<a[i][j];std::cout<<"]";}
  std::cout<<"]}"<<std::endl;
}

int main() {
  auto start=std::chrono::steady_clock::now();
  auto elapsed=[&](){return std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();};
  int initial;std::cin>>initial;
  for(int k=0;k<initial;k++){Matrix a(n,std::vector<int>(n));for(auto& row:a)for(auto& value:row)std::cin>>value;add_class(a);}
  int order[total_spins],root[total_spins],products[m][total_spins];
  uint64_t rejects[total_spins]={};
  for(int code=0;code<total_spins;code++){
    int x[n];x[0]=1;root[code]=0;for(int i=1;i<n;i++){x[i]=(code>>(i-1)&1)?-1:1;root[code]+=x[i];}
    int bit=0;for(int i=1;i<n;i++)for(int j=i+1;j<n;j++,bit++)products[bit][code]=x[i]*x[j];
    order[code]=code;
  }
  constexpr int half=14,side_count=1<<half;
  std::vector<int16_t> left(side_count*total_spins),right(side_count*total_spins);
  for(int side=0;side<2;side++){
    auto& table=side?right:left;
    for(int spin=0;spin<total_spins;spin++)for(int bit=0;bit<half;bit++)table[spin]+=products[side*half+bit][spin];
    int previous=0;
    for(int step=1;step<side_count;step++){
      int gray=step^(step>>1),bit=__builtin_ctz(unsigned(step)),new_sign=(gray>>bit&1)?-1:1;
      for(int spin=0;spin<total_spins;spin++)table[gray*total_spins+spin]=table[previous*total_spins+spin]+2*new_sign*products[side*half+bit][spin];
      previous=gray;
    }
  }
  uint64_t feasible=0;int minimum=10000;uint64_t energy_checks=0;
  for(uint32_t mask=0;mask<=all_bits;mask++) {
    int cap=0;bool pass=true;const int16_t* l=&left[(mask&(side_count-1))*total_spins];const int16_t* rr=&right[(mask>>half)*total_spins];
    for(int k=0;k<total_spins;k++){
      int spin=order[k];int e=root[spin]+l[spin]+rr[spin];energy_checks++;
      cap=std::max(cap,std::abs(e));if(cap>target){rejects[spin]++;pass=false;break;}
    }
    if(pass){
      feasible++;minimum=std::min(minimum,cap);auto it=orbit.find(mask);
      if(it==orbit.end()){
        Matrix a(n,std::vector<int>(n,1));for(int i=0;i<n;i++)a[i][i]=0;
        int bit=0;for(int i=1;i<n;i++)for(int j=i+1;j<n;j++,bit++)if(mask>>bit&1)a[i][j]=a[j][i]=-1;
        add_class(a);it=orbit.find(mask);
      }
      counts[it->second]++;
    }
    if((mask&((1<<20)-1))==((1<<20)-1)){
      std::sort(order,order+total_spins,[&](int a,int b){return rejects[a]>rejects[b];});
      std::fill(rejects,rejects+total_spins,0);
      if((mask&((1<<24)-1))==((1<<24)-1))std::cout<<"{\"kind\":\"progress\",\"checked\":"<<uint64_t(mask)+1<<",\"feasible\":"<<feasible<<",\"classes\":"<<representatives.size()<<",\"elapsed\":"<<elapsed()<<"}"<<std::endl;
    }
  }
  for(size_t i=0;i<counts.size();i++)if(counts[i]!=orbit_sizes[i]){std::cerr<<"class orbit/count disagreement\n";return 4;}
  std::cout<<"{\"kind\":\"complete\",\"order\":9,\"root_signings_checked\":268435456,\"projective_spins\":256,\"energy_checks\":"<<energy_checks<<",\"global_minimum\":"<<minimum<<",\"minimizer_count\":"<<feasible<<",\"class_counts\":[";
  for(size_t i=0;i<counts.size();i++)std::cout<<(i?",":"")<<counts[i];
  std::cout<<"],\"elapsed\":"<<elapsed()<<"}"<<std::endl;
}
