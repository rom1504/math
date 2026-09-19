// Exact p=id subfamily; reused audited integer profiles optimize every d.
#define main frozen_profile_program_main
#include "twisted_chiral_search_frozen_2026_09_18.cpp"
#undef main

int main() {
  if(!(std::cin>>r)||r<3||r>10)return 2;N=1<<(r-1);
  for(int i=0;i<r;i++)for(int j=0;j<r;j++)std::cin>>A[i][j];
  for(int x=0;x<N;x++) {
    spins[x][0]=1;for(int i=1;i<r;i++)spins[x][i]=(x>>(i-1)&1)?-1:1;
    qa[x]=0;for(int i=0;i<r;i++)for(int j=i+1;j<r;j++)qa[x]+=A[i][j]*spins[x][i]*spins[x][j];
    gray[x]=x^(x>>1);if(x)flipped[x]=1+__builtin_ctz(unsigned(x));
  }
  auto start=std::chrono::steady_clock::now();
  auto elapsed=[&](){return std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();};
  std::vector<int>p(r),s(r,1),best_s;std::iota(p.begin(),p.end(),0);Eval best;
  std::array<uint64_t,100> histogram{};
  for(int mask=0;mask<N;mask++) {
    for(int i=1;i<r;i++)s[i]=(mask>>(i-1)&1)?-1:1;
    Eval e=evaluate(p,s);histogram[e.cap]++;
    if(e.cap<best.cap||(e.cap==best.cap&&e.score<best.score)) {
      best=e;best_s=s;emit("best",p,s,e,mask+1,elapsed());
    }
  }
  emit("complete_switching_subfamily",p,best_s,best,N,elapsed());
  std::cout<<"{\"kind\":\"switching_histogram\",\"switching_count\":"<<N<<",\"matching_count_per_switch\":"<<2*N<<",\"minimum_cap_histogram\":{";
  bool first=true;for(int k=0;k<100;k++)if(histogram[k]){std::cout<<(first?"":",")<<"\""<<k<<"\":"<<histogram[k];first=false;}
  std::cout<<"}}"<<std::endl;
}
