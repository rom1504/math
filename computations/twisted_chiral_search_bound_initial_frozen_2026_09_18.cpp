#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <random>
#include <string>
#include <unordered_set>
#include <vector>

// Input: r, then hollow symmetric signing A. B_ij=s_i*s_j*A[p_i][p_j].
// Profiles use x_0=y_0=+1; the chiral action (x,y)->(y,-x)
// together with global negation meets every energy-absolute-value orbit.
// All energies and caps are integer. Fitness is only a heuristic tie break.
struct Eval { int cap=10000, d=0, top=10000, near=10000; double score=10000; };
int r, N; int A[12][12], B[12][12], spins[2048][12], qa[2048];
int lo[2048], hi[2048], order_t[2048], gray[2048], flipped[2048];
uint64_t calls=0;
std::array<uint64_t,1000> cap_hist{};
uint64_t rooted_B_count=0, symmetry_count=0, full_B_count=0;
uint64_t bound_calls=0,bound_pairs=0;
bool bound_mode=false;
int mask_words;
std::vector<uint64_t> distance_ge;

void prepare_distance_masks() {
  mask_words=(2*N+63)/64;distance_ge.assign(N*(r+2)*mask_words,0);
  for(int t=0;t<N;t++)for(int d=0;d<2*N;d++) {
    int dist=__builtin_popcount(unsigned(d^(t<<1)));
    for(int k=0;k<=dist;k++)distance_ge[(t*(r+2)+k)*mask_words+d/64]|=uint64_t(1)<<(d%64);
  }
}

// Exact feasibility test for all matching d at a fixed threshold, without
// computing unnecessary remaining spin energies after every d is excluded.
bool feasible_matching(const std::vector<int>& p,const std::vector<int>& s,int target) {
  ++bound_calls;
  for(int i=0;i<r;i++)for(int j=0;j<r;j++)B[i][j]=s[i]*s[j]*A[p[i]][p[j]];
  std::vector<uint64_t> valid(mask_words,~uint64_t(0));
  if(2*N<64)valid[0]=(uint64_t(1)<<(2*N))-1;
  std::fill(lo,lo+N,10000);std::fill(hi,hi+N,-10000);
  for(int x=0;x<N;x++) {
    int field[12]={},bilinear=0;
    for(int i=0;i<r;i++)for(int j=0;j<r;j++)field[j]+=spins[x][i]*B[i][j];
    for(int j=0;j<r;j++)bilinear+=field[j];
    for(int k=0;k<N;k++) {
      ++bound_pairs;int y=gray[k];if(k){int j=flipped[k];bilinear+=2*spins[y][j]*field[j];}
      int t=x^y,e=qa[x]-qa[y]+bilinear;
      if(e>target+r||e<-target-r)return false;
      bool upper=e>target-r&&e>hi[t],lower=e<-target+r&&e<lo[t];
      if(!upper&&!lower)continue;
      int cutoff;
      if(upper){hi[t]=e;cutoff=(r+e-target)/2;}
      else {lo[t]=e;cutoff=(r+target+e)/2+1;}
      const uint64_t* allowed=&distance_ge[(t*(r+2)+cutoff)*mask_words];
      uint64_t nonempty=0;
      for(int w=0;w<mask_words;w++){valid[w]&=upper?allowed[w]:~allowed[w];nonempty|=valid[w];}
      if(!nonempty)return false;
    }
  }
  return true;
}

struct Symmetry { std::vector<int> edge_map; uint64_t flip; };
struct Rooted { std::vector<int> p, root_s; };

uint64_t full_code(const std::vector<int>& p,const std::vector<int>& s) {
  uint64_t code=0;int bit=0;
  for(int i=0;i<r;i++)for(int j=i+1;j<r;j++,bit++)
    if(s[i]*s[j]*A[p[i]][p[j]]<0)code|=uint64_t(1)<<bit;
  return code;
}

// Canonical orbit representatives under signed automorphisms/antiautomorphisms
// of the fixed A. For an antiautomorphism also exchange the two parent halves.
// Return orbit size, or zero for a noncanonical B. This weights histograms
// back to the complete B orbit, independently of stabilizer size.
int canonical_weight(uint64_t code,const std::vector<Symmetry>& group) {
  int stabilizer=0, edges=r*(r-1)/2;
  for(const auto& h:group) {
    int cmp=0;
    for(int bit=edges-1;bit>=0;bit--) {
      int old=(code>>bit)&1;
      int changed=((code>>h.edge_map[bit])^(h.flip>>bit))&1;
      if(changed!=old){cmp=changed-old;break;}
    }
    if(cmp<0)return 0;
    if(cmp==0)stabilizer++;
  }
  if(!stabilizer||group.size()%stabilizer){std::cerr<<"bad stabilizer\n";std::exit(3);}
  return group.size()/stabilizer;
}

void profile(const std::vector<int>& p, const std::vector<int>& s) {
  for(int i=0;i<r;i++) for(int j=0;j<r;j++) B[i][j]=s[i]*s[j]*A[p[i]][p[j]];
  std::fill(lo,lo+N,10000); std::fill(hi,hi+N,-10000);
  for(int x=0;x<N;x++) {
    int field[12]={}; int bilinear=0;
    for(int i=0;i<r;i++) for(int j=0;j<r;j++) field[j]+=spins[x][i]*B[i][j];
    for(int j=0;j<r;j++) bilinear+=field[j];
    for(int k=0;k<N;k++) {
      int y=gray[k];
      if(k) { int j=flipped[k]; bilinear+=2*spins[y][j]*field[j]; }
      int t=x^y, e=qa[x]-qa[y]+bilinear;
      lo[t]=std::min(lo[t],e); hi[t]=std::max(hi[t],e);
    }
  }
}

Eval evaluate(const std::vector<int>& p, const std::vector<int>& s) {
  ++calls; profile(p,s);
  // Large unavoidable profile width first: cheap candidate-d pruning.
  std::iota(order_t,order_t+N,0);
  std::sort(order_t,order_t+N,[](int a,int b){
    return hi[a]-lo[a]>hi[b]-lo[b];
  });
  Eval best;
  for(int d=0;d<2*N;d++) {
    int cap=0, top=0, near=0;
    bool pruned=false;
    for(int k=0;k<N;k++) {
      int t=order_t[k], shift=r-2*__builtin_popcount(unsigned(d^(t<<1)));
      cap=std::max(cap,std::max(std::abs(lo[t]+shift),std::abs(hi[t]+shift)));
      if(cap>best.cap) {pruned=true;break;}
    }
    if(pruned) continue;
    for(int t=0;t<N;t++) {
      int shift=r-2*__builtin_popcount(unsigned(d^(t<<1)));
      int l=std::abs(lo[t]+shift),h=std::abs(hi[t]+shift);
      top+=(l==cap)+(h==cap); near+=(l>=cap-2)+(h>=cap-2);
    }
    double score=cap+0.08*std::log(1.0+top)+0.008*std::log(1.0+near);
    if(cap<best.cap||(cap==best.cap&&score<best.score)) best={cap,d,top,near,score};
  }
  return best;
}

void emit(const std::string& kind, const std::vector<int>& p,
          const std::vector<int>& s,const Eval& e,uint64_t count,double elapsed) {
  std::cout<<"{\"kind\":\""<<kind<<"\",\"r\":"<<r<<",\"calls\":"<<calls
   <<",\"count\":"<<count<<",\"elapsed\":"<<elapsed<<",\"cap\":"<<e.cap
   <<",\"score\":"<<e.score<<",\"top_profiles\":"<<e.top<<",\"near_profiles\":"<<e.near<<",\"p\":[";
  for(int i=0;i<r;i++)std::cout<<(i?",":"")<<p[i];
  std::cout<<"],\"s\":[";for(int i=0;i<r;i++)std::cout<<(i?",":"")<<s[i];
  std::cout<<"],\"d\":[";for(int i=0;i<r;i++)std::cout<<(i?",":"")<<((e.d>>i&1)?-1:1);
  std::cout<<"]";
  if(kind=="complete_family"){
    std::cout<<",\"rooted_B_count\":"<<rooted_B_count<<",\"signed_aut_antiaut_group_size\":"<<symmetry_count
      <<",\"full_B_count\":"<<full_B_count;
    std::cout<<",\"bound_calls\":"<<bound_calls<<",\"bound_pair_visits\":"<<bound_pairs;
    std::cout<<(bound_mode?",\"improving_witness_cap_histogram\":{":",\"best_matching_cap_histogram\":{");
    bool first=true;for(int k=0;k<1000;k++)if(cap_hist[k]){
      std::cout<<(first?"":",")<<"\""<<k<<"\":"<<cap_hist[k];first=false;
    }
    std::cout<<"}";
  }
  std::cout<<"}"<<std::endl;
}

int main(int argc,char**argv) {
  if(argc<2){std::cerr<<"mode exhaustive|search [seconds] [seed]\n";return 2;}
  std::string mode=argv[1]; double seconds=argc>2?std::stod(argv[2]):60;
  uint64_t seed=argc>3?std::stoull(argv[3]):20260918;
  if(!(std::cin>>r)||r<2||r>11)return 2; N=1<<(r-1);
  for(int i=0;i<r;i++)for(int j=0;j<r;j++)std::cin>>A[i][j];
  for(int x=0;x<N;x++) {
    spins[x][0]=1; for(int i=1;i<r;i++)spins[x][i]=(x>>(i-1)&1)?-1:1;
    qa[x]=0;for(int i=0;i<r;i++)for(int j=i+1;j<r;j++)qa[x]+=A[i][j]*spins[x][i]*spins[x][j];
    gray[x]=x^(x>>1);if(x)flipped[x]=1+__builtin_ctz(unsigned(x));
  }
  auto start=std::chrono::steady_clock::now();
  auto elapsed=[&](){return std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();};
  std::vector<int> p(r),s(r,1),bp,bs;std::iota(p.begin(),p.end(),0);
  Eval best; uint64_t count=0;
  auto consider=[&](const Eval& e){if(e.cap<best.cap||(e.cap==best.cap&&e.score<best.score)){best=e;bp=p;bs=s;emit("best",bp,bs,best,count,elapsed());}};
  if(mode=="quotient"||mode=="bound") {
    bound_mode=mode=="bound";
    if(bound_mode){prepare_distance_masks();auto e=evaluate(p,s);consider(e);}
    int edge_index[12][12],bit=0;
    for(int i=0;i<r;i++)for(int j=i+1;j<r;j++,bit++)edge_index[i][j]=edge_index[j][i]=bit;
    uint64_t rootmask=0;for(int i=1;i<r;i++)for(int j=i+1;j<r;j++)rootmask|=uint64_t(1)<<edge_index[i][j];
    std::vector<int> orig_s(r,1);for(int i=1;i<r;i++)orig_s[i]=A[0][i];
    uint64_t rootbase=full_code(p,orig_s);
    std::unordered_set<uint64_t> roots_seen;std::vector<Rooted> roots;std::vector<Symmetry> group;
    do {
      std::vector<int> rs(r,1);for(int i=1;i<r;i++)rs[i]=A[p[0]][p[i]];
      uint64_t code=full_code(p,rs);
      if(roots_seen.insert(code).second)roots.push_back({p,rs});
      int polarity=(code==rootbase)?1:((code==(rootbase^rootmask))?-1:0);
      if(polarity) {
        std::vector<int> hs(r,1);for(int i=1;i<r;i++)hs[i]=polarity*A[0][i]*A[p[0]][p[i]];
        Symmetry h;h.flip=0;
        for(int i=0;i<r;i++)for(int j=i+1;j<r;j++) {
          h.edge_map.push_back(edge_index[p[i]][p[j]]);
          if(hs[i]*hs[j]<0)h.flip|=uint64_t(1)<<edge_index[i][j];
          if(hs[i]*hs[j]*A[p[i]][p[j]]!=polarity*A[i][j])return 4;
        }
        group.push_back(h);
      }
    }while(std::next_permutation(p.begin(),p.end()));
    rooted_B_count=roots.size();symmetry_count=group.size();full_B_count=rooted_B_count*N;
    std::cout<<"{\"kind\":\"orbit_inventory\",\"r\":"<<r<<",\"rooted_B_count\":"<<rooted_B_count
      <<",\"group_size\":"<<symmetry_count<<",\"full_B_count\":"<<full_B_count<<",\"elapsed\":"<<elapsed()<<"}"<<std::endl;
    uint64_t weight_sum=0;
    for(const auto& root:roots) {
      p=root.p;
      for(int mask=0;mask<N;mask++) {
        s=root.root_s;for(int i=1;i<r;i++)if(mask>>(i-1)&1)s[i]=-s[i];
        uint64_t code=full_code(p,s);int weight=canonical_weight(code,group);if(!weight)continue;
        count++;weight_sum+=weight;
        if(bound_mode&&!feasible_matching(p,s,best.cap-2))continue;
        auto e=evaluate(p,s);cap_hist[e.cap]+=weight;consider(e);
      }
    }
    if(weight_sum!=full_B_count){std::cerr<<"orbit weight mismatch "<<weight_sum<<" != "<<full_B_count<<"\n";return 5;}
    emit("complete_family",bp,bs,best,count,elapsed());
  } else if(mode=="exhaustive") {
    if(r>8)return 2;
    std::unordered_set<uint64_t> seen;
    do {
      for(int mask=0;mask<N;mask++) {
        for(int i=1;i<r;i++)s[i]=(mask>>(i-1)&1)?-1:1;
        uint64_t code=0;int bit=0;
        for(int i=0;i<r;i++)for(int j=i+1;j<r;j++,bit++)
          if(s[i]*s[j]*A[p[i]][p[j]]<0)code|=uint64_t(1)<<bit;
        if(!seen.insert(code).second)continue;
        count++;auto e=evaluate(p,s);cap_hist[e.cap]++;consider(e);
      }
    } while(std::next_permutation(p.begin(),p.end()));
    full_B_count=count;emit("complete_family",bp,bs,best,count,elapsed());
  } else {
    std::mt19937_64 rng(seed);std::uniform_real_distribution<double> unif(0,1);
    Eval current=evaluate(p,s);consider(current);uint64_t stall=0;
    while(elapsed()<seconds) {
      auto oldp=p,olds=s;
      if(stall>3000) {
        std::shuffle(p.begin(),p.end(),rng);for(int i=1;i<r;i++)s[i]=(rng()&1)?-1:1;stall=0;
      } else if(rng()%3) {int i=rng()%r,j=rng()%r;while(j==i)j=rng()%r;std::swap(p[i],p[j]);}
      else {int i=1+rng()%(r-1);s[i]=-s[i];}
      count++;auto e=evaluate(p,s);consider(e);
      double temp=0.04+0.45*std::pow(1.0-double(count%2500)/2500,2);
      if(e.score<=current.score||unif(rng)<std::exp((current.score-e.score)/temp)) {current=e;}
      else {p=oldp;s=olds;}
      if(e.cap>best.cap)stall++;else stall=0;
      if(count%2500==0){p=bp;s=bs;current=best; if((count/2500)%4==0){std::shuffle(p.begin(),p.end(),rng);current=evaluate(p,s);}}
    }
    emit("search_finished",bp,bs,best,count,elapsed());
  }
}
