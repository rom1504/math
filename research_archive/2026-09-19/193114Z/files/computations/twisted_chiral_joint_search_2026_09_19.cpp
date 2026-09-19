// Bounded heuristic over A,p,s, with EXACT optimization over every matching d.
// Reuses the audited integer profile evaluator; no heuristic result is a lower bound.
#define main frozen_profile_program_main
#include "twisted_chiral_search_frozen_2026_09_18.cpp"
#undef main

struct JointState {
  int a[12][12]{};
  std::vector<int> p,s;
  Eval e;
};

void refresh_qa() {
  for(int x=0;x<N;x++) {
    qa[x]=0;
    for(int i=0;i<r;i++)for(int j=i+1;j<r;j++)qa[x]+=A[i][j]*spins[x][i]*spins[x][j];
  }
}
int child_cap() { int q=0;for(int x=0;x<N;x++)q=std::max(q,std::abs(qa[x]));return q; }
JointState save(const std::vector<int>&p,const std::vector<int>&s,Eval e) {
  JointState z;std::copy(&A[0][0],&A[0][0]+144,&z.a[0][0]);z.p=p;z.s=s;z.e=e;return z;
}
void restore(const JointState&z,std::vector<int>&p,std::vector<int>&s) {
  std::copy(&z.a[0][0],&z.a[0][0]+144,&A[0][0]);p=z.p;s=z.s;refresh_qa();
}
void joint_emit(const std::string&kind,const JointState&z,uint64_t proposals,
                uint64_t accepts,uint64_t rejected_child,int restarts,double elapsed) {
  std::cout<<"{\"kind\":\""<<kind<<"\",\"r\":"<<r<<",\"calls\":"<<calls
    <<",\"proposals\":"<<proposals<<",\"accepted\":"<<accepts<<",\"rejected_child_cap\":"<<rejected_child
    <<",\"restarts\":"<<restarts<<",\"elapsed\":"<<elapsed<<",\"cap\":"<<z.e.cap
    <<",\"score\":"<<z.e.score<<",\"top_profiles\":"<<z.e.top<<",\"near_profiles\":"<<z.e.near<<",\"p\":[";
  for(int i=0;i<r;i++)std::cout<<(i?",":"")<<z.p[i];
  std::cout<<"],\"s\":[";for(int i=0;i<r;i++)std::cout<<(i?",":"")<<z.s[i];
  std::cout<<"],\"d\":[";for(int i=0;i<r;i++)std::cout<<(i?",":"")<<((z.e.d>>i&1)?-1:1);
  std::cout<<"],\"searched_child_matrix\":[";
  for(int i=0;i<r;i++){std::cout<<(i?",":"")<<"[";for(int j=0;j<r;j++)std::cout<<(j?",":"")<<z.a[i][j];std::cout<<"]";}
  std::cout<<"]}"<<std::endl;
}

int main(int argc,char**argv) {
  const double seconds=argc>1?std::stod(argv[1]):3600;
  const uint64_t seed=argc>2?std::stoull(argv[2]):20260919;
  const int child_limit=argc>3?std::stoi(argv[3]):17;
  const int target=argc>4?std::stoi(argv[4]):38;
  if(!(std::cin>>r)||r!=10)return 2;N=1<<(r-1);
  for(int i=0;i<r;i++)for(int j=0;j<r;j++)std::cin>>A[i][j];
  for(int x=0;x<N;x++) {
    spins[x][0]=1;for(int i=1;i<r;i++)spins[x][i]=(x>>(i-1)&1)?-1:1;
    gray[x]=x^(x>>1);if(x)flipped[x]=1+__builtin_ctz(unsigned(x));
  }
  refresh_qa();if(child_cap()>child_limit)return 3;
  std::vector<int> p(r),s(r,1);std::iota(p.begin(),p.end(),0);
  std::vector<std::pair<int,int>> edges;
  for(int i=0;i<r;i++)for(int j=i+1;j<r;j++)edges.emplace_back(i,j);
  auto start=std::chrono::steady_clock::now();
  auto elapsed=[&](){return std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();};
  std::mt19937_64 rng(seed);std::uniform_real_distribution<double> uniform(0,1);
  Eval current=evaluate(p,s);JointState initial=save(p,s,current),best=initial;
  std::vector<JointState> pool{initial};
  uint64_t proposals=0,accepts=0,rejected_child=0;int restarts=0;
  std::array<uint64_t,100> evaluated_caps{},child_caps{};
  joint_emit("best",best,proposals,accepts,rejected_child,restarts,elapsed());
  double next_progress=60;
  const int cycle=6000;
  while(elapsed()<seconds&&best.e.cap>target) {
    // Alternate a cold local phase and a hot phase. Keep several low-cap seeds
    // rather than presuming the conference basin represents every good child.
    if(proposals&&proposals%cycle==0) {
      restarts++;
      const auto&restart_seed=(restarts%4==0)?initial:pool[rng()%pool.size()];
      restore(restart_seed,p,s);current=restart_seed.e;
    }
    JointState before=save(p,s,current);
    int move=rng()%100;
    if(move<72) {
      int changes=(rng()%5==0)?2:1;
      for(int c=0;c<changes;c++) {
        auto [i,j]=edges[rng()%edges.size()];int old=A[i][j];A[i][j]=A[j][i]=-old;
        for(int x=0;x<N;x++)qa[x]-=2*old*spins[x][i]*spins[x][j];
      }
    } else if(move<90) {
      int i=rng()%r,j=rng()%(r-1);if(j>=i)j++;std::swap(p[i],p[j]);
    } else {int i=1+rng()%(r-1);s[i]=-s[i];}
    proposals++;
    int cq=child_cap();
    if(cq>child_limit){rejected_child++;restore(before,p,s);continue;}
    Eval candidate=evaluate(p,s);evaluated_caps[candidate.cap]++;child_caps[cq]++;
    double phase=double(proposals%cycle)/cycle;
    double hot=(restarts%3==0)?2.6:((restarts%3==1)?1.0:0.32);
    double temperature=hot*std::pow(0.025/hot,phase);
    if(candidate.score<=current.score||uniform(rng)<std::exp((current.score-candidate.score)/temperature)) {
      current=candidate;accepts++;
      if(candidate.cap<best.e.cap||(candidate.cap==best.e.cap&&candidate.score<best.e.score)) {
        best=save(p,s,candidate);joint_emit("best",best,proposals,accepts,rejected_child,restarts,elapsed());
      }
      // Seed diversity is heuristic only. No coverage claim relies on this pool.
      if(candidate.cap<=42&&rng()%128==0) {
        if(pool.size()<64)pool.push_back(save(p,s,candidate));
        else pool[1+rng()%(pool.size()-1)]=save(p,s,candidate);
      }
    } else restore(before,p,s);
    if(elapsed()>=next_progress) {
      std::cout<<"{\"kind\":\"progress\",\"elapsed\":"<<elapsed()<<",\"proposals\":"<<proposals
       <<",\"calls\":"<<calls<<",\"accepted\":"<<accepts<<",\"best_cap\":"<<best.e.cap
       <<",\"current_cap\":"<<current.cap<<",\"restarts\":"<<restarts<<",\"pool_size\":"<<pool.size()<<"}"<<std::endl;
      next_progress+=60;
    }
  }
  joint_emit("joint_search_finished",best,proposals,accepts,rejected_child,restarts,elapsed());
  std::cout<<"{\"kind\":\"search_statistics\",\"target\":"<<target<<",\"child_limit\":"<<child_limit
    <<",\"random_seed\":"<<seed<<",\"seconds_budget\":"<<seconds<<",\"pool_size\":"<<pool.size()<<",\"evaluated_cap_histogram\":{";
  bool first=true;for(int k=0;k<100;k++)if(evaluated_caps[k]){std::cout<<(first?"":",")<<"\""<<k<<"\":"<<evaluated_caps[k];first=false;}
  std::cout<<"},\"evaluated_child_cap_histogram\":{";first=true;
  for(int k=0;k<100;k++)if(child_caps[k]){std::cout<<(first?"":",")<<"\""<<k<<"\":"<<child_caps[k];first=false;}
  std::cout<<"},\"evidence\":\"bounded heuristic only; failure is not an exclusion certificate\"}"<<std::endl;
}
