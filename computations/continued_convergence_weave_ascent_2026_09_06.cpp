// Reproducible lower witnesses only. Input: N runs sweeps seed, symmetric
// hollow matrix, initial-count, and initial Boolean vectors. No cap claim.
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <random>
#include <vector>
int main() {
  int n,runs,sweeps; unsigned seed;
  if(!(std::cin>>n>>runs>>sweeps>>seed)) return 2;
  std::vector<int8_t>a(n*n); int z;
  for(auto &v:a){std::cin>>z;v=z;}
  int ni; std::cin>>ni;
  std::vector<std::vector<int>> init(ni,std::vector<int>(n));
  for(auto &v:init)for(auto &w:v)std::cin>>w;
  std::mt19937 gen(seed); std::uniform_real_distribution<double> u(0,1);
  std::vector<int>x(n),h(n),bestx; int best=0,bestsign=1;
  for(int run=0;run<runs;++run){
    int sg=1-2*(run%2);
    if(run<2*ni)x=init[run/2];else for(auto &v:x)v=2*(gen()%2)-1;
    int e=0; for(int i=0;i<n;++i){h[i]=0;for(int j=0;j<n;++j)h[i]+=sg*a[i*n+j]*x[j];e+=x[i]*h[i];} e/=2;
    auto record=[&](){if(e>best){best=e;bestx=x;bestsign=sg;}};
    record();
    if(run<2*ni){bool change=true;while(change){change=false;for(int i=0;i<n;++i)if(x[i]*h[i]<0){int old=x[i];e-=2*old*h[i];x[i]=-old;for(int j=0;j<n;++j)h[j]-=2*old*sg*a[j*n+i];change=true;}}record();}
    for(int sw=0;sw<sweeps;++sw){
      double temp=1.5*std::sqrt(n)*std::pow(.002,double(sw)/std::max(1,sweeps-1));
      for(int p=0;p<n;++p){int i=gen()%n,old=x[i],de=-2*old*h[i];
        if(de>=0||u(gen)<std::exp(de/temp)){e+=de;x[i]=-old;for(int j=0;j<n;++j)h[j]-=2*old*sg*a[j*n+i];record();}
      }
    }
    bool change=true; while(change){change=false;for(int i=0;i<n;++i)if(x[i]*h[i]<0){int old=x[i];e-=2*old*h[i];x[i]=-old;for(int j=0;j<n;++j)h[j]-=2*old*sg*a[j*n+i];change=true;}}
    if(e>best){best=e;bestx=x;bestsign=sg;}
  }
  long long exact=0;for(int i=0;i<n;++i)for(int j=0;j<n;++j)exact+=bestx[i]*int(a[i*n+j])*bestx[j];
  if(exact!=2LL*best*bestsign)return 3;
  std::cout<<"{\"lower_Q\":"<<best<<",\"orientation\":"<<bestsign<<",\"witness\":[";
  for(int i=0;i<n;++i)std::cout<<(i?",":"")<<bestx[i];
  std::cout<<"]}\n";
}
