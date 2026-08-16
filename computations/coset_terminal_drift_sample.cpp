// Deterministic-seed steepest one-edge descent for terminal coset traps.
// This is exploratory unless a reported signing is then checked exhaustively.
#include <algorithm>
#include <cstdint>
#include <iostream>
#include <map>
#include <random>
#include <set>
#include <utility>
#include <vector>

struct Eval { int q; int b; uint64_t outward; };

int main(int argc, char** argv) {
  int n = argc > 1 ? std::stoi(argv[1]) : 10;
  int trials = argc > 2 ? std::stoi(argv[2]) : 100000;
  int N = n * (n - 1) / 2;
  if (N >= 64) return 2;
  uint64_t full = (uint64_t(1) << N) - 1;
  std::vector<std::pair<int,int>> edges;
  for (int i=0;i<n;++i) for (int j=i+1;j<n;++j) edges.push_back({i,j});
  std::vector<uint64_t> cuts;
  for (uint64_t s=0;s<(uint64_t(1)<<(n-1));++s) {
    uint64_t p=0;
    for (int e=0;e<N;++e) {
      auto [i,j]=edges[e];
      int xi = i==0 ? 0 : ((s>>(i-1))&1);
      int xj = j==0 ? 0 : ((s>>(j-1))&1);
      if (xi^xj) p |= uint64_t(1)<<e;
    }
    cuts.push_back(p);
  }
  auto eval = [&](uint64_t a) {
    int Q=0;
    for (auto p:cuts) Q=std::max(Q,std::abs(N-2*__builtin_popcountll(a^p)));
    uint64_t covered=0;
    for (auto p:cuts) {
      int h=N-2*__builtin_popcountll(a^p);
      if (h>=Q-2) covered |= a^p;
      if (-h>=Q-2) covered |= full^(a^p);
    }
    uint64_t outward=full^covered;
    return Eval{Q, __builtin_popcountll(outward), outward};
  };
  std::mt19937_64 rng(0x6c6472696674ULL+n);
  std::map<int,uint64_t> terminal;
  std::map<int,uint64_t> sampled;
  std::map<int,uint64_t> first_terminal;
  std::map<int,std::set<int>> bvals;
  for (int t=0;t<trials;++t) {
    uint64_t a=rng()&full;
    while (true) {
      Eval x=eval(a);
      ++sampled[x.q];
      bvals[x.q].insert(x.b);
      if (!x.outward) {
        ++terminal[x.q];
        if (!first_terminal.count(x.q)) first_terminal[x.q]=a;
        break;
      }
      int pick=std::uniform_int_distribution<int>(0,x.b-1)(rng);
      uint64_t z=x.outward;
      while (pick--) z&=z-1;
      a ^= z&-z;
    }
  }
  std::cout << "n="<<n<<" trials="<<trials<<"\nlocal_min_Q:";
  for(auto [q,c]:terminal) std::cout<<' '<<q<<':'<<c<<"@0x"<<std::hex<<first_terminal[q]<<std::dec;
  std::cout<<"\nsampled_Q:";
  for(auto [q,c]:sampled) std::cout<<' '<<q<<':'<<c;
  std::cout<<"\nb_values:";
  for(auto const& [q,bs]:bvals) {
    std::cout<<' '<<q<<'=';
    for(int b:bs) std::cout<<b<<',';
  }
  std::cout<<'\n';
}
