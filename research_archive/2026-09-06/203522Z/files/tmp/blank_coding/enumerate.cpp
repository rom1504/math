#include <algorithm>
#include <bit>
#include <cstdint>
#include <iostream>
#include <vector>

int main(int argc, char **argv) {
  int n = argc > 1 ? std::stoi(argv[1]) : 8;
  int m = (n - 1) * (n - 2) / 2;
  std::vector<uint64_t> smask;
  std::vector<int> base;
  for (uint64_t xmask = 0; xmask < (uint64_t(1) << (n - 1)); ++xmask) {
    uint64_t s = 0;
    int k = 0;
    int q = 0;
    for (int j = 1; j < n; ++j) q += (xmask >> (j - 1) & 1) ? -1 : 1;
    for (int i = 1; i < n; ++i) for (int j = i + 1; j < n; ++j, ++k) {
      int si = (xmask >> (i - 1) & 1) ? -1 : 1;
      int sj = (xmask >> (j - 1) & 1) ? -1 : 1;
      int p = si * sj;
      q += p;
      if (p < 0) s |= uint64_t(1) << k;
    }
    smask.push_back(s);
    base.push_back(q);
  }
  int best = 1 << 30;
  uint64_t count = 0, bestmask = 0;
  for (uint64_t a = 0; a < (uint64_t(1) << m); ++a) {
    int aw = std::popcount(a), Q = 0;
    for (size_t x = 0; x < smask.size(); ++x) {
      int q = base[x] - 2 * aw + 4 * std::popcount(a & smask[x]);
      Q = std::max(Q, std::abs(q));
      if (Q > best) break;
    }
    if (Q < best) { best = Q; count = 1; bestmask = a; }
    else if (Q == best) ++count;
  }
  std::cout << "n=" << n << " M=" << best << " labeled-gauge-count=" << count
            << " example-mask=" << bestmask << '\n';
}
