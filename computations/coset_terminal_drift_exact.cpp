// Exact augmented-cut-code coset census for small graph orders.  The root is
// gauge-fixed at vertex n-1 and quotiented by global edge complement.  Output
// records distance layers, outward degree, and nondeep one-edge dead ends.
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <map>
#include <set>
#include <tuple>
#include <vector>

int main(int argc, char** argv) {
  int n = argc > 1 ? std::stoi(argv[1]) : 8;
  int m = (n - 1) * (n - 2) / 2;
  int N = n * (n - 1) / 2;
  uint32_t full = (uint32_t(1) << m) - 1;
  uint32_t count = uint32_t(1) << (m - 1);

  std::vector<std::pair<int, int>> edges;
  for (int i = 0; i < n - 1; ++i)
    for (int j = i + 1; j < n - 1; ++j) edges.push_back({i, j});

  std::vector<uint32_t> spin_masks;
  std::vector<int> star_terms;
  for (uint32_t s = 0; s < (uint32_t(1) << (n - 1)); ++s) {
    uint32_t p = 0;
    for (int k = 0; k < m; ++k) {
      auto [i, j] = edges[k];
      if (((s >> i) ^ (s >> j)) & 1) p |= uint32_t(1) << k;
    }
    spin_masks.push_back(p);
    star_terms.push_back((n - 1) - 2 * __builtin_popcount(s));
  }

  std::vector<uint8_t> radius(count);
  int rho = 0;
  auto canon = [full](uint32_t q) { return std::min(q, full ^ q); };
  for (uint32_t q = 0; q < count; ++q) {
    int best = 0;
    for (size_t s = 0; s < spin_masks.size(); ++s) {
      int internal = m - 2 * __builtin_popcount(q ^ spin_masks[s]);
      best = std::max(best, std::abs(internal + star_terms[s]));
    }
    int r = (N - best) / 2;
    radius[q] = r;
    rho = std::max(rho, r);
  }

  std::vector<uint32_t> gens;
  for (int k = 0; k < m; ++k) gens.push_back(uint32_t(1) << k);
  for (int i = 0; i < n - 1; ++i) {
    uint32_t g = 0;
    for (int k = 0; k < m; ++k)
      if (edges[k].first == i || edges[k].second == i) g |= uint32_t(1) << k;
    gens.push_back(g);
  }

  struct Summary {
    uint64_t count = 0, dead = 0;
    int min_b = 1000000, max_b = -1, min_same = 1000000, max_same = -1;
    std::set<int> bvals;
  };
  std::map<int, Summary> layers;
  for (uint32_t q = 0; q < count; ++q) {
    int r = radius[q], b = 0, same = 0;
    for (uint32_t g : gens) {
      int rr = radius[canon(q ^ g)];
      b += rr == r + 1;
      same += rr == r;
    }
    auto& x = layers[r];
    ++x.count;
    x.dead += (b == 0 && r < rho);
    x.min_b = std::min(x.min_b, b); x.max_b = std::max(x.max_b, b);
    x.min_same = std::min(x.min_same, same); x.max_same = std::max(x.max_same, same);
    x.bvals.insert(b);
  }

  std::cout << "n=" << n << " cosets=" << count << " N=" << N
            << " rho=" << rho << " M=" << N - 2 * rho << "\n";
  for (auto const& [r, x] : layers) {
    double z = double(N - 2 * r) / (n * std::sqrt(double(n)));
    std::cout << "layer=" << r << " Q=" << N - 2 * r << " z=" << z
              << " count=" << x.count << " b=" << x.min_b << ".." << x.max_b
              << " same=" << x.min_same << ".." << x.max_same
              << " nondeep_dead=" << x.dead << " bvals=";
    for (int b : x.bvals) std::cout << b << ',';
    std::cout << "\n";
  }
}
