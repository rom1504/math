#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <map>
#include <unordered_set>
#include <utility>
#include <vector>

// Exhaustive switching-class audit for signed K_n.  Gauge-fix the last
// vertex so all incident edges are +1; the remaining C(n-1,2) signs give
// one representative of every switching class.
int main(int argc, char **argv) {
  const int n = argc > 1 ? std::atoi(argv[1]) : 8;
  const int m = n - 1;
  const int states = 1 << m;
  std::vector<int> eu, ev;
  for (int i = 0; i < m; ++i)
    for (int j = i + 1; j < m; ++j) {
      eu.push_back(i);
      ev.push_back(j);
    }
  const int E = static_cast<int>(eu.size());
  if (E >= 63) {
    std::cerr << "too many free edges\n";
    return 2;
  }
  const uint64_t classes = uint64_t{1} << E;

  // q[x] includes the + star to the gauge vertex and all internal edges.
  std::vector<int16_t> q(states);
  std::vector<std::vector<int8_t>> chi(E, std::vector<int8_t>(states));
  for (int x = 0; x < states; ++x) {
    int star = 0;
    for (int i = 0; i < m; ++i)
      star += ((x >> i) & 1) ? -1 : 1;
    int val = star;
    for (int e = 0; e < E; ++e) {
      const int si = ((x >> eu[e]) & 1) ? -1 : 1;
      const int sj = ((x >> ev[e]) & 1) ? -1 : 1;
      chi[e][x] = static_cast<int8_t>(si * sj);
      val += si * sj;
    }
    q[x] = static_cast<int16_t>(val);
  }

  int best = 1 << 30;
  int best_sum = 1 << 30;
  uint64_t optimal = 0, balanced = 0;
  std::vector<uint64_t> objective_masks;
  uint64_t range_optimal = 0;
  std::vector<uint64_t> range_masks;
  std::map<std::pair<int, int>, uint64_t> profiles;
  std::map<std::pair<int, int>, uint64_t> range_profiles;
  uint64_t prev_gray = 0;
  std::vector<int8_t> sign(E, 1);

  for (uint64_t step = 0; step < classes; ++step) {
    const uint64_t gray = step ^ (step >> 1);
    if (step) {
      const uint64_t diff = gray ^ prev_gray;
      const int e = __builtin_ctzll(diff);
      const int delta = -2 * sign[e];
      for (int x = 0; x < states; ++x)
        q[x] = static_cast<int16_t>(q[x] + delta * chi[e][x]);
      sign[e] = static_cast<int8_t>(-sign[e]);
    }
    prev_gray = gray;

    int p = q[0], mn = q[0];
    for (int x = 1; x < states; ++x) {
      p = std::max(p, static_cast<int>(q[x]));
      mn = std::min(mn, static_cast<int>(q[x]));
    }
    const int qq = -mn;
    const int obj = std::max(p, qq);
    const int endpoint_sum = p + qq;
    if (endpoint_sum < best_sum) {
      best_sum = endpoint_sum;
      range_optimal = 0;
      range_profiles.clear();
      range_masks.clear();
    }
    if (endpoint_sum == best_sum) {
      ++range_optimal;
      ++range_profiles[{p, qq}];
      range_masks.push_back(gray);
    }
    if (obj < best) {
      best = obj;
      optimal = balanced = 0;
      profiles.clear();
      objective_masks.clear();
    }
    if (obj == best) {
      ++optimal;
      ++profiles[{p, qq}];
      objective_masks.push_back(gray);
      if (std::abs(p - qq) <= 2) ++balanced;
    }
  }

  std::cout << "n=" << n << " classes=" << classes << " F=" << best
            << " optimal=" << optimal << " balanced=" << balanced << "\n";
  for (const auto &[pq, count] : profiles)
    std::cout << "profile P=" << pq.first << " Q=" << pq.second
              << " count=" << count << "\n";
  std::cout << "minimum endpoint sum=" << best_sum
            << " range-optimal=" << range_optimal << "\n";
  for (const auto &[pq, count] : range_profiles)
    std::cout << "range profile P=" << pq.first << " Q=" << pq.second
              << " count=" << count << "\n";
  std::unordered_set<uint64_t> rem(range_masks.begin(), range_masks.end());
  int components = 0, antipodal_components = 0, largest = 0;
  while (!rem.empty()) {
    ++components;
    std::vector<uint64_t> stack{*rem.begin()};
    rem.erase(stack.back());
    int sz = 0;
    bool has_antipode = false;
    while (!stack.empty()) {
      uint64_t u = stack.back();
      stack.pop_back();
      ++sz;
      if (std::find(range_masks.begin(), range_masks.end(),
                    u ^ (classes - 1)) != range_masks.end())
        has_antipode = true;
      for (int e = 0; e < E; ++e) {
        auto it = rem.find(u ^ (uint64_t{1} << e));
        if (it != rem.end()) {
          stack.push_back(*it);
          rem.erase(it);
        }
      }
    }
    largest = std::max(largest, sz);
    if (has_antipode) ++antipodal_components;
  }
  std::cout << "range-minimizer flip graph components=" << components
            << " largest=" << largest
            << " antipodal-components=" << antipodal_components << "\n";

  // Connectivity of the deepest-hole plateau in the true quotient Cayley
  // graph.  Besides toggling an internal edge, toggling a gauge-star edge
  // toggles the full star at that reduced vertex.
  std::vector<uint64_t> generators;
  for (int e = 0; e < E; ++e) generators.push_back(uint64_t{1} << e);
  for (int i = 0; i < m; ++i) {
    uint64_t star = 0;
    for (int e = 0; e < E; ++e)
      if (eu[e] == i || ev[e] == i) star ^= uint64_t{1} << e;
    generators.push_back(star);
  }
  std::unordered_set<uint64_t> deep_rem(objective_masks.begin(),
                                         objective_masks.end());
  int deep_components = 0, deep_largest = 0, deep_self_antipodal = 0;
  while (!deep_rem.empty()) {
    ++deep_components;
    std::vector<uint64_t> stack{*deep_rem.begin()}, comp;
    deep_rem.erase(stack.back());
    while (!stack.empty()) {
      uint64_t u = stack.back();
      stack.pop_back();
      comp.push_back(u);
      for (uint64_t s : generators) {
        auto it = deep_rem.find(u ^ s);
        if (it != deep_rem.end()) {
          stack.push_back(*it);
          deep_rem.erase(it);
        }
      }
    }
    deep_largest = std::max(deep_largest, static_cast<int>(comp.size()));
    std::unordered_set<uint64_t> comp_set(comp.begin(), comp.end());
    bool self = false;
    for (uint64_t u : comp)
      if (comp_set.count(u ^ (classes - 1))) {
        self = true;
        break;
      }
    if (self) ++deep_self_antipodal;
  }
  std::cout << "deep plateau quotient components=" << deep_components
            << " largest=" << deep_largest
            << " self-antipodal=" << deep_self_antipodal << "\n";
}
