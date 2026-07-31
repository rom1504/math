// Exact finite audit: the strong cyclic two-fiber conference ansatz has no
// solution for k=4 (fiber size 33, conference order 66).
//
// Compile and run, without using /tmp:
//   g++ -O3 -std=c++17 computations/audit_k4_strong_cyclic_nonexistence.cpp -o /home/math/quadra/tmp/audit_k4_strong_cyclic_nonexistence
//   /home/math/quadra/tmp/audit_k4_strong_cyclic_nonexistence
//
// Floating point is used only to propose a separating integer vector.  Every
// rejection is accepted only after the integer inequality z^T(65I-AA^T)z<0
// is evaluated exactly in int64 arithmetic.  The final bridge search is an
// exhaustive integer backtrack.

#include <algorithm>
#include <array>
#include <cassert>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <set>
#include <vector>

namespace {

constexpr int s = 33;
constexpr int k = 4;
constexpr int conference_multiplier = 65;
constexpr double pi = 3.141592653589793238462643383279502884;

using HalfSet = std::array<int, 6>;

std::array<int, s> make_a(const HalfSet& negative_half) {
  std::array<int, s> a{};
  a.fill(1);
  a[0] = 0;
  for (int d : negative_half) a[d] = a[s - d] = -1;
  return a;
}

// Exact quadratic form z^T(65I-AA^T)z = 65||z||^2-||A^Tz||^2.
std::int64_t exact_psd_form(const std::array<int, s>& a,
                            const std::array<std::int64_t, s>& z) {
  std::int64_t z_norm = 0;
  std::array<std::int64_t, s> at_z{};
  for (int i = 0; i < s; ++i) {
    z_norm += z[i] * z[i];
    for (int j = 0; j < s; ++j) {
      // The choice of circulant orientation is immaterial to the norm.
      at_z[j] += a[(i - j + s) % s] * z[i];
    }
  }
  std::int64_t image_norm = 0;
  for (auto value : at_z) image_norm += value * value;
  return conference_multiplier * z_norm - image_norm;
}

bool has_exact_psd_separator(const std::array<int, s>& a) {
  // Search the real Fourier modes for a candidate negative eigendirection,
  // then certify it with a rounded integer vector.  If a coarse scale ever
  // failed, larger scales are tried; only the exact integer sign matters.
  for (int frequency = 1; frequency <= (s - 1) / 2; ++frequency) {
    double real = 0.0, imag = 0.0;
    for (int index = 0; index < s; ++index) {
      const double angle = 2.0 * pi * frequency * index / s;
      real += a[index] * std::cos(angle);
      imag -= a[index] * std::sin(angle);
    }
    if (conference_multiplier - (real * real + imag * imag) >= -0.25) {
      continue;
    }
    for (std::int64_t scale : {1000LL, 10000LL, 100000LL, 1000000LL}) {
      for (int phase = 0; phase < 2; ++phase) {
        std::array<std::int64_t, s> z{};
        for (int index = 0; index < s; ++index) {
          const double angle = 2.0 * pi * frequency * index / s;
          const double coordinate = phase == 0 ? std::cos(angle) : std::sin(angle);
          z[index] = std::llround(scale * coordinate);
        }
        if (exact_psd_form(a, z) < 0) return true;
      }
    }
  }
  return false;
}

std::set<std::vector<int>> multiplier_orbit(const HalfSet& representative) {
  std::set<int> full;
  for (int d : representative) full.insert(d), full.insert(s - d);
  std::set<std::vector<int>> orbit;
  for (int unit = 1; unit < s; ++unit) {
    if (std::gcd(unit, s) != 1) continue;
    std::vector<int> image;
    for (int d : full) image.push_back((unit * d) % s);
    std::sort(image.begin(), image.end());
    orbit.insert(image);
  }
  return orbit;
}

std::vector<int> full_negative_set(const HalfSet& half) {
  std::vector<int> result;
  for (int d : half) result.push_back(d), result.push_back(s - d);
  std::sort(result.begin(), result.end());
  return result;
}

struct BridgeSearch {
  std::array<int, 17> target{};
  std::array<int, 17> count{};
  std::array<int, 3> residue_count{};
  std::array<int, 16> selected{};
  std::uint64_t nodes = 0;
  bool found = false;

  explicit BridgeSearch(const std::array<int, s>& a) {
    for (int shift = 1; shift <= 16; ++shift) {
      int correlation = 0;
      for (int i = 0; i < s; ++i) correlation += a[i] * a[(i + shift) % s];
      // If c has 16 negative entries, corr_c(h)=-31+4N_R(h), while
      // corr_c(h)=-corr_a(h).
      assert((31 - correlation) % 4 == 0);
      target[shift] = (31 - correlation) / 4;
    }
  }

  static int circular_distance(int left, int right) {
    int distance = std::abs(left - right);
    return std::min(distance, s - distance);
  }

  bool dfs(int depth, int lower) {
    ++nodes;
    if (depth == 16) {
      for (int distance = 1; distance <= 16; ++distance) {
        if (count[distance] != target[distance]) return false;
      }
      found = true;
      return true;
    }
    const int remaining = 16 - depth;
    if (s - lower < remaining) return false;
    const std::array<int, 3> residue_target{4, 6, 6};
    for (int residue = 0; residue < 3; ++residue) {
      const int need = residue_target[residue] - residue_count[residue];
      if (need < 0) return false;
      int available = 0;
      for (int value = lower; value < s; ++value) available += value % 3 == residue;
      if (available < need) return false;
    }
    for (int distance = 1; distance <= 16; ++distance) {
      const int need = target[distance] - count[distance];
      if (need < 0 || need > depth * remaining + remaining * (remaining - 1) / 2) {
        return false;
      }
    }
    for (int value = lower; value <= s - remaining; ++value) {
      const int residue = value % 3;
      if (residue_count[residue] == residue_target[residue]) continue;
      std::array<int, 17> added{};
      bool admissible = true;
      for (int index = 0; index < depth; ++index) {
        const int distance = circular_distance(value, selected[index]);
        if (++added[distance] + count[distance] > target[distance]) {
          admissible = false;
          break;
        }
      }
      if (!admissible) continue;
      selected[depth] = value;
      ++residue_count[residue];
      for (int distance = 1; distance <= 16; ++distance) count[distance] += added[distance];
      if (dfs(depth + 1, value + 1)) return true;
      for (int distance = 1; distance <= 16; ++distance) count[distance] -= added[distance];
      --residue_count[residue];
    }
    return false;
  }
};

}  // namespace

int main() {
  const auto started = std::chrono::steady_clock::now();
  std::vector<HalfSet> survivors;
  std::uint64_t candidates = 0, exact_separators = 0;
  for (int a = 1; a <= 11; ++a)
    for (int b = a + 1; b <= 12; ++b)
      for (int c = b + 1; c <= 13; ++c)
        for (int d = c + 1; d <= 14; ++d)
          for (int e = d + 1; e <= 15; ++e)
            for (int f = e + 1; f <= 16; ++f) {
              const HalfSet half{a, b, c, d, e, f};
              ++candidates;
              if (has_exact_psd_separator(make_a(half))) {
                ++exact_separators;
              } else {
                survivors.push_back(half);
              }
            }
  assert(candidates == 8008);
  assert(exact_separators == 7998);
  assert(survivors.size() == 10);

  const HalfSet representative{1, 3, 4, 5, 12, 15};
  const auto orbit = multiplier_orbit(representative);
  assert(orbit.size() == 10);
  std::set<std::vector<int>> survivor_sets;
  for (const auto& half : survivors) survivor_sets.insert(full_negative_set(half));
  assert(survivor_sets == orbit);

  BridgeSearch bridge(make_a(representative));
  // For any bridge R, summing its difference counts over multiples of three
  // forces its three residue-class sizes to be a permutation of (4,6,6).
  // Translation puts the class of size four at residue zero; a further
  // translation by a multiple of three puts an omitted point at zero.  Thus
  // the following search over subsets of {1,...,32} is exhaustive.
  const bool exists = bridge.dfs(0, 1);
  assert(!exists);

  const double seconds = std::chrono::duration<double>(
      std::chrono::steady_clock::now() - started).count();
  std::cout << "{\n"
            << "  \"classification\": \"certified finite nonexistence for the strong cyclic ansatz only\",\n"
            << "  \"k\": 4,\n"
            << "  \"fiber_size\": 33,\n"
            << "  \"conference_order\": 66,\n"
            << "  \"symmetric_internal_candidates\": " << candidates << ",\n"
            << "  \"exact_integer_psd_separators\": " << exact_separators << ",\n"
            << "  \"spectral_survivors\": " << survivors.size() << ",\n"
            << "  \"survivor_multiplier_orbits\": 1,\n"
            << "  \"representative_negative_half\": [1,3,4,5,12,15],\n"
            << "  \"bridge_backtrack_nodes\": " << bridge.nodes << ",\n"
            << "  \"bridge_exists\": false,\n"
            << "  \"elapsed_seconds_informational\": " << seconds << "\n"
            << "}\n";
}
