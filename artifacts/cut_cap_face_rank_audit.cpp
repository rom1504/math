#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <map>
#include <tuple>
#include <vector>

static int real_rank(std::vector<std::vector<double>> a) {
  if (a.empty()) return 0;
  const int rows = static_cast<int>(a.size());
  const int cols = static_cast<int>(a[0].size());
  int rank = 0;
  for (int col = 0; col < cols && rank < rows; ++col) {
    int pivot = rank;
    while (pivot < rows && std::abs(a[pivot][col]) < 1e-9) ++pivot;
    if (pivot == rows) continue;
    std::swap(a[rank], a[pivot]);
    const double scale = a[rank][col];
    for (int j = col; j < cols; ++j) a[rank][j] /= scale;
    for (int i = 0; i < rows; ++i) {
      if (i == rank || std::abs(a[i][col]) < 1e-9) continue;
      const double factor = a[i][col];
      for (int j = col; j < cols; ++j)
        a[i][j] -= factor * a[rank][j];
    }
    ++rank;
  }
  return rank;
}

int main(int argc, char **argv) {
  const int n = argc > 1 ? std::atoi(argv[1]) : 8;
  const int m = n - 1;
  const int states = 1 << m;
  std::vector<std::pair<int, int>> free_edges;
  for (int i = 0; i < m; ++i)
    for (int j = i + 1; j < m; ++j) free_edges.push_back({i, j});
  const int free_count = static_cast<int>(free_edges.size());
  if (free_count >= 63) return 2;
  const uint64_t classes = uint64_t{1} << free_count;

  auto matrix_from_mask = [&](uint64_t mask) {
    std::vector<std::vector<int>> a(n, std::vector<int>(n));
    for (int i = 0; i < m; ++i) a[i][m] = a[m][i] = 1;
    for (int e = 0; e < free_count; ++e) {
      auto [i, j] = free_edges[e];
      const int sign = ((mask >> e) & 1) ? -1 : 1;
      a[i][j] = a[j][i] = sign;
    }
    return a;
  };

  auto spins = [&](int state) {
    std::vector<int> x(n, 1);
    for (int i = 0; i < m; ++i)
      if ((state >> i) & 1) x[i] = -1;
    return x;
  };

  auto energies = [&](const std::vector<std::vector<int>> &a) {
    std::vector<int> e(states);
    for (int state = 0; state < states; ++state) {
      auto x = spins(state);
      int value = 0;
      for (int i = 0; i < n; ++i)
        for (int j = i + 1; j < n; ++j)
          value += a[i][j] * x[i] * x[j];
      e[state] = value;
    }
    return e;
  };

  int best_width = 1 << 30;
  std::vector<uint64_t> optimal;
  std::vector<int16_t> q(states);
  std::vector<std::vector<int8_t>> chi(
      free_count, std::vector<int8_t>(states));
  for (int state = 0; state < states; ++state) {
    auto x = spins(state);
    int value = 0;
    for (int i = 0; i < n; ++i)
      for (int j = i + 1; j < n; ++j) value += x[i] * x[j];
    q[state] = static_cast<int16_t>(value);
    for (int f = 0; f < free_count; ++f) {
      auto [i, j] = free_edges[f];
      chi[f][state] = static_cast<int8_t>(x[i] * x[j]);
    }
  }
  uint64_t previous_gray = 0;
  std::vector<int8_t> free_sign(free_count, 1);
  for (uint64_t step = 0; step < classes; ++step) {
    const uint64_t gray = step ^ (step >> 1);
    if (step) {
      const uint64_t diff = gray ^ previous_gray;
      const int f = __builtin_ctzll(diff);
      const int delta = -2 * free_sign[f];
      for (int state = 0; state < states; ++state)
        q[state] = static_cast<int16_t>(
            q[state] + delta * chi[f][state]);
      free_sign[f] = static_cast<int8_t>(-free_sign[f]);
    }
    previous_gray = gray;
    const auto [lo, hi] = std::minmax_element(q.begin(), q.end());
    const int width = (*hi - *lo) / 2;
    if (width < best_width) {
      best_width = width;
      optimal.clear();
    }
    if (width == best_width) optimal.push_back(gray);
  }

  // Key: |midpoint|, smaller side, larger side, internal dimension,
  // active-facet rank, number of active cut vectors.
  std::map<std::tuple<int, int, int, int, int, int>, uint64_t> profile;
  int min_rank_defect = 1 << 30, max_rank_defect = 0;
  uint64_t endpoint_pairs = 0;

  for (uint64_t mask : optimal) {
    auto a = matrix_from_mask(mask);
    auto e = energies(a);
    const int top_value = *std::max_element(e.begin(), e.end());
    const int bottom_value = *std::min_element(e.begin(), e.end());
    std::vector<int> tops, bottoms;
    for (int state = 0; state < states; ++state) {
      if (e[state] == top_value) tops.push_back(state);
      if (e[state] == bottom_value) bottoms.push_back(state);
    }

    for (int ts : tops) {
      auto xt = spins(ts);
      std::vector<std::vector<int>> switched = a;
      for (int i = 0; i < n; ++i)
        for (int j = i + 1; j < n; ++j)
          switched[i][j] = switched[j][i] =
              a[i][j] * xt[i] * xt[j];

      for (int bs : bottoms) {
        ++endpoint_pairs;
        auto xb = spins(bs);
        std::vector<int> side(n);
        int usize = 0;
        for (int i = 0; i < n; ++i) {
          side[i] = xt[i] * xb[i] == -1;
          usize += side[i];
        }
        const int vsize = n - usize;

        std::vector<std::pair<int, int>> internal;
        for (int i = 0; i < n; ++i)
          for (int j = i + 1; j < n; ++j)
            if (side[i] == side[j]) internal.push_back({i, j});
        const int internal_dim = static_cast<int>(internal.size());

        std::vector<std::vector<double>> rows;
        for (int zs = 0; zs < states; ++zs) {
          auto z = spins(zs);
          int cross_energy = 0;
          for (int i = 0; i < n; ++i)
            for (int j = i + 1; j < n; ++j)
              if (side[i] != side[j])
                cross_energy += switched[i][j] * z[i] * z[j];
          if (std::abs(cross_energy) != best_width) continue;
          std::vector<double> row(internal_dim);
          for (int f = 0; f < internal_dim; ++f) {
            auto [i, j] = internal[f];
            row[f] = z[i] != z[j] ? 1.0 : 0.0;
          }
          rows.push_back(std::move(row));
        }
        const int rank = real_rank(rows);
        const int defect = internal_dim - rank;
        min_rank_defect = std::min(min_rank_defect, defect);
        max_rank_defect = std::max(max_rank_defect, defect);
        const int midpoint = (top_value + bottom_value) / 2;
        profile[{std::abs(midpoint), std::min(usize, vsize),
                 std::max(usize, vsize), internal_dim, rank,
                 static_cast<int>(rows.size())}]++;
      }
    }
  }

  std::cout << "n=" << n << " classes=" << classes
            << " W=" << best_width
            << " minimizers=" << optimal.size()
            << " endpoint-pairs=" << endpoint_pairs
            << " rank-defect-range=[" << min_rank_defect << ","
            << max_rank_defect << "]\n";
  for (const auto &[key, count] : profile) {
    auto [mid, us, vs, dim, rank, active] = key;
    std::cout << "|d|=" << mid << " sides=" << us << "+" << vs
              << " internal-dim=" << dim << " active-rank=" << rank
              << " active-cuts=" << active << " count=" << count << "\n";
  }
}
