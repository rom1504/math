#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <limits>
#include <vector>

static int legendre(int a, int q) {
    a %= q;
    if (a < 0) a += q;
    if (a == 0) return 0;
    int r = 1, b = a, e = (q - 1) / 2;
    while (e) {
        if (e & 1) r = int((1LL * r * b) % q);
        b = int((1LL * b * b) % q);
        e >>= 1;
    }
    return r == 1 ? 1 : -1;
}

static int energy(const std::vector<std::vector<int>>& A,
                  const std::vector<int>& s,
                  const std::vector<int>& x) {
    int h = 0;
    for (int i = 0; i < (int)s.size(); ++i)
        for (int j = i + 1; j < (int)s.size(); ++j)
            h += A[s[i]][s[j]] * x[i] * x[j];
    return h;
}

struct Extremizers {
    int p = std::numeric_limits<int>::min();
    int lo = std::numeric_limits<int>::max();
    std::vector<std::vector<int>> pos, neg;
};

static Extremizers extremizers(const std::vector<std::vector<int>>& A,
                               const std::vector<int>& s) {
    const int m = int(s.size());
    Extremizers e;
    for (uint64_t mask = 0; mask < (uint64_t(1) << (m - 1)); ++mask) {
        std::vector<int> x(m, 1);
        for (int i = 1; i < m; ++i)
            if ((mask >> (i - 1)) & 1ULL) x[i] = -1;
        int h = energy(A, s, x);
        if (h > e.p) {
            e.p = h;
            e.pos.clear();
        }
        if (h == e.p) e.pos.push_back(x);
        if (h < e.lo) {
            e.lo = h;
            e.neg.clear();
        }
        if (h == e.lo) e.neg.push_back(x);
    }
    return e;
}

static int separation(const std::vector<std::vector<int>>& A,
                      const std::vector<int>& s,
                      const std::vector<int>& t,
                      const std::vector<int>& u,
                      const std::vector<int>& v) {
    int lp = 0, lm = 0;
    for (int j : t) {
        int ap = 0, am = 0;
        for (int i = 0; i < (int)s.size(); ++i) {
            ap += A[s[i]][j] * (u[i] + v[i]);
            am += A[s[i]][j] * (u[i] - v[i]);
        }
        lp += std::abs(ap);
        lm += std::abs(am);
    }
    return std::max(lp, lm);
}

int main(int argc, char** argv) {
    const int q = argc > 1 ? std::stoi(argv[1]) : 13;
    const int N = q + 1;
    const int m = argc > 2 ? std::stoi(argv[2]) : N / 2;
    std::vector<std::vector<int>> A(N, std::vector<int>(N));
    for (int j = 0; j < q; ++j) A[q][j] = A[j][q] = 1;
    for (int i = 0; i < q; ++i)
        for (int j = 0; j < q; ++j)
            if (i != j) A[i][j] = legendre(i - j, q);

    std::vector<int> all(N);
    for (int i = 0; i < N; ++i) all[i] = i;
    Extremizers efull = extremizers(A, all);
    const int K = std::max(efull.p, -efull.lo);
    const double alpha = double(m) / N;
    const double target = 2.0 * (1.0 - std::pow(alpha, 1.5)) * K;

    int best_r = -1, worst_r = std::numeric_limits<int>::max();
    int best_centered = -1, centered_count = 0;
    int best_balanced2 = -1, balanced2_count = 0, sufficient = 0;
    int min_mchild = std::numeric_limits<int>::max();
    uint64_t count = 0;
    std::vector<int> choose(N, 0);
    std::fill(choose.begin(), choose.begin() + m, 1);
    do {
        std::vector<int> s, t;
        for (int i = 0; i < N; ++i)
            (choose[i] ? s : t).push_back(i);
        Extremizers e = extremizers(A, s);
        int rmax = -1;
        for (const auto& u : e.pos)
            for (const auto& v : e.neg)
                rmax = std::max(rmax, separation(A, s, t, u, v));
        best_r = std::max(best_r, rmax);
        worst_r = std::min(worst_r, rmax);
        min_mchild = std::min(min_mchild, std::max(e.p, -e.lo));
        if (e.p == -e.lo) {
            ++centered_count;
            best_centered = std::max(best_centered, rmax);
        }
        if (std::abs(e.p + e.lo) <= 2) {
            ++balanced2_count;
            best_balanced2 = std::max(best_balanced2, rmax);
        }
        if (rmax + 1e-12 >= target) ++sufficient;
        ++count;
    } while (std::prev_permutation(choose.begin(), choose.end()));

    std::cout << "q " << q << " N " << N << " m " << m
              << " K " << K << " full_P " << efull.p
              << " full_Q " << -efull.lo << " target " << target
              << " subsets " << count
              << " min_child_M " << min_mchild
              << " best_R " << best_r
              << " worst_R " << worst_r
              << " sufficient " << sufficient
              << " centered " << centered_count
              << " best_centered_R " << best_centered
              << " balanced2 " << balanced2_count
              << " best_balanced2_R " << best_balanced2 << "\n";
}
