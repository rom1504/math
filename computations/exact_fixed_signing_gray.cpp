#include <algorithm>
#include <cstdint>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

// Exhaustively evaluate one fixed symmetric zero-diagonal signing.  Vertex 0
// is fixed to +1 and the remaining spins are traversed in Gray-code order.
// Input: n followed by the n*n integer matrix entries on standard input.
// Output: a compact JSON record on standard output.
int main(int argc, char** argv) {
    std::ios::sync_with_stdio(false);
    std::cin.tie(nullptr);

    bool collect = false;
    if (argc == 2 && std::string(argv[1]) == "--collect-extremizers") {
        collect = true;
    } else if (argc != 1) {
        std::cerr << "usage: exact_fixed_signing_gray [--collect-extremizers]\n";
        return 2;
    }
    int n;
    if (!(std::cin >> n) || n < 1 || n > 63) {
        std::cerr << "expected 1 <= n <= 63\n";
        return 2;
    }
    std::vector<std::vector<int>> a(n, std::vector<int>(n));
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            if (!(std::cin >> a[i][j])) {
                std::cerr << "incomplete matrix\n";
                return 2;
            }
        }
    }
    for (int i = 0; i < n; ++i) {
        if (a[i][i] != 0) {
            std::cerr << "nonzero diagonal\n";
            return 2;
        }
        for (int j = i + 1; j < n; ++j) {
            if (a[i][j] != a[j][i] || (a[i][j] != 1 && a[i][j] != -1)) {
                std::cerr << "matrix is not a symmetric signing\n";
                return 2;
            }
        }
    }

    std::vector<int> x(n, 1);
    std::vector<int> field(n, 0);
    std::int64_t energy = 0;
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) field[i] += a[i][j];
        for (int j = i + 1; j < n; ++j) energy += a[i][j];
    }

    std::int64_t minimum = energy, maximum = energy;
    std::uint64_t argmin = 0, argmax = 0;
    std::uint64_t min_count = 1, max_count = 1;
    std::vector<std::uint64_t> minimizers{0}, maximizers{0};
    const std::uint64_t total = std::uint64_t{1} << (n - 1);
    for (std::uint64_t step = 1; step < total; ++step) {
        const int vertex = 1 + __builtin_ctzll(step);
        const int old = x[vertex];
        energy -= 2LL * old * field[vertex];
        for (int j = 0; j < n; ++j) {
            if (j != vertex) field[j] -= 2 * old * a[j][vertex];
        }
        x[vertex] = -old;
        const std::uint64_t gray = step ^ (step >> 1);
        if (energy < minimum) {
            minimum = energy;
            argmin = gray;
            min_count = 1;
            if (collect) minimizers.assign(1, gray);
        } else if (energy == minimum) {
            ++min_count;
            if (collect) minimizers.push_back(gray);
        }
        if (energy > maximum) {
            maximum = energy;
            argmax = gray;
            max_count = 1;
            if (collect) maximizers.assign(1, gray);
        } else if (energy == maximum) {
            ++max_count;
            if (collect) maximizers.push_back(gray);
        }
    }
    const std::int64_t cap = std::max(-minimum, maximum);
    std::cout << "{\"n\":" << n
              << ",\"projective_assignments\":" << total
              << ",\"min_energy\":" << minimum
              << ",\"max_energy\":" << maximum
              << ",\"cap\":" << cap
              << ",\"argmin_gray\":" << argmin
              << ",\"argmax_gray\":" << argmax
              << ",\"min_count\":" << min_count
              << ",\"max_count\":" << max_count;
    if (collect) {
        std::cout << ",\"minimizer_gray_codes\":[";
        for (std::size_t i = 0; i < minimizers.size(); ++i) {
            if (i) std::cout << ',';
            std::cout << minimizers[i];
        }
        std::cout << "],\"maximizer_gray_codes\":[";
        for (std::size_t i = 0; i < maximizers.size(); ++i) {
            if (i) std::cout << ',';
            std::cout << maximizers[i];
        }
        std::cout << ']';
    }
    std::cout << "}\n";
    return 0;
}
