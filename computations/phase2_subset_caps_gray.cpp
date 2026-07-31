#include <algorithm>
#include <cstdint>
#include <iostream>
#include <vector>

// Batch exact-cap evaluator for principal submatrices.  Input is
//
//   N K R
//   N*N entries of a symmetric signing
//   R rows of K distinct vertex indices
//
// and output is one line "cap min_energy max_energy" per subset.  The first
// listed child spin is fixed to +1 and the remaining spins are traversed in
// Gray-code order, so every projective Boolean state is checked.
int main() {
    std::ios::sync_with_stdio(false);
    std::cin.tie(nullptr);

    int parent_n, child_n;
    std::uint64_t records;
    if (!(std::cin >> parent_n >> child_n >> records) || parent_n < 1 ||
        child_n < 1 || child_n > parent_n || child_n > 63) {
        std::cerr << "invalid N K R header\n";
        return 2;
    }
    std::vector<std::vector<int>> parent(parent_n, std::vector<int>(parent_n));
    for (int i = 0; i < parent_n; ++i) {
        for (int j = 0; j < parent_n; ++j) {
            if (!(std::cin >> parent[i][j])) {
                std::cerr << "incomplete parent matrix\n";
                return 2;
            }
        }
    }

    const std::uint64_t states = std::uint64_t{1} << (child_n - 1);
    for (std::uint64_t record = 0; record < records; ++record) {
        std::vector<int> vertices(child_n);
        std::vector<unsigned char> used(parent_n, 0);
        for (int i = 0; i < child_n; ++i) {
            if (!(std::cin >> vertices[i]) || vertices[i] < 0 ||
                vertices[i] >= parent_n || used[vertices[i]]) {
                std::cerr << "invalid subset at record " << record << "\n";
                return 2;
            }
            used[vertices[i]] = 1;
        }

        std::vector<std::vector<int>> a(child_n, std::vector<int>(child_n));
        for (int i = 0; i < child_n; ++i) {
            for (int j = 0; j < child_n; ++j) {
                a[i][j] = parent[vertices[i]][vertices[j]];
            }
        }
        std::vector<int> x(child_n, 1), field(child_n, 0);
        std::int64_t energy = 0;
        for (int i = 0; i < child_n; ++i) {
            for (int j = 0; j < child_n; ++j) field[i] += a[i][j];
            for (int j = i + 1; j < child_n; ++j) energy += a[i][j];
        }
        std::int64_t minimum = energy, maximum = energy;
        for (std::uint64_t step = 1; step < states; ++step) {
            const int vertex = 1 + __builtin_ctzll(step);
            const int old = x[vertex];
            energy -= 2LL * old * field[vertex];
            for (int j = 0; j < child_n; ++j) {
                if (j != vertex) field[j] -= 2 * old * a[j][vertex];
            }
            x[vertex] = -old;
            minimum = std::min(minimum, energy);
            maximum = std::max(maximum, energy);
        }
        std::cout << std::max(-minimum, maximum) << ' ' << minimum << ' '
                  << maximum << '\n';
    }
    return 0;
}
