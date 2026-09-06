// Reproducible lower witnesses only. No cap optimality or convergence claim.
#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <random>
#include <vector>

int main(int argc, char** argv) {
    const int seed_order = argc > 1 ? std::atoi(argv[1]) : 2;
    const int power = argc > 2 ? std::atoi(argv[2]) : 3;
    const double retention = argc > 3 ? std::atof(argv[3]) : .8;
    const int restarts = argc > 4 ? std::atoi(argv[4]) : 30;
    const unsigned rng_seed = argc > 5 ? std::strtoul(argv[5], nullptr, 10) : 20260906;
    if ((seed_order != 2 && seed_order != 6) || power < 0 || power > 5 ||
        retention <= 0 || retention > 1 || restarts < 1) return 2;
    int outer = 1;
    for (int j = 0; j < power; ++j) outer *= 4;
    int full_order = outer * seed_order;
    int n = std::max(1, int(std::floor(retention * full_order)));
    std::vector<int> seed(seed_order * seed_order, 1);
    if (seed_order == 2) seed = {1, 1, 1, -1};
    else {
        int bit = 0;
        for (int i = 1; i < 6; ++i) for (int j = i + 1; j < 6; ++j)
            seed[i * 6 + j] = seed[j * 6 + i] = 1 - 2 * ((220 >> bit++) & 1);
        // Zero trace preserves the exact hollow seed cap Q=5.
        for (int i = 0; i < 6; ++i) seed[i * 6 + i] = i < 3 ? 1 : -1;
    }
    auto entry = [&](int i, int j) {
        int a = i / seed_order, b = j / seed_order;
        int value = seed[(i % seed_order) * seed_order + j % seed_order];
        for (int d = 0; d < power; ++d) {
            value *= (a % 4 == b % 4) ? -1 : 1;
            a /= 4; b /= 4;
        }
        return value;
    };
    std::mt19937 rng(rng_seed);
    std::uniform_real_distribution<double> uniform(0, 1);
    std::vector<int> selected(full_order);
    std::iota(selected.begin(), selected.end(), 0);
    std::shuffle(selected.begin(), selected.end(), rng);
    selected.resize(n);
    std::sort(selected.begin(), selected.end());
    std::vector<int> matrix(n * n), spin(n), field(n), best_spin;
    for (int i = 0; i < n; ++i) for (int j = 0; j < n; ++j)
        matrix[i * n + j] = i == j ? 0 : entry(selected[i], selected[j]);
    long long best = -1;
    int best_sign = 0;
    for (int sign : {-1, 1}) for (int run = 0; run < restarts; ++run) {
        for (int &x : spin) x = 2 * (rng() % 2) - 1;
        long long energy = 0;
        for (int i = 0; i < n; ++i) {
            field[i] = 0;
            for (int j = 0; j < n; ++j) field[i] += sign * matrix[i * n + j] * spin[j];
            energy += spin[i] * field[i];
        }
        energy /= 2;
        for (int sweep = 0; sweep < 100; ++sweep) {
            double temperature = 1.5 * std::sqrt(n) * std::pow(.003, sweep / 99.0);
            for (int step = 0; step < n; ++step) {
                int i = rng() % n, old = spin[i], delta = -2 * old * field[i];
                if (delta >= 0 || uniform(rng) < std::exp(delta / temperature)) {
                    spin[i] = -old;
                    energy += delta;
                    for (int j = 0; j < n; ++j) field[j] -= 2 * old * sign * matrix[j * n + i];
                }
            }
        }
        bool changed = true;
        while (changed) {
            changed = false;
            for (int i = 0; i < n; ++i) if (spin[i] * field[i] < 0) {
                int old = spin[i];
                energy -= 2 * old * field[i];
                spin[i] = -old;
                for (int j = 0; j < n; ++j) field[j] -= 2 * old * sign * matrix[j * n + i];
                changed = true;
            }
        }
        long long replay = 0;
        for (int i = 0; i < n; ++i) for (int j = 0; j < i; ++j)
            replay += sign * matrix[i * n + j] * spin[i] * spin[j];
        assert(replay == energy);
        if (energy > best) { best = energy; best_spin = spin; best_sign = sign; }
    }
    std::uint64_t spin_hash = 1469598103934665603ULL;
    for (int x : best_spin) { spin_hash ^= std::uint64_t(x + 1); spin_hash *= 1099511628211ULL; }
    std::cout << std::setprecision(16)
              << "{\"seed_order\":" << seed_order << ",\"outer_order\":" << outer
              << ",\"retention\":" << retention << ",\"order\":" << n
              << ",\"restarts_per_sign\":" << restarts << ",\"rng_seed\":" << rng_seed
              << ",\"witness_absolute_energy\":" << best
              << ",\"normalized_lower_witness\":" << best / std::pow(n, 1.5)
              << ",\"objective_sign\":" << best_sign << ",\"spin_hash\":\"" << spin_hash
              << "\",\"scope\":\"lower witness, not exact cap\"}" << std::endl;
}
