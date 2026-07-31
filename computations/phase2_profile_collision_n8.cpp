#include <array>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <unordered_map>
#include <vector>

// Exhaustively seek two order-8 signings with the same complete histograms of
// switching/permutation/global-negation classes on restrictions of orders
// 4,5,6 but different Boolean caps.  Root edges are gauged to +1 and one of a
// globally-negative pair is enumerated, leaving 2^20 exact cases.

int signing_order = 8;

struct Dsu {
    std::vector<int> p;
    explicit Dsu(int n) : p(n) { std::iota(p.begin(), p.end(), 0); }
    int find(int x) { return p[x] == x ? x : p[x] = find(p[x]); }
    void join(int a, int b) {
        a = find(a); b = find(b);
        if (a != b) p[std::max(a,b)] = std::min(a,b);
    }
};

int edge_index(int n, int i, int j) {
    if (i > j) std::swap(i,j);
    int bit = 0;
    for (int a = 1; a < n; ++a)
        for (int b = a + 1; b < n; ++b, ++bit)
            if (a == i && b == j) return bit;
    return -1;
}

int sign_from_code(std::uint32_t code, int n, int i, int j) {
    if (i == j) return 0;
    if (i == 0 || j == 0) return 1;
    return (code & (std::uint32_t{1} << edge_index(n,i,j))) ? -1 : 1;
}

std::uint32_t permute_and_gauge(std::uint32_t code, int n,
                                const std::vector<int>& permutation) {
    std::uint32_t result = 0;
    int bit = 0;
    for (int i = 1; i < n; ++i) {
        for (int j = i + 1; j < n; ++j, ++bit) {
            int value = sign_from_code(code,n,permutation[i],permutation[j]) *
                        sign_from_code(code,n,permutation[0],permutation[i]) *
                        sign_from_code(code,n,permutation[0],permutation[j]);
            if (value == -1) result |= std::uint32_t{1} << bit;
        }
    }
    return result;
}

std::vector<int> class_map(int n) {
    const int bits = (n-1)*(n-2)/2;
    const int count = 1 << bits;
    Dsu dsu(count);
    for (int code = 0; code < count; ++code) {
        dsu.join(code, code ^ (count-1));
        for (int swap_at = 0; swap_at < n-1; ++swap_at) {
            std::vector<int> permutation(n);
            std::iota(permutation.begin(), permutation.end(), 0);
            std::swap(permutation[swap_at], permutation[swap_at+1]);
            dsu.join(code, permute_and_gauge(code,n,permutation));
        }
    }
    std::unordered_map<int,int> labels;
    std::vector<int> result(count);
    for (int code = 0; code < count; ++code) {
        int root = dsu.find(code);
        auto [it, inserted] = labels.emplace(root, static_cast<int>(labels.size()));
        result[code] = it->second;
    }
    std::cerr << "restriction order " << n << " classes " << labels.size() << "\n";
    return result;
}

std::vector<std::vector<int>> subsets(int n, int k) {
    std::vector<std::vector<int>> result;
    std::vector<int> chosen(k);
    auto visit = [&](auto&& self, int start, int depth) -> void {
        if (depth == k) { result.push_back(chosen); return; }
        for (int value = start; value <= n-(k-depth); ++value) {
            chosen[depth] = value;
            self(self, value+1, depth+1);
        }
    };
    visit(visit,0,0);
    return result;
}

std::uint32_t restriction_code(std::uint32_t code,
                               const std::vector<int>& vertices) {
    std::uint32_t result = 0;
    int bit = 0;
    for (int i = 1; i < static_cast<int>(vertices.size()); ++i) {
        for (int j = i+1; j < static_cast<int>(vertices.size()); ++j, ++bit) {
            int value = sign_from_code(code,signing_order,vertices[i],vertices[j]) *
                        sign_from_code(code,signing_order,vertices[0],vertices[i]) *
                        sign_from_code(code,signing_order,vertices[0],vertices[j]);
            if (value == -1) result |= std::uint32_t{1} << bit;
        }
    }
    return result;
}

int boolean_cap(std::uint32_t code) {
    int result = 0;
    const int spin_count = 1 << (signing_order-1);
    const int nonroot_edges = (signing_order-1)*(signing_order-2)/2;
    for (int spin_code = 0; spin_code < spin_count; ++spin_code) {
        int root_term = 0;
        std::uint32_t product_negative = 0;
        int bit = 0;
        for (int i = 1; i < signing_order; ++i) {
            int xi = (spin_code & (1 << (i-1))) ? -1 : 1;
            root_term += xi;
            for (int j = i+1; j < signing_order; ++j, ++bit) {
                int xj = (spin_code & (1 << (j-1))) ? -1 : 1;
                if (xi*xj == -1) product_negative |= std::uint32_t{1} << bit;
            }
        }
        int energy = root_term + nonroot_edges - 2*__builtin_popcount(code ^ product_negative);
        result = std::max(result, std::abs(energy));
    }
    return result;
}

using Key = std::array<std::uint8_t,16>;
struct KeyHash {
    std::size_t operator()(const Key& key) const {
        std::size_t value = 1469598103934665603ULL;
        for (auto item : key) { value ^= item; value *= 1099511628211ULL; }
        return value;
    }
};
struct Record { int cap; std::uint32_t code; };

int main(int argc, char** argv) {
    std::uint32_t requested = 0;
    if (argc == 3 && std::string(argv[1]) == "--order9-sample") {
        signing_order = 9;
        requested = static_cast<std::uint32_t>(std::stoul(argv[2]));
    } else if (argc != 1) {
        std::cerr << "usage: phase2_profile_collision_n8 [--order9-sample COUNT]\n";
        return 2;
    }
    auto classes4 = class_map(4), classes5 = class_map(5), classes6 = class_map(6);
    auto subsets4 = subsets(signing_order,4), subsets5 = subsets(signing_order,5), subsets6 = subsets(signing_order,6);
    std::unordered_map<Key,Record,KeyHash> seen;
    seen.reserve(1 << 19);
    const int signing_bits = (signing_order-1)*(signing_order-2)/2;
    const std::uint32_t universe = std::uint32_t{1} << (signing_bits-1);
    const std::uint32_t total = requested ? std::min(requested, universe) : universe;
    const std::uint32_t multiplier = 747796405u; // odd: a permutation modulo a power of two
    const std::uint32_t increment = 2891336453u;
    for (std::uint32_t index = 0; index < total; ++index) {
        const std::uint32_t code = requested ?
            (multiplier*index + increment) & (universe-1) : index;
        Key key{};
        for (const auto& subset : subsets4) ++key[classes4[restriction_code(code,subset)]];
        for (const auto& subset : subsets5) ++key[2 + classes5[restriction_code(code,subset)]];
        for (const auto& subset : subsets6) ++key[6 + classes6[restriction_code(code,subset)]];
        int cap = boolean_cap(code);
        auto [it, inserted] = seen.emplace(key, Record{cap,code});
        if (!inserted && it->second.cap != cap) {
            std::cout << "{\"classification\":\"exact counterexample after exhaustive-state traversal to the reported code\",";
            std::cout << "\"n\":" << signing_order << ",\"profiles\":[4,5,6],\"first_code\":" << it->second.code;
            std::cout << ",\"first_cap\":" << it->second.cap << ",\"second_code\":" << code;
            std::cout << ",\"second_cap\":" << cap << ",\"codes_checked\":" << index+1;
            std::cout << ",\"profile\":[";
            for (int i=0;i<16;++i) { if(i) std::cout << ','; std::cout << int(key[i]); }
            std::cout << "]}\n";
            return 0;
        }
        if (index && index % 100000 == 0) std::cerr << "checked " << index << "/" << total << "\n";
    }
    std::cout << "{\"classification\":\"" << (requested ? "deterministic collision sample; no collision found" : "exhaustive finite no-collision result") << "\",\"n\":" << signing_order << ',';
    std::cout << "\"profiles\":[4,5,6],\"signings_checked\":" << total;
    std::cout << ",\"distinct_profile_states\":" << seen.size() << "}\n";
}
