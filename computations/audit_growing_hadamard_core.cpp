#include <algorithm>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <string>
#include <vector>

using Matrix = std::vector<std::vector<int>>;

static Matrix kron(const Matrix &a, const Matrix &b) {
  const int m = static_cast<int>(a.size());
  const int n = static_cast<int>(b.size());
  Matrix out(m * n, std::vector<int>(m * n));
  for (int i = 0; i < m; ++i)
    for (int j = 0; j < m; ++j)
      for (int u = 0; u < n; ++u)
        for (int v = 0; v < n; ++v)
          out[i * n + u][j * n + v] = a[i][j] * b[u][v];
  return out;
}

static Matrix sylvester(int k) {
  Matrix h{{1}};
  const Matrix h2{{1, 1}, {1, -1}};
  while (static_cast<int>(h.size()) < k) h = kron(h, h2);
  return h;
}

static std::map<int, std::uint64_t> exact_histogram(const Matrix &a) {
  const int n = static_cast<int>(a.size());
  if (n >= 63) throw std::runtime_error("enumerator requires n<63");
  std::vector<int> x(n, -1), field(n, 0);
  long long energy = 0;
  for (int i = 0; i < n; ++i) {
    for (int j = i + 1; j < n; ++j) energy += a[i][j];
    for (int j = 0; j < n; ++j)
      if (j != i) field[i] += a[i][j] * x[j];
  }
  std::map<int, std::uint64_t> hist;
  ++hist[static_cast<int>(energy)];
  const std::uint64_t total = std::uint64_t{1} << n;
  for (std::uint64_t step = 1; step < total; ++step) {
    const int p = __builtin_ctzll(step);
    const int old = x[p];
    energy += -2LL * old * field[p];
    x[p] = -old;
    for (int j = 0; j < n; ++j)
      if (j != p) field[j] += -2 * old * a[j][p];
    ++hist[static_cast<int>(energy)];
  }
  return hist;
}

static long double log_cosh(long double z) {
  const long double a = std::fabs(z);
  return a + std::log1p(std::exp(-2 * a)) - std::log(2.0L);
}

static long double pressure(const std::map<int, std::uint64_t> &hist,
                            int n, long double t) {
  long double mx = -INFINITY;
  for (const auto &[e, c] : hist)
    mx = std::max(mx, std::log(static_cast<long double>(c)) + log_cosh(t * e));
  long double sum = 0;
  for (const auto &[e, c] : hist)
    sum += std::exp(std::log(static_cast<long double>(c)) + log_cosh(t * e) - mx);
  return mx + std::log(sum) - n * std::log(2.0L);
}

static int cap(const std::map<int, std::uint64_t> &hist) {
  int ans = 0;
  for (const auto &[e, c] : hist) ans = std::max(ans, std::abs(e));
  return ans;
}

static void print_hist(std::ostream &out,
                       const std::map<int, std::uint64_t> &hist) {
  out << "{";
  bool first = true;
  for (const auto &[e, c] : hist) {
    if (!first) out << ",";
    first = false;
    out << "\"" << e << "\":" << c;
  }
  out << "}";
}

int main(int argc, char **argv) {
  const Matrix a3{{0, -1, -1}, {-1, 0, -1}, {-1, -1, 0}};
  const Matrix a4{{0, 1, -1, -1},
                  {1, 0, -1, -1},
                  {-1, -1, 0, -1},
                  {-1, -1, -1, 0}};
  const std::vector<long double> betas{1, 2, 4, 8, 16, 24, 32};
  const Matrix regular4{{-1, 1, 1, 1},
                        {1, -1, 1, 1},
                        {1, 1, -1, 1},
                        {1, 1, 1, -1}};
  struct Job { int r; int k; Matrix a; Matrix t; std::string outer; };
  const std::vector<Job> jobs{
      {3, 2, a3, sylvester(2), "sylvester"},
      {3, 4, a3, sylvester(4), "sylvester"},
      {3, 8, a3, sylvester(8), "sylvester"},
      {3, 8, a3, kron(regular4, sylvester(2)), "regular4_tensor_h2"},
      {4, 2, a4, sylvester(2), "sylvester"},
      {4, 4, a4, sylvester(4), "sylvester"}};

  std::ofstream output_file;
  if (argc > 1) output_file.open(argv[1]);
  std::ostream &out = argc > 1 ? output_file : std::cout;

  out << std::setprecision(18);
  out << "{\n  \"schema\": \"exact Gray-code histograms; D=phi_(T tensor A)(beta/sqrt(kr))-k phi_A(beta/sqrt(r))\",\n";
  out << "  \"cases\": [\n";
  bool first_job = true;
  for (const auto &job : jobs) {
    if (!first_job) out << ",\n";
    first_job = false;
    Matrix core = kron(job.t, job.a);
    const auto child_hist = exact_histogram(job.a);
    const auto core_hist = exact_histogram(core);
    out << "    {\"r\":" << job.r << ",\"k\":" << job.k
              << ",\"N\":" << job.r * job.k
              << ",\"outer\":\"" << job.outer << "\""
              << ",\"child_cap\":" << cap(child_hist)
              << ",\"core_cap\":" << cap(core_hist)
              << ",\"child_histogram\":";
    print_hist(out, child_hist);
    out << ",\"core_histogram\":";
    print_hist(out, core_hist);
    out << ",\"beta_grid\":[";
    bool first_beta = true;
    for (long double beta : betas) {
      if (!first_beta) out << ",";
      first_beta = false;
      const long double pc = pressure(child_hist, job.r, beta / std::sqrt((long double)job.r));
      const long double pk = pressure(core_hist, job.r * job.k,
                                      beta / std::sqrt((long double)(job.r * job.k)));
      const long double defect = pk - job.k * pc;
      out << "{\"beta\":" << beta << ",\"child_pressure\":" << pc
                << ",\"core_pressure\":" << pk << ",\"D\":" << defect
                << ",\"D_over_kr\":" << defect / (job.k * job.r) << "}";
    }
    out << "]}";
  }
  out << "\n  ]\n}\n";
}
