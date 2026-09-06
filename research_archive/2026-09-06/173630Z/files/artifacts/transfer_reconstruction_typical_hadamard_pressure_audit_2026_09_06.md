# Independent audit: typical retained-Hadamard profiles and actual annealed pressure

Date: 2026-09-06. The finite constants and the annealed-pressure reduction
in `transfer_director_typical_hadamard_spectra_2026_09_06.md` survive
independent reconstruction. The hard-tail corollary remains a separate
statement; this note audits the exponential partition lower bound only.

## 1. The smooth-test constants

Let `f` be even and `C^4`, with every derivative through order four
bounded in absolute value by `M`. Each normalized input direction into
one spectrum coordinate has magnitude `k^-1/2`. Rademacher and standard
Gaussian inputs agree through their third moments, and their fourth
absolute moments sum to four. Taylor's remainder therefore pays

```math
k\frac{(1+3)M}{24k^2}=\frac{M}{6k}
```

for the mean. For the product test `f(z_1)f(z_2)`, its fourth derivative
in an input direction `u` is bounded by
`M^2(|u_1|+|u_2|)^4<=16M^2/k^2`. The product expectation error is thus
`8M^2/(3k)`. Both individual means are bounded by `M`; replacing their
product costs at most `M^2/(3k)`. The covariance error is at most
`3M^2/k`.

For the Gaussian comparison pair of correlation `R_ij`, the centered
even Hermite expansion gives nonnegative covariance at most
`R_ij^2 Var(f(G))<=M^2R_ij^2`, including correlations `+-1` by continuity.
For `R=H_T^TH_T/k`, exact row orthogonality gives

```math
R^2=(m/k)R,\qquad\operatorname{tr}R=m,\qquad
\sum_{ij}R_{ij}^2=m^2/k.
```

Averaging every pair covariance, including diagonal pairs, proves exactly
`Var(Y_f)<=4M^2/k`. No independence or smallness of the individual
off-diagonal correlations is assumed.

## 2. Uniform W2 and fixed quantizers

The spectrum's second moment is exactly one for every spin. Fix a large
`R` and a smooth even function `0<=f_R<=z^2` which equals `z^2` on
`[-R,R]` and vanishes outside `[-R-1,R+1]`. Its Gaussian integral can be
made arbitrarily close to one. The finite mean/variance estimates imply
that its empirical integral is close with uniformly high probability.
Subtracting from the exact second moment bounds the empirical tail energy
outside `[-R-1,R+1]` by the same small error.

On that fixed compact interval, finitely many smooth approximations of
bin indicators control the compact empirical law. Their boundaries may
be chosen with zero Gaussian mass. Transport to bin representatives and
then between the resulting finite distributions gives the compact `W_2`
bound. Choose the truncation and test functions first, and only then let
`k` grow. Every probability estimate is uniform in the Hadamard order,
basis and selector, proving uniform `W_2` convergence in probability.

For a fixed clipped finite quantizer, its bin-frequency functionals are
continuous at the Gaussian law under this convergence. Its squared
quantization error is also continuous there: its discontinuities have
Gaussian mass zero and it has at most quadratic growth, while `W_2`
convergence supplies the required tail integrability. Thus one may choose
a slowly vanishing histogram tolerance and obtain `(1-o(1))2^k` good
spins in every fibre, uniformly. This is a fraction of all row spins,
not a claim that the random spectrum coordinates are independent.

## 3. Applying the actual endpoint count

Fix the finite quantizer and a bound `epsilon` for the mean squared
endpoint error. The good row counts obey the hypothesis of
`transfer_reconstruction_homogeneous_edge_count_2026_09_06.md`. Its
lower-bound error is uniform over these row counts: the repair radius
depends only on their total discrepancy, the alphabet is fixed, and the
logarithm of the positive kernel has a fixed finite oscillation.

Condition the diagonal spectrum coordinate in each fibre to have
magnitude at most `R_0>1`. Its probability is at least `1-R_0^-2` by
the exact row norm; the product probability costs `O(m)`. The diagonal
defect equals `(1-S_ii)^2v_i(i)^2`, so its exponential factor is at least
`exp(-2tR_0^2m)`. The factor two here is correct. After the chosen
coordinate is removed, each remaining multiset is independently uniform
on its off-diagonal endpoints and its counts have changed by only one.

Let `a` be the vector of actual off-diagonal endpoint magnitudes and `b`
its quantization. For the folded kernel,

```math
\log K_t(a,b)=-t(a^2+b^2)+\log\cosh(2tab),\qquad
|\partial_a\log K_t(a,b)|\le2t(|a|+|b|).
```

Integrating along each endpoint segment and summing over edges yields
the explicit global bound

```math
\left|\sum_{i<j}\log K_t(a_{ij},a_{ji})
-\sum_{i<j}\log K_t(b_{ij},b_{ji})\right|
\le2t(\|a\|_2+\|b\|_2)\|a-b\|_2
\le2t(2+\sqrt\epsilon)\sqrt\epsilon\,m^2.
```

Here `||a||_2<=m`, `||a-b||_2<=sqrt(epsilon)m`, and
`||b||_2<=(1+sqrt(epsilon))m`. The estimate is uniform over all endpoint
permutations. It does not require a global bounded derivative or any
favorable cancellation between signed edge contributions.

The finite-color count gives edge pressure `-F_(K_t)(mu_Q)/2`, where
`mu_Q` is the quantized magnitude law. Relative-sign folding identifies
this with `Phi_t(nu_Q)` for its symmetrized source law. Refining the fixed
quantizer only after the order limit gives `nu_Q ->N(0,1)` in `W_2`,
hence `Phi_t(nu_Q)->g_t(1)`. Multiplying by the product number of good
row spins contributes `mk log2+o(m^2)`. Therefore the actual annealed
partition, before Finner, has lower rate `p log2+g_t(1)` as claimed.

All limit orders are controlled: fixed `p,t`, then fixed truncation and
finite mesh, then matrix order, then truncation/mesh refinement. A growing
alphabet, vanishing kernel minimum, or order-dependent temperature is not
silently imported into this annealed theorem.

## 4. Bounded reproducible checks

The exact endpoint-count and repair checker is

```sh
.venv/bin/python computations/transfer_reconstruction_homogeneous_edge_count_exact_checks_2026_09_06.py
```

It checks positive rational kernels, exact row-multinomial denominators,
small partition sums, and every endpoint repair in its three-vertex case.

The retained-Hadamard checker is

```sh
.venv/bin/python computations/transfer_reconstruction_typical_hadamard_exact_checks_2026_09_06.py
```

It checks five fixed selectors in orders four, eight and twelve. All 124
row-spin energies and the five Gram trace-square identities are exact
integer equalities. Fifteen cosine-test cases separately check means,
variances, Gaussian comparisons and pair covariance identities; these
transcendental checks are explicitly numerical, not exact certificates.
Both checkers passed. The uniform theorems rest on Sections 1--3, not on
the finite examples.

The pressure floor is an obstruction to making this ensemble's first
moment exponentially small. It remains compatible with domination by
rare realizations or correlated spin clusters. A second-moment or other
correlated-spin argument is a genuinely different next obligation.
