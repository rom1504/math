# Restricted random weaving: an exact soft counting reduction

This is an original-signing upper-construction candidate, not a convergence
theorem. No upper cap is claimed by the numerical lower-witness searches.

## 1. Full weaving and the matching obstruction

Let `m` be an even Hadamard order. Choose arbitrary Hadamard matrices
`H_i` of order `m`, indexed by `i in [m]`, and an arbitrary symmetric
full sign matrix `S` of order `m`. The matrix

\[
B_{(i,a),(j,b)}=S_{ij}H_i(a,j)H_j(b,i)
\tag{1}
\]

is a symmetric Hadamard matrix of order `m^2`. Indeed the inner sum over
`b` in an entry of `B^2` is `m delta_{ij}` by column orthogonality; the
remaining sum uses row orthogonality in `H_i`.

Every fixed-point-free involution `pi` of `[m]` produces Boolean
eigenvectors of **both** signs. Set

\[
x_{i,a}=u_iH_i(a,\pi(i)),\qquad
u_{\pi(i)}=\sigma S_{i,\pi(i)}u_i,
\quad \sigma\in\{-1,1\}.
\]

Then `Bx=sigma m x`. Thus arbitrarily random choices of all `H_i` and
`S` cannot lower the full weave's normalized absolute quadratic cap.
Deleting the final diagonal changes the energy by at most `m^2/2`.
This proof concerns the explicit rank-one block weave (1), not every
construction called weaving in the literature.

## 2. Restriction and the exact sufficient inequality

Choose `T_i subset [m]` of common size `k`, retain only those rows and
columns in the `i`th fibre, and call the resulting full sign matrix `B_T`.
Its order is `N=mk`. For block spins `x_i in {+1,-1}^k`, put

\[
h_i=H_i[T_i,:]^Tx_i,\qquad
E(x)=\sum_{i,j}S_{ij}h_i(j)h_j(i).
\]

Then

\[
\sum_i\|h_i\|^2=m^2k,\qquad E(x)=x^TB_Tx,
\]

and, for either `sigma`,

\[
D_\sigma(x):=\sum_{i,j}
  (h_i(j)-\sigma S_{ij}h_j(i))^2
=2(m^2k-\sigma E(x)).
\tag{2}
\]

In particular, if for some fixed `p in (0,1)`, `k/m -> p`, and `eta>0`,

\[
\min_{x,\sigma}D_\sigma(x)
\ge 2m^2k\{1-\sqrt{k/m}+\eta\},
\tag{3}
\]

then the hollow signing obtained from `B_T` has normalized cap at most
`1/2 - eta/(2 sqrt(p)) + o(1)`. The diagonal payment is at most `N/2`.
This would be a strict asymptotic upper-bound improvement, not by itself
a proof of convergence to any constant.

Typical random restriction, existence of some useful restriction, and
the best Boolean witness of a chosen restriction are separate questions.
Matching eigenvectors from Section 1 are not all possible witnesses.

## 3. Soft Finner bound for arbitrary row spectra

Here is an exact near-saturation counting reduction, including nonplateaued
and remote Boolean states. Fix a base Hadamard `H_0`. Independently choose
uniform column permutations `P_i`, set `H_i=H_0P_i`, and independently choose
the off-diagonal signs of `S`. The row selectors may first be held fixed.

For a fixed block spin `x`, let
`a=|H_0[T,:]^Tx|` be its nonnegative absolute spectrum of length `m`.
For `t>0` define the positive-semidefinite kernel

\[
K_t(a,b)=\tfrac12\left(
 e^{-t(a-b)^2/k}+e^{-t(a+b)^2/k}\right).
\tag{4}
\]

Positive semidefiniteness follows, for example, by expanding
`exp(-t(a^2+b^2)/k) cosh(2tab/k)` in nonnegative even powers.
For each coordinate value `v` occurring in `a`, remove one occurrence
to obtain `a minus v`, a list of length `m-1`. Put

\[
L_t(a)=\max_{v\in a}
\left\{\frac{\operatorname{per}
 K_t[a\setminus v]}{(m-1)!}\right\}^{1/2},
\qquad
Z_T(t)=\sum_{x\in\{\pm1\}^k}L_t(|H_0[T,:]^Tx|).
\tag{5}
\]

The permanent uses repeated rows/columns when values repeat. Its quotient
by `(m-1)!` is the expectation of the product of kernel entries under
a uniform coordinate permutation.

**Theorem.** For any fixed selectors `T_1,...,T_m` and any `gamma>=0`,

\[
\Pr\{\exists x,\sigma:\ D_\sigma(x)\le2\gamma m^2k\}
\le 2e^{t\gamma m^2}\prod_{i=1}^m Z_{T_i}(t).
\tag{6}
\]

If the selectors are also iid from a common distribution, the right-hand
side averages to

\[
2e^{t\gamma m^2}\{\mathbb E_T Z_T(t)\}^{m}.
\tag{7}
\]

### Proof

Fix all block spins, and condition on the absolute diagonal coordinate
chosen by the random column permutation in each row spectrum. Discard
the nonnegative diagonal terms of `D_sigma`. Each undirected off-diagonal
edge has an independent random sign, and its exponential weight after
averaging that sign is exactly (4).

The remaining absolute entries of each row are a uniform permutation of
a fixed multiset. Let `V_i` be the indicator tensor of its distinct
permutations and let `n_i` be their number. The conditional exponential
moment is the complete-graph tensor contraction of the `V_i/n_i`, with
the kernel `K_t` along each edge.

Factor the PSD kernel as `K_t^{1/2} K_t^{1/2}`. Absorb one half into each
endpoint tensor. The generalized Cauchy--Schwarz/Finner inequality for a
graph tensor contraction bounds its absolute value by the product of the
Euclidean norms of the transformed vertex tensors. For a fixed reference
ordering `a'` of the row multiset, the squared norm is

\[
\sum_{b,c\text{ distinct permutations of }a'}
  \prod_jK_t(b_j,c_j)
=n_i^2\frac{\operatorname{per}K_t[a']}{(m-1)!}.
\]

After division by `n_i`, its norm is the square root in (5). Maximizing
over the conditioned diagonal value bounds the unconditional moment.
Therefore

\[
\mathbb E e^{-tD_\sigma(x)/(2k)}
\le\prod_i L_t(|H_0[T_i,:]^Tx_i|).
\]

Markov's inequality, summation over all spins, and the two energy signs
give (6). Independence of the selectors gives (7).

## 4. The first unproved asymptotic obligation

For `k/m -> p`, a sufficient one-row statement for (3) is: there exist
fixed `t>0`, `eta>0` and `epsilon>0` such that

\[
\limsup_{m\to\infty}\frac1m\log\mathbb E_T Z_T(t)
\le -t(1-\sqrt p+\eta)-\epsilon.
\tag{8}
\]

Then (7) tends to zero exponentially in `m^2`. No proof of (8) is claimed.
This is a substantially smaller one-row spectral counting problem, but
it still includes every Boolean spin and every absolute spectral profile.

At `t -> infinity`, the permanent formula reduces to the exact profile
entropy factor: if the retained off-diagonal multiset has multiplicities
`c_a` and `r` nonzero entries, the squared factor is
`2^{-r} prod_a c_a!/(m-1)!`. At fixed positive defect one cannot simply
replace arbitrary profiles by exactly bent or plateaued profiles.

Potapov's exact bent and plateaued counts are relevant to special classes:
[Upper bounds on the numbers of binary plateaued and bent functions,
arXiv:2303.16547v3](https://arxiv.org/html/2303.16547v3).
Theorem 2 has `log_2 #bent <= (11/32+o(1))m`. Its exact divisibility and
algebraic-degree arguments do not automatically supply a stable count
of approximately flat spectra.

## 5. Reproducible finite falsifiers

`computations/continued_convergence_restricted_weave_2026_09_06.py`
constructs the signing and checks the energy and squared-defect identities
in integer arithmetic. Its C++ companion searches only for lower witnesses.
The saved matrices are reproducible from the row selectors, permutations,
seed signs and RNG seed in the result JSON.

`computations/continued_convergence_weave_profile_2026_09_06.py` groups all
one-row spins by absolute spectral profile and evaluates (5) by a positive
dynamic program for a random permutation. No cancellation-prone permanent
formula is used. A positive exponent in (6) is inconclusive; a negative
finite exponent gives only the corresponding finite probabilistic bound.

## 6. Gaussian typical profiles obstruct this certificate for moderate retention

This is a lower bound on the one-row partition sum, **not** an evaluation
of its extensive pressure and not an obstruction to the actual ensemble.
It applies to any sequence of real Hadamard matrices. Let `T` be uniform
among the `k`-subsets, `x` uniform among its sign vectors, and `k/m -> p>0`.
Normalize the absolute spectrum by `sqrt(k)`. For every fixed `t>0`,

\[
\liminf_{m\to\infty}\frac1m\log\mathbb E_T Z_T(t)
\ge p\log2-\frac12 F(t),
\tag{9}
\]

where `rho in (0,1)` is determined by `2t=rho/(1-rho^2)` and

\[
F(t)=2t(1-\rho)-\tfrac12\log(1-\rho^2).
\tag{10}
\]

### Typical-profile law and tail control

The empirical law of `|H_0[T,:]^T x|/sqrt(k)` converges in probability,
in the quadratic-Wasserstein topology, to the half-normal law `mu=|G|`.
For a single column this is the ordinary Rademacher central limit theorem.
For two distinct columns, split the sampled rows according to whether
the column signs agree. Orthogonality makes the two populations equally
large; the sampled split is hypergeometric, concentrates at `k/2`, and
the two independent Rademacher sums give a joint limit of independent
standard Gaussians. This argument is uniform over pairs of columns.
Consequently the variance of every bounded continuous empirical test
function tends to zero. Finally each normalized coordinate has fourth
moment `3-2/k`, so
`E integral_{a>R} a^2 d mu_m <= 3/R^2`. This gives the claimed tail
control and quadratic-Wasserstein convergence in probability.

Here and below a diagonal coordinate of bounded normalized size can be
removed: one exists by the exact empirical second moment `1`. Removing
it does not change the limiting law or quadratic moment.

For completeness, unbounded amplitudes cause no hidden uniform-entry
assumption in the permanent lower bound. Write `a_R=min(a,R)` and
`e_a=(a-R)_+`. For every `epsilon>0`, the elementary squared-sum bound
gives

\[
K_t(a,b)\ge
 e^{-2t(1+1/\epsilon)(e_a^2+e_b^2)}
 K_{t(1+\epsilon)}(a_R,b_R),
\tag{11}
\]

where in this normalized notation
`K_t(a,b)=exp(-t(a^2+b^2)) cosh(2tab)`.
Every permutation product charges the tail factor exactly twice, hence
the normalized logarithm of the permanent pays at most
`4t(1+1/epsilon) integral e_a^2 d mu_m`.
After clamping, a fixed fine finite partition and the lower method of
types for permutation contingency tables apply on a compact interval.
One takes `m -> infinity`, then the bin width to zero and `R -> infinity`,
and finally `epsilon -> 0`. Thus every coupling `pi` of `mu` with itself
having finite entropy supplies the lower bound

\[
\liminf\frac1m\log
 \frac{\operatorname{per}K_t[a\setminus v]}{(m-1)!}
\ge \int\log K_t\,d\pi-D(\pi\Vert\mu\otimes\mu)
\tag{12}
\]

along any sequence of such typical profiles. In (12), clamping the
coupling cannot increase its relative entropy; its kernel integral
converges by the quadratic tail control. The fixed-bin method of types
uses contingency tables with the actual empirical margins, approximating
the chosen limiting coupling; rounding changes only `o(m)` entries.

### An explicit Gaussian coupling

Let `pi_rho` be the law of `(|G_1|,|G_2|)` for correlated standard
Gaussians with correlation `rho`. At `2t=rho/(1-rho^2)`, its density is

\[
\frac{d\pi_\rho}{d(\mu\otimes\mu)}(a,b)
 =(1-\rho^2)^{-1/2}
 e^{t(1-\rho)(a^2+b^2)} K_t(a,b).
\]

Substitution in (12) gives `-F(t)`, since both second moments are `1`.
In fact this is the optimizing coupling: the displayed density gives a
dual certificate with a sum of two one-coordinate potentials, and the
remaining discrepancy for any coupling is its nonnegative relative
entropy with respect to `pi_rho`. Only the lower bound is needed here.
The square root in `L_t` costs a factor `1/2`; a set of typical rows has
probability tending to one. Since
`E_T Z_T(t)=2^k E_{T,x} L_t(a)`, this proves (9).

The half-bound first-moment exponent is therefore at least

\[
g_p(t)=p\log2-\tfrac12 F(t)+t(1-\sqrt p).
\]

Using `F'(t)=2(1-rho)`, its minimum occurs at `rho=sqrt(p)`, equivalently
`t=sqrt(p)/(2(1-p))`, and is exactly

\[
\inf_{t>0}g_p(t)=p\log2+\tfrac14\log(1-p).
\tag{13}
\]

For `p=7/8` this is `0.08664339757...>0`; therefore (8) is impossible
at that retention fraction for any fixed `t`, even with `eta=0`.
For `p=15/16` it is `-0.04332169878...`; the typical Gaussian class does
not falsify the criterion there. A negative value is **not** an upper
bound on the actual partition sum: exponentially rare, structured rows
could dominate it. Precisely those rows remain the obligation.

## 7. Annealing across the square root is universally too costly

This obstruction applies to a tempting simplification when randomizing
the base Hadamard. It does not obstruct the square-root quantity in (5).
Take any order-m Hadamard `F`, multiply its input coordinates by iid
random signs `D`, and use a uniform k-selector `T`. Extra randomization
of the Hadamard is allowed. For a fixed retained input spin, let

\[
P_{F,D,T}(t)=\max_v
 \frac{\operatorname{per}K_t[a\setminus v]}{(m-1)!},
\qquad a=|F D(1_Tx)|/\sqrt k.
\]

Here the kernel is in normalized notation. Input sign invariance gives
the one-row sum `2^k E sqrt(P)`. Jensen would replace this by

\[
J_m(t)=2^k\sqrt{\mathbb E P_{F,D,T}(t)}.
\tag{14}
\]

Write `delta=1-p`. For every fixed `t>0`, the actual matching event gives

\[
\liminf\frac1m\log J_m(t)
\ge \frac p2\log2-\frac12F(\delta t).
\tag{15}
\]

To see this, fix one row `r` of `F` and require the input signs on `T`
to agree with that row. This event has probability exactly `2^{-k}`,
independently of `T`. The `r`th output coordinate is `sqrt(k)` and may
be removed in the maximum defining `P`. Every other coordinate is the
normalized sum of a uniform k-subset from a balanced sign column.
For two such distinct coordinates the four sign-pair populations have
size `m/4`, by orthogonality. The univariate and bivariate finite-population
central limit theorems therefore give a typical half-normal profile of
variance `delta`. The retained squared norm is exactly `m-k`.
The fourth-moment bound and the clamping proof from Section 6 apply as
well. Thus on a set of selectors of probability tending to one,
`P >= exp{-m F(delta t)-o(m)}`. Multiplication by the probability
`2^{-k}` and then the square root prove (15), uniformly in the choice
of the underlying Hadamard.

At the half-cap target the Jensen first-moment exponent is therefore
at least

\[
\frac p2\log2-\frac12F(\delta t)+t(1-\sqrt p).
\]

Its infimum over fixed `t>0` is

\[
\frac p2\log2+
 \frac14\log\!\left(1-\frac{p}{(1+\sqrt p)^2}\right)>0
\quad (0<p<1).
\tag{16}
\]

Indeed set `rho=sqrt(p)/(1+sqrt(p))`; the optimum satisfies
`2 delta t=rho/(1-rho^2)` and the linear terms cancel. Positivity follows
from `log(1-u)>=-u/(1-u)` and
`u/(1-u)=p/(1+2sqrt(p))`: (16) is strictly larger than
`p(log(2)/2-1/4)>0`.
As `p -> 1`, (16) tends to
`log(2)/2+log(3/4)/4=0.2746530721...`.

Consequently Jensen across the square root cannot prove the half cap
at **any** fixed retention fraction, even if the underlying Hadamard
ensemble is otherwise ideal. The matching event is rare enough to be
harmless before Jensen, but the square root outside its probability
forfeits half of the input-spin entropy cancellation. This statement
concerns uniform selectors; conditioning selectors changes the expression
and requires a separate analysis.
