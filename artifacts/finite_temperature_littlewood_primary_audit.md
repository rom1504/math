# Finite-temperature and quadratic-Littlewood primary-source audit

Status: **primary-literature audit, 2026-08-14**.  The mappings below are
exact; the final verdict is negative.  No checked theorem supplies a
thermodynamic limit for the minimum over deterministic signings, and no checked
tensor/product theorem has `1+o(1)` loss in the regime relevant here: fixed
degree two, complete squarefree support, and number of variables tending to
infinity.

This note complements the internal derivations in
`soft_cap_composition_audit.md`, `finite_temperature_scalar_no_go.md`, and
`unimodular_walsh_sidon_literature.md`.  It does not claim that the two proposed
abstractions are useless; it records exactly what the closest published
theorems do and do not prove.

## 1. Exact finite-temperature object versus the SK object

For a signing `A` of `K_n`, put

$$
H_A(x)=\sum_{i<j}a_{ij}x_ix_j,
\qquad
Z_A(\beta)=\sum_{x\in\{\pm1\}^n}
 \exp\!\left(\frac{\beta H_A(x)}{\sqrt n}\right).
$$

The proposed soft absolute pressure is

$$
\Phi_n(\beta)=\frac1n\min_A
 \log\sum_x2\cosh\!\left(\frac{\beta H_A(x)}{\sqrt n}\right)
=\frac1n\min_A\log\bigl(Z_A(\beta)+Z_{-A}(\beta)\bigr).
\tag{1.1}
$$

Thus this is not the usual quenched SK pressure.  The differences are:

1. the disorder is minimized rather than sampled and averaged;
2. the same deterministic signing is used inside the logarithm;
3. `A` and `-A` are symmetrized to detect the absolute ground state.

The last difference is harmless for typical symmetric random disorder, since

$$
\max\{\log Z_A,\log Z_{-A}\}
\le \log(Z_A+Z_{-A})
\le \max\{\log Z_A,\log Z_{-A}\}+\log2.
\tag{1.2}
$$

The first difference is decisive: a theorem for
`E_A log Z_A` or for a typical `A` gives no control of `min_A log Z_A`.

### What the rigorous SK theorems control

Guerra and Toninelli,
[“The Thermodynamic Limit in Mean Field Spin Glass
Models”](https://arxiv.org/abs/cond-mat/0204280), prove convergence of the
**quenched averaged** free energy for Gaussian SK and related Gaussian
mean-field models.  Their interpolation joins a size-`N` Gaussian Hamiltonian
to two independent Gaussian subsystem Hamiltonians.  Gaussian covariance
interpolation gives the sign needed for subadditivity.

That mechanism is not an interpolation of (1.1).  In the SK interpolation,
coefficients can be multiplied by square roots of block proportions while
remaining Gaussian.  A deterministic coefficient constrained to be `+1` or
`-1` cannot be rescaled and remain admissible.  At the project normalization,
a split `N=m+n` changes a child's scaled inverse temperature from `beta` to

$$
\beta\sqrt{m/N}
\quad\hbox{or}\quad
\beta\sqrt{n/N},
\tag{1.3}
$$

which is exactly the contraction already derived internally.  Moreover,
quenched expectation does not commute with the adversarial minimum.

Carmona and Hu,
[“Universality in Sherrington--Kirkpatrick's Spin Glass
Model”](https://arxiv.org/abs/math/0403359), show that the limiting quenched
free energy is unchanged when Gaussian disorder is replaced by centered,
unit-variance independent disorder with a finite third moment.  Their theorem
therefore covers Bernoulli disorder (up to the paper's ordered-pair versus
upper-triangle normalization convention).  It proves that a **uniform random**
signing has the usual SK limit.  It does not address the extreme left endpoint
among all `2^(n choose 2)` signings.

Parisi and Rizzo's
[2007](https://arxiv.org/abs/0706.1180) and
[2008](https://arxiv.org/abs/0811.1524) papers study disorder large deviations
of the SK free energy using replicas.  They are physics calculations for
Gaussian disorder, not a rigorous large-deviation theorem for Bernoulli
signings, and they do not analyze the symmetrized minimum (1.1).  They cannot
be imported as a proof.

A recent deterministic result of Fan, Misiakiewicz, Wang, and Wen,
[“Dynamical mean-field limit and replica-symmetric free energy for the
orthogonally-invariant SK model”](https://arxiv.org/abs/2607.10102), proves a
free-energy limit under rapid mixing and, for Ising spins, the explicit
high-temperature condition `||X||_op < 1/2`, with deterministic extensions
under delocalization hypotheses.  The exact interaction matrix for (1.1) is

$$
X=\frac{\beta A}{\sqrt n},
\qquad
\frac12x^{\mathsf T}Xx=\frac{\beta H_A(x)}{\sqrt n}.
$$

Every signing satisfies

$$
\|A\|_{\rm op}
\ge \frac{\|A\|_{\rm F}}{\sqrt n}
=\sqrt{n-1},
$$

so the theorem's condition can hold only if
`beta sqrt(1-1/n) < 1/2`.  The ground-state squeeze requires fixed `beta` as
large as desired.  This theorem therefore cannot supply the required family of
limits, even if its delocalization hypotheses were verified for minimizers.

### The missing disorder-scale theorem

Let `A_n` be uniform over all edge signings and define

$$
Y_{n,\beta}(A)=\frac1n\log
 \sum_x2\cosh\!\left(\frac{\beta H_A(x)}{\sqrt n}\right),
\qquad E_n=\binom n2.
$$

Then the desired finite-temperature value is exactly the left endpoint of the
finite support:

$$
\Phi_n(\beta)=\min\operatorname{supp}Y_{n,\beta}.
\tag{1.4}
$$

Equivalently, if

$$
N_{n,\beta}(s)=
\#\{A:Y_{n,\beta}(A)\le s\},
$$

then

$$
\Pr(Y_{n,\beta}\le s)=2^{-E_n}N_{n,\beta}(s),
\qquad
\Phi_n(\beta)\le s\iff N_{n,\beta}(s)\ge1.
\tag{1.5}
$$

Thus the external theorem that would genuinely remove an obligation is not a
usual SK thermodynamic limit.  It would be a **support-sensitive disorder LDP**
at speed `E_n = Theta(n^2)`, or equivalently convergence of

$$
\Sigma_{n,\beta}(s)=E_n^{-1}\log N_{n,\beta}(s)
\tag{1.6}
$$

where `log 0=-infinity`, with enough uniformity to locate the transition from
no support to nonempty support.  For example, a unique `s_beta` for which the
limiting entropy is `-infinity` below `s_beta` and finite above `s_beta`, with
uniform one-sided bounds, would imply convergence of `Phi_n(beta)`.
Typical SK universality works at the quenched value and does not provide
(1.6).  No rigorous Bernoulli-disorder theorem of this endpoint strength was
found in the checked sources.

## 2. Exact quadratic-Littlewood formulation

Let

$$
\Lambda_{2,n}=\{x_ix_j:1\le i<j\le n\},
\qquad L_n=|\Lambda_{2,n}|=\binom n2.
$$

Our quantity is exactly the minimum uniform norm of a complete-support real
quadratic Littlewood polynomial:

$$
M_n=\min_{a_{ij}\in\{\pm1\}}
\left\|\sum_{i<j}a_{ij}x_ix_j\right\|_{L_\infty(\{\pm1\}^n)}.
\tag{2.1}
$$

Equivalently, for the equal-modulus Sidon constant

$$
U_n=\max_{a_{ij}\in\{\pm1\}}
\frac{L_n}{\left\|\sum_{i<j}a_{ij}x_ix_j\right\|_\infty},
$$

we have `U_n=L_n/M_n`.  The original convergence question is exactly
convergence of `U_n/sqrt(n)`.

This is narrower than all of the following nearby objects:

- the ordinary Sidon constant, which optimizes over arbitrary coefficient
  magnitudes;
- a bilinear Littlewood form, which has independent left and right Boolean
  inputs;
- a complex polynomial on the full torus;
- a quadratic phase `(-1)^p` over `F_2`, whose Fourier coefficients are not
  the coefficients of (2.1).

## 3. What the closest Littlewood/BH theorems control

### Boolean Bohnenblust--Hille

Defant, Mastyło, and Pérez,
[“On the Fourier spectrum of functions on Boolean
cubes”](https://arxiv.org/abs/1706.03670), prove a dimension-free Boolean
Bohnenblust--Hille inequality.  At degree two it has the form

$$
\left(\sum_{i<j}|a_{ij}|^{4/3}\right)^{3/4}
\le B^{\rm Bool}_2
\left\|\sum_{i<j}a_{ij}x_ix_j\right\|_\infty.
\tag{3.1}
$$

For flat complete support this yields only

$$
M_n\ge \frac{\binom n2^{3/4}}{B^{\rm Bool}_2}
=\left(\frac{2^{-3/4}}{B^{\rm Bool}_2}+o(1)\right)n^{3/2}.
\tag{3.2}
$$

It is a within-order lower bound.  It gives neither a relation between
`M_m` and `M_n` nor convergence of the normalized minima.

The 2026 support-sensitive extension by Defant, Galicer, Mansilla, Mastyło,
and Muro,
[“Support-Sensitive Bohnenblust--Hille Inequalities and Local Invariants on
Hamming Schemes”](https://arxiv.org/abs/2607.05594), compares Sidon,
unconditional, and Gordon--Lewis constants up to factors depending on the
interaction order.  For the Boolean squarefree level, the interaction order
is two.  The comparison factor is dimension-free but fixed; it is not
`1+o(1)` as `n` tends to infinity and it does not preserve the nonlinear flat
coefficient slice.

Caro Montoya, Núñez Alarcón, and Serrano Rodríguez,
[“Asymptotic contractivity of the Bohnenblust--Hille inequality for
polynomials with few interacting
variables”](https://arxiv.org/abs/2607.20847), prove

$$
1\le K_{m,M}
\le A_M^{M/m}m^{(M^2-1)/(2m)},
\qquad K_{m,M}\longrightarrow1
$$

for **fixed support size `M` and homogeneous degree `m -> infinity`**.  Our
regime is homogeneous degree `m=2`, support size two per monomial, and
ambient dimension `n -> infinity`.  The theorem's `1+o(1)` therefore occurs
in the wrong asymptotic parameter and supplies no estimate tending to one in
our problem.

### Ordinary and equal-modulus Sidon estimates

Defant, Galicer, Mansilla, Mastyło, and Muro,
[“Asymptotic insights for projection, Gordon--Lewis and Sidon constants in
Boolean cube function spaces”](https://arxiv.org/abs/2302.00233), determine
orders and some projection-constant limits for Boolean polynomial spaces.
For the weight-two space they control the ordinary Sidon invariant only up to
degree-dependent constants.  Their exact projection formula concerns the
average absolute all-positive reproducing kernel, not the minimum in (2.1).
No theorem there proves

$$
\operatorname{Sid}(\Lambda_{2,n})=(1+o(1))U_n
$$

or a cross-order inequality for `U_n`.

Volberg,
[“An estimate of Sidon constant for complex polynomials with unimodular
coefficients”](https://arxiv.org/abs/2205.04936), does treat equal-modulus
coefficients.  Its domain is the complex polytorus, its support consists of
all homogeneous monomials (including repetitions), and its main asymptotic is
in the degree.  Since

$$
\sup_{z\in\mathbb T^n}|P(z)|
\ge \sup_{x\in\{\pm1\}^n}|P(x)|,
$$

a lower bound for the torus norm does not lower-bound our Boolean norm.  The
paper gives no complete-squarefree, fixed-degree, cross-`n` theorem.

Astashkin and Lykov,
[“Random unconditional convergence of Rademacher chaos in `L_infinity` and
sharp estimates for discrepancy of weighted graphs and
hypergraphs”](https://arxiv.org/abs/2412.20107), are the closest source to
(2.1): their weighted-graph theorem directly includes the minimum over edge
signs.  It identifies the correct `Theta(n^(3/2))` order through row
`ell_2` functionals, but only up to fixed universal comparison constants.
It asserts neither a ratio tending to one nor a recurrence between orders.

### Quadratic phases are a different Fourier object

Becker, Slote, Volberg, and Zhang,
[“Fourier growth of degree 2
polynomials”](https://arxiv.org/abs/2412.10842), prove the sharp estimate

$$
\sum_{|S|=k}|\widehat f(S)|
\le \sqrt{\frac2\pi}\,k^{-1/2}(1+\sqrt2)^k
\tag{3.3}
$$

for `f(y)=(-1)^{p(y)}` with `deg_F2 p <= 2`.

There is an exact encoding of an edge signing into such a phase.  Set

$$
b_{ij}=\frac{1-a_{ij}}2\in\mathbb F_2,
\qquad
p_A(y)=\sum_{i<j}b_{ij}y_iy_j,
\qquad f_A=(-1)^{p_A}.
$$

Then the edge sign is the constant second multiplicative derivative

$$
f_A(y)f_A(y+e_i)f_A(y+e_j)f_A(y+e_i+e_j)=a_{ij}.
\tag{3.4}
$$

But (3.3) controls the Walsh coefficients of the **phase** `f_A`, whereas
(2.1) controls the supremum of the real polynomial formed from the phase's
second derivatives.  Every signing has an encoding (3.4), so the universal
bound (3.3) does not distinguish small-cap signings.  No inequality in the
paper converts its level weights to `M(A)` with a sharp leading constant or
relates different orders.

## 4. Tensor/product audit

There are genuine `1+o(1)` results nearby, but none remains `1+o(1)` after
mapping to (2.1).

Pellegrino and Raposo,
[“Constants of the Kahane--Salem--Zygmund inequality asymptotically bounded
by 1”](https://arxiv.org/abs/2006.12892), construct unimodular `d`-linear
forms on

$$
\ell_\infty^{n_1}\times\cdots\times\ell_\infty^{n_d}
$$

whose KSZ constant is below `1+epsilon` when all side lengths are large.
Albuquerque and Rezende,
[“Asymptotic estimates for unimodular multilinear forms with small norms on
sequence spaces”](https://arxiv.org/abs/1710.09711), determine the order of
the corresponding rectangular multilinear minimum.

At degree two these results control

$$
\min_{a_{ij}=\pm1}\max_{x,y\in\{\pm1\}^n}
\left|\sum_{i,j}a_{ij}x_i y_j\right|,
\tag{4.1}
$$

with independent left and right spins.  Our norm puts the same spin on both
sides, imposes symmetry and zero diagonal, and counts each unordered edge
once.  Polarization/decoupling from (4.1) to (2.1) has a fixed leading-order
loss.  Consequently the `1+o(1)` KSZ constant is lost before reaching our
quantity.

The usual algebraic products also leave the class:

- multiplying two quadratic Littlewood polynomials produces degree four;
- a tensor product of their coefficient arrays is supported on
  `edge x edge`, not on all edges of a larger complete graph;
- a direct sum preserves degree two but omits every bridge edge;
- bilinearizing the bridge restores independent left/right spins and hence
  the same fixed polarization loss.

The sharp bilinear-extremizer analysis of Pellegrino and Teixeira,
[“Towards sharp Bohnenblust--Hille
constants”](https://arxiv.org/abs/1604.07595), is likewise for multilinear
forms with independent inputs.  It does not furnish a same-spin tensor law.

**Conclusion of the tensor search.** Among the checked primary sources, no
tensor, product, Sidon, BH, KSZ, or RUC theorem gives

$$
M_{mn}\le (1+o(1))\,\mathcal C(M_m,M_n)
$$

for a complete-support quadratic signing, nor any additive composition with
an `o((m+n)^(3/2))` leading loss.  Every available `1+o(1)` statement either
takes degree to infinity, uses independent input channels, moves to the
complex torus, or leaves complete squarefree support.

## 5. Research judgment

The finite-temperature reduction itself is valid, but the closest rigorous
SK theorems average over disorder and therefore address a different
quantifier.  The exact missing probabilistic object is the disorder-counting
entropy (1.6), not the usual quenched pressure.  Establishing a sharp
`Theta(n^2)` Bernoulli-disorder endpoint LDP would be genuinely new and would
remove a real convergence obligation.

The Littlewood formulation is also exact and useful as search vocabulary,
but current BH/Sidon theory supplies order estimates or fixed-factor
comparisons.  It contains no asymptotically lossless cross-order mechanism
for the flat complete quadratic slice.  In particular, no checked
tensor/product inequality has `1+o(1)` loss in the required regime.
