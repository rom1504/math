# Biased Gram--Schmidt rounding and cumulative-precision augmentation

2026-09-17 paper-combination campaign. **Sections 3--7 independently audited
PASS by the campaign's BH/mechanism agent, including the primary sharp-GS
proof and exact 9/8 constant.** No claim of convergence or improvement to the
reported global cap interval is made.

**Applicability correction (18:34 UTC).** The complete near-level code
of every bounded positive-cap quadratic sequence necessarily has
entropy at least `(1/(4c)+o(1)) eta log(1/eta)`. Therefore **finite
unweighted full-code entropy slope is impossible**, not an unproved
optimizer hypothesis. The direct entropy-slope consequences quoted
below are formally valid but have no such examples. Likewise the
Section 9 finite `e s/eta` regime with `e_0>0` is impossible. Weighted
factorization with `e_0=0` and low-response-center Hamming-cover
criteria are not ruled out. See the exact theorem in
[the near-level entropy correction](paper_discrepancy_nearlevel_entropy_floor_2026_09_17.md).
The finite GS, information, response, and transfer lemmas remain valid.

## 1. Outcome, normalization, and archive collision

Write

```math
H_A(x)=\sum_{i<j}A_{ij}x_ix_j,\qquad
Q(A)=\max_{x\in\{\pm1\}^n}|H_A(x)|,\qquad
M_n=\min_A Q(A).
```

All matrices minimized over are hollow symmetric with every off-diagonal
entry in `{-1,1}`. Thus `x^T A x=2H_A(x)`. Let
`c_*=liminf_n M_n/n^(3/2)`.

The primary-paper upgrade is that the original biased Gram--Schmidt walk
has subGaussian variance parameter **1**, not 40. More importantly, a
new use of its feature representation combines *all ancestral near-level
factorizations for the same query*. The reciprocal variances add. The
existing archive used only one level's representation at a time.

Here is the resulting actual extension theorem. For any sequence of full
signings with `Q(A_n)/n^(3/2)->c<infinity`, define the absolute superlevel
code and its factorization complexity by

```math
E_n(\eta)=\{x:Q(A_n)-|H_{A_n}(x)|\le\eta n^{3/2}\},
\qquad e(\eta)=\limsup_n\frac{\gamma_{\rm fac}(E_n(\eta))^2}{n}.
```

Here `gamma_fac` is the **matrix factorization norm** usually denoted
`gamma_2`: a sign-row matrix factors as `LV`, with column norms of `V`
at most one and row norms of `L` at most `gamma_fac`. It is **not**
Talagrand's chaining functional with the same customary notation.
Let `h(u)=-u log u-(1-u)log(1-u)` and
`H(u)=h(min(u,1/2))`. Put

```math
K=\limsup_{\eta\downarrow0}\frac{e(\eta)H(e(\eta))}{\eta}.
```

**Cumulative-precision augmentation theorem.** If `K<infinity`, then for
every `tau>9K/8` and every sufficiently small fixed `epsilon>0`, there
exist actual full-sign extensions `B_(n+q)`, `q=floor(epsilon n)`, with
old principal block exactly `A_n`, such that

```math
\limsup_n\frac{Q(B_{n+q})-Q(A_n)}{n^{3/2}}
\le \tau\epsilon+\epsilon^{3/2}.                 \tag{1}
```

Consequently, every sequence realizing the global liminf satisfies

```math
\boxed{\displaystyle
\limsup_{\eta\downarrow0}\frac{e(\eta)H(e(\eta))}{\eta}
\ge\frac43c_* .}                                \tag{2}
```

The exact prior collision is
[`principle_invent_2026_09_07_hierarchy_complexity_dilution.md`](principle_invent_2026_09_07_hierarchy_complexity_dilution.md),
especially Sections 2 and 7, independently audited in the corresponding
`principle_construct_2026_09_07_hierarchy_complexity_audit.md`. Its coefficient
is `20 phi^5=221.803398874989...`. Merely importing variance parameter 1
would change that to `phi^5/2=5.54508497187474...`. Cumulative precision
further reduces it to `9/8=1.125`. Thus the mechanism changes the range
of quantitative complexity slopes permitting profitable dilution; it is
not a new name for the old one-level construction. The small-energy
*exponent* and the unproved need to understand actual minimizers remain
unchanged.

## 2. Primary sources and what was reconstructed

1. N. Bansal, D. Dadush, S. Garg, S. Lovett, *The Gram--Schmidt Walk: A
   Cure for the Banaszczyk Blues*, Theory of Computing 15(21), 2019,
   [published primary paper](https://theoryofcomputing.org/articles/v015a021/v015a021.pdf),
   [arXiv source record](https://arxiv.org/abs/1708.01079).
2. C. Harshaw, F. Savje, D. Spielman, P. Zhang, *Balancing Covariates in
   Randomized Experiments with the Gram--Schmidt Walk Design*, JASA 119,
   2024, [paper and supplement](https://arxiv.org/abs/1911.03071),
   [publication](https://doi.org/10.1080/01621459.2023.2285474).
   The relevant proof is supplement S3.5, Lemmas S3.5--S3.6 and Theorem
   6.6*. Supplement S8.1 explicitly confirms arbitrary initial bias.

The 2019 paper is CC-BY; the derivations below use our notation and retain
attribution. Local primary PDFs are under
`tmp/paper_portfolio_2026_09_17/discrepancy/`. The initial mapping was
frozen there in `first_mapping.md` before old discrepancy artifacts were
consulted.

A newer arXiv note, Bednorz--Godlewski
[arXiv:2404.03534](https://arxiv.org/abs/2404.03534), retains an alternate
random compensator involving sums of absolute orthogonal-block inner
products. It was inspected but is **not an imported premise**: its v1
HTML contains missing or mistyped inequalities, and the 2024 JASA
supplement already supplies the sharp constant needed here.

## 3. The walk, biased means, and the original 40 proof

Let `V=[v_1 ... v_m]`, with `||v_i||<=1`, and start at any
`b in [-1,1]^m`. At a state `z`, retain the active coordinates
`|z_i|<1`, choose the last active pivot `p`, and choose `u_p=1`, with
`u` supported on active coordinates, minimizing `||Vu||`. Then

```math
Vu=\Pi_{\operatorname{span}\{v_i:i\ne p\text{ active}\}^{\perp}}v_p.
```

Move by either endpoint `delta_-<0<delta_+` of the maximal feasible
line segment in the cube, with probabilities making `E[delta|past]=0`.
At least one coordinate freezes per step. Thus the output `s` has every
coordinate in `{-1,1}`, preserves already integral coordinates, and
`E s=b`. The coordinate martingale, not symmetry, supplies unbiasedness.

For completeness, the original subGaussian proof is not merely invoked.
Fix a test direction `theta`, and partition the walk into pivot phases.
The nonpivot spans decrease. Their successive orthogonal differences
decompose the changing residual of a pivot. For every time interval in
one phase, the sum of the step sizes is a difference of two values of the
pivot coordinate and has absolute value at most two. Summing by parts
against those orthogonal residual increments, then using Cauchy--Schwarz,
bounds that interval's discrepancy by `2||P_phase theta||`. Distinct
phase subspaces are orthogonal.

Call a time good when the currently exposed part of the phase has
`theta`-projection norm at most `1/8`. Once bad, the rest of the phase
stays bad. Each bad phase has `||P_phase theta||>1/8`, so its bad suffix
contributes at most `16||P_phase theta||^2`. All bad times together cost
at most `16||theta||^2` deterministically.

For good times let `Y` be the accumulated discrepancy. If `R_t` is the
active span and `W_t` the nonpivot span, use the beginning/end potentials

```math
\Phi_t=\|\Pi_{W_t}\theta\|^2
 +(1-z_p^2)\|\Pi_{R_t\cap W_t^\perp}\theta\|^2
```

on good times; omit the second term on bad times. Projection nesting
shows that a transition between steps never increases this potential
before the next random move. Set `Z=Y+4Phi`. Its one-step increment is
bounded above by

```math
X=\delta a-4\delta(\delta+2z_p)b_\theta^2,
\quad |a|\le b_\theta\le1/8.
```

Feasibility gives `|delta|<=2`, `|delta+2z_p|<=2`, and
`|delta(delta+2z_p)|<=1`. Hence `|X|<=1`,
`E[X|past]=-4 E[delta^2|past] b_theta^2`, and
`X^2<=4 delta^2 b_theta^2`. The elementary inequality
`exp(t)<=1+t+t^2` for `t<=1` therefore makes `exp(Z)` a
supermartingale after the free projection drops. Its initial value is
`exp(4||theta||^2)`. Adding the bad-time cost yields

```math
\mathbb E e^{\langle\theta,V(s-b)\rangle}
\le e^{20\|\theta\|^2}.
```

This reproduces the variance parameter 40 and identifies exactly the
phase geometry improved in the later proof.

## 4. The sharp biased Bellman proof and its random anisotropy

For a centered sign of mean `x in [-1,1]`, write

```math
\psi_x(t)=\log(\cosh t+x\sinh t)-xt.
```

Direct differentiation gives `psi_x'(0)=0` and
`0<=psi_x''(t)<=1`; consequently `|psi_x'(t)|<=|t|` for either sign
of `t`. Define a scalar function on positive semidefinite 2-by-2 matrices

```math
F_x\!\begin{pmatrix}a&t\\t&b\end{pmatrix}
 =\exp\{\psi_x(t)-ab/2\}.
```

Its negative logarithmic matrix gradient, multiplied by two, is

```math
\begin{pmatrix}b&-\psi_x'(t)\\-\psi_x'(t)&a\end{pmatrix}\succeq0,
```

because `ab>=t^2>=psi_x'(t)^2`. Thus adding a positive semidefinite
matrix can only decrease `F_x`.

Now fix a pivot phase. The orthonormal vectors exposed in that phase
have coefficients `alpha_r=<w_r,v_p>` and `beta_r=<w_r,theta>`.
Newly exposed directions add the PSD matrices
`(alpha_r,beta_r)(alpha_r,beta_r)^T` to their 2-by-2 Gram matrix.
Conditioned at any intermediate time, the *remaining total pivot move*
is a centered sign of mean the current pivot coordinate: it has the two
values `1-z_p` and `-1-z_p`. This is true for arbitrary initial bias.

Condition backwards within the phase. The preceding PSD monotonicity
allows the Gram contributions of future exposures to be dropped one
block at a time from the conditional exponential expectation. At the
initial block the same inequality compares to `F_x(0)=1`. Therefore

```math
\mathbb E\left[
 \exp\left\{D_p(\theta)
 -\tfrac12\|P_pv_p\|^2\|P_p\theta\|^2\right\}
 \mid\text{phase-start history}\right]\le1.        \tag{3}
```

Here `D_p` is the phase's scalar discrepancy and `P_p` is the random
projector onto its exposed orthogonal subspace. Padding finished phases
with zero moves makes all backwards inductions finite and deterministic
in length. Iterated conditional expectations over successive phases give

```math
\mathbb E\exp\left\{
\langle\theta,V(s-b)\rangle-\tfrac12\theta^T\Gamma\theta
\right\}\le1,
\qquad
\Gamma=\sum_p\|P_pv_p\|^2P_p\preceq P_{\operatorname{ran}V}.
                                                               \tag{4}
```

The projectors are mutually orthogonal and `||v_p||<=1`; hence

```math
\boxed{\mathbb E e^{\langle\theta,V(s-b)\rangle}
 \le e^{\|\theta\|^2/2}.}                         \tag{5}
```

Differentiating (4) at zero also gives
`Cov(V(s-b)) <= E Gamma`. But **one may not replace Gamma by E Gamma
inside (4)**. The random compensator is a potential route to finer
certificates, not an already proved deterministic covariance proxy.

## 5. A finite biased, all-query completion lemma

Let nested sign-row codes be `E_0 superset ... superset E_J`, with
factorizations `E_j=L_j V_j`, column norms at most one and row norms at
most `sqrt(n g_j)`. Let `a_j>=0` and `theta>0` obey
`theta+sum_j a_j^2<=1`. Feed the walk the vectors

```math
w_i=(\sqrt\theta\,e_i,\ a_0(V_0)_{:i},\ldots,a_J(V_J)_{:i}).
```

For `x in E_j`, all levels `k<=j` represent the **same vector x**.
For arbitrary real coefficients `c_k` with sum one, the test vector with
blocks `c_k(L_k)_x/a_k` has inner product `z dot x` with the discrepancy
vector. Optimizing its squared norm gives

```math
\min_{\sum c_k=1}\ n\sum_{k\le j}c_k^2g_k/a_k^2
 =\frac{n}{\displaystyle\sum_{k\le j}a_k^2/g_k}.   \tag{6}
```

Zero-weight blocks are omitted. This is the cumulative precision absent
from the archived argument; all test vectors remain deterministic.

In particular, given required precisions `b_j>0`, set

```math
B_j=\max_{0\le k\le j}b_k,\quad B_{-1}=0,
\qquad a_j^2=g_j(B_j-B_{j-1}).                     \tag{7}
```

If `sum_j g_j(B_j-B_(j-1))<1`, an identity weight can fill part of the
remaining budget. A query in `E_j` then has subGaussian variance at most
`n/B_j<=n/b_j`.

This remains true for a prescribed fractional mean matrix
`M in [-1,1]^(n by q)`: run independent walks from its `q` columns,
producing an actual sign matrix `C` with `E C=M`. For every old word
`x in E_j` and new word `y in {+-1}^q`,

```math
\mathbb E\exp\{\lambda x^T(C-M)y\}
\le\exp\{\lambda^2 qn/(2B_j)\}.                  \tag{8}
```

No restriction on `x,y` is hidden. Outside the codes, the identity
feature gives variance `qn/theta`. A union bound can charge each query
its own deterministic slack, including the offsets `x^T M y` and the
old/new scalar energies. Thus (8) is a biased sign-completion theorem;
the extension theorem below uses `M=0`. “Full sign” means all entries
are signs; no claim that the output law assigns positive mass to every
sign matrix is needed or made.

## 6. Proof of the 9/8 resource bound

Assume `K<infinity`; choose `K<kappa_0<kappa` and a small fixed
`eta_0>0` such that

```math
e(\eta)H(e(\eta))\le\kappa_0\eta
\quad(0<\eta\le\eta_0).
```

Choose `eta_0` smaller if necessary. Let the smooth majorant `g(eta)`
be the unique solution in `(0,1/2)` of

```math
g(\eta)h(g(\eta))=\kappa\eta.                    \tag{9}
```

Then `e(eta)<g(eta)`. Each fixed level therefore has the needed
factorization for all sufficiently large `n`. A code with
`gamma_fac^2<=gn` has VC dimension at most `gn`: factor the labels on a
shattered coordinate set and average squared signed-column sums. Sauer's
bound gives `log|E|<=n h(g)+o(n)` for `g<1/2`.

Fix a hierarchy ratio `R>1` and slack `zeta>0`. For small fixed
`epsilon`, put `a=tau epsilon`, use levels
`eta_j=eta_0 R^(-j)` down to `eta_J` of order `epsilon^2`, and set
`g_j=g(eta_j)`. For `j<J`, request precision

```math
b_j=\frac{2(1+\zeta)\epsilon
       [h(g_j)+\epsilon\log2]}{(\eta_j/R+a)^2}.
                                                               \tag{10}
```

At the terminal level replace the denominator by `a^2`. Allocate weights
by (7). With these weights, a two-sided union bound over the code,
all `2^q` new words, and the fixed finite hierarchy gives, simultaneously,

```math
\frac{|x^TCy|}{n^{3/2}}
\le \eta_j/R+\tau\epsilon+o_n(1)
\quad(x\in E_j),                                  \tag{11}
```

and `tau epsilon+o_n(1)` on the last code. The strict factor `1+zeta`
absorbs the finite union constants and every `o(n)` entropy term. More
explicitly, use threshold
`sqrt((2qn/B_j) log[4(J+1)|E_j|2^q])`. The limiting threshold is smaller
than the right side of (11) by `1/sqrt(1+zeta)`.

It remains to prove the weight budget. Define the smooth nonterminal
function `b(eta)` by (10). A running maximum satisfies

```math
B_j-B_{j-1}\le(b_j-b_{j-1})_+.
```

Since `g` increases, summing the positive variations yields

```math
\sum_{j<J}g_j(B_j-B_{j-1})
\le g(\eta_0)b(\eta_0)
 +\int_0^{\eta_0}g(\eta)(-b'(\eta))_+\,d\eta.     \tag{12}
```

The terminal contribution is bounded by `g_J b_J=o(1)`: (9) gives a
part of order `eta_J/epsilon`, and the new-word entropy part is of order
`g_J`. The boundary term at `eta_0` is `O(epsilon)`.

Separate the entropy and new-word terms of `b`, using
`(u+v)_+<=u_++v_+`. The new-word part contributes `o(1)`: since
`h(g)>=2(log2)g` on `[0,1/2]`, equation (9) gives
`g(eta)<=sqrt(kappa eta/(2log2))`, and its integral is

```math
O\left(\epsilon^2\int_0^{\eta_0}
\frac{\sqrt\eta}{(\eta/R+a)^3}\,d\eta\right)
=O_{R,\tau,\kappa}(\sqrt\epsilon).
```

For the principal entropy term define

```math
\rho(g)=\frac{g h'(g)}{h(g)+g h'(g)}\longrightarrow\frac12
\quad(g\downarrow0).
```

Differentiating (9) shows
`g(eta) d[h(g(eta))]/deta=kappa rho(g(eta))`. Therefore the principal
integral in (12) is exactly

```math
2(1+\zeta)\epsilon\kappa
\int_0^{\eta_0}
\left[
\frac{2\eta/R}{(\eta/R+a)^3}
-\frac{\rho(g(\eta))}{(\eta/R+a)^2}
\right]_+d\eta.                                  \tag{13}
```

Substitute `eta=Ra u`. Since `rho>=0`, the scaled integrand is dominated
by `2u/(1+u)^3`, an integrable function. Dominated convergence gives

```math
\lim_{\epsilon\downarrow0}(13)
=\frac{2(1+\zeta)\kappa R}{\tau}
\int_0^\infty
\left[\frac{2u}{(1+u)^3}
 -\frac1{2(1+u)^2}\right]_+du
=\frac{9(1+\zeta)\kappa R}{8\tau}.              \tag{14}
```

Indeed the bracket is positive exactly for `u>1/3`, and its integral
is `9/16`. This verifies both the factor two from the subGaussian tail
and the factor `R` from the shell deficit.

If `tau>9K/8`, first choose `kappa>K`, then `R>1` and `zeta>0` close
enough to one and zero, respectively, so the right side of (14) is
strictly less than one. Choose a fixed identity weight `theta>0` smaller
than the remaining resource slack. For every sufficiently small fixed
`epsilon`, (12)--(14) and the terminal estimates yield total feature
weight at most `1-theta`.

## 7. Paying all responses and the order of limits

For `x in E_j\E_(j+1)`, its old absolute-energy deficit exceeds
`eta_(j+1)=eta_j/R`, which pays the first term in (11). For the terminal
code no old deficit is needed. Outside `E_0`, the identity feature and
a union bound over **all** old/new words give normalized bridge bound

```math
\sqrt{2\epsilon(1+\epsilon)\log2/\theta}+o_n(1),
```

which is below the fixed `eta_0` once `epsilon` is sufficiently small.
Choose a full-sign new principal block of order `q` with cap at most
`q^(3/2)`; elementary independent-edge signing supplies one for all large
`q`. For every parent word,

```math
|H_{B}(x,y)|\le |H_A(x)|+|x^TCy|+Q(\text{new block}).
```

The absolute bridge is fully paid, so independent reversal of the new
spins causes no missing polarity case. Every old edge is untouched,
every new edge is filled by a sign, and (1) follows.

The quantifier order is important. Choose `tau,kappa,R,zeta,theta,eta_0`
first. Then choose small **fixed** `epsilon` and its finite hierarchy.
Only then let `n` tend to infinity. The strict factorization majorants
and all finite-level entropy estimates apply simultaneously. Last let
`epsilon` decrease. No claim at `n`-dependent vanishing energy windows is
used.

For a liminf-realizing sequence, if `K<4c_*/3`, choose
`9K/8<tau<3c_*/2`. Formula (1) gives a parent subsequence with normalized
cap at most

```math
\frac{c_*+\tau\epsilon+\epsilon^{3/2}}
{(1+\epsilon)^{3/2}}<c_*
```

for sufficiently small fixed `epsilon`, contradicting the global liminf.
This proves (2). It is a stronger quantitative obstruction to simple
augmentation, not a proof that favorable augmentation holds for minimizers.

## 8. Direct edge rounding: useful certificate and unavoidable scope limits

Let `m=binom(n,2)`, `q_x=(x_i x_j)_(i<j)`, and fractional edge vector
`b in [-1,1]^m`. For any positive definite matrix `Q` on edge space with
`Q_ee<=1`, apply (5) to the columns of `Q^(1/2)`. This gives a sign
rounding `s` with `E s=b` and

```math
\mathbb E e^{\langle t,s-b\rangle}
\le e^{t^TQ^{-1}t/2}.                            \tag{15}
```

For a proposed cap `L>=Q(b)`, if

```math
\sum_{x:x_1=1}\sum_{\sigma=\pm1}
\exp\left\{-\frac{(L-\sigma\langle b,q_x\rangle)^2}
{2q_x^TQ^{-1}q_x}\right\}<1,                     \tag{16}
```

then one full-sign output has cap at most `L`. All Boolean queries and
both polarities remain. This is GS plus an exact slack-aware union bound,
not separately claimed as a new convergence theorem.

Three limits prevent overinterpreting it.

- Boolean orthogonality gives `E_x q_x q_x^T=I_m`. Thus the average
  proxy in (16) is `tr(Q^-1)>=m` under `Q_ee<=1`; deterministic
  preconditioning cannot lower *every* query's proxy below independent
  rounding scale. It can move budget toward a smaller exposed family.
- Arbitrary bias does not permit an all-lambda subGaussian proxy bounded
  by a universal constant times true covariance. In one dimension let
  `P(s=-1)=p`, `b=1-2p`. At
  `lambda=-2log(1/p)/(1+b)`, the negative-sign event gives
  `log E exp(lambda(s-b))>=log(1/p)`. Any subGaussian proxy must therefore
  be at least `(1+b)^2/[2log(1/p)]`, while the variance is `4p(1-p)`.
  Their ratio diverges as `p` decreases. A Bernstein/local-Laplace or
  rare-event term is necessary for genuinely covariance-scaled bias.
- If `V=Q^(1/2)`, the sharper random phase proxy in (4) yields useful
  tails on events where its directional quadratic form is small. But
  turning its *average* into a uniform MGF proxy is invalid. The missing
  input is a joint, all-relevant-query control of that random geometry.

## 9. Actual-entropy refinement and a nonzero inner complexity

**Independent audit PASS (BH/mechanism agent).** The factorization-only
profile above deliberately paid for all codewords through Sauer's bound.
The cumulative representation also gives a stronger type of statement
when the actual code entropy is smaller. Define

```math
s(\eta)=\limsup_n n^{-1}\log|E_n(\eta)|,\qquad
e_0=\lim_{\eta\downarrow0}e(\eta),\qquad
K_{es}=\limsup_{\eta\downarrow0}\frac{e(\eta)s(\eta)}\eta.
```

Both profiles are nondecreasing in `eta`, so `e_0` exists. Suppose
`K_es<infinity`. Then the actual extension conclusion (1) holds for every
positive `tau` such that

```math
\boxed{\frac{2K_{es}}\tau+
       \frac{2(\log2)e_0}{\tau^2}<1.}             \tag{17}
```

The second term is essential. If the innermost code has positive
factorization complexity density but only subexponentially many words,
the old-word entropy vanishes, but maximizing over the new spins still
costs a positive amount. Dropping this term is not justified.

**Proof.** Choose `kappa>K_es`, a small `eta_0`, a ratio `R>1`, and
positive tail slack `zeta`. For a fixed small `epsilon`, use the same
finite hierarchy down to `eta_J` of order `epsilon^2`. Choose monotone
strict upper bounds `g_j>e(eta_j)` and `s_j>s(eta_j)` with
`g_j s_j<=kappa eta_j`; the trivial bounds `g_j<=1` and `s_j<=log2`
may be used exactly at their endpoints. To see that this choice has no
hidden regularity requirement, first bound the original products by
`kappa_0 eta` with `K_es<kappa_0<kappa`, then add to both profiles a
common sufficiently small positive number, at most `epsilon^3`, and
truncate at the trivial endpoints. There are only finitely many levels.
Monotonicity is preserved and all strict positive gaps needed for the
limsup bounds remain. This perturbation tends to zero with `epsilon`.

Use (7) with

```math
b_j=\frac{2(1+\zeta)\epsilon(s_j+\epsilon\log2)}
{(\eta_j/R+\tau\epsilon)^2}\quad(j<J),
\qquad
b_J=\frac{2(1+\zeta)\epsilon(s_J+\epsilon\log2)}
{(\tau\epsilon)^2}.
```

Separate these requested precisions as `b_j=d_j+y_j`, for the old-word
and new-word entropies, respectively. The sequence `y_j` increases.
For `1<=j<J`, monotonicity `s_j<=s_(j-1)` gives

```math
(d_j-d_{j-1})_+
\le2(1+\zeta)\epsilon s_j
\left[(\eta_j/R+a)^{-2}-(\eta_{j-1}/R+a)^{-2}\right],
\qquad a=\tau\epsilon.
```

Multiply by `g_j`, use `g_j s_j<=kappa eta_j`, and bound the resulting
Riemann sum by the integral with `eta_j<=eta` on its interval. Apart
from the `O(epsilon)` initial and terminal old-entropy terms, the total
old-entropy resource is at most

```math
\frac{4(1+\zeta)\epsilon\kappa}{R}
\int_0^{\eta_0}\frac{\eta}{(\eta/R+a)^3}\,d\eta
\le\frac{2(1+\zeta)\kappa R}{\tau}.             \tag{18}
```

The terminal old-entropy term is bounded by
`g_J d_J<=2(1+zeta)kappa eta_J/(tau^2 epsilon)=O(epsilon)`.
For the new entropy, charge the **increment** `g_J(y_J-y_(J-1))` at
the terminal level, not the whole `g_J y_J` a second time. Fix `delta>0`.
For levels with `eta_j>=delta`, their new-entropy resource is
`O(epsilon^2/delta^2)`, since `g_j<=1`. Below `delta`, every `g_j` is
at most `e(delta)+o_epsilon(1)`, and the increments telescope to at most
`y_J=2(1+zeta)log2/tau^2`. Hence

```math
\limsup_{\epsilon\downarrow0}\sum_jg_j(y_j-y_{j-1})
\le\frac{2(1+\zeta)(\log2)e_0}{\tau^2},          \tag{19}
```

where `y_(-1)=0`; send `delta` to zero last. Combining positive
variations uses
`(d_j-d_(j-1)+y_j-y_(j-1))_+ <= (d_j-d_(j-1))_++y_j-y_(j-1)`.
Equations (18)--(19) bound the running-maximum resource from (7).
Under (17), choose `kappa,R,zeta` close enough to `K_es,1,0` to leave
a positive identity budget. All tail, all-query, and limit-order steps
of Sections 6--7 then apply verbatim, proving (1).

For any liminf-realizing family, (17) cannot hold at a `tau<3c_*/2`.
By continuity at that endpoint this proves the necessary condition

```math
\boxed{\frac{4K_{es}}{3c_*}
       +\frac{8(\log2)e_0}{9c_*^2}\ge1.}          \tag{20}
```

In particular, if `e_0=0`, then `K_es>=3c_*/4`. If `K_es=0`, then
`e_0>=9c_*^2/(8log2)`. This can be stronger than the factorization-only
test: for example, abstract profiles `e(eta)=eta^(1/4)` and
`s(eta)=eta^(3/4+delta)` with `delta>0` have `K_es=0` and `e_0=0`,
while `e(eta)H(e(eta))/eta` diverges. Such profiles are now ruled out
for liminf signings without asserting that any actual signing realizes
the profiles.

## 10. Reproduction and current status

The script
`computations/paper_discrepancy_2026_09_17_precision.py` independently
checks the Bellman gradient determinant, harmonic representation identity,
the exact `9/16` integral, finite-grid resource bounds, convergence of
the smooth entropy integral toward `9/8`, and the one-coordinate bias
obstruction. The additional script
`computations/paper_discrepancy_2026_09_17_exact_walk.py` exhaustively
branches 100 small biased GS instances, checking their means, sharp MGF,
and random phase compensator. It found no counterexample to these
imported statements. A search for failure after replacing the random
compensator by its expectation was inconclusive, not a proof of that
replacement. Both scripts use the repository `.venv` and no ignored
mathematical inputs. Numerical checks are corroboration; the proof is
Sections 3--9.

Sections 5--7 received an independent mathematical PASS for the
cumulative-precision allocation, integral constant, and limit-order
argument. Section 9 also received an independent PASS, including the
nonzero-`e_0` terminal bookkeeping. No steering/ledger changes or commits
are made by this agent.

## 11. Counterexample-first search for a usable quadratic profile upper bound

The extension tests are not automatically applicable to all bounded-cap
matrices. There is an exact obstruction to hoping that bounded cap alone
forces `e_0=0`, or even `e_0<1`.

Let `n=4^a`, `H=(J_4-2I_4)^(tensor a)`, and
`A=H-(-1)^a I`. The existing
[`Boolean-eigenbasis audit`](resumed_bound_audit_boolean_eigenbasis_hamming_rigidity_2026_09_06.md)
establishes `Q(A)=(n sqrt(n)+n)/2`. Tensor products of the four Walsh
vectors give an orthogonal Boolean eigenbasis of `H`; every vector in
this basis has absolute-energy gap either zero or `n`. Thus every
fixed positive `eta` eventually puts the whole basis inside `E_n(eta)`.
If a sign-row code contains an isotropic law `E xx^T=I`, every
factorization `G=LV` has

```math
n=\|\operatorname{diag}(\sqrt\mu)G\|_*
\le\|\operatorname{diag}(\sqrt\mu)L\|_F\|V\|_F
\le\gamma_{\rm fac}(G)\sqrt n.
```

Together with the identity factorization this gives
`gamma_fac(G)^2=n`. Therefore this bounded-cap family has exactly
`e(eta)=e_0=1` for every positive `eta`.

A *genuinely* sub-`1/2` asymptotic cap might behave differently. Its naive
finite analogue is false: the archived order-6, order-12, and order-14
signings have exact caps `5,18,21`, respectively, and rational isotropic
laws supported on their exact absolute grounds. These laws were replayed
independently in the present campaign. The normalized caps are about
`0.340207`, `0.433013`, and `0.400892`. The first two orders have uniform
isotropic ground laws; the stored order-14 rational certificate has 92
atoms. Each therefore has maximally large factorization norm on its
ground code despite its finite subhalf cap.

The precise archive collision is
[`Failed isotropic-ground shortcut`](principle_director_ground_isotropy_failed_probe_2026_09_07.md)
and
[`finite exact-minimizer LP certificates`](transfer_adversary_exact_minimizer_isotropy_finite_2026_09_06.md).
The present replay does not turn finite data into an asymptotic
counterexample, and imports rather than reruns the global optimality
certificates at those orders.

All Paley conference matrices have isotropic absolute-ground laws by
signed pair-transitive symmetry and a signed permutation carrying `A`
to `-A`, as proved in
[`the Paley ground-law audit`](transfer_fresh_cavity_and_isotropic_ground_law_2026_09_06.md).
This assertion does not require Boolean eigenvectors. Consequently an
asymptotic subhalf cap bound on any infinite Paley subfamily would falsify
the proposed asymptotic isotropy-rigidity lemma. The square-field family
instead approaches `1/2`, and no strict-asymptotic-subhalf Paley theorem
is imported here.

The reproducible probe
`computations/paper_discrepancy_2026_09_17_isotropy.py` independently
replays the rational finite laws and uses the existing exact integer
Gray-code evaluator for prime-field Paley examples through order 30.
Neither a universal useful upper profile nor an asymptotic isotropic
subhalf counterexample has been proved by this check.

## 12. Quantitative covariance meaning of maximal factorization complexity

**Derived finite theorem; independently audited by the Bernoulli track.** For a nonempty
sign code `E subset {+-1}^n`, define

```math
g(E)=\gamma_{\rm fac}(E)^2/n,\qquad
D(E)=\frac1n\min_{\mu\text{ on }E}
\|\mathbb E_\mu xx^T-I\|_* .
```

The star denotes nuclear norm, not the Boolean cap. Then

```math
\boxed{\left(1-\frac{D(E)}2\right)^2
\le g(E)\le
\frac{1+\sqrt{1-D(E)^2/4}}2.}                    \tag{21}
```

In particular `g(E)=1` if and only if `E` admits an isotropic law.
Quantitatively, `g(E)>=1-epsilon` with `0<=epsilon<=1/2` implies
`D(E)<=4sqrt(epsilon(1-epsilon))`; conversely `D(E)=o(1)` forces
`g(E)=1-o(1)`. Thus near-maximal factorization complexity is equivalent
to *some* probability law on the full code having asymptotically isotropic
covariance in normalized nuclear norm. This is an existence statement
about a law, not a claim about the uniform or Gibbs law.

For the lower inequality take a law with covariance `Sigma`. For every
factorization `G=LV`, nuclear/Frobenius duality gives

```math
\operatorname{tr}\sqrt\Sigma
=\|\operatorname{diag}(\sqrt\mu)G\|_*
\le\gamma_{\rm fac}(E)\sqrt n.
```

Since `tr Sigma=n` and `sqrt(lambda)>=min(lambda,1)`,
`tr sqrt Sigma>=n-||Sigma-I||_*/2`. Minimize the latter defect to obtain
the lower bound in (21).

For the upper inequality use convex separation in nuclear norm. There
is a symmetric matrix `C` with `||C||op<=1` such that

```math
\min_{x\in E}\operatorname{tr}C(xx^T-I)=nD(E).
```

This is the ordinary finite-dimensional minimax dual of the distance
from `I` to `conv{xx^T}`. Hollow `C` and divide by two:
`B=(C-diag C)/2` has zero diagonal, `||B||op<=1`, and
`min_E x^T Bx>=nD(E)/2`. More generally, if any such `B` has
`min_E x^T Bx>=delta n`, then the Gram metric `M=I+tB`, `0<t<1`,
has diagonal one and supplies a valid factorization. Spectral calculus
and the scalar chord inequality

```math
\frac1{1+tz}\le\frac{1-tz}{1-t^2}\quad(-1\le z\le1)
```

give

```math
\frac{\gamma_{\rm fac}(E)^2}{n}
\le\max_{x\in E}\frac{x^T(I+tB)^{-1}x}{n}
\le\frac{1-t\delta}{1-t^2}.
```

Optimizing at `t=delta/(1+sqrt(1-delta^2))`, with endpoint limits when
needed, yields `(1+sqrt(1-delta^2))/2`. Set `delta=D(E)/2` to obtain
the upper bound in (21). Notice that `0<=D(E)<=2`, since the compared
positive matrices both have trace `n`.

This gives a concrete covariance-sensitive target for the paper-combination
campaign: a uniformly positive normalized nuclear separation from isotropy
would force an inner factorization gap. Exact finite failure of isotropy,
or an unnormalized separation supported on a mesoscopic planted block,
does not provide that scale.

## 13. Exact failure of averaging the random GS compensator

**Independent exact-arithmetic audit: PASS by the Bernoulli track.**

The warning in Section 3 is now an exact counterexample, not merely an
absence of a proof. The random-phase certificate

```math
\mathbb E\exp\{\langle\theta,V(s-b)\rangle
-\tfrac12\theta^T\Gamma\theta\}\le1
```

does **not** imply, and the GS walk does not in general satisfy,

```math
\mathbb E e^{\langle\theta,V(s-b)\rangle}
\le e^{\theta^T(\mathbb E\Gamma)\theta/2}.          \tag{22, false}
```

Take the unit-column rational matrix, biased start, and test vector

```math
V=\begin{pmatrix}-1&1/3&0\\0&2/3&3/5\\0&-2/3&-4/5\end{pmatrix},
\quad b=(99/100,1/2,-4/5),\quad\theta=(8,-3,3).
```

Use the largest live index as pivot, retaining it until it freezes.
Exact rational enumeration gives eight paths. The phase compensator is

```math
\Gamma_0=\frac1{2025}
\begin{pmatrix}457&392&-392\\392&1045&-784\\-392&-784&1045\end{pmatrix}
\quad\text{with probability }150/169,
```

and is `I` on the remaining paths. Consequently

```math
\mathbb E\Gamma=\frac1{4563}
\begin{pmatrix}1427&784&-784\\784&2603&-1568\\-784&-1568&2603\end{pmatrix},
\qquad \tfrac12\theta^T\mathbb E\Gamma\theta=\frac{45571}{4563}.
```

The output `s=(-1,-1,-1)` occurs on one path with probability `13/14000`
and has `theta^T V(s-b)=469/25`. This single path implies

```math
\mathbb E e^{\langle\theta,V(s-b)\rangle
-\theta^T(\mathbb E\Gamma)\theta/2}
\ge\frac{13}{14000}\exp\left(\frac{1000772}{114075}\right)>1.
```

The final inequality needs no numerical transcendental assertion:
`1000772/114075>8` and
`sum_{j=0}^{11}8^j/j!=412782941/155925>14000/13`.
For scale only, the full logarithmic violation is approximately `1.79208`,
whereas the valid random-compensated log moment is `-4.96115`.

The reproduction
`computations/paper_discrepancy_2026_09_17_compensator.py`
computes all endpoint moves, probabilities, orthogonal projectors, and
the counterexample inequalities with exact rational arithmetic. This
supersedes the inconclusive moderate-bias numerical search in Section 10.
It does not contradict `Cov(V(s-b)) <= E Gamma`, which is a second-moment
conclusion. It specifically blocks importing the covariance average as
an all-lambda subGaussian metric in the fixed-law Bernoulli composition.

## 14. Basis-cover stress test: entangled ground states defeat every linear-size cover

The new sparse-in-basis full-sign column theorem in
[`the Bernoulli artifact, Section 15`](paper_bernoulli_2026_09_17.md)
has been independently audited here. Its sufficient extension slope is
`sqrt(2/pi)/sqrt(k)+2k K_r`, where
`K_r=limsup_(eta down to0) r(eta)h(r(eta))/eta` and the full absolute
near-level code must be covered by the antipodal rows of a tensor
Hadamard basis. In particular `r(eta)=O(eta)` gives `K_r=0`, because
`r h(r)~r^2 log(1/r)`. The factorization norm of the full center family
is maximal, but this does not imply the required cover statement.

Here is an exact obstruction to applying that theorem to the familiar
near-half tensor-Hadamard example. It rules out **every** linear-size
center family, not only one natural basis.

Let `H_4=J_4-2I_4`, `H_p=H_4^{tensor a}`, `p=4^a`, and
`H_n=H_4 tensor H_p`, `n=4p`. The matrix `H_p` has a full orthogonal
Boolean eigenbasis. Select `d>=p/2` basis vectors `u_1,...,u_d` sharing
one eigenvalue `lambda in {+-sqrt(p)}`. For every ordered pair form

```math
x_{ij}=(u_i,-u_i,u_j,-u_j)\in\{\pm1\}^{4p}.
```

Their four blocks sum to zero, so
`H_n x_ij=-2lambda x_ij`. After hollowing,
`A_n=H_n-(-1)^(a+1)I`, these are exact absolute ground states or have
absolute deficit exactly `n`. Thus all `d^2` words belong to
`E_n(eta)` eventually for each fixed `eta>0`.

Consider any sign center `f=(f_1,f_2,f_3,f_4)` and put

```math
A_i=\langle u_i,f_1-f_2\rangle,\qquad
B_j=\langle u_j,f_3-f_4\rangle.
```

Then `|A_i|,|B_j|<=2p`, and orthogonality gives
`sum_i A_i^2<=4p^2`, `sum_j B_j^2<=4p^2`. If
`d_H(x_ij,f)<=(1/4-delta)n`, `0<delta<1/4`, then

```math
A_i+B_j\ge(2+8\delta)p.
```

Since either summand is at most `2p`, both must be at least `8delta p`.
There are at most `1/(16delta^2)` choices for either index. Therefore
one center covers at most `1/(256delta^4)` of the selected words.
Any `L_n=O(n)` centers cover only `O(n/delta^4)` of the `d^2=Omega(n^2)`
words. Consequently every such center family satisfies, for every
fixed positive energy window,

```math
\liminf_n\frac1n\max_{x\in E_n(\eta)}
             \min_{f\in F_n}d_H(x,f)\ge\frac14.       \tag{23}
```

This includes any antipodal Hadamard basis, even one chosen adaptively
from the signing. For the canonical tensor eigenbasis the explicit
words with `i!=j` already have exact distance `n/4` from their nearest
center. The reproduction
`computations/paper_discrepancy_2026_09_17_basis_cover.py`
checks the eigenvector identity, near-level deficit, and canonical
distance at finite orders. The asymptotic arbitrary-center conclusion
uses the displayed exact Bessel/counting argument, not those checks.

This is a scope test, not a counterexample to the new cover theorem.
It shows why isotropic laws and maximal center factorization cannot
substitute for control of the **full** near-extreme landscape.

## 15. Stronger stress test: no fixed repeated-block mode covers the full ground landscape

**New finite-moment argument; independently audited PASS by the Bernoulli
track, including its improved exceptional-edge count.** The extension
of the Bernoulli track's basis-column theorem permits subexponential
inner center codes, not merely a linear number of orthogonal centers.
Section 14 by itself does not exclude those larger codes. The following
argument excludes, on the same near-half family, even the repeated-block
mode family with its inner code enlarged to the **entire** Boolean cube.

Let `n=4^a` and `H_n=(J_4-2I_4)^{tensor a}`. Construct a sign-valued random
Boolean eigenvector `X_a` recursively. Start with an independent unbiased
sign `X_0`. At each step choose uniformly one of the three perfect matchings
of four block labels. In one matched pair place `(U,-U)`, and in the other
place `(V,-V)`, where `U,V` are independent copies of `X_(a-1)`. Orient
each pair by its increasing labels; symmetry of the child law makes this
choice immaterial. The four blocks sum to zero, so inductively

```math
H_n X_a=(-2)^a X_a.                              \tag{24}
```

After hollowing, every output has deficit zero or `n` from the absolute
cap. Hence its entire support lies in `E_n(eta)` eventually, for every
fixed `eta>0`. This law is used only as an averaging witness; it need
not be isotropic, and no such claim is made.

### 15.1 Uniform second and fourth moments

Index coordinates by words of length `a` over `{0,1,2,3}`. If distinct
indices `i,j` differ in `d(i,j)` positions, recursion gives

```math
\mathbb E X_iX_j=(-1/3)^{d(i,j)}.                 \tag{25}
```

In particular every off-diagonal pair correlation has absolute value
at most `1/3`. For four **distinct** indices,

```math
\left|\mathbb E X_iX_jX_lX_m\right|\le11/27,       \tag{26}
```

except when the indices form a coordinate-axis four-line: all digits
agree except one digit, which takes all four possible values. There are
exactly `a n/4` such unordered four-sets, and their fourth moment is one.

Here is the complete induction for (26), splitting indices by their
first digit. If all lie in one block, recurse. For block multiplicities
`3+1`, only the matching joining those two blocks contributes, giving
absolute value at most `1/3`. For `2+2`, that matching has probability
`1/3` and contributes at most one; the other two matchings contribute
products of off-diagonal child correlations, each at most `1/9`.
The total is at most `1/3+2/27=11/27`. For `2+1+1`, only the matching
joining the two singletons contributes. The remaining doubled block
has distinct child indices, so its pair moment is at most `1/3`;
the total is at most `1/9`.

For `1+1+1+1`, the moment is the average over three pairings of products
of two child pair correlations. If all child indices agree, this is
the exceptional four-line and the moment is one. Otherwise at most
one pairing can have both child pairs equal; it contributes at most
one, and the others at most `1/9`. If no pairing has both pairs equal,
a product may instead have just one equal pair and absolute value
`1/3`, so the total is at most `1/3`, still below `11/27`. These cases
exhaust the possibilities. Odd child moments vanish by antipodality.

### 15.2 Arbitrary signed matchings cannot have a common mode

Let `M` be **any** matching of `p` disjoint coordinate pairs, carrying
arbitrary fixed signs `sigma_e`. Set

```math
S_M(X)=\sum_{e=\{i,j\}\in M}\sigma_e X_iX_j,
\qquad c_4=11/27.
```

The diagonal terms in `E S_M^2` contribute `p`. Each pair of distinct
matching edges uses four distinct coordinates, so (26) applies. An
edge can belong to at most one coordinate-axis four-line. Consequently
at most `p/2` unordered edge pairs can form exceptional four-lines.
This sharpening of the initial global four-line count was supplied
by the independent Bernoulli-track auditor. Thus

```math
\boxed{\mathbb E S_M(X)^2
\le c_4 p^2+2(1-c_4)p=\frac{11p^2+32p}{27}.}   \tag{27}
```

There is therefore an output `x` with `|S_M(x)|` at most the square root
of the right side. To make all `p` signed pair products equal to either
`+1` or `-1`, at least `(p-|S_M(x)|)/2` coordinates must be changed.
The disjointness of the pairs is essential and is part of the hypothesis.

For a perfect matching, `p=n/2`, this proves an asymptotic distance of
at least

```math
\frac14(1-\sqrt{11/27})n-o(n)
```

from the union of the two common-product mode cubes, uniformly over the
choice of matching and its signs. Thus neither a subexponential inner
code nor even the full inner cube fixes the cover failure.

More generally fix `k>=2`, let `n=kp`, and consider any family of modes
`h_l tensor g`, where each `h_l` is a fixed sign vector in dimension
`k` and `g` ranges over any subset of `{+-1}^p`. Arbitrary fixed coordinate
permutations and switches are permitted. Pair corresponding coordinates
of the first two blocks. Every mode has a single common signed pair
product `h_(l,1)h_(l,2)`. Applying (27) to this partial matching gives

```math
\liminf_n\frac1n\max_{x\in E_n(\eta)}
  \min_{f\text{ in the mode family}}d_H(x,f)
\ge\frac{1-\sqrt{11/27}}{2k}>0                 \tag{28}
```

for every fixed `eta>0` and every fixed admissible `k>=2`, along these
orders. An `o(n)` leftover set of coordinates has no effect on the
conclusion. This falsifies the cover hypothesis on the near-half family
for **any fixed repeated-block mode geometry**, while leaving the new
conditional extension theorem untouched.

The rational moment recursion and the finite eigenvector checks are
reproduced by
`computations/paper_discrepancy_2026_09_17_recursive_ground.py`.

## 16. Exact finite-minimizer common-mode diagnostics

The stored even-order exact-minimizer examples were independently
enumerated at the level of their full projective absolute-ground codes.
Global optimality labels are imported from the prior solver certificates;
the cap and all code calculations below were replayed using integer
arithmetic. Let

```math
c_4(E)=\max_{i,j,l,m\text{ distinct}}
       \left|\mathbb E_{x\text{ uniform on }E}x_ix_jx_lx_m\right|.
```

Antipodal symmetrization leaves these moments unchanged. The results are:

| order | projective ground words | exact `c_4(E)` | common signed perfect-matching mode possible? |
|---|---:|---:|---|
| 4 | 2 | 1 | yes |
| 6 | 12 | 1/3 | no |
| 8, each of two classes | 8 | 1/2 | no |
| 10 | 40 | 3/10 | no |
| 12 | 20 | 3/5 | no |
| 14 | 156 | 3/13 | no |

For every listed order at least six, the edge-character columns
`(x_i x_j)_(x in E)` are pairwise non-proportional even up to a minus
sign. If a common signed matching mode existed, all its edge-character
columns would be proportional. Thus the obstruction excludes every
signed perfect matching at once, even if the inner code is the full
cube. It also excludes exact fixed-`k` mode covers whenever the two
selected blocks contain at least two matched pairs.

The fourth moments give a quantitative finite certificate. For any
signed matching with `p` edges,

```math
\mathbb E S_M^2\le c_4(E)p^2+(1-c_4(E))p.
```

Consequently some ground word requires at least
`(p-sqrt(c_4(E)p^2+(1-c_4(E))p))/2` coordinate changes to reach either
common-product mode. One may round this bound upward to an integer,
because the number of required changes is integral. In particular the
order-14 example requires at least two changes for every signed perfect
matching. These finite observations do not imply an asymptotic lower
cover radius for minimizers. Their role is to reject an unsupported
inference from isotropy or small cap to a common repeated-block pattern.

Reproduction: `computations/paper_discrepancy_2026_09_17_mode_probe.py`.

## 17. A random-involution ground law defeats every fixed finite union of mode frames

**New application of an explicit Maiorana--McFarland-type sign law;
independently audited PASS by the Bernoulli track.** This strengthens Section 15 from one or
two mode geometries to every fixed finite union, with asymptotically
vanishing matching correlations. The classical construction template
`u dot pi(v)+g(v)` is discussed in Carlet--Danielsen--Parker--Sole,
[*Self-dual bent functions*, Section 4.1.2](https://www.codetables.de/larsed/sdbent.pdf).
The Fourier kernel used below is the alternating kernel, **not** the
ordinary dot-product kernel in that paper's definition of self-duality.
The exact eigenvector identity is proved here, so no stronger imported
self-duality classification is assumed.

### 17.1 Full-sign eigenvectors from random involutions

Let `p=2^a>=8`, `n=p^2`, and index physical coordinates by
`(u,v) in F_2^a times F_2^a`. Choose uniformly a fixed-point-free
involution `pi` of the `p` labels. Independently choose a fair bit on
each of its two-cycles, giving a function `g` with `g(pi(v))=g(v)`.
Define

```math
X_{u,v}=(-1)^{u\cdot\pi(v)+g(v)}.
```

For the alternating Fourier matrix
`F_(s,t),(u,v)=(-1)^(s dot v+t dot u)`, direct summation over `u` gives

```math
(FX)_{s,t}
=p(-1)^{s\cdot\pi(t)+g(\pi(t))}=pX_{s,t}.       \tag{29}
```

Set `Q(u,v)=u dot v+sum_i(u_i+v_i)` over `F_2` and
`D_(u,v)=(-1)^Q(u,v)`. After regrouping the binary digits,

```math
(J_4-2I_4)^{\otimes a}=(-1)^a DFD.
```

Thus `Y=DX` is a full-sign eigenvector with eigenvalue `(-1)^a p`.
For the hollow matrix `A_n=H_n-(-1)^a I`, every such vector has absolute
energy deficit exactly `n` from its cap. Therefore the law's full support
lies in `E_n(eta)` eventually for each fixed `eta>0`.

### 17.2 Four-coordinate bound, with exactly localized exceptions

For distinct physical coordinates, pair moments have absolute value
at most `1/(p-1)`: if their `v` labels agree, average a nontrivial
character over the `p-1` possible partners; if they differ, the fair
cycle bits cancel only when those two labels are paired by `pi`.

For four distinct physical coordinates define

```math
c_p=\frac{2p-1}{(p-1)(p-3)}<1.
```

Their fourth moment has absolute value at most `c_p`, except when all
four `v` labels agree and their four `u` labels have xor zero. In that
exceptional case it equals one. Fixed multiplication by `D` changes
only the signs of moments, not these absolute bounds.

Here is the complete classification by `v` multiplicities:

- `4`: the moment is the average of the character with frequency the
  xor of the four `u` labels. It equals one when that xor is zero;
  otherwise its absolute value is `1/(p-1)`.
- `3+1`: the two distinct labels must be paired, so the bound is
  `1/(p-1)`.
- `2+2`: let the nonzero xor frequencies in the two columns be `d,e`.
  The labels are paired with probability `1/(p-1)`. Otherwise their
  distinct partners are uniform ordered elements of the set `U` of
  size `p-2` excluding those labels. The character numerator is
  `sum_U chi_d sum_U chi_e - sum_U chi_(d+e)`, of magnitude at most
  `4+(p-2)=p+2`. The resulting bound is
  `1/(p-1)+(p+2)/[(p-1)(p-3)]=c_p`.
- `2+1+1`: the two singleton labels must pair, with probability
  `1/(p-1)`. The partner of the doubled label is then uniform outside
  three specified labels. Its nontrivial character average is at most
  `3/(p-3)`, giving `3/[(p-1)(p-3)]`.
- `1+1+1+1`: all four labels must pair internally, an event of probability
  `3/[(p-1)(p-3)]`.

### 17.3 A uniform matching theorem and simultaneous frame exclusion

For any signed matching `M` of `m` physical-coordinate pairs, put
`S_M=sum_(ij in M) sigma_ij Y_iY_j`. Exceptional pairs of its edges must
lie in the same `v` column and have the same nonzero xor difference in
their `u` endpoints. Group such within-column edges by `(v,d)`.
Each group has at most `p/2` edges; therefore the total number of
unordered exceptional edge pairs is at most `m(p-2)/4`. Consequently

```math
\boxed{\mathbb E S_M^2
\le c_p m^2+(1-c_p)\frac{mp}{2}.}               \tag{30}
```

In particular for a matching between two blocks of a fixed-`k` mode
geometry, `m=n/k+O(1)`,

```math
\frac{\mathbb E S_M^2}{m^2}
\le c_p+(1-c_p)\frac{k}{2p}+o(1/p)=O_k(1/p).
```

For any fixed number of independently specified mode geometries,
Markov's inequality and a finite union bound give a single near-ground
word with `|S_(M_r)|=o(m_r)` for all their matchings simultaneously.
The matchings, switches, and permutations can depend on `n` and on the
old signing; the estimates are uniform over those choices. The number
of differing coordinates from any mode in the `r`-th geometry is at
least `(m_r-|S_(M_r)|)/2`. Thus, writing `k_max=max_r k_r`,

```math
\liminf_n\frac1n\max_{x\in E_n(\eta)}
 \min_{f\in\bigcup_r F_{r,n}}d_H(x,f)
\ge\frac1{2k_{\max}}                            \tag{31}
```

for every fixed `eta>0`. Each inner dictionary may be the entire cube;
its cardinality is irrelevant to this obstruction. Leftovers `o(n)`
are harmless. This excludes every fixed finite union of the proposed
repeated-block mode dictionaries on this actual near-half family. It
does not exclude more general column laws or establish a property of
unknown asymptotic minimizers.

The reproduction `computations/paper_discrepancy_2026_09_17_involution_ground.py`
enumerates all involution/cycle-bit choices at `p=8`, checks exact
eigenvectors and pair moments, tests the classified fourth moments,
and replays (30) for signed matchings with integer arithmetic.

## 18. Growing low-effective-mode dictionaries also fail on the full Hadamard landscape

**New consequence of Section 17; independently audited PASS by the
Bernoulli track.** The
localization track has strengthened the positive extension theorem to
balanced deterministic mode counts, permitting `k=k_n` to grow and
covering a larger dictionary than pure mode words. The effective mode
count of a Boolean vector is

```math
R_k(z)=\frac1n\left(\sum_{a=1}^k\|P_a z\|_2\right)^2,
```

where `P_a` are the mutually orthogonal Hadamard-mode projections of
rank `n/k`. The whole code `R_k(z)<=u_n=o(k)` is now allowed. A
single matching obstruction would not exclude it. Nevertheless the
random-involution law gives a stronger energy-spreading certificate
for **every switched and permuted Sylvester mode frame** when
`k=o(sqrt(n))`.

Continue with `n=p^2=4^a`, and let `k` be a power of two dividing `n`.
The grouping, coordinate permutation, and switches are arbitrary and
may depend on `n`. For `Y` from Section 17 set

```math
w_b=\frac{\|P_bY\|_2^2}{n},\qquad
D_k(Y)=k\sum_b w_b^2-1\ge0.
```

Then

```math
\boxed{\mathbb E D_k(Y)
\le(k-1)\left(c_p+\frac{1-c_p}{p}\right)=O(k/p).} \tag{32}
```

To see this index block positions by the additive group of order `k`.
For every nonzero xor difference `d`, pair block positions `b` and
`b+d` within each of the `n/k` coordinate groups. These edges form a
perfect matching of the `n` physical coordinates. Incorporate the
fixed coordinate switches in its signs and denote its statistic by
`U_d`. For the Sylvester row indexed by `b`, write
`T_b=k||P_bY||^2`. Expansion and finite-group Parseval give the exact
identities

```math
T_b=n+2\sum_{d\ne0}\chi_b(d)U_d,
\qquad \sum_bT_b^2=kn^2+4k\sum_{d\ne0}U_d^2.
```

Apply (30), with matching size `n/2`, to each of the `k-1` statistics.
This proves (32), uniformly over the frame.

The mode-amplitude vector is close to uniform whenever `D_k` is small:

```math
\sum_b\left(\sqrt{w_b}-\frac1{\sqrt k}\right)^2
\le k\sum_b(w_b-1/k)^2=D_k(Y).
```

For **every** Boolean vector `z` with `R_k(z)<=u`, Cauchy--Schwarz
within the orthogonal modes and then in the mode index yield

```math
\frac{|\langle Y,z\rangle|}{n}
\le\sum_b\sqrt{w_b}\frac{\|P_bz\|_2}{\sqrt n}
\le\sqrt{u/k}+\sqrt{D_k(Y)}.                    \tag{33}
```

Now let `k_n->infinity`, `k_n=o(p)=o(sqrt(n))`, and `u_n=o(k_n)`.
Equation (32) allows a choice of near-ground word `Y` with `D_k(Y)=o(1)`.
Equation (33), uniform over the entire low-effective-mode code, gives

```math
\liminf_n\frac1n\max_{x\in E_n(\eta)}
 \min_{z:R_{k_n}(z)\le u_n}d_H(x,z)\ge\frac12   \tag{34}
```

for every fixed `eta>0`. A finite union of such frames has the same
obstruction: use Markov's inequality and a finite union bound in (32)
to choose one word with small `D` for all frames simultaneously.
More generally a growing number of frames is permitted whenever the
sum of their bounds in (32) tends to zero.

This includes the whole `k=o(n^(1/4))` range used by the localization
track's stronger Gaussian-comparison transfer, and a larger range for
the exact sign construction. It does **not** address the remaining
`sqrt(n)`-to-`o(n)` mode range. No conclusion about unknown minimizing
sequences is inferred from this explicit near-half family.

The script in Section 17 also checks the exact Parseval identity and
the expected defect bound over its complete finite probability law.

## 19. Uniform entropy from small response to an arbitrary isotropic sign law

**Status update (18:16 UTC).** The nuclear-norm necessity and finite
bound below remain valid. The director has since strengthened the
asymptotic entropy conclusion to `Psi(mu^2)~mu^2 log(1/mu^2)` using a
response-adapted, flat-diagonal Gaussian covariance; this agent has
independently reconstructed that proof. The BH/mechanism agent has also
removed the finite-support restriction from the realization paragraph
below: a sampled mode list simultaneously controls its covariance and
all matrix responses when `k^2=o(n)`. Thus the older sufficient
atom-count condition is not a necessary assumption of the current
combined theorem. These stronger results are recorded in their owners'
canonical artifacts; the argument below is retained as a separate
spectral obstruction and an auditable intermediate theorem.

**New response-to-nuclear-rank argument; independent audit requested.**
This removes the orthogonal-dictionary premise from the entropy step.
It applies uniformly to every isotropic sign-pattern law, including
nonuniform finite laws found by the director's exact linear programs.
The entropy conclusion does not depend on the number of atoms, their
smallest mass, or a selected factorization of the sign code.

**Subsequent improvement.** The director has now obtained the stronger
support-independent entropy bound
`inf_delta [mu^2/(2delta)+h(Q(1/sqrt(delta)))]`
by a trace-norm net of signal covariance and a regularized weighted-frame
Gaussian noise covariance with flat diagonal. This was independently
reconstructed here and passed. Its small-`mu` scale is
`mu^2 log(1/mu^2)`, improving (39). The finite nuclear inequality (36)
and its geometric rank consequence remain useful independently; the
entropy estimate below is retained as the earlier valid mechanism,
not presented as the strongest campaign conclusion.

### 19.1 A finite support-independent entropy theorem

Let `n=kp`, identify a Boolean vector with a `k by p` sign matrix `X`,
and let `nu` be any probability law on `{+-1}^k` such that
`E_nu vv^T=I_k`. Centering of `nu` is unnecessary. Define

```math
\Phi_\nu(X)=\frac{\mathbb E_\nu\|v^TX\|_2}{\sqrt n},
\qquad F_\nu(\mu)=\{X\in\{\pm1\}^{k\times p}:
                                  \Phi_\nu(X)\le\mu\}.
```

For every integer `1<=r<=k`, every `delta>0` and `w>0`, put
`t=r/k` and

```math
e=\left(\frac{\mu}{2\sqrt t}+2\delta\right)^2+w.
```

Whenever `e<=1/2`, the following finite bound holds:

```math
\boxed{\log|F_\nu(\mu)|
\le k^2\log(1+2/\delta)
 +\frac{nt}{2}\log(1+1/w)+n h(e).}               \tag{35}
```

Empty codes cause no issue. The proof has three elementary components.

First let `S=XX^T/n`, so `tr S=1`. Isotropy and Cauchy--Schwarz give

```math
\operatorname{tr}\sqrt S
=\mathbb E_\nu v^T\sqrt S\,v
\le\sqrt k\,\mathbb E_\nu\sqrt{v^TSv}
=\sqrt k\,\Phi_\nu(X).
```

Equivalently, the nuclear norm satisfies

```math
\boxed{\|X\|_*^2\le nk\,\Phi_\nu(X)^2.}          \tag{36}
```

Second, for any decreasing nonnegative sequence `a_i` with sum `L`,

```math
\sum_{i>r}a_i^2\le\frac{L^2}{4r}.               \tag{37}
```

Indeed, setting `b=a_(r+1)` gives tail at most
`b sum_(i>r)a_i<=b(L-rb)<=L^2/(4r)`. The factor `1/4` is sharp:
take exactly `2r` equal nonzero terms. Applying (37) to singular values
and using (36), the top-`r` left singular space approximates `X` with
squared Frobenius error at most `n mu^2/(4t)`.

Third, take an operator-norm `delta`-net of the orthogonal group `O(k)`
with at most `(1+2/delta)^(k^2)` points. Such a net follows directly by
disjoint-ball volume comparison in the `k^2`-dimensional normed space
of real matrices; its net points may be chosen orthogonal. There is
**no extra `log k` factor**, since the volume comparison uses operator
norm rather than a Frobenius-radius bound. An eigenbasis net point
within `delta` changes the top-`r` projector by operator norm at most
`2delta`. Its `rp`-dimensional image in the vectorized matrix space
therefore approximates `X` with squared error at most
`n(mu/(2sqrt(t))+2delta)^2`.

Label each codeword by such a net point, breaking ties deterministically.
Conditional on the label, observe its orthogonal projection plus
independent Gaussian noise of variance `w/t` on that `rp`-dimensional
space. The noise contributes total mean squared error `nw`, and
Gaussian maximum entropy gives conditional mutual information at most
`(nt/2)log(1+1/w)`. Coordinatewise nearest-sign prediction has average
error probability at most `e`; hence the remaining conditional entropy
is at most `n h(e)`. Add the label entropy to obtain (35).

### 19.2 Precise growing-dimension and energy-profile consequences

Suppose `k_n->infinity`, `k_n^2=o(n)`, and `nu_n` is **any** isotropic
sign law in dimension `k_n`, allowed to vary arbitrarily with `n` and
with the old signing. For every fixed `mu>0`, (35) implies

```math
\limsup_n\frac1n\log|F_{\nu_n}(\mu)|
\le\mathcal E(\mu),
\quad
\mathcal E(\mu)=\inf_{\substack{0<t\le1,\ w>0\\
                              \mu^2/(4t)+w<1/2}}
\left\{\frac t2\log(1+1/w)
                  +h\left(\frac{\mu^2}{4t}+w\right)\right\}. \tag{38}
```

For clarity, the order of limits is: fix `mu,t,w` with strict error
margin; choose a sufficiently small fixed `delta`; set
`r_n=ceil(t k_n)`; send `n` to infinity; then let `delta` decrease to
zero. The net label cost vanishes because `k_n^2/n->0`. Infimize only
after this argument. This is uniform in the sign law. For arbitrary
orders, use `p=floor(n/k)` and append the fewer-than-`k` leftover signs;
their entropy is `o(n)` and the same conclusion holds on the main array.

For `0<mu<=e^(-2)`, put `L=log(1/mu)`, `t=mu/sqrt(2)`, and
`w=mu/(2sqrt(2)L)`. Then

```math
\mathcal E(\mu)
\le\frac{\mu}{2\sqrt2}
          \log\left(1+\frac{2\sqrt2L}{\mu}\right)
 +h\left(\frac{\mu}{2\sqrt2}(1+1/L)\right)
=\left(\frac1{\sqrt2}+o(1)\right)
                         \mu\log(1/\mu).          \tag{39}
```

Thus, if the **full** absolute near-level code satisfies
`limsup_n max_(x in E_n(eta)) Phi_(nu_n)(X)<=mu(eta)` with
`mu(eta)->0`, its entropy profile is bounded by (38)--(39). Strict
majorants of `mu(eta)` handle the limsup exactly; they are sent down
only after the fixed-level order limit. In particular, writing

```math
L_\mu=\limsup_{\eta\downarrow0}
              \frac{\mu(\eta)\log(1/\mu(\eta))}{\eta},
```

the actual entropy slope obeys `K_s<=L_mu/sqrt(2)` whenever `L_mu` is
finite. This entropy implication is unconditional once the displayed
response bound is supplied.

To use it for sign augmentation, one still needs an **actual balanced
column realization** of the laws `nu_n`: aggregate subGaussian proxy
`(q+o(n))I` and uniform mean absolute response at most
`q sqrt(n) Phi_(nu_n)(X)+o(n^(3/2))`. Under those explicitly separate
construction certificates, the director's direct entropy theorem yields
every extension slope `tau>L_mu/(2sqrt(2))`. A liminf-realizing family
satisfying all those hypotheses must therefore have
`L_mu>=3sqrt(2)c_*`. No unrestricted realization of arbitrary enormous
supports is being assumed in this last step.

For a concrete sufficient realization condition, if `nu_n` has `R_n`
atoms and `k_n R_n=o(n)`, round its probabilities to integer mode counts
`q_r` summing to `q` with `|q_r-q nu_n(r)|<1`. The aggregate covariance
error is at most `k_n R_n I`; the uniform mean rounding loss is at most
`R_n sqrt(k_n n)=o(n^(3/2))`. Fair independent scalar signs on the
inner coordinates turn each sign pattern into actual full-sign bridge
columns. Thus finite-template laws grown sufficiently slowly satisfy
the required construction conditions. The entropy theorem itself has
no atom-count restriction.

## 20. Finite-state correction: a universal scalar-response entropy conjecture is false

The proposed bound

```math
\log|\{x\in\{\pm1\}^N:
      \mathbb E_\nu|x\cdot V|\le\delta\sqrt N\}|
\le C\delta^2N\log(1/\delta)+O(\log N)             \tag{40, false}
```

cannot hold uniformly over isotropic sign laws with absolute constants.
This is already false for a Walsh basis law. The example is a finite-state
correction to the additive term in (40), **not** an exponentially large
small-response code.

Let `N=2^(2r)` and identify coordinates with `F_2^r times F_2^r`.
Let `nu` be uniform on the antipodal Walsh characters of this group;
it has mean zero and covariance `I_N`. For every `r by r` binary matrix
`T`, define its graph subspace and a sign word

```math
S_T=\{(u,Tu):u\in\mathbb F_2^r\},\qquad
f_T=1-2\,1_{S_T}.
```

There are exactly `2^(r^2)` distinct graph subspaces and hence that many
distinct words. The normalized Fourier transform of `1_(S_T)` is
`2^(-r)` on its annihilator and zero elsewhere. Therefore, for `r>=1`,

```math
\mathbb E_\nu|f_T\cdot V|
=\|\widehat f_T\|_1
=3-2^{2-r}\le3.                                \tag{41}
```

Use `delta=3/sqrt(N)` and let `r` grow. The right side of (40) would be
`O(log N)`, but the left side is at least

```math
r^2\log2=\frac{(\log N)^2}{4\log2}.
```

Thus at least a squared-logarithmic additive/structural term is needed
in such a universal scalar conjecture. The classical fact that affine
subspace indicators have small Fourier `ell_1` norm underlies the
construction; see also Green--Sanders,
[*Boolean functions with small spectral norm*](https://arxiv.org/abs/math/0605524).
All numerical constants and the contradiction here follow from the
displayed elementary Fourier calculation, not an imported asymptotic
counting theorem.

Section 21 below rules out `exp(cN)` codewords with `delta->0` by proving
a corrected universal double-logarithmic bound. In particular this
example does not contradict the director's matrix-array entropy theorem:
its covariance-label cost vanishes only in the stated `k^2=o(n)`
aspect-ratio regime, which is absent when a single vector is treated
as an `N by 1` matrix.

The exact Walsh calculation for all binary graph maps through `r=3`
is reproduced in
`computations/paper_discrepancy_2026_09_17_scalar_entropy_counterexample.py`.

### 20.1 A squared-logarithmic additive repair alone also fails

The same construction can be amplified to refute the replacement of
`O(log N)` in (40) by an absolute `O((log N)^2)` term. Let `t=2^s`,
`N=t 4^r`, and index the coordinates by
`F_2^s times F_2^r times F_2^r`. Independently for each outer coordinate
`a`, choose an arbitrary binary `r by r` matrix `T_a`, and put

```math
U=\bigcup_{a\in\mathbb F_2^s}
       \{a\}\times\{(u,T_a u):u\in\mathbb F_2^r\},
\qquad f=1-2\,1_U.
```

The `t` slabs are disjoint, so there are exactly `2^(t r^2)` different
sign words. Each slab indicator is a product of two affine-subspace
indicators. Its normalized Fourier `ell_1` norm is one. The triangle
inequality consequently gives

```math
\mathbb E_\nu|f\cdot V|
=\|\widehat f\|_1\le1+2t,\qquad
\log|\mathcal F|=t r^2\log2.                    \tag{42}
```

Take `t=2^(floor((log_2 r)/2))`, so `t` is between `sqrt(r)/2` and
`sqrt(r)`, and set `delta=(1+2t)/sqrt(N)`. Then `delta->0`,

```math
\delta^2 N\log(1/\delta)=O(r^2),\qquad
(\log N)^2=O(r^2),\qquad
\log|\mathcal F|=\Omega(r^{5/2}).
```

Thus the failure is not repaired by a bare squared-logarithmic additive
term either. A valid universal bound, if one exists, must couple its
finite-state correction to the unnormalized response or use additional
structure. This remains a **subexponential** family in `N`. The replay
script also enumerates every two-slab choice at `r=2` and verifies
the exact word count and the response upper bound.

## 21. A support-free scalar response entropy theorem by importance encoding

The scalar question has a positive corrected answer, with one extra
logarithm. This result has **no tensor aspect-ratio assumption** and
applies to an arbitrary isotropic sign law, regardless of its support.

Let `nu` be any probability law on `{-1,1}^N` with
`E_nu hh^T=I_N`, and let `C` be a nonempty subset of the Boolean cube
such that

```math
a_x:=\mathbb E_\nu|h\cdot x|\le A=\delta\sqrt N
\quad\hbox{for every }x\in C.
```

Necessarily `A>=1`: since `|h dot x|<=N`, isotropy gives
`N=E(h dot x)^2<=N a_x`. For every integer `m` with
`epsilon_m=exp(-m/(2A^2))<=1/2`,

```math
\boxed{\displaystyle
\log|C|\le
 m\log(1+2/\delta^2)+N h(\epsilon_m).}          \tag{43}
```

In particular, if `0<delta^2<=1/2` and
`m=ceil(2A^2 log(1/delta^2))`, then

```math
\frac1N\log|C|
\le \delta^2\bigl(2\log(1/\delta^2)+1\bigr)
                    \log(1+2/\delta^2)+h(\delta^2)
=\bigl(8+o(1)\bigr)\delta^2\log^2(1/\delta).
                                                        \tag{44}
```

Thus for **every** sequence of isotropic sign laws, a uniformly
vanishing normalized mean response forces subexponential code size.
The Walsh-affine construction of Section 20 shows that the squared
logarithmic order in (44) cannot uniformly be replaced by a single
logarithm with only a dimension-logarithmic remainder. At
`delta=3/sqrt(N)` it gives order `(log N)^2` words in logarithmic scale,
matching the scale of (44), though not its numerical constant.

### 21.1 An exact response-adapted encoding law

Use as a common reference law `R=nu times Uniform{-1,1}` on pairs
`(h,sigma)`. For each word `x`, define a probability law `Q_x` by its
density relative to `R`:

```math
q_x(h,\sigma)
=1-\frac{a_x}{A}+\frac{2(\sigma\,h\cdot x)_+}{A}.
                                                        \tag{45}
```

It is nonnegative because `a_x<=A`. Averaging over `sigma` shows that
its mass is one. Its signed-vector mean is

```math
\mathbb E_{Q_x}\sigma h
=\frac1A\mathbb E_\nu h(h\cdot x)=\frac{x}{A}.
                                                        \tag{46}
```

The zero-mean reference contribution vanishes; this identity uses only
isotropy, not that `nu` itself is centered. A direct second-moment
calculation gives

```math
\mathbb E_R q_x^2
=1+\frac{2N-a_x^2}{A^2}
\le1+\frac2{\delta^2}.
```

Indeed, putting `b=2(sigma h dot x)_+`, one has `E_R b=a_x` and
`E_R b^2=2N`. Jensen's inequality under `Q_x` therefore gives

```math
D(Q_x\Vert R)=\mathbb E_{Q_x}\log q_x
\le\log\mathbb E_{Q_x}q_x
=\log\mathbb E_Rq_x^2
\le\log(1+2/\delta^2).                           \tag{47}
```

Zero density values are harmless and contribute zero to relative
entropy. The construction is an information-theoretic encoding, not
a proposed physical bridge law depending on the unknown query.

### 21.2 Information and coordinatewise reconstruction

Let `X` be uniform on `C`. Conditional on `X=x`, draw independent
`Z_1,...,Z_m` from `Q_x`, retaining both `h` and `sigma` in each sample.
Comparison with the common reference `R^m` gives

```math
I(X;Z_1,\ldots,Z_m)
\le\mathbb E_X D(Q_X^m\Vert R^m)
\le m\log(1+2/\delta^2).                         \tag{48}
```

For coordinate `i`, predict `X_i` by the majority of the `m` signs
`sigma_j h_{j,i}`; ties may be broken in any fixed way. Conditional on
`X=x`, these signs are independent and have mean `x_i/A` by (46).
Hoeffding's elementary bounded-variable inequality yields error at
most `exp(-m/(2A^2))`. Consequently,

```math
H(X\mid Z_1,\ldots,Z_m)
\le\sum_{i=1}^N h\bigl(\Pr\{\widehat X_i\ne X_i\}\bigr)
\le N h(\epsilon_m).
```

Adding this to (48) proves (43). For (44), the chosen integer `m`
has error at most `delta^2`; also
`m<=A^2(2log(1/delta^2)+1)` because `A^2>=1`. This absorbs integer
rounding without an additive dimension-dependent remainder.

All density normalizations, signed means, second moments, and the
relative-entropy comparison are replayed exactly (except logarithms)
for Walsh examples in
`computations/paper_discrepancy_2026_09_17_scalar_entropy_encoding.py`.
