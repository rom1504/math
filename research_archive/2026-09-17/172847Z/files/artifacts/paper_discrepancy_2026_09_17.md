# Biased Gram--Schmidt rounding and cumulative-precision augmentation

2026-09-17 paper-combination campaign. **Sections 3--7 independently audited
PASS by the campaign's BH/mechanism agent, including the primary sharp-GS
proof and exact 9/8 constant.** No claim of convergence or improvement to the
reported global cap interval is made.

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
