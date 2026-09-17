# Wave 48 main-agent audit notes

## Fixed-density complement product wall

The Wave 48A argument was independently reconstructed against (10.1022),
(10.1186), and (10.1190).

For `a=2p_2-1>=0`, the definitions give exactly

```math
L_S(d)=2X_d(S)+aE_d.
```

On `theta Gamma_S>(1-theta)q`, one has
`L_S(d)=q+Gamma_S>q/theta`.  Since `E_d<=q`, `theta<=1`, and
`p_2<=p_1^2`,

```math
X_d(S)>{1\over2}(theta^{-1}-a)q
\ge (1-p_1^2)q.
```

Thus the unconditional high-slack mass is exactly `alpha_d beta_d` and is
bounded by the fixed-slice upper tail (10.1022).  Its three exponents reduce,
up to fixed constants, to

```math
q^2/R_2(d),\qquad n,\qquad n^{3/4},
```

using `q>=c_0n^(3/2)` and `||A||_op<=sqrt(2q)`.

In fact the same calculation with only `Gamma_S>=0` gives the stronger
direct row--mass wall on every compact fixed-density window

```math
\alpha_d\le 3(n+1)\exp\left\{-c\min\left(
{q^2\over R_2(d)},n,n^{3/4}\right)\right\}.
```

Indeed `S in I_d` and `|E_d|<=q` imply
`X_d(S)>=(1-|a|)q/2=min(p_2,1-p_2)q>=delta_0q`. Hence every cut with
`R_2(d)=O(n^(9/4-c))` has `-log alpha_d=Omega(n^(3/4))`, with no scalar
optimality or high-slack hypothesis.  This by itself rules out the required
project-row/project-saving-mass pair for all fixed-density complement
columns.  The product argument below additionally explains why using
arithmetic high-slack retention to manufacture the row bound cannot evade
the wall.

The Cantelli lower bound has no hidden exponential loss.  The positive
margin is a nonzero rational with denominator dividing `n(n-1)`, hence at
least `1/[n(n-1)]`.  Principal monotonicity gives `Q(A[T])<=q`, so
`|L_S(d^U)|<=Q(B_S)<=q+2Q(A[T])<=3q<=3n(n-1)`.  Therefore each high-slack
incidence contributes at least a fixed multiple of `n^(-8)` to the survival
factor.  Consequently

```math
R_2(d) <= C[n^2+n^(3/2)(-log beta_d+O(log n))]
```

for the same scalar optimizer.  If `-log beta_d>=n^(3/4)` the product wall
is immediate; otherwise this row bound is `O(n^(9/4))`, and every exponent
in (10.1022) is `Omega(n^(3/4))`.  Hence

```math
-log(alpha_d beta_d)>=c n^(3/4).
```

All inequality directions, normalizations, and uniformity requirements check
out.  In particular, this contradicts two simultaneous
`O(n^(3/4-c))` logarithmic costs for every fixed `c>0`.  It does not apply to
the bare arbitrary-cut event, whose threshold includes the restriction gap
and tolerance rather than the much stronger complement condition `L_S>=q`.

## Immediate strategic consequence

The structured fixed-density arithmetic complement implementation is
retired, not merely demoted.  A successor cannot use constant-scale slack to
regularize the same scalar column while retaining project-saving mass.  The
bare arbitrary-cut lemma remains the sharp sufficient target; direct
profile-conditioned counting and hard exceptional-center/coarea mechanisms
must be compared after the other Wave 48 audits.

## Row-truncated coarea sharpening

The Wave 48B harmonic calculation and checker were independently rerun. On
the `m`-slice, the covariance of the centered coordinates on the sum-zero
subspace is `p(1-p)n/(n-1)`. If `f=1_F`, `a=E f`, and
`mu_i=Pr(i in S|F)`, orthogonal projection therefore gives exactly

```math
W_1={ (n-1)a^2\over p(1-p)n}\sum_i(\mu_i-p)^2.
```

The hypersimplex constraint `0<=mu_i<=1`, `sum mu_i=np` makes the last sum
at most `np(1-p)`, so `W_1<=(n-1)a^2`. Writing
`delta=1-lambda_1`, `g=lambda_1-lambda_2`, and `kappa=g/delta`, the spectral
decomposition then yields

```math
{B\over\delta a}\ge1+\kappa-(1+\kappa n)a.
```

The inequality remains valid when the right side is negative. Weighting it
by `1_Ca^2/E[1_Ca^2]` verifies, for every selector-independent truncation
`C`,

```math
{E[1_CaB]\over\delta E[1_Ca^2]}\le1-\eta
\quad\Longrightarrow\quad
\max_{z\in C}a_z\ge{\kappa+\eta\over1+\kappa n}.
```

Thus a non-strict coarea bound (`eta=0`) already extracts degree
`Theta(1/n)` whenever `kappa` is bounded below. The symbolic formula for
`ell=m-s` checks and gives `kappa=1-o(1)` for fixed selector density and
`s=o(n)`. This rigorously removes the exponentially small strict-margin
obligation from the live coarea target.

The denominator `E[1_Ca^2]>0` is essential. For `A=J-I`, positive ground
degree at density `p>1/2` forces spin imbalance at least `(2p-1)n`, while
`R_2(z)=n+(n-2)(2k-n)^2=Omega(n^3)`. Hence every ground lift lies outside
the project cap. This is a scalable signing obstruction to a law-free mass
claim, but not an exact-minimizer obstruction.

### FKN upgrade

There is a further dimension-free sharpening. If the aggregate non-strict
inequality holds and `E[1_Ca^2]>0`, some `z in C` satisfies
`B_z<=delta a_z`. For this individual Boolean family the exact spectral
identity cancels level one and gives

```math
\sum_{j\ge2}(\lambda_1-\lambda_j)W_j\le\delta a_z^2,
\qquad
W_{>1}\le a_z^2/\kappa.
```

On a fixed-density slice and at a scale with `kappa>=kappa_0>0`, the
balanced-slice FKN theorem (Filmus, *FKN theorem for the multislice, with
applications*, arXiv:1809.03089, Theorem 1) supplies a Boolean dictator `g`
with

```math
Pr(f_z\ne g)\le C_{p,\kappa_0}a_z^2.
```

The possible means of a Boolean function depending on one coordinate are
`0,1,p,1-p`. If `g=0`, the left side equals `a_z`, forcing
`a_z>=1/C`. For the other three possibilities,
`Pr(f_z ne g)>=|a_z-Eg|`; when `a_z` is below half the fixed slice-balance
constant this is bounded below by a fixed positive number and again
contradicts `Ca_z^2` for sufficiently small `a_z`. Therefore

```math
a_z\ge c(p_0,p_1,\kappa_0)>0.
```

Subject to the cited theorem's fixed-density and sufficiently-large-`n`
hypotheses, non-strict truncated coarea plus nonzero mass thus extracts a
constant-degree center, not merely degree `Theta(1/n)`. The remaining live
obligations are still the existence of a project-row favorable incidence
and the non-strict boundary inequality; FKN does not supply either one.

### Project-row extension audit

For a child ground `y` on `S`, uniform completion on `T=S^c` gives exactly

```math
E_w R_2(y,w)=||A[:,S]y||_2^2+|T|(n-1).
```

The full completion energy has orthogonal constant, linear, and quadratic
Walsh parts, so

```math
Q(A[S])^2+4||A[T,S]y||_2^2+2|T|(|T|-1)\le q_n^2.
```

These verify the sufficient project-mass package: one child ground with
internal Gram `O(R_*)` and principal shortfall
`q_n-Q(A[S])=O(R_*/q_n)` has a completion of row at most `O(R_*)`. At
`R_*=n^(9/4-c)` the shortfall target is `O(n^(3/4-c))`.

Known one-spin stability gives only
`||A[S]y||^2<=(m-1)Q(A[S])=O(n^(5/2))`; the minimizer operator estimate gives
the same scale for the full column. Both miss the project cap by the factor
`n^(1/4+c)`. Thus the extension calculation is an exact conditional bridge,
not a proof of nonzero project mass.

## Cross-block tight-decomposition audit

The Wave 48C checker was independently rerun and passed.

For `R subset S`, optimizing a parent completion in two stages gives exactly

```math
k_R(y_R)=\max_{y_S|_R=y_R}\{k_S(y_S)+c_S(y_S)\}-c_R(y_R),
\qquad
s_R(y_R)=\min_{y_S|_R=y_R}s_S(y_S).
```

Thus active block supports are projections of the one global ground
relation, although block-dependent dual laws need not align.

If maximal `m`-selectors have nonnegative weights of total `Lambda` and
constant edge load `r`, double counting one parent ground gives
`sum_S lambda_S c_S(d[S])=r q_n`. Failure of every tight pair lowers each
term by at least four, so `r q_n>Lambda(q_*-4)` proves a tight pair. The
normalization and modulo-four gap check.

The universal `m=3` theorem also checks. For an oriented parent ground,
write `P=sum_e b_e=q_n/2` and `N=binom(n,2)`. Orthogonality gives
`P>=sqrt(N)`, hence the positive-edge graph has
`(N+P)/2>n^2/4` edges for `n>2`. Mantel supplies an all-positive triangle,
whose inherited energy is the universal triangle norm six. The `m=2` case
is immediate from a positive edge.

For a block fractional dual, lift every child state to a completion attaining
`k_S`. If `t=E s`, the expected total signed edge sum is `(q_n-t)/2`.
With `K=binom(n,2)-binom(m,2)` outside edges,

```math
\sum_{e\in S}(E M_e)_+\ge((q_n-t)/2-K)_+,
```

and minimizing the dual objective verifies

```math
\gamma_{frac}(S)\ge(q_n-2K)_+.
```

Hence the zero-dual antipodal mechanism is impossible for
`n-m=o(sqrt(n))`, but a positive fractional/integral gap remains possible.

The exact A8 audit confirms that this caveat is real at finite scale. Its two
bad complementary maximal four-blocks have disjoint zero-dual ground faces;
the common ground-face minimax is `1/2`. A four-edge pressure certificate
lowers all old grounds from 20 to 12, while eight slack-eight states rise
from 12 to 20. Thus ground-face cross-block pressure alone cannot produce an
integral descent; a valid successor must control positive slack layers.

The exact order-ten MILP sample (twelve deterministic, distinct minimizers)
was reproducible but is correctly labeled nonexhaustive. All samples retain
an all-size common parent.
