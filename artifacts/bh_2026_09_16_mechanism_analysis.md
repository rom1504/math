# BH mechanism audit: exact joint-field consequences and boundaries

Date: 2026-09-16. Independent mechanism pass, then comparison with the
project's existing lower-bound and high-order-chaos artifacts. This note
contains proved finite identities, proved scoped obstructions, and clearly
marked numerical diagnostics. It does not improve the current
`.4333221116640807` lower bound or prove convergence of the original minima.

## 1. Paper reconstruction and the actual point of contact

Primary source: Paata Ivanisvili, [Polynomial growth of Bohnenblust--Hille
constants on the Hamming cube, arXiv:2609.12427v1](https://arxiv.org/html/2609.12427v1).
The entire paper, including proofs and the flat-level influence corollary,
was read. Its mechanism is as follows. A weighted square function sums
the degree-r coefficient `ell_(2r/(r+1))` norms with squared weights
`r^-10`. Random coordinate bipartitions retain balanced portions of high
levels. Mixed-norm interpolation, hypercontractivity of whole coefficient
rows, and the weights together contract the high-level part. The supremum
is taken over a restriction-closed class at fixed ambient dimension.
Low degrees through eleven are handled separately; their Markov extraction
cost and removal of weights give polynomial degree growth. Level-flat
magnitudes then identify every coordinate's influence and yield the stated
polynomial influence lower bound.

The useful question here is not to substitute degree two into that
high-degree theorem. Freezing coordinates of a full sign quadratic creates
an inhomogeneous quadratic with a random linear field. Thus the most
literal possible new mechanism is a restriction-closed certificate that
retains BOTH its linear and quadratic contributions. The rest of this
note investigates that possibility quantitatively.

Notation:

```
Q_A(x)=sum_(i<j) a_ij x_i x_j,
W(A)=(max_x Q_A(x)-min_x Q_A(x))/2,
M(f)=max_x |f(x)|,
f(x)=c+Q_A(x)+h dot x,
H=sum_i |h_i|.
```

The local comparison files were read only after the independent initial
derivations. `high_order_flat_chaos.md` already contains the hafnian/parity
power obstruction, so that observation is not claimed as new. The current
half-range constant is justified in
`decisive_audit_certified_minimum_width_lower_2026_09_07.md`, using the
marked-response proof reconstructed in
`decisive_audit_fresh_full_lower_chain_2026_09_07.md`.

## 2. A full-sign, natural-field-scale obstruction

An additive claim such as

```
M(c+Q_A+h dot x)^2 >= H^2 + positive_quadratic_coefficient_term
```

is false even with every off-diagonal coefficient a sign, `c=0`, and
every `h_i` of order `sqrt(n)`. This is stronger than the elementary
two-variable example with a freely chosen constant.

### 2.1 Exact Paley plateau

Let n be an odd prime with `n=1 mod 4`, and let
`a_ij=chi(i-j)`, where chi is its quadratic character, extended by
`chi(0)=0`. The elementary character-sum identities give

```
A 1=0,             A^2=n I-J,             ||A||op=sqrt(n).
```

For a Boolean x put `s=sum_i x_i`. Orthogonally removing the constant
component yields

```
|Q_A(x)| <= sqrt(n)/2 * (n-s^2/n).
```

Consequently, for `h=t sqrt(n) 1`,

```
M(Q_A+h dot x)/n^(3/2)
 <= max_(0<=u<=1) [(1-u^2)/2+t u]
 = (1+t^2)/2,                   0<=t<=1,
 = t,                          t>=1.
```

For `t>=1`, equality holds because `x=1` has `Q_A(1)=0` and linear
value `t n^(3/2)`. In particular, at t=1 the full quadratic term costs
NOTHING beyond `H=n^(3/2)`.

This does not model an independently Gaussian external field, and does
not prove such a plateau occurs with substantial probability under an
actual random restriction of a near-minimizing whole signing. Its exact
scope is the failure of a uniform field-independent joint certificate.

### 2.2 The entire upper parabola is asymptotically attained

The preceding spectral curve is not merely an upper bound with unknown
slack. Let p be an odd prime, choose a nonsquare d modulo p, and index
an order `n=p^2` signing by pairs `(a,b)` in `F_p^2`. Define

```
A_((a,b),(c,e))=chi_p((a-c)^2-d(b-e)^2).
```

The argument of chi is zero only on the diagonal. This is the quadratic
character of `F_(p^2)` expressed through the field norm, so again
`A1=0`, `A^2=nI-J`, and `||A||op=p`.

There is also a direct finite block calculation. Group coordinates by b.
A diagonal block has every off-diagonal entry +1, hence row sum p-1.
Each off-diagonal block has row sum -1, by
`sum_a chi_p(a^2-D)=-1` for nonzero D. Thus A acts on row-constant
vectors through `p I_p-J_p`. In particular every centered row-constant
vector is a +p eigenvector.

Make k of the p whole rows negative and the others positive. Its
magnetization density and exact normalized energy are

```
u=1-2k/p,           Q_A(x)/n^(3/2)=(1-u^2)/2.
```

For each `0<=t<=1`, choose the grid point u nearest t. Its distance
is at most 1/p, and its perturbed energy is

```
(1-u^2)/2+t u=(1+t^2)/2-(u-t)^2/2.
```

Therefore these FULL signings satisfy

```
(1+t^2)/2-1/(2p^2)
 <= M(Q_A+t sqrt(n) sum_i x_i)/n^(3/2)
 <= (1+t^2)/2,                         0<=t<=1.
```

For t>=1 the value is exactly t. This provides an asymptotically sharp
parabolic obstruction on the whole natural field scale. For example,
the stronger matrix-specific proposal `M(Q+h)^2>=M(Q)^2+H^2` fails
asymptotically for every fixed `0<t<1`: its proposed lower curve is
`sqrt(1/4+t^2)`, strictly above `(1+t^2)/2`.

## 3. Exact low-field completion: a genuine joint value bound

The obstruction above does NOT prevent all additive gain. For every
subset T of coordinates, with no condition on its size or field mass,

```
M(c+Q_A+h dot x) >= H-sum_(i in T)|h_i|+W(A_T).          (1)
```

This is a finite, deterministic inequality for arbitrary real c and h
and arbitrary real symmetric hollow A. No random independence is used.

Proof. Switch coordinates so that h is nonnegative. Global reversal
gives the exact parity identity

```
M(c+Q_A+h dot x)=max_x (|c+Q_A(x)|+|h dot x|).
```

Fix every coordinate outside T at +1. For z on T, the even expression
is `u+v`, where `u=c'+Q_(A_T)(z)` and v is linear in z. Write
`H_out=sum_(i notin T)h_i` and `w=h_T dot z`. The two choices z,-z
have expressions

```
|u+v|+|H_out+w|,          |u-v|+|H_out-w|.
```

Their maximum is at least their average, which is at least
`|u|+H_out`, by the triangle inequality separately for both pairs.
Maximize over z. Since
`max_z |c'+Q_(A_T)(z)|>=W(A_T)`, (1) follows.

In particular, fields with many very small coordinates admit a strictly
positive completion gain over H. Uniform natural-scale fields need not,
as Section 2 shows.

## 4. Genuine random restrictions: dependence-safe bootstrap

Assume a universal asymptotic half-range lower bound

```
W(B)>= (w-o(1)) k^(3/2)
```

for every complete signing B of order k. Split a signing of order n into
a fixed set I of `alpha n+o(n)` coordinates and its complement J.
Freeze y on J as independent fair signs. The restricted field is

```
h_i=sum_(j in J) a_ij y_j,         i in I.
```

Every h_i separately has exactly the simple random-walk distribution
with |J| steps. Distinct h_i can be arbitrarily dependent; no joint
Gaussian approximation is required below.

For fixed t>=0 choose the ACTUAL subset
`T(y)={i in I: |h_i|<=t sqrt(|J|)}`. Let

```
p(t)=P(|G|<=t)=erf(t/sqrt(2)),      kappa=sqrt(2/pi).
```

Applying (1) for every y and averaging gives the lower certificate

```
M(Q_A)/n^(3/2)
 >= kappa alpha sqrt(1-alpha) exp(-t^2/2)
       +w [alpha p(t)]^(3/2)-o(1).                    (2)
```

The same holds for W(A): apply (1) to `c+Q_A` for every c and then
minimize c. Here are the limit details:

* The retained field mass uses ONLY the scalar central limit theorem
  and uniform integrability from its exact second moment.
* `E |T|/n -> alpha p(t)` follows from the marginal law and linearity
  of expectation, without concentration of |T|.
* Jensen supplies `E |T|^(3/2)>=(E |T|)^(3/2)` even with fully
  correlated rows.
* The half-range bound applies to each realized principal signing A_T.
  For each fixed epsilon choose a finite threshold K above which its
  error is at most epsilon. The finitely many orders below K contribute
  only an O(K^(3/2)) uniform correction. Divide by n^(3/2), then send
  epsilon to zero. Selection of T does not invalidate this uniformity.

Thus this is an actual value implication, not merely an informal analogy
with the paper's restriction closure.

## 5. A rigorous ceiling below the existing lower constant

The bootstrap (2), even iterated, cannot improve the current w. Define
its critical value by

```
w_crit=sup_(0<alpha<1,t>=0)
  [kappa alpha sqrt(1-alpha) exp(-t^2/2)]
       /[1-(alpha p(t))^(3/2)].                        (3)
```

Numerical optimization gives approximately

```
w_crit=0.3142961430858171,
alpha=0.7011962438,       t=0.4050494263.
```

These last decimals are diagnostics, not an interval-certified optimum.
No numerical optimization is needed for the decisive rigorous statement

```
w_crit <= kappa/2 = 1/sqrt(2pi)
       =0.3989422804014327 <0.4333221116640807.          (4)
```

To prove it, a Gaussian square is contained in the disk of radius
sqrt(2)t, so

```
p(t)^2 <= 1-exp(-t^2).
```

It remains to show, for all alpha,p in [0,1],

```
2 alpha sqrt((1-alpha)(1-p^2))+(alpha p)^(3/2)<=1.       (5)
```

Put `u=1-p^2`. Concavity gives
`p^(3/2)=(1-u)^(3/4)<=1-3u/4`. For alpha>0, complete the square:

```
2 alpha sqrt((1-alpha)u)-(3/4)alpha^(3/2)u
 <=(4/3)sqrt(alpha)(1-alpha).
```

Hence the left side of (5) is at most
`sqrt(alpha)(4-alpha)/3`, which is increasing on [0,1] and equals one
at alpha=1. The alpha=0 case is immediate. This proves (4).

If w is at least kappa/2, (2)'s displayed right side is never larger
than w: replace its linear term by
`(kappa/2)[1-(alpha p)^(3/2)]`. In particular the known half-range
constant is a fixed-point upper barrier for this entire bootstrap.

## 6. Soft reserves and the precise coefficient-norm boundary

There is an exact continuum extension of (1). Let `0<=r_i<=1`, let
`D_r=diag(r_i)`, and first switch coordinates so h>=0. Then

```
M(c+Q_A+h dot x)
 >=sum_i h_i(1-r_i)+W(D_r A D_r).                     (6)
```

Indeed put `m_i=1-r_i` and use the two cube points `m+D_r z` and
`m-D_r z`. The parity identity remains valid on the continuous cube,
whose supremum equals the vertex supremum by multilinearity. The same
two triangle inequalities as in Section 3 leave
`h dot m+|c+Q_A(m)+Q_(D_r A D_r)(z)|`. Maximization gives (6).

This WOULD expose more information if the weighted residual width could
be bounded using its actual signing and the field that generated r.
But if the residual is bounded only by a BH-type coefficient norm,
soft reserves provably add nothing. For any p>=1 and positive C, the
resulting objective is

```
Phi(r)=sum_i h_i(1-r_i)
       +(1/C)[sum_(i<j)|a_ij r_i r_j|^p]^(1/p).
```

Fix every coordinate except r_i. Its nonlinear term has the form
`(b+a r_i^p)^(1/p)`, a norm of an affine two-coordinate vector, hence
is convex in r_i. Adding the linear term preserves convexity. Replacing
each r_i successively by one endpoint never decreases Phi. Therefore
its global maximum on `[0,1]^n` is attained at a hard subset.

This applies in particular to p=4/3, and is exact at finite order. A
genuine degree-two inhomogeneous BH bound may be minimized over its
constant coefficient to give the required residual width bound. The
statement does not assume an unavailable sharp BH constant.

Scope: this rules out a particular natural soft-reserve completion using
ONLY coefficient-norm residual data. It does not rule out every joint
square-function inequality, all optimizer-specific arguments, or
residual widths exploiting correlations with the original random field.

## 7. Reproducibility and status

Source: `computations/bh_mechanism_2026_09_16_joint_fields.py`.

```
.venv/bin/python computations/bh_mechanism_2026_09_16_joint_fields.py \
  --output /home/math/quadra/tmp/bh_2026_09_16/bh_mechanism_joint_fields_replay.json
```

The replay checks exact integer Paley identities at orders
5,13,17,29,37,101; exhaustive field curves through order17; square-order
Paley row witnesses at 9,25,49,121,289; 1,984 hard-reserve and 1,600
soft-reserve finite inequalities; coordinate convexity diagnostics; and
the scalar optimization. All checks PASS. Random tests supplement the
proofs and are not their justification.

The paper's flat-level influence corollary applied directly to the
normalized original quadratic supplies no new structure: variance is of
order 1/n and each influence of order 1/n^2. Applying it to nonlinear
selectors or powers would require a separate proof of their level-flat
coefficient magnitudes, which is not true automatically. This audit
therefore supplies exact useful boundaries and a genuine but exhausted
restriction bootstrap, not a new original-value theorem.

## 8. Actual balanced construction: weighted-port extension

**Scope first.** The following is a new deterministic relative
core--tail inequality. Its application extends the sparse-active theorem
to arbitrarily many nonconstant fibres, but only when a sparse core
carries almost all variance. Consequently total variance is o(N), and
the existing global operator estimate already makes the absolute energy
o(N^(3/2)). This is NOT progress on the leading upper bound or the dense
mixed-profile gap. It may be useful where a relative variance coefficient,
rather than an additive cap error, is required.

The construction and previously proved frame estimates were reconstructed
from `principle_synthesis_2026_09_07_balanced_bulk_sparse_active_bound.md`
(including all of Section 6),
`principle_director_joint_packet_audit_2026_09_07.md`, and
`principle_invent_2026_09_07_stratified_marked_selector.md`.

### 8.1 Exact theorem for the actual reciprocal sign bulk

Let m physical fibres each have q coordinates. At fibre i let T_i be a
q-by-m real matrix with balanced columns; omitted/self/masked columns
may be zero. Let S_ij=S_ji be arbitrary signs on the remaining edges,
and form the actual block matrix

```
W_ij=S_ij T_i(:,j) T_j(:,i)^T,        W_ii=0.
```

Write `P_i=I-J/q`, `y_i=P_i x_i`, and define directed ports
`h_i(j)=T_i(:,j)^T y_i`. Assume that for some `0<=K<=M`, and a fixed
set C of at most s fibres,

```
||T_i||op^2<=M,       ||T_i(:,C)||op^2<=K       for EVERY i.    (7)
```

The letter M here is a column-map squared operator bound, not the
original minimum M_n. Define

```
V_C=sum_(i in C)||y_i||^2,       V_L=sum_(i notin C)||y_i||^2.
```

Then, simultaneously for either energy sign,

```
|H_W(x)| <= M(V_C+V_L)/2
          -(sqrt((M-K)V_C)-sqrt(K V_L))_+^2/2.          (8)
```

This holds for arbitrary real x, in particular every physical Boolean
word, every outer seed, and every symmetric mask. If (7) holds uniformly
over all sets of at most s fibres, C may be selected AFTER seeing x.
No independence or entropy estimate is used in the deterministic theorem.

### 8.2 The weighted square function and exact minimization

For arbitrary positive weights lambda_i, define

```
S_lambda(x)=sum_i lambda_i^(-1) sum_j lambda_j h_i(j)^2.
```

Weighted Young on each reciprocal pair gives

```
|H_W(x)| <= S_lambda(x)/2.                             (9)
```

This is a literal value inequality for the same matrix; the weights do
not change or discard its cross interaction. Take lambda_i=t on C and
lambda_i=1 elsewhere, where t>=1. For every row,

```
sum_j lambda_j h_i(j)^2
 =||h_i||^2+(t-1)||h_i restricted to C||^2
 <=[M+(t-1)K]||y_i||^2.
```

Consequently

```
|H_W(x)| <= [M+(t-1)K](V_C/t+V_L)/2.                 (10)
```

When `K V_L>0`, its minimizer over t>=1 is

```
t=max(1,sqrt((M-K)V_C/(K V_L))).
```

Expanding (10) at that value gives exactly (8). Zero-variance or K=0
cases follow directly or by limits. Thus the theorem pays the cross
term jointly, instead of adding separately certified phase caps.

This is an adaptation of weighted-square-function bookkeeping prompted
by the new paper. It is an independently proved port-space theorem,
NOT a direct application of that paper's Fourier-degree theorem.

### 8.3 Application to the existing conditioned actual selector

The already constructed operator-controlled balanced repair and its
conditioned restricted-Gram event give, uniformly across all physical
fibres and every column set of size at most s,

```
M=m+o(m),      K=q+o(m),      q/m -> p=31/32,
s log^4(m)/m ->0.
```

The global bound follows also directly from the full Hadamard norm
sqrt(m), orthogonal centering, and the o(sqrt(m)) operator repair.
The restricted bound requires that operator-controlled repair; arbitrary
individual balancing flips are not silently substituted for it.

Put `N=mq`, `V=V_C+V_L`, and `theta=V_C/V` when V>0. For theta>=p,
(8) gives the actual uniform RELATIVE estimate

```
|H_W(x)|/(sqrt(N)V)
 <= B_p(theta)+o(1),
B_p(theta)=[1-(sqrt((1-p)theta)-sqrt(p(1-theta)))^2]
                                      /(2sqrt(p)).     (11)
```

For theta<=p the displayed positive-part version is simply the old
global coefficient `1/(2sqrt(p))`. The errors in (11) multiply V itself,
even when V is microscopic; no division by an uncontrolled additive
error has occurred. B_p is decreasing on [p,1], and
`B_p(1)=sqrt(p)/2`, recovering the sparse-active endpoint.

At theta=999/1000 and p=31/32 an exact rational square-root enclosure
proves

```
B_(31/32)(999/1000) < 497237/1000000
                   = .497237 < .498.
```

The unrounded value is about `.49723661580914424`. Thus the SAME actual
conditioned ensemble eventually has coefficient .498 uniformly on every
word for which some set of at most s fibres carries at least 99.9% of
the variance. There is no restriction on the number, densities, spectral
coherence, or relative scales of nonconstant fibres in the remaining tail.
The numeric threshold for coefficient .499 is approximately
theta=.9980469950928259; that optimized threshold is only diagnostic.

The limitation is mathematically immediate:

```
V_C<=s q,     theta>=.999   =>   V<=s q/.999
                                      =o(N/log^4(m)).
```

Hence this theorem does not control a previously uncovered macroscopic
variance sector. It only strengthens relative control of a microscopic
mixed-tail sector, relevant if a later composition needs the variance
coefficient itself.

### 8.4 Why multilevel weights do not help with only a step RIP input

Suppose the ONLY row-square information retained is the global cap M
and the uniform restricted cap K on every set of at most s columns.
For arbitrary positive lambda let a be its largest entry and b its
(s+1)-st largest entry. A layer-cake bound gives

```
sum_j lambda_j h_i(j)^2
 =integral_0^infty sum_(j:lambda_j>u) h_i(j)^2 du
 <=[M b+K(a-b)]||y_i||^2.                             (12)
```

For fixed a,b, at most s indices have weights above b. Replacing all
those weights by a, and all other weights by b, can only DECREASE
`sum_i ||y_i||^2/lambda_i`, while the coefficient in (12) stays fixed.
Therefore the best certificate obtainable through (12) is already a
two-level potential. Its high-weight set consists of the s largest
variance fibres, and its exact optimization is (10).

In particular, assigning many polynomial or logarithmic weights to
minority-density strata does not itself produce a new contraction from
these inputs. A stronger actual column-profile estimate, or a jointly
paid packet/port constraint beyond the step RIP bounds, is needed for
the remaining dense mixed branch. This is a limit of the specified
certificate, not a counterexample to such a stronger theorem.

### 8.5 Replay

`computations/bh_mechanism_2026_09_16_weighted_ports.py` reuses the
existing exact Walsh builder. It checks 1,920 weighted inequalities
on actual individually balanced sign-frame weaves at m=8,16,32,64,
with each frame's true finite global/restricted Gram norms; it also
checks the analytic minimization in 500 numerical instances. The
asymptotic frame estimates are inherited from the existing proved
conditioned-selector construction, not extrapolated from these small
examples. Integer square roots give the exact rational phase certificate.
All checks PASS.

## 9. Adjacent multiscale attempt: a sharp temperature-mixture obstruction

The exact one-hole packet curve from
`principle_construct_2026_09_07_rare_packet_temperature_curve.md` has a
useful extension. For every dephased Hadamard H_L, L>=4, and EVERY
c>0, not only c in [1,L],

```
rho(c)=min_(nonempty proper S)
 [sum_j min(C_j(S)^2/c,1)-1]/|S|
 =min(1/c,2/L,L/c-1).                                 (13)
```

The old proof gives the same lower bound without a restriction on c:
if DC is not clipped, Parseval gives at least L/c-1; if it is clipped,
the non-DC estimate gives at least min(1/c,2/L). The old equality
packets also cover all new temperatures. A half-row packet attains 2/L
for c<=L/2; a complement-of-one packet attains 1/c on [L/2,L-1]; a
singleton attains L/c-1 for c>=L-1, including c>L.

For the one-hole node k=L-1, put `ell=log(1/r)` and `tau=t/ell`.
The leading optimized rare-row exponent consequently equals

```
Gamma_t(1-2r)/(r ell)
 = max(-4 tau/L, -2k/L^2, k/L-4 tau)+o(1).             (14)
```

The maximum of these affine functions is convex in tau. Therefore any
fixed finite common nonnegative temperature mixture with weights w_j summing to one
satisfies, at leading order and for EACH density r,

```
sum_j w_j Gamma_(t_j)(1-2r)
 >= Gamma_(sum_j w_j t_j)(1-2r)+o(r ell).              (15)
```

The Markov variance term is linear in temperature, so it has the same
mean on both sides. Thus taking common convex combinations of the
existing row-exponent certificates cannot outperform using their mean
temperature. In particular it cannot enlarge the already sharp
less-than-four logarithmic density bandwidth merely by a weighted
average or a Holder interpolation of those same inputs.

This is another specified-certificate obstruction, not a proof that a
new jointly conditioned multiscale estimate is impossible. It identifies
why the paper's degree-weight contraction has no immediate analog here:
the available scalar temperature potential has the wrong convexity for
the proposed interpolation gain, and the genuinely stronger physical
joint constraint is still missing.

The weighted-port replay includes exact packet checks at L4,L8 with
positive c both below one and above L, and numerical Jensen regression.
The proof, not these tests, establishes the all-c statement.

Canonical generated replay outputs are preserved as
`computations/results/bh_mechanism_2026_09_16_joint_fields_replay.json`
and `computations/results/bh_mechanism_2026_09_16_weighted_ports_replay.json`.

## 10. Stronger column-profile inputs: exact constraints and a frame-uniform obstruction

### 10.1 A rank/trace floor valid for every actual balanced sign frame

Let T have q rows and D active balanced sign columns. Thus each column
has squared norm q and `rank(T)<=q-1`. For `2<=s<=D` there exists an
s-column set J with

```
||T(:,J)||op^2
 >=q+(s-1)q(D-q+1)/[(D-1)(q-1)].                     (16)
```

Proof. Put `C=T^T T-q I_D`. It has zero diagonal, zero trace, and
`C>=-q I`. The rank/trace inequality gives

```
||C||F^2=||T^T T||F^2-Dq^2
 >=Dq^2[D/(q-1)-1].
```

The average squared Frobenius norm of an s-principal submatrix is
`s(s-1)||C||F^2/[D(D-1)]`. For any such submatrix, its eigenvalues lie
between -q and `t=lambda_max(C_J)`. Summing
`lambda^2<=(t-q)lambda+qt` and using its zero trace yields
`||C_J||F^2<=q s t`. Choose a submatrix at least as large as the average
to obtain (16).

For the actual construction `D=m-1`, `q/m->p`, and `s/m->delta>0`,
this forces

```
K(s)>=q+(1-p)delta m+o(m).
```

Therefore a positive-density extension retaining the exact same
`q+o(m)` matching-scale restricted Gram bound is impossible even for
actual balanced sign columns. This does NOT rule out a weaker linear-
support bound with a still-useful coefficient; it specifies a mandatory
loss that such a new theorem must pay.

### 10.2 Sylvester one-hole selectors rule out broad frame-uniform RIP

There is a stronger obstruction if a proposed improvement, like the
existing RIP proof, is required uniformly conditional on EVERY full
Hadamard frame and its prescribed groups. Fix a legitimate full frame
`H_L tensor H_M`, with L a fixed power of two and `M=2^d`. Index
physical coordinates by `(a,g)` with a in F_2^(log L) and g in F_2^d.
The group at g contains its L choices of a; independently omit one
uniform a from every group. Retain fraction `p=1-1/L`.

Fix a_0 and partition the g-space into affine cosets of a fixed
t-dimensional subspace, of size `K=2^t`. A coset is entirely retained
on the slice a=a_0 with probability p^K, and these events are independent
over its M/K cosets. Fix epsilon in (0,1), and choose K to be the largest
power of two at most `(1-epsilon)log(M)/(-log p)`. Then K=Theta(log M)
and

```
P(no such retained coset)<=exp[-(M/K)p^K]
                         <=exp[-M^epsilon/K] ->0.       (17)
```

On an entirely retained affine coset U of size K choose a nontrivial
character of U, extended by zero outside U and normalized in L2. It
has physical sum zero. Its Walsh transform is supported on exactly
`s=m/K=O(m/log m)` columns, on a coset of the annihilator of U, and
does not contain DC. If u is that normalized coefficient vector, then

```
||P H_selected u||^2=m,       ||u||=1,       |supp(u)|=s.
```

Thus the selected centered frame has restricted squared norm EXACTLY m
on s columns. An operator-controlled exact balanced repair differing by
o(sqrt(m)) retains a lower bound `m-o(m)` on that restricted norm.
The obstruction survives conditioning on an event of probability bounded
away from zero: (17)'s exceptional probability remains o(1) after paying
the bounded density factor.

In particular, on this allowed full frame a uniform restricted-column
certificate on ANY fixed positive fraction of columns has to pay
`K>=m-o(m)`. Its direct reciprocal-port coefficient is at least
`1/(2sqrt(p))-o(1)>1/2`, so this frame-uniform RIP extension alone
cannot give the desired strict-subhalf phase.

This uses the actual one-hole selector and actual balanced repair; it is
not a fractional Gram counterexample. Its scope is nevertheless precise:
the full frame was fixed to an allowed Sylvester tensor. This refutes a
uniform-in-full-frame extension of the old restricted-Gram input. It does
NOT prove that the typical independent-child recursive frame law has
these affine obstructions. An estimate using that extra law, rather
than discarding it before the selector argument, remains a distinct
possible route.

No quantitative improvement of a macroscopic actual cap was obtained in
the bounded adjacent pass. Degree-weight contraction alone supplies
neither that missing frame-specific estimate nor the joint packet/port
budget needed to close the dense mixed branch.

## 11. Root-proposed ellipse barrier, independently verified

The root agent found the following stronger discriminating consequence
of the square-Paley family. Its derivation was independently checked in
this mechanism pass. It rules out every norm-only quadratic joint-weight
bootstrap of a natural type, not just the coefficient-one additive bound
or the particular reserve certificate above.

Suppose there is a uniform joint inequality for full sign quadratics,
arbitrary constants, and arbitrary linear fields (in particular the
fields needed by random restrictions),

```
M(c+Q_A+h dot x)^2 >= a ||h||_1^2+b n^3-o(n^3),
a,b>=0.                                               (18)
```

The error must be uniform on that class (or have the corresponding
uniform sequence formulation); a bound valid only for a typical
independently Gaussian field is not being assumed here. Applying (18)
to the square-Paley sequence and the deterministic uniform field in
Section 2 gives necessary constraints

```
0<=a<=1,
b<=1/4                           if 0<=a<=1/2,
b<=a-a^2                        if 1/2<=a<=1.         (19)
```

Indeed t>=1 forces a<=1, while minimizing
`(1+t^2)^2/4-a t^2` over 0<=t<=1 gives the displayed b bounds.
For a>=1/2 the minimum occurs at `t^2=2a-1`. No assertion that this
uniform field is typical under a random restriction is used: it tests
the claimed uniform joint inequality itself.

It suffices to assume uniformity on every bounded natural-scale class
`||h||_1<=C n^(3/2)`, for each fixed C; all the Paley tests use fixed t.
A theorem restricted to a smaller prescribed field-distribution class
need not satisfy these tests and is outside this barrier.

Freeze a complement of size `(1-theta)n` as independent signs and
apply (18) on theta*n free coordinates. The marginal random-walk law
and Jensen give the prospective parent lower certificate

```
M_parent^2/n^3 >= max_(0<=theta<=1)
 [a s theta^2(1-theta)+b theta^3]-o(1),
s=2/pi.                                               (20)
```

Let B=w^2 be an existing universal half-range squared lower constant.
To isolate bootstrapping rather than an independently stronger zero-
field theorem, impose b<=B. Then (19)--(20) cannot improve w whenever

```
B >= B_*=(s/3)(1-s/3)=.16717495361482146...,
w >= sqrt(B_*)=.4088703383895944...,                    (21)
```

in particular at the current w=.4333221116640807.

Here is the complete optimization proof for B in [B_*,1/4]; larger B
is irrelevant to the present comparison. Put

```
a_+=(1+sqrt(1-4B))/2,       a_*=1-s/3.
```

Since B>=B_*, `a_+<=a_*`. For fixed a and b, set K=a s. The maximum
of the polynomial in (20) equals b when K<=3b, and otherwise equals
`4K^3/[27(K-b)^2]`.

* If a<=a_+, then K<=3B and b<=B. Replacing b by B makes the
  polynomial at most B throughout [0,1].
* If a_+<a<=a_*, replace b by its larger permitted value a(1-a).
  This value is at most B, and K<=3a(1-a), so the polynomial's maximum
  is at most a(1-a)<=B.
* If a>a_*, make the same replacement b=a(1-a). The maximum is
  `4a s^3/[27(a+s-1)^2]`. This is strictly decreasing for a>a_*:
  its derivative has sign `s-1-a<0`. Its value at a_* is B_*, hence
  the maximum is at most B_*<=B.

This proves the barrier. It does not rule out field-distribution-sensitive
certificates, nonlinear joint functionals, optimizer-specific hypotheses,
or an ellipse whose zero-field coefficient b independently improves B.
The field replay includes 20,000 numerical checks of this analytic
optimization, preserving the threshold and qualification.
