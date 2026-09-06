# Wave 23 memo: minimax boosting for selector codebooks

## Status

The minimax identities, greedy theorem, fixed-slice obstruction, and finite
`A_9` certificates below are **proved**.  The finite certificates are checked
with exact rational arithmetic by `tmp/minimax_boosting_r23.py`.  The result
is a positive abstract codebook theorem, but neither fixed-slice
concentration nor presently known exact-minimizer structure supplies its
weak-learner premise at the required scale.

## 1. The exact weak-learner game

Fix an order-`n` signing `A`, let `X=Omega_{n,m}` be the selector slice, and
let `H=D_n` be the full oriented projective cuts.  Write

```math
q_S=Q(A[S]),\qquad
\ell(S,d)=q_S-c_A(S,d).
```

In the ordered normalization used in the ledger,

```math
0\le \ell(S,d)\le L_m:=2m(m-1).
```

For a tolerance `t`, put

```math
G_t(d)=\{S:\ell(S,d)\le t\}.
```

Define the finite zero-sum game value

```math
\boxed{
\alpha_t
=\min_{w\in\Delta(X)}\max_{d\in H}w(G_t(d))
=\max_{\nu\in\Delta(H)}\min_{S\in X}
\nu\{d:\ell(S,d)\le t\}.}
\tag{B.1}
```

The equality is von Neumann minimax applied to the zero-one incidence
matrix.  Therefore the exact dual weak-learner condition `alpha_t>=alpha`
has two equivalent forms:

1. for every probability weighting `w` of the selectors, some full cut is
   `t`-good on at least `alpha` of the weighted mass;
2. there is one randomized common prior on full cuts which is `t`-good for
   every selector with probability at least `alpha`.

This is a genuine duality, not an exchange of `max` and expectation by
assertion.  Choosing one exact child ground for each selector and averaging
uniformly over all full cuts gives only the universal baseline
`alpha_0>=2^{-m}`.  It leads to a linear logarithmic codebook size and uses
no minimality of `A`.

## 2. Greedy boosting and the exact cardinality bound

Assume `alpha_t>=alpha>0`.  Start with all selectors uncovered.  At a
nonempty residual set `R`, apply the weak-learner condition to the uniform
law on `R` and choose a cut covering at least an `alpha` fraction of `R`.
After `K` choices,

```math
U_m(R_K)\le(1-\alpha)^K.
```

Consequently, for every `epsilon in (0,1)`, there is a codebook `C` with

```math
\boxed{
|C|\le
\max\left\{1,
\left\lceil\frac{\log(1/\epsilon)}{-\log(1-\alpha)}\right\rceil
\right\},
\qquad
\mathbb E_{U_m}\min_{d\in C}\ell(S,d)
\le t+L_m\epsilon.}
\tag{B.2}
```

For `alpha=1`, one word covers every selector and the maximum with one handles
the extended-value denominator.  In particular, if `alpha>=e^{-R}`, then

```math
\boxed{
\log|C|\le R+\log(1+\log(1/\epsilon)).}
\tag{B.3}
```

Indeed, `-log(1-alpha)>=alpha`, so the ceiling in (B.2) is at most
`1+e^R log(1/epsilon)`.

For a requested average deficit `Delta<L_m`, take

```math
t\le\frac\Delta2,
\qquad
\epsilon=\frac{\Delta}{2L_m}.
```

Then (B.2) gives average deficit at most `Delta`.  Thus (10.714) would follow
from the concrete weak-learning statement

```math
\boxed{
\alpha_{O(n^{3/2-c})}
\ge\exp\{-O(n^{1/2-2c})\}.}
\tag{B.4}
```

The extra term in (B.3) is only `O(log log n)`.

### The average-only capped version

Condition (B.1) asks to cover even a selector weighting concentrated on one
bad tail.  For the average target one can stop at residual mass `epsilon`.
Let

```math
W_\epsilon=
\left\{w:X\to[0,1/\epsilon]:\mathbb E_{U_m}w=1\right\}
```

and define

```math
\boxed{
\alpha_{t,\epsilon}^{cap}
=\min_{w\in W_\epsilon}\max_d
\mathbb E_{U_m}[w(S)1_{G_t(d)}]
=\max_{\nu\in\Delta(H)}\min_{w\in W_\epsilon}
\mathbb E_{U_m}[w(S)p_\nu(S)],}
\tag{B.5}
```

where `p_nu(S)=nu{d:ell(S,d)<=t}`.  This is again finite minimax.  The inner
minimum on the right is exactly the mean of the lowest `epsilon`-fraction of
the values `p_nu(S)`, with fractional mass at the boundary when necessary.
Thus the dual asks for a common prior with a controlled **lower coverage
tail**, rather than pointwise coverage.

Whenever the uncovered set has mass at least `epsilon`, its conditional
density is at most `1/epsilon`; hence the same greedy proof uses
`alpha_{t,epsilon}^{cap}` and stops with residual mass below `epsilon`.
Equations (B.2)--(B.3) remain valid.  This capped formulation is the precise
average-tail weak learner relevant to (10.712).

It is also necessary at the exponential scale.  If a `K`-word codebook has
`t`-uncovered uniform mass at most `eta epsilon`, where `0<=eta<1`, then for
every `w in W_epsilon` that uncovered set has weighted mass at most `eta`.
The union bound in the reverse direction gives

```math
\sum_{d\in C}\mathbb E_{U_m}[w1_{G_t(d)}]\ge1-\eta,
```

and therefore

```math
\boxed{
\alpha_{t,\epsilon}^{cap}\ge\frac{1-\eta}{K}.}
\tag{B.5a}
```

The same statement holds for any restricted hypothesis class containing
the codebook.  Thus capped weak learning and average threshold codebooks are
equivalent up to the logarithmic greedy overhead and a constant slack in the
exceptional mass; (B.5) is not merely an artifact of the proof algorithm.

### Restricting the learner to row-regular parent grounds

The support-dependent fixed-slice refinement permits a substantially larger
codebook if every selected word is a regular exact parent ground.  Define

```math
H_B^{gr}(A)=
\left\{d\in\mathcal D_n:
\langle A,d\rangle=q_n,
\ R_\infty(d):=\max_i r_i(d)\le B
\right\}.
```

Every parent-ground row field is nonnegative, since flipping coordinate `i`
changes its payoff by `-4r_i`.  Thus this `R_infty` is also the absolute row
maximum needed in the slice mgf.  Replace `H` by `H_B^{gr}(A)` in (B.1) and
(B.5), and denote the resulting values by

```math
\alpha_t^{gr}(B),
\qquad
\alpha_{t,\epsilon}^{gr,cap}(B).
\tag{B.5b}
```

The value is defined to be zero when this restricted ground class is empty.

Finite minimax and the greedy proof are unchanged.  In particular, a clean
exact sufficient condition for the support-dependent restriction theorem is

```math
\boxed{
\begin{aligned}
B&=O(n^{3/4-c}),\\
t&=O(n^{3/2-c}),\\
\epsilon&\asymp n^{-1/2-c},\\
\alpha_{t,\epsilon}^{gr,cap}(B)
&\ge\exp\{-O(n^{3/4-c})\}.
\end{aligned}}
\tag{B.5c}
```

Indeed, boosting gives a parent-ground codebook with average deficit
`O(n^{3/2-c})` and

```math
\log|C|=O(n^{3/4-c}).
```

For such grounds, `sum_i r_i=q_n` and
`R_2<=Bq_n=O(n^{9/4-c})`.  Taking `lambda=eta n^{-3/4}` in the
support-dependent entropy inequality gives the error terms

```math
n^{3/4}\log|C|,
\quad
n^{-3/4}R_2,
\quad
n^{-3/4}n^2,
```

which are respectively `O(n^{3/2-c})`, `O(n^{3/2-c})`, and
`O(n^{5/4})`.  The operator denominator is uniformly positive for small
fixed `eta`, and `lambda B=O(n^{-c})`.  Thus the exponents in the
support-dependent refinement are consistent, including the restriction
`c<1/4` needed to absorb `n^{5/4}`.

This formulation also isolates two separate missing facts.  There need not
be any known theorem making `H_B^{gr}(A)` nonempty, and even nonemptiness
does not give a weak learner for an adversarial selector weighting.  The
identity `sum_i r_i=q_n` controls the number of heavy rows within one ground,
not its maximum, and it gives no way to change a heavy ground while
preserving its selector coverage.

## 3. Why expected-loss weak learning is circular

There is a tempting real-valued alternative: require for every selector law
`w` a full cut with small weighted mean deficit.  Its exact value exposes the
problem.  Put

```math
r_{ij}(w)=\Pr_{S\sim w}\{i,j\in S\},
\qquad R_w=(r_{ij}(w)),
```

and extend `Q` homogeneously to real zero-diagonal matrices.  Linearity of a
fixed cut payoff gives

```math
\boxed{
\min_d\mathbb E_w\ell(S,d)
=\mathbb E_wQ(A[S])-Q(A\circ R_w).}
\tag{B.6}
```

For `w=U_m`, every off-diagonal inclusion probability is
`p_2=(m)_2/(n)_2`, so

```math
\boxed{
\min_d\mathbb E_{U_m}\ell(S,d)
=\mathbb E_{U_m}Q(A[S])-p_2q_n.}
\tag{B.7}
```

The right side is exactly the optimized-restriction excess that the route is
trying to prove power-saving.  Therefore a mean-deficit weak learner with
the desired error already assumes the conclusion at its initial uniform
weighting.  Point-mass selector laws have zero gap in (B.6), but that gives
no control of mixtures.  Global minimality of `A` among complete signings
does not bound `Q(A circ R_w)`, which is the norm of a weighted matrix.

## 4. Fixed-slice concentration gives a converse, not a learner

Retain the notation of (10.691)--(10.692):

```math
p=m/n,\qquad p_2=(m)_2/(n)_2,\qquad
\varepsilon_{sl}=p^2-p_2,
```

the conditioning defect `chi_{n,m}`, and the uniform variance proxy
`Vbar_{n,m}`.  A safe common mgf domain is

```math
0<\lambda\le\lambda_0
:=\min\left\{
\frac3{8p(n-1)},
\frac1{2\sqrt2\lVert A\rVert_{op}}
\right\}.
```

If `ell(S,d)<=t`, then

```math
c_A(S,d)\ge q_m-t,
\qquad
p_2\langle A,d\rangle\le p_2q_n.
```

Chernoff's inequality and (10.692) therefore give, uniformly in `d`,

```math
\boxed{
U_m(G_t(d))\le e^{-\Psi_{n,m}(t)},}
\tag{B.8}
```

where

```math
\boxed{
\Psi_{n,m}(t)=
\left[
\sup_{0<\lambda\le\lambda_0}
\left\{
\lambda\big[q_m-t-(p_2+\varepsilon_{sl})q_n\big]
-\lambda^2\overline V_{n,m}
\right\}
-\chi_{n,m}
\right]_+.}
\tag{B.9}
```

Taking `w=U_m` in either weak-learner game shows

```math
\boxed{
\alpha_t,\ \alpha_{t,\epsilon}^{cap}
\le e^{-\Psi_{n,m}(t)}.}
\tag{B.10}
```

There is also a direct codebook converse.  If a codebook has average deficit
at most `Delta`, then for every `T>Delta`, Markov and the union bound give

```math
1-\frac\Delta T
\le U_m\{\min_{d\in C}\ell(S,d)\le T\}
\le |C|e^{-\Psi_{n,m}(T)}.
```

Hence

```math
\boxed{
\log|C|
\ge\Psi_{n,m}(T)+\log(1-\Delta/T),
\qquad
\log|C|\ge\Psi_{n,m}(2\Delta)-\log2.}
\tag{B.11}
```

This proves precisely why the Wave 21 mgf cannot itself construct the Wave
22 codebook: it controls the size of each fixed cut's good-selector set from
above.

The same audit matches the new regular-ground exponent.  If
`B=O(n^{3/4-c})`, then on `H_B^{gr}(A)` the variance proxy improves to

```math
O(Bq_n+n^2)=O(n^{9/4-c}+n^2).
```

Under a leading gap `q_m-p_2q_n>=g n^{3/2}`, take
`lambda=theta n^{-3/4}` with small fixed `theta`.  The favorable Chernoff
term is `Theta(n^{3/4})`, whereas the row-variance terms are
`O(n^{3/4-c})+O(n^{1/2})`.  Consequently

```math
\boxed{
\max_{d\in H_B^{gr}(A)}U_m(G_t(d))
\le\exp\{-\Omega(n^{3/4})\}}
\tag{B.11a}
```

for `t=o(n^{3/2})`.  Hence a regular-ground codebook would require
`log|C|=Omega(n^{3/4})`, just beyond the allowed
`O(n^{3/4-c})`.  The relaxed exponent is therefore correctly matched to the
support-dependent mgf; it still cannot cross a pre-existing leading
restriction gap.

At a fixed density `m/n -> rho`, suppose along some sequence

```math
q_m-p_2q_n\ge g n^{3/2}
```

for a constant `g>0`.  If `t=o(n^{3/2})`, choose
`lambda=theta/n` with a sufficiently small fixed `theta>0` in (B.9).
Since `Vbar=O_rho(n^{5/2})`, `epsilon_sl q_n=O(n^{1/2})`, and
`chi=O(log n)`, this gives

```math
\boxed{\Psi_{n,m}(t)=\Omega_\rho(\sqrt n).}
\tag{B.12}
```

Thus every codebook with subleading average deficit has
`log|C|=Omega(sqrt n)` and cannot satisfy (10.714), whose logarithm is
`o(sqrt n)`.  Equivalently, the weak-learner exponent `R` in (B.4) must be
at least order `sqrt n`.

There is a rigorous density range where this obstruction is unconditional.
Using the current bounds

```math
q_m\ge(a-o(1))m^{3/2},\quad
q_n\le(1+o(1))n^{3/2},\quad
a=0.672986728862\ldots,
```

the leading gap is positive whenever

```math
0<\rho<a^2=0.452911137224\ldots.
```

Therefore the codebook target is impossible at every such fixed density.
This does not touch the strategically relevant range `rho>=1/2`, where the
known global constants do not force a positive gap.  In that range (B.8)
remains a diagnostic: proving the large weak coverage in (B.4) must exploit
new overlap of minimizer sections, and cannot be read off from fixed-slice
concentration.

## 5. Exact `A_9`, `m=8` audit

The finite game is nontrivial but compressible.  In the ordered ledger
normalization, exhaustive enumeration of all 512 full oriented cuts gives:

| tolerance `t` | game value `alpha_t` | minimum deterministic cover |
|---:|---:|---:|
| `0` | `4/13` | `4` cuts |
| `4` | `5/8` | `2` cuts |
| `8` | `1` | `1` cut |

At `t=0`, the dual selector weighting is

```math
w=(1,2,1,1,2,2,2,1,1)/13.
```

Every realizable exact-ground incidence pattern has weight at most `4/13`.
An exact primal prior, listed and checked in the verifier, gives every
selector coverage at least `4/13`, proving equality.  At `t=4`, a rational
primal and the dual weighting

```math
(1,1,1,0,2,1,1,1,0)/8
```

prove `alpha_4=5/8`.  The integral cover sizes are certified by explicit
covers and exhaustive exclusion of smaller ones.

These constants explain the good one-deletion behavior already seen in
(10.715)--(10.717), but they do not yield a fixed-density asymptotic lower
bound on `alpha_t`.

All 25 oriented parent grounds of `A_9` have maximum row field `4`, `6`, or
`8`; their multiplicities are respectively `2`, `16`, and `7`.  Restricting
the hypothesis class gives the exact tradeoff

| row cap `B` | `alpha_0^{gr}` | `alpha_4^{gr}` | `alpha_8^{gr}` |
|---:|---:|---:|---:|
| `4` | `0` | `0` | `1` |
| `6` | `2/7` | `1/2` | `1` |
| `8` | `4/13` | `5/8` | `1` |

At cap four, the two very regular grounds miss some selectors entirely at
tolerances zero and four, so a point-mass selector adversary makes the game
value zero.  Cap six restores a positive learner but weakens its exact
coverage.  Thus the greedy theorem does not manufacture regularity: a new
structural result must jointly guarantee row-regular parent grounds and
their coverage of every nonexceptional selector weighting.  These displayed
values are for the uncapped game.  In particular, at cap four and tolerance
four only one of the nine selectors is missed, so a capped average-tail game
may discard it when its allowed exceptional mass is at least `1/9`.  The
table is a finite tradeoff diagnostic, not an asymptotic obstruction to
(B.5c).

## Updated route conclusion

The exact missing combinatorial lemma can be stated as the capped
weak-learner bound (B.4)--(B.5): every selector weighting that has not yet
fallen into the allowed exceptional tail must admit one cut covering
`exp{-O(n^{1/2-2c})}` of its mass at power-saving deficit.  Minimax turns
this into a common-prior lower-tail overlap statement, and greedy boosting
then constructs the desired codebook with no hidden linear loss.

Neither available ingredient supplies that premise.  Expected-loss weak
learning is circular by (B.7); fixed-slice mgf is a converse by
(B.8)--(B.11); and exact child extendability gives only `2^{-m}` coverage.
The support-dependent variant (B.5c) relaxes the required exponent from
`n^{1/2-2c}` to `n^{3/4-c}`, but introduces a genuine row-regularity burden;
the `A_9` cap audit shows that regularity and coverage cannot simply be
separated.  A positive continuation must prove minimizer-specific overlap
across a large family of principal sections, optionally inside a regular
parent-ground class.  The theorem sharpens, but does not replace, the
parent-overlap target (10.727).
