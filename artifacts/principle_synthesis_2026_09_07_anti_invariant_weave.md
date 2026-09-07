# A paired-row weave inside the skew anti-invariant lift family

2026-09-07. Status: PROVED by the analytic construction and a completed
directed scalar certificate. All twenty point bounds and all ten density
chords passed; the construction and scalar formulas were independently
reconstructed. The certified family upper is below `0.498795641`, hence
strictly below `499/1000`. This is not an improvement to the smaller original
all-signing upper endpoint.

## 1. Prospective original-problem implication and the bounded advance

The initial cross-order target was an actual-sign lift for selected original
near-minimizers, at two independent fixed multipliers, with summable error.
For multiplier two, the skew lift

```math
L(A,C)=\begin{pmatrix}A&C\\-C&-A\end{pmatrix},
\qquad C^{\mathsf T}=-C,
```

has only its matching edges missing and obeys the exact identity

```math
Q(L)=2\max_{u,v\text{ disjoint signed shores}}
\bigl(|u^{\mathsf T}Av|+|u^{\mathsf T}Cv|\bigr).
```

Thus a seed-sensitive joint orientation bound, not cancellation on a cut,
would be needed. The archive proves a spectral-certificate half-floor for
this family but explicitly leaves an actual Boolean half-floor open.

The present construction resolves that structural hard step. Its conclusion is

```math
\limsup_{n\to\infty}\;
\min_{A,C}\frac{Q(L(A,C))}{(2n)^{3/2}}<\frac{499}{1000}<\frac12,
```

where `A` is hollow symmetric full sign and `C` is hollow skew full sign.
It therefore rules out an actual Boolean half-floor for this exact
family. It does NOT prove the seed-sensitive orientation bound, preserve a
given seed's cap, or establish the original limit.

## 2. Fix the tournament first; do not condition independent incidences

Let `m=2q` and choose an almost-regular tournament on the macro vertices:
every outdegree is `q` or `q-1`. For example, remove one vertex from a regular
tournament on `m+1` vertices. Write its incidence signs as `eta_i(j)`, with

```math
\eta_i(j)\eta_j(i)=-1\quad(i\ne j).
```

Choose the diagonal sign `eta_i(i)` so that each entire row has exactly `q`
plus signs and `q` minus signs. No probability or entropy payment is made
for this tournament: it is fixed before choosing the row bases.

In each macro fibre independently take two order-`q` Hadamard bases `H_i^+`
and `H_i^-` from the existing fixed-depth recursive ensemble, and form

```math
H_i=\begin{pmatrix}H_i^+&H_i^-\\H_i^+&-H_i^-\end{pmatrix}.
```

Map the first `q` output columns uniformly to the positions with `eta_i=+1`
and the second `q` columns uniformly to the positions with `eta_i=-1`.
These two output permutations are independent. Multiplying output columns
by additional signs, if desired, does not change the following relation.
The involution `pi_i` exchanging the paired rows satisfies

```math
H_i(\pi_i(a),j)=\eta_i(j)H_i(a,j).
```

Retain `k=2r` rows in each fibre, using the SAME selector of `r` coordinates
in both paired halves. Put `p_m=k/m=r/q`; eventually `p_m -> p=24/25`.
For independently fair symmetric macro signs `S_ij`, define only the
macro-off-diagonal entries

```math
W_{(i,a),(j,b)}=S_{ij}H_i(a,j)H_j(b,i),\qquad i\ne j.
```

All these entries are signs. Applying all paired-row involutions changes
every such entry's sign, because the two incidence characters multiply to
minus one. The absent macro-diagonal blocks are addressed in Section 6.

## 3. The joint off-diagonal defect identity avoids a leading deletion cost

For block spins `x_i`, let `h_i=H_i[T_i,:]^T x_i` and let
`E_off=x^T W_off x`. Row orthogonality gives
`sum_i ||h_i||^2=m^2 k`. For either polarity,

```math
D_\sigma^{\rm off}
 :=\sum_{i\ne j}(h_i(j)-\sigma S_{ij}h_j(i))^2
 =2\left(m^2k-\sum_i h_i(i)^2-\sigma E_{\rm off}\right).
```

Consequently `sigma E_off>=m^2k(1-gamma)` implies

```math
D_\sigma^{\rm off}\le2\gamma m^2k-2\sum_i h_i(i)^2
\le2\gamma m^2k.
```

This is a direct high-energy implication for the partial weave. It is NOT
triangle deletion of the old rank-one macro-diagonal blocks: that triangle
bound would cost order `m k^2`, a leading amount when `k~m`.

Averaging the independent sign of each undirected macro edge in
`exp[-t D_sigma^off/(2k)]` gives the same folded Gaussian kernel `K_t` as in
the archived all-spin weave proof. Each edge appears twice. The row-orbit
measure is now a PRODUCT of the permutation groups on the two prescribed
column classes. In its Gram tensor, the squared row norm is therefore
exactly

```math
P_t(v_i^+)P_t(v_i^-),\qquad v_i^\pm=h_i^\pm/\sqrt k,
```

before removing the one diagonal coordinate. Removing that coordinate from
one of the two groups costs at most `sqrt(2q)` in the row norm, using the
archived permanent deletion bound. Graph Cauchy--Schwarz applies to the
same PSD kernel on every edge, irrespective of which endpoint groups label
that edge. Its product row bound is thus `L_t(v_i^+)L_t(v_i^-)` times a
polynomial factor. There is no independent-incidence conditioning argument.

## 4. Exact paired-spin enumeration and one uniform finite depth

For a retained row pair, its two spins have either a nonzero sum or a
nonzero difference, but not both. If `a` of the `r=k/2` pairs use the sum
channel, write `theta=a/r`. There are exactly

```math
2^r\binom ra
```

such block spin assignments. The two child input vectors have respective
nonzero densities `p_m theta` and `p_m(1-theta)`, and common nonzero
amplitude `sqrt(2/p_m)`. Their second moments are `2theta` and
`2(1-theta)`; the sum is exactly two.

For fixed `p,t`, write

```math
\mu_s=(1-s)\delta_0+\frac s2(\delta_{-\sqrt{2/p}}+
                                  \delta_{\sqrt{2/p}}),
\qquad 0\le s\le p,
```

and let `h(theta)` denote binary entropy with natural logarithms. The
fixed-depth child row theorem and Stirling's formula bound the expected
all-spin row factor, per macro dimension, by

```math
R_{t,p}(\theta)
 =\frac p2\,[\log2+h(\theta)]
  +\frac12\left[E_t(\mu_{p\theta})+
                       E_t(\mu_{p(1-\theta)})\right].       (1)
```

Equivalently the last two terms are
`E_(2t theta)(nu_(p theta))` and
`E_(2t(1-theta))(nu_(p(1-theta)))`, but the fixed-amplitude formulation is
regular at zero and is the one used in the certificate.

One fixed recursive depth must work for ALL `theta`. This follows from the
archived explicit Bellman stopping estimate, not from convexity of any
finite iterate. Both child source moments are at most two. Their initial
budget `Phi-G` is at most `-g_t(2)`. Thus the cutoff `C`, the compact-set
positive drift constant `kappa`, and then a depth can all be chosen
uniformly over this source family. The resulting error in (1) can be made
arbitrarily small before `m` tends to infinity. Fixed-depth alphabet counts
are polynomial and terminal orbital errors are `o(m)`, uniformly over the
terminal Hadamards. The actual `p_m` converges to `p`; the fixed common
amplitude and these compact estimates handle that harmless change.

The independent fibres and the union over all full spins and both polarities
then yield an actual partial signing with asymptotic normalized cap at most

```math
\frac{t+\sup_{0\le\theta\le1}R_{t,p}(\theta)}{2t\sqrt p}.
                                                               (2)
```

As usual, first choose a strict margin in (2), next a sufficiently large
FINITE depth, and only then let the terminal orders tend to infinity.

## 5. Reducing the entire scalar continuum to twenty certified point bounds

Fix `p=24/25`, `t=7/2`, and `a^2=2/p=25/12`. The exact archived ternary
reproduction theorem applies to every `mu_s`, `0<s<1`: an optimal symmetric
reproduction law has support `{0,+az,-az}` for `0<=z<=1`. Eliminating its
mixture weight gives

```math
E_t(\mu_s)=\sup_{0<\lambda\le t,\ 0\le z\le1}
 \left[c_t(\lambda)-\lambda s a^2+
       \max_{0\le w\le1}\{s\log(1+wB)-\log(1+wA)\}\right],
```

where

```math
c_t(\lambda)=\tfrac14\log\frac{\lambda(2t-\lambda)}{t^2},
\quad A=e^{a^2\lambda z^2}-1,
\quad B=\cosh(2a^2\lambda z)-1,
\quad w_* =\operatorname{clip}_{[0,1]}
             \frac{sB-A}{(1-s)AB}.
```

Since every source is centered and bounded by `a`, it is
`a^2`-subgaussian. Entropy duality gives
`I(X;L)>=E(E[X|L]^2)/(2a^2)`. Therefore
`J_lambda=lambda s a^2` for `lambda<=1/(2a^2)=6/25`.
The low-precision contribution is bounded by the explicit Gaussian value
`g_t(s a^2)`. The remaining CLOSED rectangle
`[6/25,7/2] times [0,1]` is covered by the directed-interval verifier
`computations/principle_synthesis_2026_09_07_anti_weave_interval.py`.

For every density `s_j=p j/20`, `j=1,...,20`, the certified rational bound is

```text
j:       1    2    3    4    5    6    7    8    9   10
1000 e: -224 -378 -465 -529 -581 -623 -660 -691 -720 -745
j:      11   12   13   14   15   16   17   18   19   20
1000 e: -768 -789 -808 -826 -843 -859 -874 -888 -901 -832
```

The actual upper target is `e_j+1/2000`; at zero the exact value is zero.
The targets were selected by a floating diagnostic and then verified by
complete outward-rounded rectangle covers. The twenty covers checked
100,466 boxes and accepted 50,243 leaves, with no numerical pruning.

The continuum interpolation has a separate analytic justification.
Rate-distortion duality writes `J_lambda(mu)` as an infimum over
reproduction laws of a LINEAR functional of the source `mu`. Thus `J_lambda`
is concave in the source, and `E_t=sup_lambda(c_t-J_lambda)` is convex.
Because `mu_s` is affine in `s` at the FIXED amplitude `a`, each interval
between adjacent density points is bounded above by its endpoint chord.
No convexity under simultaneous density/normalizing-amplitude changes is
asserted.

On each `theta` interval `[j/20,(j+1)/20]`, the two upper chords in (1)
give an affine function `b theta+d`. It remains to bound

```math
\frac p2 h(\theta)+b\theta+d.
```

The maximum is at an endpoint if the derivative has one sign throughout;
otherwise its unrestricted maximum is
`(p/2) log(1+exp(2b/p))+d`. These are directed interval computations at
rational inputs. Symmetry covers `theta>=1/2`. The verifier checks (2)
strictly below `499/1000` on all ten half-intervals.

The canonical complete rational result is
`computations/results/principle_synthesis_2026_09_07_anti_weave_interval.json`.
Its largest displayed chord cap is approximately `0.4987956402157831`;
the exact rational upper endpoint in that file is the certificate. The
independent entropy-tangent implementation gives an upper below
`0.498795642` without relying on the critical-point branch formula.
Replay and exploratory failure records are in
`principle_synthesis_2026_09_07_anti_weave_verification.md`.

The generalized envelope derivatives, clipped-weight intervals,
low-precision branch, density convexity, and product-orbit construction
were independently checked in
`principle_construct_2026_09_07_skew_lift_and_paired_weave_audit.md`.

## 6. Completing every nonmatching edge and filling all orders

Each absent macro-diagonal fibre has even order `k`. Fill it with an
independent small skew anti-invariant partial signing, relative to the same
paired-row involution. Such a filler with cap `O(k^(3/2))` exists by an
elementary random-sign union bound: in its order-`k=2r` lift, each test has
independent coefficients in `{0,+2,-2}` and total squared coefficient sum
at most `k^2/2`. A threshold `k sqrt((k+2)log2)` makes the union failure
probability strictly less than one. This leaves only the `k/2` paired
matching edges absent in that fibre.

The sum of filler caps is `O(m k^(3/2))=O(N^(5/4))=o(N^(3/2))` when
`k~m`, where `N=mk`. The resulting partial signing has EVERY nonmatching
edge filled and is negated by the global involution. Ordering one vertex
from every pair first gives exactly `L(A,C)` with the stipulated symmetric
`A` and skew `C`. Filling its remaining `N/2` matching edges, if desired,
costs only `N/2` in cap and produces an actual complete signing.

Use the archived terminal orders generated by `H_2,H_12`; their ratios are
dense, and multiplying by the fixed recursion-depth factor preserves that
property. Hence available values `N=m*2 floor(pm/2)` are ratio-dense.
Restricting to any prescribed number of vertex PAIRS preserves the exact
skew anti-invariant form. Principal cap monotonicity follows by averaging
the omitted independent Boolean coordinates, and the ratio loss tends to
one. Thus the asymptotic family bound is all-order in its half-order `n`,
not merely a subsequence construction.

## 7. Scope retained

This is an existence theorem for a newly constructed structured family.
It is not a theorem for every skew completion of a given seed, nor for
every original exact minimizer. The independent column signs still erase
the macro seed in this ensemble. Accordingly a successful subhalf
certificate removes one prospective family obstruction but does not yet
provide seed landing, a multiplier-two recurrence for `M_n`, multiplier
three, width-to-cap equivalence, or convergence of the original sequence.
