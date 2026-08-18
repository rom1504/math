# Entropy and unoriented overlap do not decide the actual-child row resource

Status: **rigorous optimizer-specific entropy bound, exact actual-minimizer
no-go, and finite-amplitude tangent theorem**.  The results concern the
actual child Gibbs priors and the canonical hybrid path of ES.21--ES.29.
No conference or generic surrogate law is used.

The positive result is that one joint, sector-oriented four-spin overlap
tensor gives the exact leading integrated row-total-correlation response
when the bridge channel is turned on.  The negative result has two
deliberately separated strengths.  Child pressure/entropy profiles fail even
for an orientation-invariant response tangent.  The complete collection of
separate-child spin/unoriented replica-overlap laws additionally fails to
label the two fixed orientation sectors.  The latter witness swaps the two
sectors and therefore is **not** an obstruction to a criterion which takes
their minimum or unordered pair.

This is a finite exact obstruction to an entropy-only or unoriented-overlap
decision of ES.42.  It is not an asymptotic lower bound and does not decide
which ES.42 branch a sequence of large optimizing children occupies.

## 1. What pressure minimality says about child entropy

For a hollow signing `A` of order `d`, put `E_d=binom(d,2)` and

```math
 F_A(t)=\log E_{x,\tau}\exp\{t\tau H_A(x)\},
 \qquad
 {d\nu_{A,t}\over dU_{d+1}}(x,\tau)
 =\exp\{t\tau H_A(x)-F_A(t)\}.                     \tag{EO.1}
```

Here both `x` and `tau` are fair under `U_(d+1)`.  Thus

```math
D(\nu_{A,t}\Vert U_{d+1})=tF_A'(t)-F_A(t).         \tag{EO.2}
```

**Proposition EO.1 (optimizer entropy ceiling).**  If `A` minimizes
`F_A(t)` over all order-`d` signings and `t>0`, then

```math
\boxed{
0\le D(\nu_{A,t}\Vert U_{d+1})
\le E_d\,t\tanh t
\le E_d t^2.}                                      \tag{EO.3}
```

The minimizing pressure itself satisfies

```math
0\le F_A(t)\le E_d\log\cosh t.                    \tag{EO.4}
```

At `t=beta/sqrt(N)`, (EO.3) is only an order-`N` entropy bound when
`d` is proportional to `N`; it supplies no entropy-density saving.

*Proof.*  The identity (EO.2) is immediate from the Gibbs density.  The
one-edge consequence AC.33 of exact sign-flip minimality is

```math
a_eE_{\nu_{A,t}}[\tau x_ux_v]\le\tanh t.
```

Summing over edges gives `F_A'(t)<=E_d tanh t`.  Jensen gives `F_A(t)>=0`,
so (EO.2) proves (EO.3).  Finally, averaging the unnormalized partition
over a uniformly random signing makes its edge factors independent and
gives exactly `(cosh t)^(E_d)`.  A minimizer is no larger than that average,
which proves the upper bound in (EO.4). `square`

The ceiling is optimizer-specific, but it is too large to control either
weighted path mass in ES.29.  More importantly, the next sections show an
information mismatch: scalar entropy forgets a relative sector orientation
which the row path can observe.

## 2. The oriented overlap tensor controls the row-correlation tangent

Fix two finite child Gibbs measures at internal temperature `t`, one
orientation `epsilon`, and the exact zero-bridge joint sector law

```math
\nu_\epsilon(s,x,y)
=\pi_s^{(\epsilon)}\mu_{A,s}(x)\mu_{D,\epsilon s}(y). \tag{EO.5}
```

For a separate bridge amplitude `u`, let `rho=tanh u` and let `p_u` be the
forward bridge likelihood relative to the fair bridge:

```math
p_u(B)=E_{\nu_\epsilon}
       \prod_{i,j}(1+\rho B_{ij}sX_iY_j).            \tag{EO.6}
```

Let `r_u`, `q_(s,u)`, and `h_u` be the canonical row product, hybrid path,
and interaction from IC.1--IC.4.  Define the two weighted resources

```math
\begin{aligned}
\mathsf T_u
 &=\lambda\int_0^\lambda
   {\operatorname {TC}(q_{s,u})\over s^2}\,ds,\\
\mathsf M_u
 &=\lambda\int_0^\lambda
   {\sum_iD((q_{s,u})_i\Vert(r_u)_i)\over s^2}\,ds.
                                                               \tag{EO.7}
\end{aligned}
```

Thus ES.28 says `J_u=mathsf T_u+mathsf M_u`.  Introduce the genuinely joint,
sector-oriented tensor

```math
\Gamma_{ik;j\ell}^{(\epsilon)}
=E_{\nu_\epsilon}[X_iX_kY_jY_\ell],
\qquad
K_\epsilon=\sum_{i<k}\sum_{j,\ell}
 (\Gamma_{ik;j\ell}^{(\epsilon)})^2.                \tag{EO.8}
```

**Theorem EO.2 (oriented-overlap row-correlation response).**  At fixed
finite children, internal `t`, and `lambda>0`, as `u` tends to zero,

```math
\boxed{
\mathsf T_u={\lambda^2\over2}u^4K_\epsilon+O(u^6),
\qquad
\mathsf M_u=O(u^8),
\qquad
\mathcal J_u={\lambda^2\over2}u^4K_\epsilon+O(u^6).} \tag{EO.9}
```

All constants are finite-system constants.  In particular, the tensor in
(EO.8), not child entropy or an unlabelled pair-overlap norm, is the exact
first statistic seen by the row-dependence component of the actual bridge
escort.

*Proof.*  Central symmetry of `Q_(ij)=sX_iY_j` gives the uniform finite-cube
expansion

```math
\log p_u(B)
=u^2\sum_{e<f}E[Q_eQ_f]B_eB_f+O(u^4),              \tag{EO.10}
```

up to an irrelevant constant.  Removing the exact one-row marginals leaves

```math
h_u(B)=u^2H_2(B)+O(u^4),
\qquad
H_2(B)=\sum_{i<k}\sum_{j,\ell}
 \Gamma_{ik;j\ell}^{(\epsilon)}B_{ij}B_{k\ell}.    \tag{EO.11}
```

The Walsh monomials in (EO.11) are orthonormal, so
`E_U H_2^2=K_epsilon`.  Also `r_u=U+O(u^2)`.  From

```math
{dq_{s,u}\over dr_u}\propto e^{-s h_u}
```

one obtains, uniformly for `0<=s<=lambda`,

```math
{dq_{s,u}\over dU}
=1+\text{one-row scores}-su^2H_2+O(u^4).            \tag{EO.12}
```

The product of the row marginals removes the one-row scores.  The quadratic
expansion of KL therefore gives

```math
\operatorname {TC}(q_{s,u})
={s^2u^4\over2}K_\epsilon+O(s^2u^6).               \tag{EO.13}
```

The leading cross-row score has zero projection onto every single row.
Consequently `(q_(s,u))_i/(r_u)_i=1+O(su^4)`, and another quadratic KL
expansion gives

```math
\sum_iD((q_{s,u})_i\Vert(r_u)_i)=O(s^2u^8).        \tag{EO.14}
```

Divide by `s^2`, integrate on the compact interval, and use ES.28. `square`

Theorem EO.2 is a strict finite-response implication: `K_epsilon` is a
child-only four-spin table and not a parent bridge landscape.  It is not
uniform in child order or sufficient at the physical large-system
amplitude without a separate cumulant theorem.

## 3. Exact optimized triangle pair: a fixed-sector obstruction

Take both child orders equal to three and let `A` be the all-positive
triangle signing.  Every order-three signing is an exact minimizer of the
augmented pressure at every `t`, because the high-temperature expansion
gives

```math
\boxed{F_A(t)=3\log\cosh t}                         \tag{EO.15}
```

independently of the sign product around the triangle.  In particular,
both `A` and `-A` are exact optimizing children.  Their augmented Gibbs
laws are carried into one another by the bijection

```math
(x,\tau)\longmapsto(x,-\tau).                       \tag{EO.16}
```

It follows that they have exactly the same:

- pressure profile at every temperature and the same entropy
  `3(t tanh t-log cosh t)` relative to the fair law;
- spin marginal, since it is proportional to `cosh(tH_A(x))`;
- joint law of every finite array of spin-replica overlaps; and
- every replica statistic invariant under reversal of all replica sector
  signs, including the usual even `tau`-overlap hierarchy.

Now compare the two actual minimizing pairs `(A,A)` and `(A,-A)`.  Equivalently,
keep the displayed signings fixed and use the two relative orientations
`epsilon=+1` and `epsilon=-1`.  Put `q=tanh t`.  In sector `s`, the ordinary
triangle partition and any edge correlation are

```math
Z_s=(\cosh t)^3(1+s q^3),
\qquad
c_s={s q+q^2\over1+s q^3}.                          \tag{EO.17}
```

For `epsilon=-1`, the shared-sector weights are `1/2`; for `epsilon=+1`,
they are proportional to `(1+s q^3)^2`.  If

```math
a_\epsilon=\sum_s\pi_s^{(\epsilon)}c_s,
\qquad
b_\epsilon=\sum_s\pi_s^{(\epsilon)}c_sc_{\epsilon s},
```

then direct simplification gives

```math
\begin{array}{c|cc}
 &a_\epsilon&b_\epsilon\\ \hline
+&{q^2\over q^4-q^2+1}&{q^2\over q^4-q^2+1}\\[4pt]
-&{q^2\over q^4+q^2+1}&-{q^2\over q^4+q^2+1}.
\end{array}                                         \tag{EO.18}
```

There are three row pairs, three diagonal column pairs, and six ordered
off-diagonal column pairs.  Hence the exact oriented-overlap norms are

```math
\boxed{
K_+={27q^4\over(q^4-q^2+1)^2},
\qquad
K_-={27q^4\over(q^4+q^2+1)^2},}                    \tag{EO.19}
```

and

```math
\boxed{
K_+-K_-
={108q^6(q^4+1)\over
 (q^4-q^2+1)^2(q^4+q^2+1)^2}>0\quad(0<q<1).}        \tag{EO.20}
```

Thus the oriented tensor distinguishes the labels of two exact optimizing
sectors which all the separate-child entropy and unoriented-overlap data
above identify.  Replacing the right child by its negative swaps `epsilon=+`
and `epsilon=-`.  Consequently the unordered pair `{K_+,K_-}`, and in
particular its minimum, is unchanged.  This witness is a no-go for a
fixed-sector criterion which forgets relative `tau` alignment, not for an
orientation-combined criterion.

## 4. The distinction persists at the physical bridge amplitude

The preceding distinction is not only a derivative obtained by varying an
unphysical bridge parameter.  Set the bridge amplitude equal to the child
temperature and take `lambda=1`.  The channel parameter is then the same
`q=tanh t` as in (EO.17).  Let `J_+(q),J_-(q)` denote the exact canonical
errors.

**Theorem EO.3 (physical actual-minimizer orientation separation).**  As
`q` tends to zero,

```math
\boxed{
\begin{aligned}
J_+(q)&=18q^8+27q^{10}-300q^{12}+O(q^{14}),\\
J_-(q)&=18q^8-27q^{10}-36q^{12}+O(q^{14}),\\
J_+(q)-J_-(q)&=54q^{10}+O(q^{12}).                  \tag{EO.21}
\end{aligned}}
```

Moreover, for the ES.28 split,

```math
\mathsf M_+(q)+\mathsf M_-(q)=O(q^{12}),            \tag{EO.22}
```

and consequently

```math
\boxed{
\mathsf T_+(q)-\mathsf T_-(q)=54q^{10}+O(q^{12})>0} \tag{EO.23}
```

for all sufficiently small positive `q`.

*Proof.*  Everything is a finite sum over `2^9` bridges.  A convenient
exact expansion avoids exponentials.  Relative to the fair bridge, put

```math
N_\epsilon(B;q)=2^{-7}\sum_{s,x,y}
 \prod_{e\in K_3}(1+qsx_e)
 \prod_{f\in K_3}(1+qs\epsilon y_f)
 \prod_{i,j}(1+qB_{ij}s x_i y_j).                  \tag{EO.24}
```

The common child normalization cancels from both inverse escorts.  If
`n_epsilon(b;q)` is the analogous one-row numerator, the exact bridge
probabilities are

```math
R_\epsilon(B)
={\prod_i n_\epsilon(B_i;q)^{-1}
  \over[\sum_b n_\epsilon(b;q)^{-1}]^3},
\qquad
Q_\epsilon(B)
={N_\epsilon(B;q)^{-1}
  \over\sum_C N_\epsilon(C;q)^{-1}}.               \tag{EO.25}
```

The eight row words split into six and two, with numerator polynomials

```math
\begin{array}{c|cc}
 &\text{six rows}&\text{two rows}\\ \hline
+&1-q^4&1+3q^4+4q^6\\
-&1-q^4&1+3q^4-4q^6.
\end{array}                                         \tag{EO.26}
```

Substitution of the polynomial (EO.24) in the finite identity

```math
J_\epsilon(q)=\sum_B R_\epsilon(B)
 \log{R_\epsilon(B)\over Q_\epsilon(B)}            \tag{EO.27}
```

and ordinary Taylor expansion gives (EO.21).  This calculation uses only
integer polynomial multiplication; in particular the coefficients in
(EO.21) are exact, not numerical fits.

For completeness, (EO.24) is `1+O(q^4)` uniformly in `B`, and the same is
true of every one-row numerator.  Hence `h_epsilon=q^4h_(4,epsilon)+O(q^6)`.
The order-four term has zero one-row projection because the `n_epsilon`
are the exact row marginals.  It follows uniformly for `0<=s<=1` that

```math
{(q_{s,\epsilon})_i\over(r_\epsilon)_i}
=1+O(sq^6).
```

The quadratic expansion of KL gives
`sum_iD((q_(s,epsilon))_i||(r_epsilon)_i)=O(s^2q^12)`.
Integration proves (EO.22).  Subtract it from ES.28 and use (EO.21) to
obtain (EO.23). `square`

## 5. An orientation-invariant entropy no-go

The loss of sector labels in the triangle example is not the only failure
of scalar child entropy.  There is also an exact orientation-invariant
witness, at the response-tangent level.

Let `A_0,A_1` be the two order-eight signings displayed in FC.20--FC.21 of
[`actual_child_flip_averaging_ceiling.md`](actual_child_flip_averaging_ceiling.md).
Their common projective absolute-energy histogram is

```math
\begin{array}{c|rrrrrr}
|H|&0&2&4&6&8&10\\ \hline
\#&12&32&32&24&20&8.
\end{array}                                         \tag{EO.28}
```

The certified exhaustive classification proves that these are the only two
cap-`10` classes and that every other order-eight signing has cap at least
`12`.  The elementary comparison FC.22 therefore proves that both are exact
thermal-pressure minimizers for every `t>=3`.  Since (EO.28) is identical,
they have the same complete pressure and entropy profiles.

Pair either child with the unique order-two minimizer `D`, and define the
oriented tensor (EO.8).  Put `q=tanh t` and, under the augmented Gibbs law of
the order-eight child, define

```math
a_{ik}(t)=E[X_iX_k],
\qquad
b_{ik}(t)=qE[\tau X_iX_k].                           \tag{EO.29}
```

**Theorem EO.4 (entropy-equivalent optimizers have different unoriented
row-TC tangents).**  For this `8+2` child pair,

```math
\boxed{
K_\epsilon(A_r,D;t)
=2\sum_{i<k}\{a_{ik}(t)^2+b_{ik}(t)^2\},
\qquad r\in\{0,1\},}                                \tag{EO.30}
```

and the right side is independent of `epsilon`.  At zero temperature its
two exact values are

```math
\boxed{
K(A_0,D;\infty)=20,
\qquad
K(A_1,D;\infty)=12.}                                \tag{EO.31}
```

Consequently, for every sufficiently large finite internal temperature
`t`, both children are exact pressure minimizers with identical pressure and
entropy, while

```math
\min_{\epsilon}\mathsf T_u(A_0,D,\epsilon)
-\min_{\epsilon}\mathsf T_u(A_1,D,\epsilon)
={\lambda^2\over2}u^4
 \{K(A_0,D;t)-K(A_1,D;t)\}+O(u^6)>0                \tag{EO.32}
```

for all sufficiently small positive bridge amplitudes `u`.  The same holds
if the unordered orientation pair is retained instead of its minimum.

*Proof.*  The order-two sector partition is independent of its sector sign,
and its only off-diagonal correlation is `epsilon s q`.  Hence the shared
sector weights are exactly those of the augmented order-eight child and are
independent of `epsilon`.  For `j=ell`, (EO.8) equals `a_(ik)`; for
`j!=ell`, it equals `epsilon b_(ik)`.  There are two ordered choices of each
kind, which proves (EO.30).

The exact ground-shell enumeration in FC.24 gives

```math
\sum_{i<k}a_{ik}(\infty)^2=3\quad(A_0),
\qquad
\sum_{i<k}a_{ik}(\infty)^2=1\quad(A_1).             \tag{EO.33}
```

On the same eight projective ground states, direct integer summation with
`tau=sign(H_A(X))` gives

```math
\sum_{i<k}b_{ik}(\infty)^2=7\quad(A_0),
\qquad
\sum_{i<k}b_{ik}(\infty)^2=5\quad(A_1).             \tag{EO.34}
```

Equations (EO.30), (EO.33), and (EO.34) prove (EO.31).  Finite Gibbs sums
are analytic in `e^(-t)`, so the strict separation persists for all
sufficiently large finite `t`.  Both matrices are then still exact pressure
minimizers.  Apply EO.2; because its leading coefficient is independent of
orientation and differs strictly, minimization over the two orientations
preserves the difference for sufficiently small `u`. `square`

The bridge amplitude in EO.32 is varied independently of the internal child
temperature.  Thus EO.4 is a decisive orientation-invariant failure of an
entropy/pressure-profile *response tangent*, not a claim that the two
physical-amplitude values at `u=t` differ.

## 6. Consequence for the ES.42 search

The triangle pair is an exact actual-minimizer falsifier to the fixed-sector
implication

```text
separate child pressure/entropy data plus unoriented replica-overlap data
    determine the labelled canonical row-TC path resource.
```

Even retaining the full such hierarchy does not repair the missing relative
sector alignment.  The two children are isomorphic for the complete declared
statistic and nevertheless have different fixed-sector physical `mathsf T`.
However, their two orientation values are merely exchanged.  No conclusion
about `min_epsilon mathsf T_epsilon` follows from the triangle.

EO.4 supplies the complementary orientation-invariant statement: even the
complete child pressure/entropy profile does not determine the minimum or
unordered row-TC response tangent.  Its two optimizing classes do have
different replica geometry, so EO.4 does not falsify a suitably oriented
overlap criterion.

Neither falsifier decides an asymptotic zero/positive-density phase.  The
precise combined conclusion is: entropy cannot decide even the unoriented
response tangent, while a separate-child unoriented overlap hierarchy cannot
label the fixed sectors.  An orientation-combined overlap theorem remains
open.

The smallest viable overlap statistic must retain cross-child `tau`-oriented
alignment.  The tensor `Gamma^(epsilon)` is the first such statistic and,
by EO.2, it controls the exact infinitesimal row-TC response with constant
`lambda^2/2`.  At physical amplitude in growing order, however, every fixed
overlap truncation remains subject to the archived nonuniform cumulant
barrier.  A useful optimizer-specific theorem must therefore either control
a growing oriented overlap hierarchy or prove deterministic synchronization
which makes that hierarchy a function of a smaller signed order parameter.
