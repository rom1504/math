# Logarithmic statistics of an automatic extremal phase

Status: rigorous general theorem draft.  It applies to every sequence with
the uniform mantissa law proved in Theorem 30.1.  It does not assert ordinary
pointwise convergence.

## Theorem PA.1 (unique logarithmic phase law)

Fix an integer `h>=2`.  Let `(q_n)` be bounded, put

```math
r(n)=floor(log_h n),\qquad t_n=n/h^(r(n)) in [1,h),
```

and suppose one continuous `L:[1,h]->R` satisfies

```math
sup_(h^r<=n<h^(r+1))|q_n-L(n/h^r)|->0.                \tag{PA.1}
```

Then, for every continuous test function `psi`,

```math
{1\over log N}\sum_(n<=N){psi(q_n)\over n}
\longrightarrow
{1\over log h}\int_1^h {psi(L(t))\over t}\,dt.       \tag{PA.2}
```

Thus the logarithmic empirical law of the extremal response is the
pushforward of `dt/(t log h)` by `L`.  It is a singleton exactly when `L` is
constant.

### Proof

On a complete scale block `h^r<=n<h^(r+1)`, uniform continuity of `psi` and
(PA.1) replace `psi(q_n)` uniformly by `psi(L(n/h^r))`.  The weighted Riemann
sum gives

```math
\sum_(n=h^r)^(h^(r+1)-1){psi(q_n)\over n}
\longrightarrow \int_1^h{psi(L(t))\over t}\,dt.      \tag{PA.3}
```

The error is `o(1)` per sufficiently late block.  There are
`log_h N+O(1)` complete blocks; the finitely many early blocks and one final
partial block contribute `O(1)` after division by `log N`.  Cesaro averaging
of the block errors proves (PA.2). `square`

## Theorem PA.2 (ordinary means retain a phase)

For fixed `s in [1,h]` and `N_R=floor(s h^R)`,

```math
{1\over N_R}\sum_(n<=N_R)q_n
\longrightarrow
C(s):={1\over s}\left\{
 {1\over h-1}\int_1^hL(t)dt+\int_1^sL(t)dt
\right\}.                                             \tag{PA.4}
```

Hence ordinary Cesaro averaging need not remove discrete scale phase;
logarithmic averaging always does.

### Proof

An unweighted Riemann sum on one complete block gives

```math
\sum_(n=h^r)^(h^(r+1)-1)q_n
=h^r\int_1^hL(t)dt+o(h^r).                            \tag{PA.5}
```

Summing the geometric block sizes through `r=R-1` contributes
`h^R int_1^hL/(h-1)+o(h^R)`.  The current partial block contributes
`h^R int_1^sL+o(h^R)`.  Division by `s h^R+o(h^R)` proves (PA.4). `square`

## Theorem PA.3 (every positive power bias retains the phase)

Under the hypotheses of Theorem PA.1, fix `alpha>0`.  For
`s in [1,h]` and `N_R=floor(s h^R)`, put

```math
C_alpha(s)={alpha\over s^alpha}\left\{
 {1\over h^alpha-1}\int_1^h t^(alpha-1)L(t)dt
 +\int_1^s t^(alpha-1)L(t)dt
\right\}.                                             \tag{PA.6}
```

Then

```math
{\sum_(n<=N_R)n^(alpha-1)q_n
 \over \sum_(n<=N_R)n^(alpha-1)}
\longrightarrow C_alpha(s).                           \tag{PA.7}
```

Moreover `C_alpha` is constant on `[1,h]` if and only if `L` is
constant.  Thus every positive power bias remembers the discrete scale
phase; the logarithmic boundary `alpha=0` in Theorem PA.1 is qualitatively
different.

### Proof

On a complete block, the weighted Riemann sum is

```math
\sum_(n=h^r)^(h^(r+1)-1)n^(alpha-1)q_n
=h^(alpha r)\int_1^h t^(alpha-1)L(t)dt+o(h^(alpha r)). \tag{PA.8}
```

The geometric sum of complete blocks through `r=R-1` and the final partial
block give the numerator

```math
h^(alpha R)\left\{
 {1\over h^alpha-1}\int_1^h t^(alpha-1)L(t)dt
 +\int_1^s t^(alpha-1)L(t)dt+o(1)\right\}.            \tag{PA.9}
```

The denominator is `(s h^R)^alpha/alpha+o(h^(alpha R))`, proving
(PA.7).  Direct differentiation of (PA.6) gives

```math
C_alpha'(s)={alpha\over s}\bigl(L(s)-C_alpha(s)\bigr).\tag{PA.10}
```

If `C_alpha` is constant, (PA.10) forces `L=C_alpha`; the converse is
immediate.  More generally the entire phase is reconstructed from the
power-mean phase by

```math
L(s)=C_alpha(s)+{s\over alpha}C_alpha'(s).             \tag{PA.11}
```

Thus positive power averaging is not even a nontrivial quotient of the
continuous phase profile. `square`

### Interpretation

The measure `dt/(t log h)` is Haar probability on one fundamental
multiplicative interval.  When `L(1)=L(h)`, as in the Walsh application, it
descends to Haar probability on the scale-phase circle.  Theorem PA.1
therefore averages all scale epochs equally.  Positive
power weights put a fixed geometric fraction of their mass on the most recent
epoch, and (PA.10) shows that no nonconstant continuous phase can disappear
under such an averaging rule.  This is a statement about the declared
sampling query, not about recovery or synchronization of the landscape.

## Corollary PA.4 (logarithmic statistics forget phase order)

The logarithmic empirical law in Theorem PA.1 is not a complete invariant
of the phase profile.  For example, with `u=log_h t`, the two continuous
profiles

```math
L_1(t)=cos(2 pi u),\qquad L_2(t)=cos(4 pi u)           \tag{PA.12}
```

have the same pushforward law under `dt/(t log h)`, but they are distinct.
For every fixed `alpha>0`, their power-mean phase functions cannot coincide,
because (PA.11) would then imply `L_1=L_2`.

Consequently a logarithmic thermodynamic statistic can exist uniquely while
forgetting information that every last-scale-sensitive continuation can
recover.  What is lost is the cyclic ordering of response values across
scale, not their empirical distribution. `square`

## Walsh consequence

The explicit prefix family in Corollary 30.2 has a nonconstant continuous
`L`, so it has a nontrivial limiting logarithmic response distribution even
though its pointwise response does not converge.  The logarithmic mean exists
canonically, whereas its ordinary Cesaro means retain the nonconstant phase
`C(s)`.  In fact Theorem PA.3 says that every positive power-weighted mean
retains a nonconstant subsequential phase.

This is a useful distinction for limit theories: compact scale phase can be
integrated out by a declared sampling law, but that does not synchronize the
underlying extremal state or prove an all-order pointwise limit.
