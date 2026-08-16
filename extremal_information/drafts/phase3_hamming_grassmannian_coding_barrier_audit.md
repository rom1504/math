# Audit of the Hamming Grassmannian coding barrier

**Verdict.** HG.1--HG.6 in the companion draft are valid with the stated
strict thresholds.  The main reduction is finite and self-contained; it does
not assume existence of an asymptotic coding-rate limit.

## 1. Threshold and floor audit

Hausdorff distance is integer-valued.  Therefore, for every real `Delta`,

```math
d_{\rm Hs}>\Delta
\quad\Longleftrightarrow\quad
d_{\rm Hs}\ge\lfloor\Delta\rfloor+1.
```

The draft puts `t=floor(Delta)`, so the associated coding distance is exactly
`d=t+1`.  The disjoint-ball radius in the Hamming bound is

```math
\left\lfloor{d-1\over2}\right\rfloor
=\left\lfloor{t\over2}\right\rfloor,
```

as used in HG.3 and HG.5.  Puncturing exactly `t` coordinates gives fibre
diameter `t`, so equal projected carriers cannot coexist in a strict `>t`
packing.  There is no missing `+1` in that quotient argument.

For presented responses, the carrier law loses `p`, not `2p`, because both
one-sided presentation errors lie in the common interval `[0,p]`.  Therefore
the carrier thresholds `s+p` and `s-p` in HG.3a are correct.

For `t=L=D-k`, Singleton gives one word and the evaluation construction in
HG.6 is interpreted as the singleton code.  For `t>L`, the chart alphabet
has no two words at distance `t+1`; HG.10 is therefore stated only for
`t<=L`.

## 2. Low-weight line exception

The exact equivalence is

```math
d_{\rm Hs}(L_v,L_w)>t
\Longleftrightarrow
\operatorname{wt}(v+w)>t
\text{ and }
\max\{\operatorname{wt}(v),\operatorname{wt}(w)\}>t.
```

It is **not** legitimate to identify the line metric pointwise with Hamming
distance.  The proof of HG.2 handles the discrepancy optimally: a packing
contains at most one low-weight representative.  Removing that representative
and adjoining zero preserves cardinality; if none exists, adjoining zero
increases cardinality by one.  This is exactly why the finite result is the
two-sided `A_2-1 <= Pi <= A_2` sandwich rather than an unjustified equality.

## 3. Asymptotic-rate audit

No limit for

```math
D^{-1}\log_2 A_2(D,\lfloor\delta D\rfloor+1)
```

is asserted.  HG.2 transfers its limsup and liminf separately.  The additive
one-point ambiguity is exponentially negligible for `0<delta<1/2`, because
the Gilbert bound is exponential there.

The strict quotient gap uses only the Hamming upper bound:

```math
\log_2\Pi_{D,1}(t_D)
\le D-DH_2(\delta/2)+o(D),
```

whereas the puncturing state space has logarithm
`D-t_D=(1-delta)D+o(D)`.  Their difference is at least
`(H_2(delta/2)-delta-o(1))D`; strict concavity gives positivity.  This does
not assume that nonlinear codes outperform linear codes.

## 4. General-chart audit

For one fixed information set, systematic form is unique.  The comparison

```math
d_{\rm Hs}(C_X,C_Y)
\le d_{\rm col}(X,Y)
```

uses the same coefficient vector `u` only as a candidate match, so it is an
upper bound and has the correct direction: a Hausdorff-separated family must
be column-separated.  The converse is false because a different coefficient
vector can improve the match.  HG.6 deliberately saturates only this
necessary column-code condition and does not claim that Reed--Solomon words
give a Grassmannian packing.

The factor `binom(D,k)` is harmless on the `D^2` scale but not always on the
fixed-`k` `D` scale.  HG.9 retains it exactly; no fixed-`k` exponent claim
drops it silently.

## 5. What is and is not resolved

The draft proves two structural facts:

1. unrestricted determination of the Grassmannian exponent contains the
   classical nonlinear binary coding-rate problem at `k=1`; and
2. the anticode quotient has a positive linear state-count excess there.

It does not determine the `k=kappa D` exponent.  In that regime the column
alphabet has size `2^(kappa D)`, and the systematic shadow meets Singleton.
Any improvement must exploit coherent same-switch recoupling rather than a
standard coding bound on column differences.

The exhaustive verifier checks the exact formulas and inequalities, not the
asymptotic entropy estimates, which follow from standard binomial bounds.
