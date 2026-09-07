# Wave 38 harmonic attack: collapse the matched covariance before estimating it

## Status

- **Verified:** after summing the nonreference vertex coordinates, the full
  matched covariance in (10.1066) is one state-space covariance
  `Cov_mu_t(g,L_t)`, where `L_t(d)` is the sum of the current conditional
  binary KL costs over vertex coordinates.  This identity keeps the omitted
  external-partition and exclusion-shape cancellation exact.
- **Verified:** the vertex weighted energy has an exact boundary formula in
  terms of endpoint binary cost and the **signed time integral** of this one
  covariance.  Consequently (10.983) and the transient-flux condition
  (10.1065) take positive parts much earlier than the parent-energy argument
  requires.
- **Verified sufficient lemma:** an approximate isotonicity estimate saying
  that `L_t` is close in `L^2(mu_t)` to a nondecreasing function of `g`
  controls the signed migration directly.  A pairwise discordance version is
  also given below.
- **Open:** no exact-minimizer estimate of the isotonic defect, endpoint
  vertex cost, or signed covariance is proved.  This is a sharper interface,
  not a proof of the harmonic package.
- **Numerical:** on `A_9,m=7,beta=2`, cancellation across vertex coordinates
  reduces the integrated adverse charge from `0.00285723` to `0.00146794`.
  The exact zero-defect monotonicity conjecture is already incompatible with
  the negative finite covariance, so only a quantitative approximate theorem
  is viable.

The checker `tmp/combined_covariance_r38_check.py` verifies the identities,
the matched fixed-size external-partition sum, and the boundary formula on
`A_4,A_6,A_8,A_9`.  Orientation, restoring correlation, and
adjacent-selector Hellinger control remain separate throughout.

## 1. Collapse from coordinate contexts to one state load

Let `V_*` be the `N=n-1` nonreference vertex coordinates in the oriented-cut
chart.  At interpolation time `t`, retain the notation from (10.978)--
(10.980): `M_(t,i)(e)` is the `i`-context law,

```math
X_{i,e}(t)=A'_{i,e}(t)
=\mathbb E_{\mu_t}[g\mid D_{-i}=e],
```

and

```math
k_{i,e}(t)
=D(\mu_t(D_i\mid e)\Vert\nu(D_i\mid e)).
```

Define the statewise total vertex KL load

```math
L_t(d)=\sum_{i\in V_*}k_{i,d_{-i}}(t).
\tag{R38.H1}
```

The context score has the same mean for every coordinate:

```math
\mathbb E_{M_{t,i}}X_{i,e}(t)
=\mathbb E_{\mu_t}g=:\psi'(t).
```

Since `k_(i,D_-i)` is measurable with respect to `D_-i`, the conditional
expectation tower gives

```math
\begin{aligned}
\sum_{i\in V_*}\operatorname{Cov}_{M_{t,i}}(X_{i,e},k_{i,e})
&=\sum_{i\in V_*}\mathbb E_{\mu_t}
  [(X_{i,D_{-i}}-\psi')k_{i,D_{-i}}]\\
&=\sum_{i\in V_*}\mathbb E_{\mu_t}
  [(g-\psi')k_{i,D_{-i}}]\\
&=\boxed{\operatorname{Cov}_{\mu_t}(g,L_t)}.
\end{aligned}
\tag{R38.H2}
```

This is an exact identity for every finite parent interpolation.  In the
matched problem it recombines the two terms of (10.1066) before any absolute
value is taken.

There is an equivalent coordinate-context formulation.  Let `Lambda_t` be
the law obtained by choosing `I` uniformly from `V_*` and then
`E sim M_(t,I)`.  Put `X=X_(I,E)` and `Y=k_(I,E)`.  Then

```math
\boxed{
\sum_{i\in V_*}\operatorname{Cov}_{M_{t,i}}(X_{i,e},k_{i,e})
=N\operatorname{Cov}_{\Lambda_t}(X,Y).
}
\tag{R38.H3}
```

Thus cancellation across different vertex coordinates is legitimate, not a
numerical coincidence produced by adding unrelated means.

## 2. The exact boundary identity and the minimal migration target

Set

```math
C_V(t)=\sum_{i\in V_*}\mathbb E_{M_{t,i}}k_{i,e}(t)
=\mathbb E_{\mu_t}L_t.
```

The one-edge identity `partial_t k=t A''`, differentiation of the moving
state law, and (R38.H2) give

```math
\boxed{
C_V'(t)
=t\mathcal E_{t,V}(g)
+\operatorname{Cov}_{\mu_t}(g,L_t),
}
\tag{R38.H4}
```

where `mathcal E_(t,V)` includes only nonreference vertex edges.  Since
`C_V(0)=0`, integration yields

```math
\boxed{
\int_0^1t\mathcal E_{t,V}(g)\,dt
=C_V(1)-\int_0^1\operatorname{Cov}_{\mu_t}(g,L_t)\,dt.
}
\tag{R38.H5}
```

Consequently the sharp combined migration condition is only

```math
\boxed{
\left[-\int_0^1\operatorname{Cov}_{\mu_t}(g,L_t)\,dt\right]_+
=O(n^{1/2-2c}).
}
\tag{R38.H6}
```

Together with `C_V(1)=O(n^(1/2-2c))`, this proves the desired vertex weighted
energy bound.  The conditions already in the ledger imply (R38.H6) through
the strict hierarchy

```math
\begin{aligned}
\left[-\int_0^1\sum_i c_i(t)\,dt\right]_+
&\le \int_0^1\left[-\sum_i c_i(t)\right]_+dt\\
&\le \int_0^1\sum_i[-c_i(t)]_+dt\\
&\le \sum_i\mathfrak T_i,
\end{aligned}
\tag{R38.H7}
```

where `c_i(t)=Cov_(M_(t,i))(X_i,k_i)` and the last line is (10.1064).
Thus (10.1065) is a useful selector-aware majorant but is not the logically
minimal migration theorem.  In particular, the proof should not discard
cross-coordinate or time cancellation unless the resulting quantity can
actually be controlled.

Equation (R38.H5) is also a warning: (R38.H6) is sharp bookkeeping, not
independent progress on the energy bound.  Given endpoint control, a large
negative signed covariance is exactly the excess vertex energy one is trying
to rule out.

## 3. A falsifiable approximate-order sufficient lemma

Let `D,D'` be independent with law `mu_t`.  Symmetrization of (R38.H2) gives

```math
\operatorname{Cov}_{\mu_t}(g,L_t)
=\frac12\mathbb E
[(g(D)-g(D'))(L_t(D)-L_t(D'))].
\tag{R38.H8}
```

Define the discordant-pair pressure

```math
\mathfrak D_t
=\frac12\mathbb E
\left[[-(g(D)-g(D'))(L_t(D)-L_t(D'))]_+\right].
```

Then

```math
[-\operatorname{Cov}_{\mu_t}(g,L_t)]_+\le\mathfrak D_t.
\tag{R38.H9}
```

Hence `integral mathfrak D_t dt=O(n^(1/2-2c))` is one combined sufficient
lemma.  It asks only that pairs whose parent likelihood score and total local
KL load have opposite order carry little weighted discrepancy.

A smoother version may be more tractable.  Let `M_up` be the nondecreasing
Borel functions on the range of `g` and define the isotonic defect

```math
R_t=\inf_{\phi\in\mathcal M_\uparrow}
\mathbb E_{\mu_t}[L_t-\phi(g)]^2.
\tag{R38.H10}
```

For every nondecreasing `phi`, independent-copy symmetrization gives
`Cov(g,phi(g))>=0`.  Cauchy--Schwarz therefore proves

```math
\boxed{
[-\operatorname{Cov}_{\mu_t}(g,L_t)]_+
\le\sqrt{\operatorname{Var}_{\mu_t}(g)R_t}.
}
\tag{R38.H11}
```

The concrete reduced sufficient lemma

```math
\boxed{
\int_0^1
\sqrt{\operatorname{Var}_{\mu_t}(g)R_t}\,dt
=O(n^{1/2-2c})
}
\tag{R38.H12}
```

therefore proves (R38.H6).  It is genuinely about the **combined** score:
neither `log B_i` nor `H'_i` is fitted separately.  A theorem could replace
the full isotonic cone by an explicit nondecreasing comparison profile.

Exact monotonicity `R_t=0` is too strong.  The negative integrated covariance
on finite exact minimizers `A_8` and `A_9` is numerical evidence against any
universal zero-defect statement.  The project-scale version (R38.H12) is
falsified only by an unbounded actual-minimizer family for which its left side
is `omega(n^(1/2-2c))`; such a falsifier would defeat this implementation,
not every possible proof of (R38.H6).

## 4. Where the matched external partitions enter

Normalize the omitted external partition by

```math
\widetilde B_i(e)=\frac{B_i(e)}{\mathbb E_\nu U}.
```

The matched identities from Wave 37 become

```math
\widetilde B_i(D_{-i})=f(D)r_i(D),
\qquad
X_{i,e}=\log\widetilde B_i(e)+\mathcal H'_{i,e}(t).
\tag{R38.H13}
```

Moreover

```math
\mathcal H'_{i,e}(t)
=\mathbb E_{\mu_t(D_i\mid e)}[-\log r_i(D)].
```

Thus the cancellation already occurs in the coordinate mean:

```math
\mathbb E_{M_{t,i}}X_{i,e}
=\mathbb E_{\mu_t}[\log\widetilde B_i-\log r_i]
=\mathbb E_{\mu_t}g.
\tag{R38.H14}
```

Multiplication by the context-measurable `k_(i,e)` and the conditional tower
then gives (R38.H2).  This is the precise algebraic reason that separately
bounding the two terms in (10.1066) can destroy the relevant cancellation.

There is also a pointwise fixed-size relation over all `n` vertices,
including the chart reference vertex:

```math
\boxed{
\sum_{i=1}^n\widetilde B_i(D_{-i})
=(n-m)f(D).
}
\tag{R38.H15}
```

It follows immediately from `sum_i r_i(D)=n-m`.  This is real coupling among
the external partitions, but it is linear rather than logarithmic.  By
itself it controls only the omission-weighted entropy

```math
\sum_i r_i\log\frac{n-m}{r_i}
=(n-m)H\left(\frac{r_i}{n-m}\right)
\le(n-m)\log n,
\tag{R38.H16}
```

which has the wrong weight and is far above the project target at fixed
density.  A proof of (R38.H12), or directly of (R38.H6), therefore needs an
additional exact-minimizer correlation between the parent score `g` and the
total local cost `L_t`; fixed selector size alone does not supply it.

The dynamic erasure identity remains compatible with this collapse.  For
`bar r_(i,t)(e)=E_(Q_(t,e)) r_i(e,D_i)`, (10.1061) gives

```math
Q_{t-1,e}(b)
=Q_{t,e}(b)\frac{r_i(e,b)}{\bar r_{i,t}(e)},
```

and hence

```math
q^{\leftarrow}_{i,e}(t)
=\log\bar r_{i,t}(e)
-\mathbb E_{Q_{t,e}}\log r_i(e,D_i),
\qquad k_{i,e}(t)\le tq^{\leftarrow}_{i,e}(t).
\tag{R38.H17}
```

One possible successor is to prove the order/discordance estimate for
`L_t` through this selector-erasure representation, while retaining the sum
over `i` until after fixed-size cancellation is used.

## 5. Finite audit

The checker uses 96-point Gauss--Legendre quadrature.  The boundary error in
(R38.H5) is below `2e-15` in the displayed cases, and the pointwise collapse,
mixture, pair-symmetrization, matched-split, and fixed-size identities agree
to at worst `1.1e-12`.

```text
case             signed integral    [-signed]_+   aggregate adverse   coordinate adverse
A4, b=.5,m=3      0.000006998        0              0                   0
A6, b=.5,m=3      0.000005609        0              0                   0
A8, b=2,m=5      -0.007441065        0.007441065    0.007441065         0.007441065
A9, b=2,m=7      -0.001467944        0.001467944    0.001467944         0.002857230
```

For the last case, the Wave 37 separate matched pieces had signed integrals
`+0.08670755` and `-0.08817549`; recombination first leaves
`-0.00146794`.  Summing coordinates before taking the adverse part then
removes another `0.00138929`.  The Wave 37 dynamic-flux majorant is about
`0.00547038`.  These are finite numerical comparisons only, but they show
that every early positive-part operation can have a material cost even on a
known exact minimizer.

## Frontier

The mathematically correct next harmonic target is the aggregate signed
lower bound (R38.H6), not separate absolute estimates for the omitted
partition and exclusion shape, and not necessarily the coordinatewise
positive-part target in (10.983).  Equations (R38.H2)--(R38.H5) reduce it to
one global object: the covariance between parent score `g` and total local
binary-KL load `L_t`.

The sharply falsifiable structured successor is approximate order:
prove (R38.H12), or a stronger usable discordance estimate (R38.H9), from the
fact that the `B_i` are matched principal-restriction external partitions
and obey (R38.H15).  A generic exact-monotonicity sign theorem is already the
wrong target.  Endpoint vertex cost, restoring comparison, orientation, and
adjacent-selector Hellinger estimates remain independent open inputs.
