# The marginal-union wall for an iid dense bridge

**Status.** Proved cross-order method-class no-go.  This note isolates the exact obstruction
to using an unbiased iid bridge together with marginal first-moment
accounting at a comparable split.  The obstruction persists even if the
certificate retains the complete energy tables of both children and uses
the two one-sided caps separately.

It does **not** say that a good bridge does not exist.  It leaves open a
correlated analysis of the violation events, a non-product bridge law, or a
bridge chosen deterministically from the full child landscapes.

## 1. Setup

For a hollow complete signing `A`, put

```math
H_A(x)=\sum_{i<j}a_{ij}x_ix_j,
\qquad Q(A)=\max_x|H_A(x)|.
```

Let `A` and `C` have orders `m` and `n`, with

```math
p=Q(A),\qquad q=Q(C),\qquad
R=(p^{2/3}+q^{2/3})^{3/2}.
\tag{1.1}
```

For child spins `x,y`, write

```math
U_{x,y}=H_A(x)+H_C(y).
\tag{1.2}
```

Since `0<2/3<1`,

```math
R\ge p+q\ge |U_{x,y}|,
\tag{1.3}
```

with strict first inequality when `p,q>0`.  Let `B` be an `m` by `n`
matrix of independent unbiased signs.  For every fixed `x,y`,

```math
x^TBy\ \buildrel d\over=\ S_{mn},
\tag{1.4}
```

where `S_d` is a sum of `d` independent Rademacher variables.

## 2. Exact marginal-union barrier

### Theorem MU.1 (energy-resolved Hoeffding ledgers have the raw entropy wall)

Let `T=R+E`, where `E>=0`, and form the fully energy-resolved Hoeffding
ledger

```math
\Phi_T(A,C)=
\sum_{x,y}\sum_{\sigma\in\{+-1\}}
\exp\left\{-{(T+\sigma U_{x,y})^2\over2mn}\right\}.
\tag{2.1}
```

Then

```math
\boxed{
\Phi_T(A,C)
\ge 2^{m+n}\exp\left\{-{T^2\over2mn}\right\}.}
\tag{2.2}
```

Consequently, the sufficient condition `Phi_T<1` forces

```math
\boxed{
T>\sqrt{2mn(m+n)\log2}.}
\tag{2.3}
```

**Proof.**  For each pair `(x,y)`, choose `sigma` opposite to the sign of
`U_(x,y)`.  By (1.3), its one-sided threshold is positive and

```math
|T+\sigma U_{x,y}|=T-|U_{x,y}|\le T.
```

Thus that single summand in (2.1) is at least
`exp(-T^2/(2mn))`.  Sum it over all `2^(m+n)` pairs.  This proves (2.2),
and (2.3) follows immediately. `square`

The point is the absolute parent cap.  Child energy can make one polarity
harder, but it makes the opposite polarity at least as easy.  No information
about the distribution, tails, overlaps, or multiplicities of the values
`U_(x,y)` can lower (2.2).

The same threshold is intrinsic to a union bound using the **exact**
marginal probabilities, rather than to Hoeffding slack.

### Theorem MU.2 (exact marginal probabilities have the same exponent)

Define

```math
\Psi_T(A,C)=\sum_{x,y}
 \Pr_B\{|U_{x,y}+x^TBy|>T\}.
\tag{2.4}
```

For every `T>=R`,

```math
\boxed{
\Psi_T(A,C)\ge2^{m+n}\Pr\{S_{mn}>T\}.}
\tag{2.5}
```

If `m,n` are comparable and `T=O((m+n)^(3/2))`, then

```math
\log\Psi_T(A,C)
\ge(m+n)\log2-{T^2\over2mn}-O(\log(m+n)).
\tag{2.6}
```

In particular, `Psi_T<1` requires

```math
\boxed{
T\ge\sqrt{2mn(m+n)\log2}-o((m+n)^{3/2}).}
\tag{2.7}
```

**Proof.**  If `U_(x,y)>=0`, the positive violation contains
`{S_(mn)>T-U_(x,y)}`; if `U_(x,y)<0`, symmetry gives the analogous negative
violation.  In either case its probability is at least
`Pr{S_(mn)>T}`, proving (2.5).

For completeness, choose the least attainable value `k>T`.  Stirling's
formula for the single binomial mass at `S_(mn)=k` gives

```math
\log\Pr\{S_{mn}=k\}
=-mn I(k/(mn))-O(\log(m+n)),
```

where

```math
I(t)={1\over2}[(1+t)\log(1+t)+(1-t)\log(1-t)].
```

At a comparable split, `k=O((m+n)^(3/2))` and `mn=Theta((m+n)^2)`, so

```math
mn I(k/(mn))={k^2\over2mn}+O(k^4/(mn)^3)
={T^2\over2mn}+O(1).
```

Insert this lower bound in (2.5) to obtain (2.6).  Equation (2.7) follows.
`square`

Thus replacing a subgaussian tail estimate by an exact binomial tail does
not repair the method.  The missing ingredient would have to exploit the
dependence among the exponentially many violation events.

## 3. Direct recurrence-defect consequence

Put `b_k=M_k^(2/3)`.  If the children are optimal, then the desired parent
radius is

```math
R=(b_m+b_n)^{3/2}.
```

Any iid-marginal-union certificate at a comparable split therefore incurs
the energy-radius defect

```math
\boxed{
E\ge
\left[\sqrt{2mn(m+n)\log2}-(b_m+b_n)^{3/2}\right]_+
-o((m+n)^{3/2}),}
\tag{3.1}
```

and the corresponding defect in the proposed `b` recurrence is at least

```math
\boxed{
\left[\{2mn(m+n)\log2\}^{1/3}-(b_m+b_n)\right]_+
-o(m+n).}
\tag{3.2}
```

For the equal split `m=n=r`, write `p=M_r`.  Then

```math
R=2\sqrt2\,p,
\qquad T\ge(2\sqrt{\log2}-o(1))r^{3/2}.
\tag{3.3}
```

The rigorous upper frontier

```math
M_r\le(1/2+o(1))r^{3/2}
\tag{3.4}
```

therefore forces

```math
\boxed{
E\ge
(2\sqrt{\log2}-\sqrt2-o(1))r^{3/2}.}
\tag{3.5}
```

The leading constant is

```math
2\sqrt{\log2}-\sqrt2=0.2508956599\ldots.
```

More directly, the certified `b`-defect is at least

```math
\boxed{
\left[(2\sqrt{\log2})^{2/3}-2^{1/3}-o(1)\right]r
= (0.1449241895\ldots-o(1))r.}
\tag{3.6}
```

In terms of the parent order `N=2r`, this is
`(0.0724620948...-o(1))N`.  It rules out every
`O(N^(1-delta))`, `delta>0`, recurrence defect within this method class.
These are floors on what this certificate can prove, not lower bounds on
the globally optimized value `M_(m+n)`.

## 4. Archive comparison and scope

This note was derived before consulting the archive.  The later comparison
found two nearby but distinct uses of a raw bridge union bound:

* `bounded_cap_near_top_tail.md` invokes the archived random-bridge estimate
  only for a near-order deletion inequality.
* Lemma CT.1 of `thin_tail_entropy_bridge_no_go.md` uses the safe cap
  `2r^(3/2)` and obtains failure exponent
  `2-2log2` for the bridge cut norm.

The increment here is the sharp target-relative statement (2.2)--(3.6).
It shows that allowing a different soft cap for every child-state pair, and
even accounting with exact marginal binomial probabilities, cannot cross
the equal-split target furnished by the `1/2` upper frontier.  The reason is
the polarity choice in the proof of MU.1, not a coarse estimate of the
typical child-energy bulk.

The precise scope is:

1. bridge entries are independent unbiased signs;
2. feasibility is certified by summing one-constraint marginal failure
   probabilities, or upper bounds for those probabilities;
3. the desired parent constraint is two-sided absolute energy.

Correlated violation accounting, deterministic or correlated bridge design,
and a joint child--bridge optimization remain outside the theorem.
