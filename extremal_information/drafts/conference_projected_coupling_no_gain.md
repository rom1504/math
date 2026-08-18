# Projected couplings cannot create a favorable conference-pressure phase

**Status.** Task-local theorem.  This abstracts the mechanism behind the
full-range row-magnitude theorem.  It is a one-sided pressure comparison,
not an assertion that the conditioned law has the same full pressure rate.

## 1. Setup

Let `A_r` be symmetric conference signings, fix an orientation, and write

```math
f_r(C)=\log\left[2^{-2r}\sum_{x,y}
\cosh\left\{{\beta\over\sqrt{2r}}
\big(H_A(x)+\epsilon H_A(y)+x^TCy\big)\right\}\right],
\tag{PC.1}
```

where

```math
0<\beta<{\sqrt2\over6}.
\tag{PC.2}
```

Let `W_r` be an iid Rademacher bridge and let `B_r` be any random real
bridge coupled to `W_r`.  Let `P_r` be a deterministic orthogonal projection
on the column space and put

```math
W_r^\circ=W_r(I-P_r),\qquad
B_r^\circ=B_r(I-P_r).
\tag{PC.3}
```

The entries of `B_r` need not be independent or exact signs in this theorem.

## 2. Projected-coupling criterion

### Theorem PC.1 (sublinear exceptional response rank is one-sidedly harmless)

Assume the following three conditions.

1. There is a fixed `delta>0` satisfying

   ```math
   {\beta(3+\delta)\over\sqrt2}<{1\over2},
   \tag{PC.4}
   ```

   and events `G_r` with `Pr(G_r^c)=o(1)` on which

   ```math
   \max\{\|W_r\|_{op},\|B_r^\circ\|_{op}\}
   \le(2+\delta)\sqrt r.
   \tag{PC.5}
   ```

   (The omitted bound for `W_r^circ` follows because `I-P_r` is a
   contraction.)

2. The projected coupling has subcritical Frobenius cost:

   ```math
   \mathbb E\|(B_r-W_r)(I-P_r)\|_F=o(r).
   \tag{PC.6}
   ```

3. The removed components have subcritical nuclear cost:

   ```math
   \mathbb E\|W_rP_r\|_*=o(r^{3/2}),
   \qquad
   \mathbb E\|B_rP_r\|_*=o(r^{3/2}).
   \tag{PC.7}
   ```

Then, with

```math
h_\beta=2\psi(\beta/\sqrt2)+{\beta^2\over4},
```

one has

```math
\boxed{
\mathbb E\left[\left(h_\beta-{f_r(B_r)\over r}\right)_+\right]
\longrightarrow0.}
\tag{PC.8}
```

Thus a favorable linear pressure phase cannot be produced solely by a
sublinear-nuclear-cost component outside a regularly coupled bulk.

**Proof.**  Choose `kappa<1/2` strictly above the left side of (PC.4).
On `G_r`, both projected parents lie in the same strict operator-temperature
ball.  The audited high-temperature pressure stability theorem gives

```math
|f_r(B_r^\circ)-f_r(W_r^\circ)|
\le {K_\kappa\beta\over\sqrt2}
\|(B_r-W_r)(I-P_r)\|_F.
\tag{PC.9}
```

Conditions (PC.6) and (PC.5) make the right comparison `o(r)` on `G_r` in
mean.  The same stability theorem, now applied
between `W_r` and `W_r^circ`, gives on the common regular event

```math
|f_r(W_r)-f_r(W_r^\circ)|
\le {K_\kappa\beta\over\sqrt{2r}}\|W_rP_r\|_*.
\tag{PC.10}
```

It remains to restore `B_rP_r`; the full endpoint need not be regular.
On `G_r`, consider pressure on the affine line from `B_r^circ` to `B_r`.
Convexity and the high-temperature covariance bound only at the base point
give the global supporting-line estimate

```math
f_r(B_r)\ge f_r(B_r^\circ)
-{K_\kappa\beta\over\sqrt{2r}}\|B_rP_r\|_*.
\tag{PC.11}
```

On `G_r`, chain (PC.9)--(PC.11) and (PC.10) to compare `f_r(B_r)` from below
with `f_r(W_r)`.  Every comparison error has expectation `o(r)` by
(PC.6)--(PC.7), while the uniform conference theorem gives

```math
\mathbb E[(h_\beta-f_r(W_r)/r)_+]\longrightarrow0.
```

On `G_r^c`, the normalized positive shortfall in (PC.8) is at most
`h_beta`, because pressure is nonnegative.  Therefore `Pr(G_r^c)=o(1)` is
enough; no uniform-integrability assumption on the irregular endpoint is
needed.  This proves (PC.8).  `square`

## 3. Rank interpretation

If `rank(P_r)=k_r`, then

```math
\|CP_r\|_*\le\sqrt{k_r}\|CP_r\|_F.
\tag{PC.13}
```

For the iid bridge,

```math
\mathbb E\|W_rP_r\|_F^2=rk_r,
\qquad
\mathbb E\|W_rP_r\|_*\le k_r\sqrt r.
\tag{PC.14}
```

The same estimate, up to a fixed factor, holds for a row law obtained by
conditioning each uniform row on an event of probability bounded below.
Consequently (PC.7) is automatic for such laws whenever `k_r=o(r)`.

The nontrivial hypothesis is (PC.5): after deleting the exceptional response
subspace, the coupling must have a sharp regular bulk.  In the row-magnitude
theorem `k_r=1`, and nested Hamming-layer coupling proves this via a
uniform-subset covariance identity and matrix Bernstein.

## 4. Consequence for the basin search

The theorem rules out a broader mechanism than one distinguished row
magnitude.  A speed-`r` favorable law cannot be explained by finitely many,
or even `o(r)` suitably controlled, exceptional column-response directions
attached to a regular coupled bulk.  A surviving construction must violate
at least one explicit condition: it must have extensive exceptional rank,
fail every projected regular coupling, or retain leading nuclear mass after
all sublinear-dimensional projections.

This is not yet an information lower bound: a nonlinear law may fail the
coupling criterion while still having only linear relative entropy.  It is
also not a theorem about arbitrary row dependence.  Its value is to replace
the vague phrase “operator irregular” by a falsifiable structural target:
**irreducible extensive irregularity after every low-rank response peel.**
