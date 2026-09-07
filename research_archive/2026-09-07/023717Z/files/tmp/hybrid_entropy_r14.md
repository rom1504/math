# Wave 14: the critical entropy route collapses to cross-only overshoot

This memo audits ledger (10.549) against the exact Kwan--Sauermann results.
The main conclusion is stronger than an entropy estimate: at every sublinear
block scale, hybrid selection is asymptotically irrelevant.  The sole
remaining issue is the norm of the cross-only mosaic.

## 1. Exact two-sided comparison

Let `C` have the original cross blocks and zero diagonal blocks.  For any
choice of local minimizers `G_i`, put

```math
K_{\mathbf G}=\bigoplus_iG_i,
\qquad B_0=\sum_iq_{m_i}.
```

The triangle and reverse-triangle inequalities give, for **every** tuple,

```math
\boxed{
\left|Q(C+K_{\mathbf G})-Q(C)\right|
\le Q(K_{\mathbf G})
\le B_0.
}
\tag{R14.1}
```

The lower bound can also be seen by evaluating `C+K_G` at a Boolean ground
state of `C` and losing at most `B_0`.  If

```math
E=\min_{\mathbf G}Q(C+K_{\mathbf G})-q_n,
\qquad d_C=Q(C)-q_n,
```

then

```math
\boxed{|E-d_C|\le B_0.}
\tag{R14.2}
```

The standard random-sign union bound recorded in the ledger is stated for
the undoubled Hamiltonian optimum `M_m`.  Since `q_m=2M_m`, its correct
form in the present `Q`-normalization is

```math
q_m\le2\sqrt{2\binom m2(m+2)\log2}
=2\sqrt{m(m-1)(m+2)\log2}.
```

Since

```math
(m-1)(m+2)\le\frac98m^2
```

(the difference after multiplying by eight is `(m-4)^2`),

```math
B_0\le\kappa\sum_im_i^{3/2}
\le\kappa n\sqrt{s},
\qquad
\kappa=\frac{3\sqrt{\log2}}{\sqrt2}<1.767,
\qquad s=\max_i m_i.
\tag{R14.3}
```

Consequently, if `s=o(n)`, then `B_0=o(n^{3/2})`, and (R14.2) proves the
exact equivalence

```math
\boxed{
E=o(n^{3/2})
\quad\Longleftrightarrow\quad
[Q(C)-q_n]_+=o(n^{3/2}).
}
\tag{R14.4}
```

Necessity follows from `E\ge Q(C)-q_n-B_0`; sufficiency follows from
`E\le[Q(C)-q_n]_++B_0`.  At the critical scale `s=Theta(sqrt n)`, the
uncertainty is only `B_0=O(n^{5/4})`.

Thus (10.549) is unnecessary for the `o(n^{3/2})` goal.  In its own notation,
choosing `t=[Q(C)-q_n]_++B_0` makes the low-margin set empty.  More sharply,
no entropy or correlated selection can repair a leading cross overshoot,
because the reverse half of (R14.1) applies to every tuple.  The additional
theorem actually required is precisely

```math
[Q(C)-q_n]_+=o(n^{3/2})
```

for cross mosaics obtained from global minimizers and sublinear-block
partitions.

## 2. What the two Kwan--Sauermann papers actually prove

### Point anti-concentration (arXiv:2312.13826)

Theorem 1.1 states: if a real polynomial `P` of degree at most two in `n`
independent Rademachers cannot be made constant by specifying outcomes of
any `m-1` variables, then

```math
\sup_x\Pr(P(\xi)=x)\le C/\sqrt m
```

for an absolute constant `C`.  This is a point-probability theorem.  The
paper explicitly lists a small-ball counterpart as an open direction.

For a cross-only signing mosaic with largest block `s`, the coefficient
support is the complete multipartite graph.  Any fixing set which makes the
restriction constant must cover every quadratic monomial, so its complement
lies in one block.  Hence the polynomial robustly depends on at least `n-s`
variables, and the theorem gives

```math
\sup_x\Pr(z^{\mathsf T}Cz=x)=O((n-s)^{-1/2}).
\tag{R14.5}
```

Even on one exact level this permits

```math
O(2^n/\sqrt n)
```

states, whose logarithm is `n log 2-O(log n)`, not `o(n)`.  For an interval
of width `W`, the energies occupy a lattice of spacing four, so a union bound
only gives

```math
2^n\min\{1,O((W+1)/\sqrt n)\}.
```

This is trivial at the relevant `W=O(n^{5/4})`.  Point
anti-concentration controls the multiplicity of a scalar value; (10.549)
needs few near-cap states (or a legitimate quotient of their random
replacement events).  These are different, essentially opposite,
statements.

### Algebraic inverse theorem (arXiv:1909.02089, v2)

Theorems 1.1 and 1.2 say the following.  Fix `r>=3` and `epsilon>0`.  If a
quadratic polynomial with coefficients bounded by one (respectively,
degree-two coefficients in a fixed finite set `S`) has an atom of size at
least

```math
C(r,\epsilon[,S])
\frac{(\log n)^{r/2}}{n^{1-2/(r+2)}},
```

then it is within coefficient `L_1` distance `epsilon n^2` (respectively,
differs in at most `epsilon n^2` coefficients) from a homogeneous quadratic
form of rank strictly less than `r`.  The contrapositive, even when its
low-rank-distance hypothesis is proved, saves only `O(log n)` from the
logarithm of an exact fibre.  The constants are for fixed `r`; the theorem
does not permit an uncontrolled choice `r=r(n)`.

Their Lemma 3.2 is a genuine small-ball result under a stronger robust-rank
hypothesis: `delta n` disjoint `r`-tuples of rows must each form a
`delta`-non-degenerate `r x n` matrix.  With coefficients bounded by one it
gives

```math
\Pr\{|P(\xi)-x|\le n^{2/(r+2)}\}
\le C(r,\delta)
\frac{(\log n)^{r/2}}{n^{1-2/(r+2)}}.
```

That hypothesis is not supplied by global minimality.  Even granting it,
covering a window of width `n^{5/4}` by these small balls yields the trivial
factor `n^{1/4}(\log n)^{r/2}`.  Neither the direct nor inverse result is a
large-deviation entropy theorem.

If one wanted a refinement below the `B_0` scale, the missing external-style
result would have to be an upper-tail complexity statement such as

```math
\log\#\{[z]:Q(C)-|z^{\mathsf T}Cz|\le O(B_0)\}=o(n),
```

or, more faithfully, a bound on the number of distinct replacement-event
profiles.  It would need the special hypothesis that `C` is cut out of a
global minimizer.  Point anti-concentration and fixed-rank inverse theory do
not imply it.

## 3. Finite audit

For the `3+6` A9 cross mosaic, exact projective enumeration gives `Q(C)=24`
and energy histogram

```text
-24:1, -20:4, -16:8, -12:16, -8:31, -4:44, 0:48,
  4:44, 8:31, 12:16, 16:8, 20:4, 24:1.
```

The numbers of 256 projective states within absolute gaps
`0,4,8,12,16,20` of the cap are respectively

```text
2, 10, 26, 58, 120, 208.
```

The exact minimum fixing certificate has size seven (the support-cover lower
bound alone is only three).  Over all `8*384=3072` local-minimizer pairs,
`B_0=q_3+q_6=16`, and every hybrid satisfies (R14.1); the best norm is 28 as
already recorded.  The Hoeffding sum (10.549) is `17.631...>1` at `t=4`,
whereas `t=B_0` makes its index set empty.  This cleanly illustrates why the
entropy calculation is weaker than the two-sided comparison for the
asymptotic target.

Nine seeded complete-multipartite mosaics of orders 8--10 were also
enumerated exactly.  Their cap multiplicities were 1--4 projective states,
and their exact fixing numbers were between `n-s` and `n-s+2`.  This is
consistent with strong finite anti-concentration, but supplies no asymptotic
entropy theorem.  All finite claims are checked by
`tmp/check_hybrid_entropy_r14.py`.

Primary sources:

- Matthew Kwan and Lisa Sauermann, [Resolution of the quadratic
  Littlewood--Offord problem](https://arxiv.org/abs/2312.13826).
- Matthew Kwan and Lisa Sauermann, [An algebraic inverse theorem for the
  quadratic Littlewood--Offord problem, and an application to Ramsey
  graphs](https://arxiv.org/abs/1909.02089).
