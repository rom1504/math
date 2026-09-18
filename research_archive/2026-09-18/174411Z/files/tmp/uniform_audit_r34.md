# Independent Wave 34 audit: uniform coset and rare-center reduction

## Verdict

The uniform-selector weakening is correct.  The max-center hard and soft
moment sandwiches are also correct.  The only material correction is to avoid
calling the hard and strict-margin soft criteria fully bidirectionally
equivalent: the strict-margin soft criterion implies a hard high-degree
center, while a hard high-degree center proves the star event directly but
need not recover the same strict soft margin at the prescribed `lambda`.

This is a genuine weakening of the Wave 33 all-selector-law/common-weight
implementation, but the one-coset part is mostly a recombination of the
already verified fixed-cut lemma (10.795) with the block-coset support.  It is
not a new signing estimate.  The exceptional-center sandwich itself appears
new in the ledger and usefully removes any need for center abundance or a
Jensen first moment.

## Uniform one-coset reduction

For each eligible coset,

```math
H_a=\bigcup_{d\in C_a}E_d(t)
```

follows exactly from (10.806), (10.831), and (10.832), including orientation.
Moreover `|C_a|<=2^k` (the root memo's `2^(k+1)` is harmless).  Hence

```math
\max_{d\in C_a}U_m(E_d(t))\ge U_m(H_a)/|C_a|.
```

At `k=O(L0/log n)`, the support loss is `O(k)=o(L0)`.  If
`U_m(H_a)>=exp{-O(TL0)}` for fixed `T<=n^eta`, `eta<c0`, then the selected
cut has coverage `exp{-O(TL0)}`.  Put `c'=c0-eta>0`; then
`TL0<=n^(3/4-c')`, while the row cap `O(n^(9/4-c0))` and tolerance
`O(n^(3/2-c0))` are stronger than the caps required at `c'`.  Thus all three
conditions of (10.795) hold.  No new landing argument is needed: exact
landing and geometric-window summability are inherited from the verified
sufficiency of (10.795).  The hypothesis must, as always, hold uniformly for
every target `(n,m)` in the fixed active ratio window.

Conversely, (10.853)--(10.854) plant any cut from (10.795) into a balanced
eligible coset, whose hit set contains the cut's favorable selector set.
Thus the `T=1` uniform one-coset formulation is exponent-equivalent to
(10.795), up to the subleading support factor.

The collision conversion also has the right direction:

```math
\Pr(J_{r,s})\le {r\choose s}M h_U^s,
\qquad h_U=\max_aU_m(H_a).
```

Here `log M=O(n log k)=O(rL0)`, `r/s<=T`, and
`s^(-1)log binom(r,s)<=log(eT)`.  Therefore a lower bound
`Pr(J)>=exp{-O(rL0)}` forces `h_U>=exp{-O(TL0)}`.  Hence it is sufficient to
prove (10.926) only for `w=U_m`; the adversarial infimum in (10.907) is not
needed for this implication.

## Exceptional-center exponents

For either a hard degree `p(z)` or soft degree `q(z)`, uniformity of `nu_2`
gives, with `N=|C_2|<=2^n`,

```math
N^{-1}u_*^s\le E_{nu_2}u(z)^s\le u_*^s.
```

The scale identities are

```math
rL0=n\log(2k0)+O(L0),\qquad sT/r\longrightarrow1,
\qquad n/(rL0)\longrightarrow0,\qquad n/s=o(TL0).
```

Thus a uniform hard moment of size `exp{-O(rL0)}` is equivalent to one hard
center of degree `exp{-O(TL0)}`, and the analogous algebraic statement holds
for the soft moment and `q_*`.  One center's atom is affordable because it
costs only `exp{-O(n)}=exp{-o(rL0)}`.

For the Laplace lower bound (10.947), however, unspecified `O` constants are
not enough.  At

```math
lambda=Lambda rL0/(D+1)
```

one needs a fixed strict margin.  If

```math
q_*\ge e^{-A TL0},\qquad A<Lambda,
```

then

```math
L_lambda(U_m)\ge N^{-1}e^{-A sT L0},
```

whose exponent coefficient tends to `A`, whereas the subtraction is
`e^{-Lambda rL0}`.  The `N^{-1}` term contributes `o(rL0)` and hence does
not consume the fixed gap.

With `d=floor(D/s)`, one has `s/D=o(1)` and therefore

```math
lambda(d+1)=(Lambda+o(1))TL0.
```

The one-column inequality then turns the strict soft condition into

```math
g_d(z_*)\ge e^{-(A+o(1))TL0}.
```

This hard degree directly yields the center-star event.  The reverse
implication to a strict-margin soft certificate at the same `lambda` is not
proved: the elementary lower bound from a hard degree pays an additional
`e^{-lambda d}`, i.e. approximately another `Lambda TL0`.  Hence the safe
wording is: **the strict-margin soft route reduces to a hard exceptional
center, and the latter is already sufficient directly**.

## Remaining identities

The distance--deficit inequality (B34.13) is correct because each child-spin
flip changes the quadratic energy by at most `4(m-1)`, and absolute value is
1-Lipschitz.  Exponentiation gives (B34.14); (B34.15) follows because a center
whose absolute deficit is within allowance is itself a favorable projective
label.  The scale `lambda/[4(m-1)]` is
`n^(-3/4+c0)` up to logarithmic factors.

The parent identities (B34.16) are exact:

```math
E_S delta_S(d)=\bar Q-p_2\sigma z^TAz
=B+G+p_2Delta_A(d),\qquad E_S\widehat\ell(S,d)=G.
```

They remain valid for selector-independent parent mixtures.  They do not
apply to an `S`-dependent conditional state law and do not rule out an
exceptional nonlinear selector tail.

The finite checker `tmp/common_capacity_r34_check.py` passes on the synthetic
one-center extremizer and on `A_6,A_8,A_9`.

Minor editorial issue: `tmp/common_capacity_r34.md` refers to root `R34.6`
as the average-center Jensen target, but root `R34.6` is the max hard-center
condition; the Jensen target is `R34.9`/`B34.7`.
