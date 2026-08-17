# Independent audit: the sharp orientation-visibility threshold

**Verdict:** PASS.  The internal-cap ceiling, biased exact-sign bridge
existence, both outer-channel envelopes, bounded-parent estimate, and claimed
`n^(3/4)` sharp exponent are correct.  No mathematical repair is needed.
The hypotheses that the two landscapes in OV.1 are even are harmless but
stronger than necessary: the proof of OV.2 only uses the global inversion of
the newly appended spins after temporarily deleting `K`.

I independently reconstructed the argument and ran
`experiments/verify_orientation_visibility_threshold.py`.  Its scalar
envelopes pass; exact enumeration gives `(Q_+,Q_-)=(11,9)` at `(n,m)=(4,2)`
and `(126,96)` at `(16,8)`.

## 1. The continuation-cap law

Delete `K` and denote the resulting caps by `R_sigma^0`.  Then

```math
\begin{aligned}
R_-^0
&=\max_{x,y}|-H(x)+x^TBy|\\
&=\max_{x,y}|-H(x)-x^TBy|\\
&=\max_{x,y}|H(x)+x^TBy|=R_+^0.
\end{aligned}
```

The first equality after the definition is the bijection `y -> -y`; the
second is multiplication inside the absolute value by `-1`.  Adding `K`
changes each cap by at most `||K||_infinity`, hence

```math
|R_+-R_-|<=2||K||_infinity.
```

For a hollow `m`-vertex signing, `||K||_infinity=Q(C)` and every energy is a
sum of `binom(m,2)` unit terms.  Thus

```math
|R_+-R_-|<=2Q(C)<=m(m-1).
```

If `Q(C)<=K_0m^alpha`, a gap `epsilon n^(3/2)` therefore forces

```math
m^alpha>=epsilon n^(3/2)/(2K_0),
```

which is exactly OV.4.  In particular the exponents are `3/4` for an
unrestricted exact-sign internal block and `1` when the appended block is
itself constrained to natural `O(m^(3/2))` cap.

## 2. Existence of the biased flat exact-sign bridge

Put `a=m/n` and choose independent exact signs with

```math
Pr(B_ij=1)=(1+a)/2,
\qquad E=B-aJ.
```

For fixed Euclidean unit vectors `u,v`, `u^TEv` is a sum of independent,
centred, uniformly bounded variables with squared coefficient sum

```math
\sum_{i,j}u_i^2v_j^2=1.
```

Hoeffding therefore gives a dimension-free subgaussian tail.  Taking
`1/4`-nets of the two unit spheres and a union bound at
`t=C sqrt(n+m)` proves, with positive probability,

```math
||E||_(2->2)<=C'(sqrt n+sqrt m)<=2C'sqrt n.
```

This proves existence of an **exact** sign matrix; no rounding of `aJ` is
being asserted.  Uniformly on Boolean inputs,

```math
|x^TEy|<=||E||sqrt(nm)=O(nsqrt m).
```

At `m=floor(n^(3/4))`, this is `O(n^(11/8))=o(n^(3/2))`.

## 3. Old-block normalization and exact signs

On the orders `n=4^j`, the regularized Walsh construction supplies a
symmetric sign matrix `mathcal H` with

```math
mathcal H^2=nI,
\qquad mathcal H1=sqrt(n)1,
\qquad tr(mathcal H)=0.
```

Deleting its diagonal gives a hollow complete signing `A`.  Since Boolean
coordinates square to one and the trace is zero,

```math
H_A(x)=x^Tmathcal Hx/2
```

exactly.  Both orientations in OV.9 are therefore complete hollow signings:
their old blocks are `+-A`, the cross block is the exact sign matrix `B`, and
their new block is the positive clique.

Let

```math
S=n^(3/2),
\qquad lambda=m^2/S,
\qquad p=(1^Tx)/n,
\qquad s=(1^Ty)/m.
```

Here `lambda<=1` and `lambda -> 1`.

## 4. Positive-orientation lower bound

At `x=y=1`, regularity, the biased bridge decomposition, and the clique
identity give

```math
P_+(1,1)
=S/2+m^2+(m^2-m)/2+1^TE1.
```

Writing `D=||E||sqrt(nm)=O(n^(11/8))`, this yields the explicit bound

```math
Q(P_+)>=S/2+3m^2/2-m/2-D=(2-o(1))S.
```

The main term is eventually positive, so passing to the absolute cap causes
no sign issue.

## 5. Both negative-orientation outer channels

Write `x=p1+x_perp`.  Symmetry and regularity make the decomposition
orthogonal for `mathcal H`, while `||mathcal H||=sqrt n`.  Consequently

```math
x^Tmathcal Hx
>=S(2p^2-1),
\qquad
x^Tmathcal Hx<=S.
```

The first inequality gives

```math
-x^Tmathcal Hx/2<=S(1/2-p^2).
```

The mean bridge and clique terms are exactly

```math
a(1^Tx)(1^Ty)=m^2ps=lambda Sps,
\qquad
H_C(y)=lambda Ss^2/2-m/2.
```

For the positive outer channel of `P_-`, discard the favourable `-m/2` and
optimize first in `p`:

```math
\begin{aligned}
P_-/S
&<=1/2-p^2+lambda ps+(lambda/2)s^2+D/S\\
&<=1/2+(lambda^2/4+lambda/2)s^2+D/S\\
&<=5/4+D/S.
\end{aligned}
```

The optimizer `p=lambda s/2` lies in `[-1,1]`, so no boundary case was
missed.  For the negative outer channel, use the second old-block bound:

```math
\begin{aligned}
-P_-/S
&<=1/2+lambda(-ps-s^2/2)+m/(2S)+D/S\\
&<=1/2+lambda(|s|-s^2/2)+m/(2S)+D/S\\
&<=1+m/(2S)+D/S.
\end{aligned}
```

Thus the absolute maximum, including both signs, satisfies

```math
Q(P_-)<=5S/4+D+m/2=(5/4+o(1))S.
```

Combining the two explicit estimates even gives

```math
Q(P_+)-Q(P_-)
>=3m^2/2-3S/4-m-2D
=(3/4-o(1))S.
```

## 6. Parent cap and sharpness scope

Every channel is uniformly bounded as follows:

```math
Q(A)<=S/2,
\qquad
\max_{x,y}|a(1^Tx)(1^Ty)|<=m^2,
\qquad
\max_{x,y}|x^TEy|<=D=o(S),
\qquad
Q(C_m)={m\choose2}<=m^2/2.
```

Triangle inequality gives `Q(P_+),Q(P_-) = O(S)`.  Since
`N=n+m=(1+o(1))n`, this is also `O(N^(3/2))`.

OV.1 forbids a fixed `Theta(n^(3/2))` orientation gap when
`m=o(n^(3/4))`; OV.3 supplies such a gap at
`m=(1+o(1))n^(3/4)`.  The exponent is therefore sharp for arbitrary
exact-sign quadratic continuations when only the **final** parent cap is
required to stay natural.  The example does not contradict the stronger
linear-width barrier for a continuation whose own internal cap is
`O(m^(3/2))`: its clique deliberately spends `Theta(m^2)` cap.

Finally, this construction is not a universal pin.  It does not force one
old configuration to optimize against every child.  It compares only the
two orientations of one regular-Hadamard child, and its proof uses a scalar
mean correlation plus a uniformly negligible flat residual.  Thus the
universal-pin no-go theorem does not subsume it.
