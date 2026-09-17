# Wave 34 Route B: the uniform soft moment is a rare-center problem

## Status

The max-center sandwich, the one-column hard/soft inequalities, the
distance--deficit implication, and the parent-slack identities below are
**proved**.  They are checked on the exact minimizers `A_6,A_8,A_9` by
`tmp/common_capacity_r34_check.py`.  I also independently audited the
uniform-selector reduction in `tmp/uniform_partial_coset_r34.md`; the audit is
recorded in Section 1.

No signing estimate, convergence proof, or asymptotic counterexample is
obtained.  The useful correction is that the exact uniform soft target does
not collapse to the Jensen first moment.  It can be driven by one exceptional
low-row center.  The common all-selector weight (10.950) and the average-center
Jensen target (R34.9 in the updated root memo) are both stronger than this
rare-center target.

## 1. Audit of the uniform-selector weakening

For an eligible row-good coset `a`, let `C_a` be its oriented support and put

```math
E_d(t)=\{S:\widehat\ell(S,d)\le t\}.
```

Equations (10.806), (10.831), and (10.832) give exactly

```math
H_a
=\{S:\max_{d\in C_a}X_d(S)\ge Y_A(S)-t\}
=\bigcup_{d\in C_a}E_d(t).                         \tag{B34.1}
```

Hence

```math
\max_{d\in C_a}U_m(E_d(t))
\ge \frac{U_m(H_a)}{|C_a|}.                         \tag{B34.2}
```

In fact (10.831) gives `|C_a|<=2^k`; the looser `2^(k+1)` in the root memo
is harmless.  At `k=O(L_0/log n)`, this costs `o(L_0)` in the logarithm.
Thus one coset with uniform hit mass `e^{-O(TL_0)}`, `T<=n^eta`, produces a
cut satisfying (10.795) at the weaker saving `c_0-eta`.

Conversely, (10.853)--(10.854) put any cut satisfying (10.795) into an
eligible balanced row-good coset, and (B34.1) then gives `H_a superset
E_d(t)`.  Thus the uniform one-coset target with `T=1` and the arbitrary-cut
lemma are exponent-equivalent up to the harmless codebook factor.  The
adversarial quantity `inf_w max_a w(H_a)` is sufficient but is not necessary.

The use of (10.904) is also in the correct direction.  For uniform selectors,

```math
\Pr(J_{r,s})
\le {r\choose s}\sum_aU_m(H_a)^s
\le {r\choose s}|\mathscr A_C|\max_aU_m(H_a)^s.     \tag{B34.3}
```

Therefore an `e^{-O(rL_0)}` *lower* bound on the collision event forces one
coset to have uniform mass `e^{-O(TL_0)}`.  It is enough to prove the
center-star or soft event only at `w=U_m`.

## 2. Exact max-center sandwich

Write `nu=nu_2`, `N=|C_2|`, and

```math
q_\lambda(z)=\mathbb E_{S\sim U_m}e^{-\lambda a_z(S)},
\qquad
L_\lambda(U_m)=\mathbb E_{z\sim\nu}q_\lambda(z)^s,
\qquad
q_*=\max_{z\in C_2}q_\lambda(z).
```

Because `nu_2` is uniform, one maximizing summand and the pointwise upper
bound give

```math
\boxed{
N^{-1}q_*^s\le L_\lambda(U_m)\le q_*^s,
\qquad
N^{-1/s}q_*\le L_\lambda(U_m)^{1/s}\le q_*.
}                                                       \tag{B34.4}
```

Here `N<=2^n`, while

```math
rL_0=n\log(2k_0)+O(L_0),
\qquad
\frac{n}{s}=o(TL_0),
\qquad
sT=(1+o(1))r.                                         \tag{B34.5}
```

Consequently the root-mass loss `N^{-1}` is `e^{-o(rL_0)}`, and its loss
after taking the `s`th root is `e^{-o(TL_0)}`.  At the project exponent, the
exact uniform soft moment is therefore equivalent to finding **one** low-row
center with

```math
q_\lambda(z)\ge e^{-O(TL_0)}.                         \tag{B34.6}
```

This is strictly weaker than the Jensen sufficient condition

```math
Z_\lambda:=\mathbb E_{z\sim\nu,S\sim U_m}
e^{-\lambda a_z(S)}\ge e^{-O(TL_0)}.                 \tag{B34.7}
```

Indeed (10.949) reads

```math
L_\lambda(U_m)
=Z_\lambda^s
 \exp\{(s-1)D_s(P_\lambda\Vert\nu)\},               \tag{B34.8}
```

where `P_lambda` has density `q_lambda/Z_lambda`.  Since
`D_s(P_lambda||nu)<=log N`, an exceptional center can recover up to
`N^(s-1)` relative to `Z_lambda^s`.  If `q_lambda` is supported at one
center, equality holds: Jensen pays the rarity of that center `s` times,
whereas the exact moment pays it once.  Thus R34.6 is a valid sufficient
lemma, but “the soft target collapses to a first moment” is too strong.

## 3. A strict-margin soft proof yields one hard high-degree center

For an integer `d>=0`, put

```math
g_d(z)=U_m\{S:a_z(S)\le d\}.
```

For every center separately, integrality of `a_z(S)` gives

```math
\boxed{
e^{-\lambda d}g_d(z)\le q_\lambda(z)
\le g_d(z)+e^{-\lambda(d+1)}[1-g_d(z)],
}
                                                               \tag{B34.9}
```

and hence

```math
\boxed{
g_d(z)\ge
\max\left\{0,
\frac{q_\lambda(z)-e^{-\lambda(d+1)}}
     {1-e^{-\lambda(d+1)}}\right\}.
}                                                               \tag{B34.10}
```

Take `d=floor(D/s)` and
`lambda=Lambda rL_0/(D+1)`.  At the project scales `d->infinity`, and

```math
\lambda(d+1)=(\Lambda+o(1))TL_0.                  \tag{B34.11}
```

Suppose, for example, that

```math
q_*\ge e^{-ATL_0},\qquad 0<A<\Lambda.
```

Equations (B34.10)--(B34.11) produce one `z_* in C_2` with
`g_d(z_*)>=e^{-(A+o(1))TL_0}`.  Conversely, such a hard degree gives the
common-star event directly:

```math
\mathbb E_{z\sim\nu}g_d(z)^s
\ge N^{-1}g_d(z_*)^s
=e^{-(A+o(1))rL_0}.                                  \tag{B34.12}
```

On this event every one of the `s` distances is at most `d`, so their sum
is at most `D`.  Thus it proves the uniform center-star event, (B34.3), one
uniformly large coset, (10.795), and convergence.

More generally, if `L_lambda(U_m)>=e^{-A rL_0}` with a fixed strict margin
`A<Lambda`, (B34.4) gives (B34.6) with constant `A+o(1)`, and (B34.10)
again extracts the hard center.  Therefore the exact strict-margin soft
implementation **reduces to** a max-center hard-degree theorem at exponent
scale; the hard theorem already has the same convergence consequence
directly.  The converse implication to the same fixed soft margin is not
asserted.  The average multiplicity
wall R34.8--R34.9 applies to the stronger Jensen route (B34.7), not to this
rare-center mechanism.

## 4. Exact distance--deficit obstruction

For an unoriented projective center define

```math
\delta_S^{\rm abs}(z)
=Q(A[S])-|z_S^{\mathsf T}A[S]z_S|.
```

One vertex flip changes a child quadratic energy by at most `4(m-1)`.
If `a_z(S)=a`, choose a favorable label `y` at that distance and follow a
shortest projective flip path.  Since
`delta_S^abs(y)<=B_(n,m)+t`, the triangle inequality gives

```math
\boxed{
a_z(S)\ge
\frac{[\delta_S^{\rm abs}(z)-B_{n,m}-t]_+}{4(m-1)}.
}                                                               \tag{B34.13}
```

Consequently every prospective rare center must satisfy

```math
\boxed{
q_\lambda(z)
\le\mathbb E_{S\sim U_m}
\exp\left\{-\frac{\lambda}{4(m-1)}
[\delta_S^{\rm abs}(z)-B_{n,m}-t]_+\right\}.
}                                                               \tag{B34.14}
```

At the project scale, `lambda/[4(m-1)]` is of order
`n^(-3/4+c_0)` up to logarithms.  Thus failure of the corresponding
positive-part selector Laplace tail for every `z in C_2` is a sharp
obstruction to the uniform soft-center route.  Conversely,

```math
q_\lambda(z)\ge
U_m\{S:\delta_S^{\rm abs}(z)\le B_{n,m}+t\},         \tag{B34.15}
```

because such a center is itself favorable and has distance zero.  The
remaining gap is exactly a nonlinear fixed-center upper tail, not a uniform
degree statement over all selectors.

## 5. Parent slack cannot improve the uniform first moment

Let `d=(sigma,z)` be any oriented parent state and put

```math
\Delta_A(d)=q_n-\sigma z^{\mathsf T}Az,
\qquad
\delta_S(d)=Q(A[S])-\sigma z_S^{\mathsf T}A[S]z_S,
```

and let

```math
\overline Q=\mathbb E_{S\sim U_m}Q(A[S]),
\qquad
G=\overline Q-p^{3/2}q_n.
```

Uniform pair inclusion gives the exact identities

```math
\boxed{
\mathbb E_{S\sim U_m}\delta_S(d)
=B_{n,m}+G+p_2\Delta_A(d),
\qquad
\mathbb E_{S\sim U_m}\widehat\ell(S,d)=G.
}                                                               \tag{B34.16}
```

They remain true after averaging under **any** selector-independent law on
parent states, including a law reweighted by parent slack or replacement
surplus.  Such reweighting cannot lower the uniform mean effective loss;
parent slack actually raises the mean signed child deficit.  This recovers
the circularity in (10.824) in the precise setting of the proposed weight
construction.

Equation (B34.16) does not obstruct the route in Sections 2--4: a nonlinear
tail can be large at one exceptional low-row center even when every
selector-independent parent mixture has the same unfavorable mean.  It does
show that a successful parent-slack construction must prove the nonlinear
positive-part Laplace tail (B34.14), or an equivalent rare-center statement;
no first-moment use of (10.945) can establish it.

## 6. Resulting frontier for this route

The common-weight certificate (10.950) should no longer be treated as the
sharp scalar target after the uniform-selector reduction.  The logically
weakest center implementation now visible is:

> Find one `z in C_2` for which a uniform `e^{-O(TL_0)}` fraction of
> selectors has a favorable completion within projective distance
> `D/s`, for some `T<=n^eta`, `eta<c_0`.

This theorem implies convergence by (B34.12) and the audited uniform-coset
reduction.  The exact soft version reduces to this hard theorem at exponent
scale when it has a strict margin over the subtraction.  A proof may exploit rare-center
Rényi amplification and need not establish near-maximal *average* favorable
multiplicity.  The sharp obstruction is failure of (B34.14) for every
low-row center.  Parent-slack and replacement laws remain potentially useful
only if they create this exceptional nonlinear tail; their uniform first
moments are exactly circular by (B34.16).
