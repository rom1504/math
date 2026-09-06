# Wave 36 memo: agreement does not by itself transfer the row bound

## Scope and conclusion

I tested the following proposed implication.  On an anchored family
`G` of `m`-selectors, choose favorable projective labels, orient them at the
anchor, and suppose their weighted conflict is `O(d)`.  Does the plurality
word, or a word at sublinear Hamming distance from it, automatically satisfy

```math
R_2(z)=\lVert Az\rVert_2^2=O(n^{9/4-c'})?
```

The generic answer is **no**: agreement controls the Hamming perturbation
cost at more than enough precision, but it controls neither (a) the row
second moment of the favorable labels nor (b) the degree-three selector
statistics needed to pass from a rare anchored family to the uniform slice.
Both gaps remain even after using all elementary ground-state stability.

I do **not** have an actual-minimizer counterexample.  I isolate below an
exact sufficient transfer lemma and two explicit, deliberately scoped
nonminimal mechanism walls.  Thus the automatic transfer statement for
exact minimizers remains an open minimizer-specific theorem; it is not a
consequence of conflict, entropy, Hamming perturbation, or generic local
stability.

Throughout, `A` is symmetric, has zero diagonal and off-diagonal entries in
`{+1,-1}`, `z` is a full spin, and

```math
L_S(z):=\lVert A[S]z_S\rVert_2^2,
\qquad R_2(z):=\lVert Az\rVert_2^2.
```

## 1. Exact uniform and anchored identities

Put `p_k=(m)_k/(n)_k`.  Expanding the three distinct indices in
`L_S(z)` gives the exact identity

```math
\boxed{
\mathbb E_{S\sim U_m}L_S(z)
=p_3R_2(z)+(p_2-p_3)n(n-1).
}
\tag{A36.1}
```

Indeed, the terms with the two inner indices equal contribute
`p_2 n(n-1)`, while every term with three distinct indices has inclusion
probability `p_3`; the corresponding full sum is `R_2(z)-n(n-1)`.

There is an equally explicit singleton-anchored version.  Gauge by `z`:

```math
W=D_zAD_z,
\qquad r=W\mathbf1,
\qquad
c_{ijk}=w_{ij}w_{ik}+w_{ij}w_{jk}+w_{ik}w_{jk}.
```

For a fixed anchor `v`, define

```math
C_v=\sum_{\substack{j<k\\j,k\ne v}}c_{vjk},
\qquad
\alpha_2=\frac{(m-1)_2}{(n-1)_2},
\qquad
\alpha_3=\frac{(m-1)_3}{(n-1)_3}.
```

Then

```math
\boxed{
\mathbb E[L_S(z)\mid v\in S]
=m(m-1)+\alpha_3\{R_2(z)-n(n-1)\}
+2(\alpha_2-\alpha_3)C_v.
}
\tag{A36.2}
```

The anchor correction itself has the exact row formula

```math
\boxed{
C_v=\frac12r_v^2+(Wr)_v-\frac32(n-1)
=\frac12r_v^2+z_v(A^2z)_v-\frac32(n-1).
}
\tag{A36.3}
```

Also `|C_v|<=3 binom(n-1,2)`.  Consequently, uniformly when
`m/n>=rho>0`,

```math
\mathbb E[L_S(z)\mid v\in S]=\alpha_3R_2(z)+O(n^2),
\qquad \alpha_3\ge c_\rho>0.
\tag{A36.4}
```

Thus a target-scale bound on the **uniform anchored** local-row mean would
immediately give the desired full row bound.  The issue is that agreement
only bounds a conditional mean on `G`, not this uniform anchored mean.

For reference, under an arbitrary law `P` on the `m`-slice one has

```math
\mathbb E_PL_S(z)
=m(m-1)+2\sum_{i<j<k}P\{i,j,k\in S\},c_{ijk}(z).
\tag{A36.5}
```

So the exact missing selector statistic is degree-three representativeness
against the signed coefficient tensor `c(z)`.

## 2. Conflict makes the Hamming perturbation affordable

Let `y^S` be the canonically oriented favorable label and let `z` be its
coordinatewise plurality.  Write

```math
e_S=d_H(z_S,y^S).
```

The anchored decoder in (10.1002) says

```math
\mathbb E_{S\sim P}e_S\le \mathcal C_{v_*}(\mathbf y),
\qquad P=U(\mathcal G).
\tag{A36.6}
```

Pointwise, principal compression and the triangle inequality give

```math
\boxed{
\left|\sqrt{L_S(z)}-\sqrt{L_S(y^S)}\right|
\le 2\lVert A\rVert_{\rm op}\sqrt{e_S}.
}
\tag{A36.7}
```

Taking the `L^2(P)` norm yields the useful averaged form

```math
\boxed{
\sqrt{\mathbb E_PL_S(z)}
\le
\sqrt{\mathbb E_PL_S(y^S)}
+2\lVert A\rVert_{\rm op}
\sqrt{\mathcal C_{v_*}(\mathbf y)}.
}
\tag{A36.8}
```

At the Wave 34 scales,

```math
d=O\!\left(\frac{Tn^{1/2-2c_0}}{(\log n)^2}\right),
\qquad T\le n^\eta,
\qquad c'=c_0-\eta,
```

and exact minimality gives `||A||_op^2<=2q_n=O(n^(3/2))`.  Therefore

```math
\lVert A\rVert_{\rm op}^2(d+1)
=O\!\left(n^{3/2}
+\frac{Tn^{2-2c_0}}{(\log n)^2}\right)
=o(n^{9/4-c'}).
\tag{A36.9}
```

So if the selected favorable labels already had conditional row mean
`O(n^(9/4-c'))`, conflict `O(d)` would preserve it.  The Hamming error is
not the obstruction.

The same estimate shows that changing the plurality in only `O(d)` global
coordinates cannot be expected to repair a genuinely high-row plurality:
for any `z'` at Hamming distance `k`,

```math
\sqrt{R_2(z')}
\ge \sqrt{R_2(z)}-2\lVert A\rVert_{\rm op}\sqrt{k}.
\tag{A36.10}
```

## 3. Generic favorable/ground stability is off scale

For one label choose its favorable orientation and put

```math
Q_S=Q(A[S]),
\qquad
\sigma (y^S)^TA[S]y^S=Q_S-\delta_S.
```

In its gauge, let `r_i^S=\sigma y_i^S(A[S]y^S)_i`.  If `F subset S` and
`a_S(F)` is the gauged boundary weight, flipping `F` gives exactly

```math
a_S(F)\ge-\delta_S/4,
\qquad r_i^S\ge-\delta_S/4.
\tag{A36.11}
```

For an exact child ground (`delta_S=0`) all row fields are nonnegative and

```math
\sum_{i\in S}r_i^S=Q_S,
\qquad
L_S(y^S)=\sum_i(r_i^S)^2
\le(m-1)Q_S.
\tag{A36.12}
```

Every principal cap satisfies `Q_S<=q_n`: extend the child ground by
independent uniform outside spins and average its oriented full energy.
Thus even for exact child grounds the generic conclusion is only

```math
L_S(y^S)\le(m-1)q_n=O(n^{5/2}),
\tag{A36.13}
```

whereas the needed scale is `n^(9/4-c')`.  It misses by
`n^(1/4+c')`.  For a merely favorable label at a fixed density bounded
away from one, the allowed deficit `B_(n,m)+t=Theta(n^(3/2))` makes the
singleton lower bound in (A36.11) asymptotically vacuous.  The universal spectral bound
`L_S(y)<=m||A||_op^2=O(n^(5/2))` gives the same wall.

Therefore a proof needs a genuinely new **row-regular favorable-label
selection** theorem, not another use of singleton or block stability.

## 4. A precise sufficient row-transfer package

Let `U_v` be the uniform anchored slice, `P=U(\mathcal G)`, and let `z` be
a plurality median of the chosen labels.  The following two estimates are a
clean sufficient package:

```math
\boxed{
\begin{aligned}
\mathbb E_P L_S(y^S)&=O(n^{9/4-c'}),
&&\text{(favorable-label row regularity)},\\
\mathbb E_{U_v}L_S(z)&\le
\mathbb E_PL_S(z)+O(n^{9/4-c'}).
&&\text{(one-sided degree-three representativeness).}
\end{aligned}
}
\tag{A36.14}
```

Together with conflict `O(d)`, (A36.8)--(A36.9) bound the right side in the
second line by `O(n^(9/4-c'))`; (A36.2)--(A36.4) then give
`R_2(z)=O(n^(9/4-c'))`.

Equivalently, the exact single missing lemma can be stated without auxiliary
conditions:

> **Anchored favorable row-transfer lemma.**  Uniformly for every relevant
> exact minimizer and target pair, if an anchored family has
> `log beta^(-1)=O(TL_0)` and admits favorable labels of conflict `O(d)`,
> then some coordinatewise Hamming median of those labels satisfies
> `E_(U_v)L_S(z)=O(n^(9/4-c'))`.

By (A36.2), this is precisely the desired weaker arbitrary-cut row clause.
The decomposition (A36.14) says what a proof of that lemma must manufacture.
Neither line follows from the present agreement hypothesis.

An actual-minimizer falsifier would be an unbounded family satisfying the
agreement hypotheses for which **every** Hamming median (including all
choices on ties and unseen coordinates) has
`R_2(z)=omega(n^(9/4-c'))`.  No such family is known.

## 5. Entropy alone does not give degree-three representativeness

Here is an explicit scoped signing example.  Take a symmetric conference
core `B` of order `ell`, so `B^2=(ell-1)I`.  Adjoin a hub set `H` of size
`h=ceil(ell^(1/4))` and put `+1` on every edge incident with a hub; call the
resulting signing `A`, of order `n=ell+h`.  Let `z=1` and fix an anchor
`v` in the core.

The hub rows give

```math
R_2(z)\ge h(n-1)^2=\Theta(n^{9/4}).
\tag{A36.15}
```

Let `G` be all anchored `m`-sets which avoid every hub.  At fixed density
`m/n->rho in (0,1)`, its relative density is

```math
\beta=
\frac{\binom{\ell-1}{m-1}}{\binom{n-1}{m-1}},
\qquad
\log\beta^{-1}
=h\log\frac1{1-\rho}+O(h^2/n)
=\Theta(n^{1/4}).
\tag{A36.16}
```

This is much smaller than the allowed `O(TL_0)` entropy loss.  Yet every
`S in G` lies inside the conference core, so

```math
L_S(z)
\le \lVert B[S]\rVert_{\rm op}^2m
\le(\ell-1)m
=O(n^2).
\tag{A36.17}
```

In contrast, (A36.2) says that its mean on the full anchored slice is
`Omega(hn^2)`.  Hence a subproject entropy deficit can condition directly
onto a very atypical value of the cubic statistic.  This example does not
claim that `z_S` is favorable and does not use an exact minimizer; it
falsifies only the proposed inference from `log beta^(-1)` to
representativeness.  The signing is still competitive in the coarse sense:
the conference spectral bound and the positive-hub cap give
`Q(A)<=ell sqrt(ell-1)+O(hn)=O(n^(3/2))`.

## 6. A nonminimal hub join makes the stability wall sharp

The following construction shows that optimal-order cap, exact ground
stability, perfect agreement, and sublinear Hamming repair still do not
control rows without exact signing minimality.

Choose any signing `B` of order `ell` with `Q(B)<=C ell^(3/2)`.  Gauge and,
if needed, negate it so that `1` is a positive absolute ground.  Thus every
gauged cut weight lies in `[0,Q(B)/2]`.  Adjoin
`h=ceil(K sqrt(ell))` positive hubs as above.  Then `1` remains an exact
absolute ground and, exactly,

```math
Q(A)=Q(B)+2h\ell+h(h-1)=O_K(n^{3/2}),
\qquad
R_2(1)\ge h(n-1)^2=\Theta_K(n^{5/2}).
\tag{A36.18}
```

For `G={S:H subset S}` and an anchor in `H`, choose `y^S=1_S`.  The labels
have zero conflict and

```math
\log\beta^{-1}=\Theta(h)=\Theta(\sqrt n)=O(L_0).
\tag{A36.19}
```

For `S=H union R`, the positive-hub part has `1` as an absolute ground and
`Q(B[R])<=Q(B)`, whence

```math
Q(A[S])-\mathbf1^TA[S]\mathbf1\le2Q(B).
\tag{A36.20}
```

If the favorable allowance is formed with this signing's own parent cap,
`
[(m/n)^(3/2)-p_2]Q(A)`, choosing the fixed constant `K` sufficiently large
makes every displayed label favorable with `t=0`.  The plurality is uniquely
`1`.  Moreover `||A||_op=O_K(n^(3/4))`, so (A36.10) shows that every word at
Hamming distance `o(n)` from it still has row square `Omega_K(n^(5/2))`.

This is **not** an actual counterexample: the joined signing is not shown to
minimize `Q` at order `n`, and the actual project notation `q_n` cannot be
replaced by its larger cap.  It is a mechanism wall showing that any proof
using only the cap, generic ground/favorable stability, agreement, and
Hamming perturbation cannot work.  Exact edge-sign minimality must force the
new row-transfer property if that property is true.

## 7. Frontier from this attack

1. The exact identities (A36.1)--(A36.4) make the weaker arbitrary-cut row
   clause equivalent to a uniform anchored local-row mean.
2. Conflict at the required `d` scale is cheap enough; it is not the row
   bottleneck.
3. Generic exact-ground or favorable-label stability stops at `n^(5/2)`.
4. `log beta^(-1)=O(TL_0)` does not make a rare anchored family
   degree-three representative; the conference-core example gives an exact
   obstruction.
5. The viable new theorem is (A36.14), or the combined anchored favorable
   row-transfer lemma.  It must use exact minimizer structure to select
   row-regular favorable labels and/or prevent the selector family from
   hiding the cubic row statistic.

The checker `tmp/agreement_row_r36_check.py` independently verifies the
uniform and anchored identities, (A36.3), the exact-ground row bound, the
finite hub-join ground claim, and the Hamming perturbation inequality.
