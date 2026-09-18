# Wave 53: low-cross child grounds reduce exactly to internal-field regularity

## Status

No child-ground box theorem for exact minimizers is proved.  The attack gives:

1. the exact random-completion expectation and centered minimization formula;
2. a sharp sufficient field-regularity lemma and its converse obstruction;
3. an audit showing that all currently available generic and
   exact-minimizer estimates stop at `O(n^(5/2))`; and
4. a scalable nonminimal construction with competitive `O(n^(3/2))` cut
   norm, zero cross energy, and `Theta(n^(5/2))` internal Gram.  Thus any
   positive theorem has to use discrete exact signing minimality beyond its
   presently extracted spectral/trace consequences.

The exact finite calculations are reproduced by
`tmp/lowcross_box_r53_check.py` and `tmp/lowcross_box_r53_check.out`.

## 1. Exact completion expectation and minimum

Let `A` be an order-`n` signing, `S` an `m`-selector, `T=S^c`, and
`k=n-m`.  Let `y` be a child ground of `A[S]`.  Orient it by `sigma` so that

```math
Q_S=Q(A[S])=\sigma y^{\mathsf T}A[S]y.
```

Put

```math
I=\lVert A[S]y\rVert_2^2,\qquad
X=\lVert A[T,S]y\rVert_2^2.
```

For uniform `w in {+-1}^T`, blockwise orthogonality gives the **Verified
exact decomposition**

```math
\begin{aligned}
\mathbb E_w\lVert A[S]y+A[S,T]w\rVert_2^2&=I+mk,\\
\mathbb E_w\lVert A[T,S]y+A[T]w\rVert_2^2&=X+k(k-1),
\end{aligned}
```

and hence

```math
\boxed{
\mathbb E_w\lVert A[:,S]y+A[:,T]w\rVert_2^2
=I+X+k(n-1).}
\tag{L53.1}
```

In particular, for the exact box quantity `mathcal V(S,y)` of (10.1237),

```math
\boxed{\mathcal V(S,y)\le I+X+k(n-1).}
\tag{L53.2}
```

There is also an exact centered formula for the minimum.  Set

```math
b=A[:,S]y,\qquad C=A[:,T],\qquad
H=A^2-(n-1)I_n,\qquad d=C^{\mathsf T}b=(A^2)[T,S]y.
```

Since `diag(C^T C)=(n-1)1`, direct expansion proves

```math
\boxed{
\mathcal V(S,y)=I+X+k(n-1)
+\min_{w\in\{\pm1\}^T}
\{2d^{\mathsf T}w+w^{\mathsf T}H[T]w\}.}
\tag{L53.3}
```

Pairing `w` and `-w`, then averaging the better member of each pair, also
gives the rigorous but generally lower-order improvement

```math
\boxed{
\mathcal V(S,y)
\le I+X+k(n-1)-2\mathbb E_w|d^{\mathsf T}w|
\le I+X+k(n-1)-\sqrt2\lVert d\rVert_2.}
\tag{L53.4}
```

The last step is the sharp `p=1` Khintchine inequality.  The exact-minimizer
operator bound only gives `||d||_2=O(n^2)`, so (L53.4) by itself cannot
remove a possible `Theta(n^(5/2))` internal term.  The quadratic minimum in
(L53.3) can nevertheless be important; it is not legitimate to identify
the true box minimum with its uniform-completion mean.

## 2. The minimal internal-field reduction

Define the oriented internal local fields

```math
r_i=\sigma y_i(A[S]y)_i,\qquad i\in S.
```

One-spin optimality of the absolute child ground gives

```math
0\le r_i\le\min\{m-1,Q_S/2\},\qquad
\sum_{i\in S}r_i=Q_S,
```

and therefore the **Verified identity and bound**

```math
\boxed{
I=\sum_{i\in S}r_i^2
\le \lVert r\rVert_\infty Q_S.}
\tag{L53.5}
```

Fix `0<c<1/4` and `R_n=n^(9/4-c)`.  Uniformly on a compact fixed-density
window, `Q_S<=q_n=O(n^(3/2))` and `k(n-1)=O(n^2)=o(R_n)`.  Equations
(L53.2) and (L53.5) prove the exact sufficient reduction

```math
\boxed{
X=O(R_n),\qquad
\lVert r\rVert_\infty=O(n^{3/4-c})
\quad\Longrightarrow\quad
\mathcal V(S,y)=O(R_n).}
\tag{L53.6}
```

Equivalently, the single condition

```math
X+Q_S\lVert r\rVert_\infty=O(R_n)
\tag{L53.7}
```

is an **Open sufficient exact-minimizer lemma** for the child-ground box
witness.

The converse identifies the obstruction without an asymptotic constant
loss.  For every `R>X+k(n-1)`, (L53.2) implies

```math
\boxed{
\mathcal V(S,y)>R
\quad\Longrightarrow\quad
\lVert r\rVert_\infty>
\frac{R-X-k(n-1)}{Q_S}.}
\tag{L53.8}
```

Thus, if the box target fails and a child ground lies in the branch
`X<=theta R_n` for a fixed `theta<1`, that ground must have a local field of
size at least

```math
(1-\theta-o(1))R_n/Q_S=\Omega(n^{3/4-c}).
```

This is the exact low-cross versus spiky-internal-field dichotomy.  It does
not produce such a low-cross child ground, nor does an existence statement
alone supply the saved population needed by the separate K-functional
route.

## 3. Why the known structures do not prove field regularity

The following implications were audited separately.

- **One-spin optimality.**  It gives only `r_i<=m-1` and hence
  `I<=(m-1)Q_S=O(n^(5/2))`.  The desired maximum-field estimate saves the
  factor `n^(1/4+c)`.
- **Exact-minimizer operator control.**  From `||A||op^2<=2q_n`,

  ```math
  I+X=\lVert A[:,S]y\rVert_2^2
  \le2q_nm=O(n^{5/2}),
  ```

  exactly the same wall.  The `tr(A^4)` bound is an average spectral budget
  and does not improve a particular optimized child ground.
- **Maximal versus arbitrary selectors.**  If `S` maximizes `Q(A[S])`, the
  port theorem (10.1271) gives `|t_ji|<=r_i`; it supplies no upper bound on a
  large `r_i`.  Its square budget is

  ```math
  kI-(m-2)X-mk\ge0,
  ```

  which has the wrong direction for (L53.6).  Moreover the box minimizer may
  use a nonmaximal selector, so the port inequalities cannot be applied to
  it without a separate selection theorem.
- **Selector averaging.**  For one fixed full word `z`, the exact identity

  ```math
  \mathbb E_{S\sim U_m}\lVert A[:,S]z_S\rVert_2^2
  =p_2R_2(z)+(p-p_2)n(n-1)
  ```

  is useful.  But the restrictions `z_S` need not be child grounds.  When a
  separately constructed low-information certificate law is available,
  (10.1087) transports this mean; it does not construct a child-ground
  incidence.  Replacing each restriction by an independently optimized
  child ground destroys the fixed-word averaging identity.
- **Averaging the fields inside one ground.**  It shows only
  `min_i r_i<=Q_S/m=O(sqrt(n))`; the box estimate needs control of the
  maximum or the squared-field sum.

Consequently no presently verified selector, trace, or minimizer identity
forces (L53.7).  The missing input would have to use exact discrete
edge-sign minimality in a new way, or exploit the full quadratic
cancellation in (L53.3) rather than the mean certificate.

## 4. A scalable low-cross spiky obstruction without exact minimality

The power loss is not an artifact of loose inequalities.  The following
standard probabilistic construction satisfies all the generic scales but is
not asserted to be an exact minimizer.

Let `m/n->p in (1/2,1)`, put `h=Theta(sqrt(n))`, and choose a core signing
`B` of order `m-h` with

```math
Q(B)=O(n^{3/2}),\qquad \lVert B\rVert_{op}=O(\sqrt n).
```

Such signings exist probabilistically.  Gauge and, if needed, negate `B` so
that `1` is a positive absolute ground.  Adjoin `h` positive hubs to form
the child block `A[S]`: all hub--hub and hub--core edges are `+1`.  Triangle
inequality shows exactly that `1_S` is an absolute child ground, because its
energy is

```math
Q(B)+2h(m-h)+h(h-1),
```

the sum of the separate absolute caps.  Every hub local field equals
`m-1`, so

```math
I\ge h(m-1)^2=\Theta(n^{5/2}),\qquad
\lVert r\rVert_\infty=m-1.
\tag{L53.9}
```

Choose the cross block `C=A[T,S]` with balanced rows and
`||C||op=O(sqrt(n))`; random balanced sign rows provide such a matrix.  Then

```math
C1_S=0,\qquad X=0.
\tag{L53.10}
```

Complete `A[T]` by any `O(n^(3/2))`-cap signing.  For every full spin, the
two diagonal block energies and Cauchy--Schwarz on the cross block give

```math
Q(A)\le Q(A[S])+Q(A[T])+2\lVert C\rVert_{op}\sqrt{mk}
=O(n^{3/2}).
```

The whole matrix also has `||A||op=O(n^(3/4))`: the positive hub--core
join has norm `Theta(sqrt(hn))=Theta(n^(3/4))`, and all other blocks are
smaller.  Thus optimal-order cut cap, the known exact-minimizer operator
scale, one-spin child optimality, and even zero cross energy coexist with
the full `n^(5/2)` internal wall.  The signing is not known to minimize
`Q` at order `n`; this construction isolates exact edge-sign minimality as
the only available hypothesis not yet used.

## 5. Exact finite audit

All child-ground incidences and all outside completions were enumerated for
the stored `A8`, `A9`, and sampled `A10` exact minimizers.  The checker
verifies (L53.1)--(L53.5) exactly.

- On `A9,m=6`, a selector maximizing the child cap has
  `Q_S=22`, `X=0`, internal fields `(3,5,5,3,3,3)`, and `I=86`.
- On sampled `A10,m=6`, a maximizing selector likewise has
  `Q_S=22`, `X=0`, `r_max=5`, and `I=86`.
- Thus `X=0` does not force equal or exceptionally small internal fields,
  even for a maximal selector of a stored exact minimizer.
- The best sampled `A10,m=6` box has `Q_S=10`, `I=30`, `X=8`, true minimum
  row `10`, and uniform-completion mean `74`.  It is nonmaximal, and the
  quadratic correction in (L53.3) is decisive.  This confirms that
  (L53.6) is sufficient rather than necessary.

These finite facts are mechanism warnings, not asymptotic falsifiers.

## 6. Remaining exact target

The cleanest low-cross continuation is:

> **Open minimizer-specific field lemma.**  For some fixed `c in (0,1/4)`,
> every required exact order-`n` minimizer and every compact fixed-density
> window admit an `m`-selector and a child ground satisfying (L53.7).

It immediately proves the child-ground box witness by (L53.1)--(L53.2).
It would be falsified as an implementation by a scalable exact-minimizer
family for which the minimum of

```math
X+Q_S\lVert r\rVert_\infty
```

over all child-ground incidences is `n^(9/4-o(1))` or larger at every usable
`c`.  Such a family would not falsify the box witness itself, because the
quadratic completion term in (L53.3) can beat the mean dramatically.

An alternative continuation should therefore attack (L53.3) directly:
prove that whenever the low-cross branch has a spike violating (L53.7), the
outside `A^2` quadratic/linear discrepancy has a completion canceling that
spike at project scale.  No known maximal-selector or trace identity gives
this cancellation.
