# Wave 27 Route 1: low-cost completion of an exact child ground

## Status

- **Proved:** An exact oriented child ground has favorable retained effective
  loss under *every* full spin completion.  There is no remaining
  loss/parent-deficit tradeoff; the problem is purely row-square completion.
- **Proved:** Exact product-rounding expectation, variance, parent-energy
  variance, and covariance formulas are given below.  They yield deterministic
  conditional-expectation, zonotope-distance, Bhatia--Davis, and one-sided
  hypercontractive completion bounds.
- **Falsified (finite, scoped):** Exact child-ground optimality alone does not
  force a favorable variance correction or zonotope membership.  Two exact
  one-deletion child-ground fibers of `A_9` have row square identically `128`
  over both completions, hence conditional variance zero.  All of
  `A_6,A_8,A_9` are invertible, so the exact zonotope-membership premise fails
  for every nonempty partial spin, ground or otherwise.
- **Open:** None of the per-selector completion results compresses the outputs
  for different selectors into one fixed cut with the coverage required by
  (10.795).  Nor is the partial-field or zonotope-distance bound below proved
  at the target scale for asymptotic exact minimizers.

The checker is `tmp/low_cost_completion_r27.py`.  It uses the full outside
cube and exact integer/rational arithmetic; in particular, it does not make
the projective-outside enumeration error of identifying `z` with `-z` while
the child spin is held fixed.

## 1. Exact child payoff makes every completion loss-favorable

Let `A` be an exact order-`n` minimizer, so `Q(A)=q_n`.  Fix an `m`-set
`S`, put `T=[n]\setminus S`, `k=n-m`, and let `(sigma,y)` be an exact
oriented child ground:

```math
\sigma y^{\mathsf T}A[S]y=Q(A[S]).
```

For any $`z\in\{\pm1\}^T`$, let `x=(y,z)`, let
`d=(sigma,x)`, and write

```math
H(z)=\langle A,d\rangle=\sigma x^{\mathsf T}Ax.
```

The child loss is exactly zero.  Substituting in (10.792), with
`B_{n,m}=(p^{3/2}-p_2)q_n`, gives the pointwise identity

```math
\boxed{
\widehat\ell(S,d)=p_2H(z)-p^{3/2}q_n.
}
\tag{R27.1}
```

Every oriented full cut satisfies `H(z)<=Q(A)=q_n`.  Since
`p_2<p^2<=p^{3/2}` for `0<m<n`, (R27.1) gives

```math
\boxed{
\widehat\ell(S,d)\le(p_2-p^{3/2})q_n=-B_{n,m}<0
\quad\hbox{for every completion }z.
}
\tag{R27.2}
```

Equivalently, the exact child payoff and the nonnegative parent deficit
jointly pay more than the full retained allowance.  Thus choosing outside
spins to reduce `R_2` cannot spoil the effective-loss condition.  This is
stronger than the conditional-mean statement (10.803).

## 2. Exact biased completion formulas

Put

```math
b=A[:,S]y,
\qquad V=A[:,T],
\qquad h=V^{\mathsf T}b,
\qquad K=V^{\mathsf T}V.
```

Every column of `V` has squared norm `n-1`, so `K_{jj}=n-1`.  Independently
round the outside coordinates with
`E z_j=u_j in [-1,1]`, and put `v_j=1-u_j^2`.  Direct expansion gives

```math
\boxed{
\Phi(u):=\mathbb E_u R_2(d)
=\lVert b+Vu\rVert_2^2+(n-1)\sum_{j\in T}v_j.
}
\tag{R27.3}
```

Define

```math
\alpha_j(u)=h_j+\sum_{l\ne j}K_{jl}u_l.
```

The centered variables `z_j-u_j` and their distinct quadratic products are
pairwise orthogonal in `L^2` of the product law.  Therefore

```math
\boxed{
\operatorname{Var}_u R_2
=4\sum_jv_j\alpha_j(u)^2
 +4\sum_{j<l}v_jv_lK_{jl}^2.
}
\tag{R27.4}
```

For completeness, put `g=A[T,S]y`, `D=A[T]`, and

```math
\beta_j(u)=g_j+\sum_{l\ne j}D_{jl}u_l.
```

The corresponding parent formulas, with the ordered quadratic-form
normalization used in the ledger, are

```math
\boxed{
\begin{aligned}
\mathbb E_uH
 &=\sigma\left(y^{\mathsf T}A[S]y
     +2g^{\mathsf T}u+2\sum_{j<l}D_{jl}u_ju_l\right),\\
\operatorname{Var}_uH
 &=4\sum_jv_j\beta_j(u)^2
   +4\sum_{j<l}v_jv_lD_{jl}^2,\\
\operatorname{Cov}_u(R_2,H)
 &=4\sigma\sum_jv_j\alpha_j(u)\beta_j(u)
   +4\sigma\sum_{j<l}v_jv_lK_{jl}D_{jl}.
\end{aligned}}
\tag{R27.5}
```

For uniform outside spins (`u=0`), these specialize to

```math
\boxed{
\begin{aligned}
\mu&:=\mathbb E R_2=\lVert A[:,S]y\rVert_2^2+k(n-1),\\
v&:=\operatorname{Var}R_2
 =4\lVert(A^2)[T,S]y\rVert_2^2
  +4\sum_{j<l\in T}(A^2)_{jl}^2,\\
\mathbb EH&=Q(A[S]),\\
\operatorname{Var}H&=4\lVert A[T,S]y\rVert_2^2+2k(k-1).
\end{aligned}}
\tag{R27.6}
```

Thus (10.804) at zero child noise is recovered, but (R27.3)--(R27.6) also
give biased cancellation, exact fluctuation, and deterministic extraction.

## 3. Four deterministic completion bounds

### 3.1 Conditional expectation and an exact greedy identity

For every `u in [-1,1]^T`, some vertex `z` in the support of its product
rounding satisfies

```math
\boxed{R_2(y,z)\le\Phi(u).}
\tag{R27.7}
```

This can be derandomized one coordinate at a time.  In the uniform case,
fix an ordering `j_1,...,j_k` and choose recursively

```math
z_{j_r}=-\operatorname{sign}\left(
h_{j_r}+\sum_{s<r}K_{j_rj_s}z_{j_s}\right),
```

with either sign at a zero.  The final row square obeys the exact identity

```math
\boxed{
R_2(y,z)=\mu-2\sum_{r=1}^k
\left|h_{j_r}+\sum_{s<r}K_{j_rj_s}z_{j_s}\right|
\le\mu.
}
\tag{R27.8}
```

This is a deterministic improvement of merely sampling the outside-uniform
channel.

### 3.2 Zonotope-distance completion

Let

```math
\delta(S,y)=\operatorname{dist}_2(-b,V[-1,1]^T).
```

Choose a box least-squares minimizer `u_*`.  Equations (R27.3) and (R27.7)
give

```math
\boxed{
\min_zR_2(y,z)
\le \delta(S,y)^2+(n-1)(k-\lVert u_*\rVert_2^2)
\le \delta(S,y)^2+k(n-1).
}
\tag{R27.9}
```

In particular, exact membership `-b in V[-1,1]^T` gives an `O(n^2)`
completion, enough for every fixed `c<1/4` in (10.795).  Membership has the
transparent equivalent form

```math
\boxed{
-b\in V[-1,1]^T
\quad\Longleftrightarrow\quad
\exists w\in\ker A:\ w_S=y, \lVert w_T\rVert_\infty\le1.
}
\tag{R27.10}
```

Thus it is a genuine kernel/zonotope hypothesis, not a consequence of child
optimality.

### 3.3 Bhatia--Davis extraction

For an exact minimizer, (10.650) bounds every full-cut row square by

```math
M=2(n-1)q_n.
```

Let `a=min_z R_2(y,z)`.  Since `(R_2-a)(M-R_2)>=0`, taking expectation
under any product rounding gives

```math
\operatorname{Var}_uR_2
\le(M-\Phi(u))(\Phi(u)-a).
```

Consequently, whenever `Phi(u)<M`,

```math
\boxed{
\min_zR_2(y,z)
\le\Phi(u)-\frac{\operatorname{Var}_uR_2}{M-\Phi(u)}.
}
\tag{R27.11}
```

If `Phi(u)=M`, boundedness forces zero variance and constant row square on
the support.  No sign or factor of two is hidden in (R27.11).

### 3.4 A one-sided degree-two fluctuation gain

For uniform outside spins, `F=R_2-mu` is a centered Rademacher polynomial
of degree at most two.  Bonami--Beckner gives
`||F||_3<=2||F||_2`.  A short one-sided argument turns this into an actual
low-row completion:

```math
\boxed{
\min_zR_2(y,z)\le\mu-\frac1{10}\sqrt v.
}
\tag{R27.12}
```

Indeed, if the negative part `F_-` were pointwise at most
`delta sqrt(v)`, then
`E F_+=E F_-<=delta sqrt(v)`.  Hence

```math
\mathbb EF_-^2\le\delta^2v,
\qquad
\mathbb EF_+^2
\le(\mathbb EF_+\,\mathbb EF_+^3)^{1/2}
\le\sqrt{8\delta}\,v.
```

With `delta=1/10`, their sum is strictly less than `v`, a contradiction.
This proves (R27.12).  It is useful only when the fiber variance is nonzero;
the `A_9` audit below shows that this qualification is essential.

## 4. Exact conditional sufficient theorem and its scope

Combine (R27.2) with any of (R27.7), (R27.9), (R27.11), or (R27.12).
For example, a directly checkable sufficient condition for the desired
per-selector completion is

```math
\boxed{
\min\left\{
\delta(S,y)^2+k(n-1),
\ \mu-\frac{v}{M-\mu},
\ \mu-\frac1{10}\sqrt v
\right\}
=O(n^{9/4-c}),
}
\tag{R27.13}
```

where the middle entry is omitted if `mu=M`.  Then there is a full cut `d`
with

```math
R_2(d)=O(n^{9/4-c}),
\qquad \widehat\ell(S,d)\le-B_{n,m}<0.
```

A simpler sufficient hypothesis is

```math
\lVert A[:,S]y\rVert_2^2=O(n^{9/4-c}),
```

because the remaining term `k(n-1)=O(n^2)` is lower order for fixed
`c<1/4`.  Exact zonotope membership is stronger still and gives `O(n^2)`.

This is a **per-selector** theorem.  The selected completion may depend on
both `S` and the chosen child ground.  Criterion (10.795), in contrast,
requires one fixed full cut to work on a stretched-exponential fraction of
the entire slice.  Neither conditional expectation nor (R27.13) bounds the
number of distinct outputs or forces collisions between different child
fibers.  Therefore this theorem does not prove (10.795), even if its row
premise were established for every selector.

## 5. Exact `A_6,A_8,A_9` obstruction audit

For `m=n-1`, exhaustive enumeration gives the following types.  Multiplicity
counts oriented projective child grounds over all deleted vertices; the two
outside signs are both retained.

| signing | multiplicity | `mu` | `v` | `min R_2` | `max R_2` |
|:---:|---:|---:|---:|---:|---:|
| `A_6` | 60 | 30 | 0 | 30 | 30 |
| `A_8` | 24 | 68 | 16 | 64 | 72 |
| `A_9` | 2 | 88 | 64 | 80 | 96 |
| `A_9` | 6 | 96 | 0 | 96 | 96 |
| `A_9` | 4 | 104 | 64 | 96 | 112 |
| `A_9` | 6 | 112 | 0 | 112 | 112 |
| `A_9` | 6 | 112 | 256 | 96 | 128 |
| `A_9` | 2 | 120 | 64 | 112 | 128 |
| `A_9` | 2 | 128 | 0 | 128 | 128 |

The final `A_9` type is the sharp finite variance obstruction: on two exact
child-ground fibers every full completion has row square `128`.  Thus no
theorem using only child exactness and a strictly positive variance rebate
can cover all fibers.  This is not an asymptotic falsifier for the
`O(n^{9/4-c})` target.

The exact determinants are

```math
\det A_6=-125,
\qquad \det A_8=729,
\qquad \det A_9=808.
```

By (R27.10), zonotope membership fails for every nonempty partial spin in
all three examples.  For the one-deletion child-ground fibers, the exact
squared distances are

| signing | squared zonotope distances (multiplicity) |
|:---:|:---|
| `A_6` | `25` (60) |
| `A_8` | `423/7` (24) |
| `A_9` | `78` (2), `88` (6), `94` (4), `96` (6), `104` (6), `110` (2), `120` (2) |

On the worst `A_9` fiber, (R27.9) gives `120+8=128`, exactly the constant
row square.  Hence even the distance completion can be sharp on a finite
exact minimizer.

## 6. Verdict

- **Proved:** favorable retained loss is automatic for every exact-child
  completion; exact biased moments and several deterministic row-completion
  bounds hold.
- **Falsified, finite and scoped:** exact child-ground status does not imply
  zonotope membership, nonzero completion variance, or a row improvement
  beyond the exact mean on every fiber.
- **Open target:** prove that some target-specific child ground (or enough
  child grounds with shared completions) satisfies the partial-field,
  distance, or variance-corrected bound in (R27.13), together with the
  fixed-cut coverage/compression demanded by (10.795).
