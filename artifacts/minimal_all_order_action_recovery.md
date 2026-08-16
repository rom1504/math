# Minimal all-order recovery for the quadratic signing problem

Date: 2026-08-16.

Status: proved implication and quantitative error audit. This note does **not**
prove the recovery hypothesis introduced below.  The selected-phase target in
this note is valid but is no longer the weakest retained optimizer-free target;
the second-checkpoint envelope weakening is in
`extremal_envelope_recovery.md`.

The final implication is

```math
\boxed{
\text{tolerance-dependent spectral purification}
+\text{ action compactness}
+\mathrm{AR}_{\min}^{\to}
\Longrightarrow
\frac{M_n}{n^{3/2}}\text{ converges}.}
```

## 1. Normalization

For a symmetric hollow signing (A) of order (n), put

```math
H_A(x)=\sum_{i<j}a_{ij}x_ix_j,
\qquad
Q(A)=\max_{x\in\{\pm1\}^n}|H_A(x)|,
\qquad
a_n=\frac{M_n}{n^{3/2}}.
```

On the uniform (n)-point probability space let (T_A=A/\sqrt n), with
the repository's action-operator convention. Then

```math
\Phi(T_A)
:=\sup_{|f|\le1}|\langle f,T_Af\rangle|
=\frac{2Q(A)}{n^{3/2}}.                                    \tag{AR.1}
```

For a (P)-operator (T), write ​\(\mathcal S_1(T)\) for its closed
one-profile: the set of laws of ((f,Tf)) with ​\(|f|\le1\). Define the
directed one-profile distance

```math
\partial_1(S,T)
=\sup_{\mu\in\mathcal S_1(S)}
  \inf_{\nu\in\mathcal S_1(T)}d_{\rm LP}(\mu,\nu).          \tag{AR.2}
```

This is deliberately not a metric. It only says that (S) has no
one-profile far outside the profile set of (T); it need not realize all
profiles of (T), and no joint (k\)-profile with (k\ge2) is retained.

## 2. The one-sided continuity estimate

The proof of quantitative action continuity is directed. If

```math
\lVert S\rVert_{2\to2},\lVert T\rVert_{2\to2}\le D,
\qquad
\partial_1(S,T)\le\delta,
```

then the same Strassen-coupling and output-truncation argument gives

```math
\boxed{
\Phi(S)
\le \Phi(T)+5D\sqrt\delta+\delta .}                         \tag{AR.3}
```

The constants can be replaced by the slightly looser action-metric form

```math
\Phi(S)\le\Phi(T)+5D\sqrt{2d_M(S,T)}+2d_M(S,T).             \tag{AR.4}
```

Only the inclusion from ​\(\mathcal S_1(S)\) toward
​\(\mathcal S_1(T)\) is used for this upper bound. Consequently, if
the realizing norm bound is (D_j), the quantitative requirement is

```math
D_j\sqrt{\delta_j}\longrightarrow0,                         \tag{AR.5}
```

not necessarily a fixed (D). With a common bound, merely
​\(\delta_j=o(1)\) suffices.

## 3. Purified extremal cluster objects

Put (L=\liminf_n a_n). The proved purification theorem has the following
two-parameter consequence. For every ​\(\eta>0\), there are exact hollow
signings (B_j^{(\eta)}) of orders (n_j\to\infty), a finite constant
​\(C_\eta\), and an action-limit object (T_\eta) such that

```math
\lVert T_{B_j^{(\eta)}}\rVert_{2\to2}\le C_\eta,
\qquad
T_{B_j^{(\eta)}}\longrightarrow T_\eta,
\qquad
2L\le\Phi(T_\eta)\le2L+\eta.                               \tag{AR.6}
```

Here the harmless reparameterization of ​\(\eta\) absorbs the
​\(O(\sqrt\varepsilon)\) purification loss. In particular, a single
operator bound uniform as ​\(\eta\downarrow0\) is not needed for the
convergence implication: first work at fixed ​\(\eta\), then send
​\(\eta\downarrow0\).

## 4. The logically weakest recovery statement

Call a set of orders ​\(\mathcal N\subset\mathbb N\) **upward
ratio-dense** if

```math
s(N):=\min\{m\in\mathcal N:m\ge N\}
\quad\text{satisfies}\quad
\frac{s(N)}N\longrightarrow1.                               \tag{AR.7}
```

Equivalently, every large (N) has an (m\in\mathcal N) with
​\(N\le m=N+o(N)\), uniformly in (N).

Let

```math
p_n:=\frac{2M_n}{n^{3/2}},
\qquad
\alpha:=\liminf_{n\to\infty}p_n=2L.
```

The literally weakest sufficient recovery statement is the following.
Choose any sequence \(\eta_\ell\downarrow0\).  For each \(\ell\), require
one cluster object \(T_\ell\) from (AR.6), an upward ratio-dense set
\(\mathcal N_\ell\), and exact signings \(C_m\),
\(m\in\mathcal N_\ell\), such that

```math
\limsup_{\substack{m\to\infty\\m\in\mathcal N_\ell}}
\Phi(T_{C_m})
\le \Phi(T_\ell)+r_\ell,
\qquad r_\ell\longrightarrow0.                            \tag{MR}
```

This is enough for convergence, but it is not a useful reduction.  Once
purification and deletion are available, (MR) is equivalent to the missing
scalar \(\limsup\le\liminf\) assertion: if convergence is already known,
exact minimizers at all orders witness (MR).  A noncircular recovery theorem
must construct its witnesses from independently controlled information in
the limit object, rather than selecting target-order minimizers.

## 5. The weakest candidate structural recovery statement

> **Directed extremal recovery ​\(\mathrm{AR}_{\min}^{\to}\).** For every
> member of some null sequence of tolerances \(\eta_\ell\downarrow0\), at
> least one purified extremal cluster \(T_{\eta_\ell}\) satisfying
> (AR.6) has an upward ratio-dense set ​\(\mathcal N_\eta\) and exact
> symmetric hollow signings (C_m), (m\in\mathcal N_\eta), such that,
> for some bounds (D_m\ge
> \max\{\lVert T_{C_m}\rVert_{2\to2},\lVert T_\eta\rVert_{2\to2}\}),
>
> ```math
> \delta_m:=\partial_1(T_{C_m},T_\eta)\longrightarrow0,
> \qquad
> D_m\sqrt{\delta_m}\longrightarrow0.                       \tag{AR.8}
> ```

This statement is weaker than the archived full AR in four ways:

1. it concerns only one selected purified liminf cluster at each fixed
   tolerance, not every signed action limit;
2. it requires only directed one-profile inclusion, not full action
   convergence;
3. it allows recovery only on an upward ratio-dense set, not at every order;
4. the operator bound may grow if the profile error shrinks faster than its
   inverse square.

No recovery is required for every signed limit, every purified cluster, or
even every positive tolerance.  Conversely, the action requirement cannot
be dropped in favor of its scalar consequence without reverting to (MR).

## 6. Convergence theorem

> **Theorem.** The proved purification theorem, bounded-operator action
> compactness, and ​\(\mathrm{AR}_{\min}^{\to}\) imply that
> ​\(M_n/n^{3/2}\) converges.

**Proof.** Fix ​\(\eta>0\), select (T_\eta) as in (AR.6), and apply
(AR.3) to the recovery sequence. Equations (AR.6) and (AR.8) give

```math
\limsup_{\substack{m\to\infty\\m\in\mathcal N_\eta}}
\frac{Q(C_m)}{m^{3/2}}
\le \frac12\Phi(T_\eta)
\le L+\frac\eta2.                                           \tag{AR.9}
```

For (N\le m), take any (N)-vertex principal submatrix (C_m[S]).
For fixed (x\in\{\pm1\}^S), extend it by independent uniform signs (Y)
on the deleted vertices. All cross terms and all terms internal to the
deleted set have mean zero, hence

```math
H_{C_m[S]}(x)=\mathbb E_Y H_{C_m}(x,Y).
```

Therefore

```math
Q(C_m[S])\le Q(C_m),
\qquad M_N\le M_m.                                          \tag{AR.10}
```

Choose (m=s(N)). Since (m/N\to1), (AR.9)--(AR.10) imply

```math
\limsup_{N\to\infty}a_N
\le L+\frac\eta2.
```

Finally send ​\(\eta\downarrow0\). The reverse inequality is the
definition of (L). ​\(\square\)

The same proof works if the recovery orders are listed increasingly as
​\(m_j\) with (m_{j+1}/m_j\to1).

## 7. Near-order transfer and exact error scales

### 7.1 Deletion

Equation (AR.10) is lossless. Thus a realization at any
​\(m=N+o(N)\) transfers to exact order (N) with only the normalization
factor

```math
\left(\frac mN\right)^{3/2}=1+o(1).                         \tag{AR.11}
```

No action-distance estimate for the deleted matrix is needed.

### 7.2 Padding

The archived random rectangular-padding inequality is

```math
M_{n+h}\le M_n+M_h+
\sqrt{2nh(n+h+2)\log2}.                                    \tag{AR.12}
```

It gives an (o(n^{3/2})) insertion defect for every (h=o(n)). Hence
near-order transfer itself is already proved in both directions; what it
does not supply is a ratio-dense sequence carrying the liminf value.

There is also a direct completion form that does not require a separately
optimized order-\(h\) block.  Fix any signing \(A\) of order \(n\), put
\(N=n+h\), and independently sign the

```math
r=nh+\binom h2
```

new edges.  A union bound over the \(2^{N-1}\) spin pairs modulo global
sign gives a deterministic completion \(B\) with

```math
Q(B)\le Q(A)+\sqrt{2r(N+1)\log2}.                          \tag{AR.13}
```

Thus \(h=o(N)\) has insertion cost \(O(N\sqrt h)=o(N^{3/2})\).
Some balancing is essential: if every old--new edge is \(+1\), comparing
an old spin vector with its negative shows \(Q(B)\ge nh\), regardless of
the two internal blocks.  Arbitrary \(o(n)\)-vertex padding is therefore
not safe.

### 7.3 Edge edits

If (A) and (B) differ on (r) unordered off-diagonal edges, then

```math
|Q(A)-Q(B)|\le2r.                                           \tag{AR.14}
```

Thus a worst-case edit repair is negligible at the target scale only when
​\(r=o(n^{3/2})\). The condition ​\(r=o(n^2)\) is not sufficient.

For (E=A-B), the more structural sufficient bounds are

```math
|Q(A)-Q(B)|
\le\frac12\max_x|x^{\mathsf T}Ex|
\le\frac12\lVert E\rVert_{\infty\to1}
\le\frac n2\lVert E\rVert_{\rm op},                        \tag{AR.15}
```

so it is enough that the same-spin or ​\(\infty\to1\) error be
​\(o(n^{3/2})\), or that ​\(\lVert E\rVert_{\rm op}=o(\sqrt n)\).
The generic Frobenius estimate

```math
|Q(A)-Q(B)|\le\frac n2\lVert E\rVert_F                     \tag{AR.16}
```

requires ​\(\lVert E\rVert_F=o(\sqrt n)\). It is often much weaker than
(AR.14), since \(r\) sign edits have
​\(\lVert E\rVert_F=2\sqrt{2r}\).

Equivalently,

```math
|\Phi(T_A)-\Phi(T_B)|
\le \frac{\lVert E\rVert_{\ell^1}}{n^{3/2}}
\le \frac{\lVert E\rVert_F}{\sqrt n}.                       \tag{AR.17}
```

### 7.4 Weighted intermediate matrices

Intermediate matrices need not be hollow or sign-valued. A bounded
diagonal contributes only (O(n)=o(n^{3/2})), so hollowness can be imposed
at the end. A weighted ​\([-1,1]\) realization (W_n) is useful if a
dependent rounding produces a final signing (A_n) with one of the error
bounds in (AR.15)--(AR.17), and, when (AR.3) rather than direct objective control is
used, with ​\(D_n\sqrt{\partial_1(T_{A_n},T)}=o(1)\).

Merely obtaining ​\(\lVert A_n\rVert_{\rm op}=O(\sqrt n)\) does not
control the rounding defect: a residual of operator norm
​\(c\sqrt n\) can change the normalized objective by a fixed constant.

## 8. Strictness and circularity audit

The scalar statement (MR) is a reformulation, not a strict reduction.
Directed one-profile recovery does discard vertex labels, reverse profile
inclusion, and every joint profile of depth at least two, so it is a strict
information quotient of full action recovery.  It is nevertheless a stronger
proposition than scalar convergence, and it has not yet been shown easier to
prove.  Its unresolved universal quantifier says that the microscopic sign
residual creates no new Boolean profile with excessive energy.

The following assumptions would make a purported proof circular:

1. choosing \(C_m\) by minimizing the target-order objective;
2. assuming the desired asymptotics of \(M_m\) to bound \(\Phi(T_{C_m})\);
3. storing the complete map \(x\mapsto H_{C_m}(x)\), or an equivalent rooted
   coset histogram, as the realization state; or
4. proving only recovery along the original sparse liminf subsequence.

The first genuine strict reduction would be a theorem deriving
\(\mathrm{AR}_{\min}^{\to}\), or merely its scalar consequence, from a
state whose size and definition are independently controlled and
demonstrably smaller than the full Boolean landscape.  Whether its **proof
difficulty** is strictly lower remains open.  The next campaign question is
whether sampling, dependent rounding, absorption, or an order-dependent
microcanonical construction can enforce the outer profile without
reconstructing all Boolean spikes.

## 9. Fixed-accuracy information content

At any fixed accuracy, the directed one-profile has a finite abstract
description, although satisfying it remains a global constraint.

Suppose \(\lVert S\rVert_{2\to2}\le C\), and round any
\(f:\Omega\to[-1,1]\) pointwise to a grid of mesh \(h\), obtaining \(g\).
Under the identity coupling,

```math
|f-g|\le h,
\qquad
\lVert S(f-g)\rVert_2\le Ch.
```

Chebyshev's inequality at threshold \(\sqrt{Ch}\), followed by the identity
coupling and the triangle inequality in \(\mathbb R^2\), gives, when
\(Ch\le1\),

```math
d_{\rm LP}\bigl(\mathcal L(f,Sf),\mathcal L(g,Sg)\bigr)
\le h+\sqrt{Ch}.                                            \tag{AR.18}
```

Thus an \(\epsilon\)-approximation needs only an alphabet of size
\(O((1+C)/\epsilon^2)\). The output laws have uniformly bounded second
moment, so truncation followed by a finite bounded-Lipschitz test-function
net gives a finite description whose size depends on \((C,\epsilon)\),
not on \(n\).

This answers one design-theoretic question only halfway. The **state** is
finite at fixed accuracy, but the finite matrix must satisfy its outer
profile condition simultaneously for every coloring by that fixed alphabet,
still \(q^n\) colorings. Ordinary design and graphon theorems prescribe
finitely many local densities; they do not automatically enforce a universal
local-global coloring profile. A valid absorption theorem would have to
control this universal quantifier without enumerating all colorings.
