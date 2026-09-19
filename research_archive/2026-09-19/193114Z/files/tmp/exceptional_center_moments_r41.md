# Wave 41: fixed-slice moments and the exceptional-center tail

## 1. Verdict

The signed principal energy of a fixed center has an exact second moment
depending only on its parent energy and row square.  Combining Bernoulli
conditioning, Hanson--Wright, and the competitive spectral bound gives
stretched-exponential concentration at exactly the required speed
`exp{-Omega(n^(3/4-c))}` for every low-row center.

This does **not** prove the exceptional-center lower tail.  The absolute
deficit also contains the random principal norm `Q(A[S])`, and neither its
mixed term with the center energy nor its nonlinear slice cumulants are
controlled by global signing minimality.  An exact `A_8` example shows that
the same parent energy, row square, first two signed-energy moments, and even
the same marginal principal-norm distribution can have different favorable
tails.  Thus a second-moment or structure-free Paley--Zygmund argument stops
at a genuine mixed-correlation obstruction.

There is nevertheless a clean nonlinear sufficient successor: a low-row
center whose `Q(A[S])`-weighted signed-energy exponential moment has positive
pressure and project-scale replica gap.  Paley--Zygmund then proves (10.967).
No current consequence of exact minimality establishes those two cumulant
bounds.

## 2. Exact fixed-slice first and second moments

Fix a spin `z`, switch the signing by it, and write

```math
b_{ij}=A_{ij}z_iz_j,
\qquad
E=z^TAz,
\qquad
R=z^TA^2z=\|Az\|_2^2.
\tag{R41.1}
```

For `S` uniform among the `m`-subsets, put

```math
C_S=z_S^TA[S]z_S
=2\sum_{i<j}b_{ij}\mathbf1_{\{i,j\in S\}},
\qquad
p_j=\frac{(m)_j}{(n)_j}.
\tag{R41.2}
```

Then

```math
\boxed{
\mathbb E C_S=p_2E,
}
\tag{R41.3}
```

and

```math
\boxed{
\mathbb E C_S^2
=p_4E^2+4(p_3-p_4)R
+2n(n-1)(p_2-2p_3+p_4).
}
\tag{R41.4}
```

Consequently

```math
\boxed{
\operatorname{Var}(C_S)
=(p_4-p_2^2)E^2+4(p_3-p_4)R
+2n(n-1)(p_2-2p_3+p_4).
}
\tag{R41.5}
```

To prove (R41.4), put `N=binom(n,2)`.  In the square of (R41.2), equal
edges have inclusion probability `p_2`, distinct edges meeting at one vertex
have probability `p_3`, and disjoint edges have probability `p_4`.  If

```math
W_1=\sum_{\substack{e<f\\|e\cap f|=1}}b_eb_f,
\qquad
W_0=\sum_{\substack{e<f\\e\cap f=\varnothing}}b_eb_f,
```

then

```math
W_1=\frac{R-n(n-1)}2,
\qquad
W_0=\frac{(E/2)^2-N}{2}-W_1.
```

Substitution into
`4[Np_2+2p_3W_1+2p_4W_0]` gives (R41.4).  No inequality or minimality is
used.

## 3. Exact higher moments and an affordable concentration upper bound

There is a complete high-moment expansion.  For an ordered `k`-tuple of
edges `mathbf e=(e_1,...,e_k)`, let `vmathbf(e)` be the number of vertices in
their union.  Then

```math
\boxed{
\mathbb E C_S^k
=2^k\sum_{\mathbf e\in E(K_n)^k}
p_{v(\mathbf e)}\prod_{r=1}^k b_{e_r}.
}
\tag{R41.6}
```

For `k>=3`, this introduces signed counts of triangles, paths, stars, and all
higher multigraph patterns; they do not collapse to `E` and `R`.  At the
moment order needed for a stretched-exponential tail, (R41.6) is therefore a
large new hierarchy.

A norm bound avoids enumerating that hierarchy.  Assume the density `m/n`
stays in a fixed compact subinterval of `(0,1)`.  Let independent
`xi_i~Ber(m/n)` and condition on `sum_i xi_i=m`.  Before conditioning,

```math
\xi^TB\xi-p^2E
=2p(B\mathbf1)^T(\xi-p\mathbf1)
+(\xi-p\mathbf1)^TB(\xi-p\mathbf1),
\tag{R41.7}
```

where `B=diag(z)A diag(z)`, `||B 1||^2=R`, and
`||B||_F^2=n(n-1)`.  Scalar subgaussian concentration and Hanson--Wright,
followed by division by
`P{Bin(n,p)=m}=Theta(n^(-1/2))`, give the audited bound

```math
\boxed{
\begin{aligned}
\Pr_{U_m}\{|C_S-p_2E|\ge u\}
\le C\sqrt n\exp\left[-c\min\left\{
\frac{(u-q_n/n)_+^2}{R+n^2},
\frac{(u-q_n/n)_+}{\|A\|_{op}}
\right\}\right].
\end{aligned}
}
\tag{R41.8}
```

Constants depend only on the fixed density window.  The harmless `q_n/n`
shift is the difference between the product center `p^2E` and the slice
center `p_2E`.

Equivalently, for `k>=C log n`,

```math
\boxed{
\|C_S-p_2E\|_{L^k(U_m)}
\le C\left[
\sqrt{k}(\sqrt R+n)+k\|A\|_{op}+\frac{q_n}{n}
\right].
}
\tag{R41.9}
```

The exact norm interpolation (10.67) gives

```math
\|A\|_{op}\le\sqrt{2Q(A)}.
\tag{R41.10}
```

For an exact minimizer, `q_n=O(n^(3/2))`; for a low-row center,
`R=O(n^2)`.  Taking `u=K n^(3/2-c)` in (R41.8), with `0<c<1/4`, yields

```math
\boxed{
\Pr_{U_m}\{|C_S-p_2E|\ge K n^{3/2-c}\}
\le\exp\{-c_K n^{3/4-c}\}
}
\tag{R41.11}
```

after increasing the fixed threshold constant.  The spectral term
`u/||A||_op` is the limiting exponent; the quadratic term is
`n^(1-2c)`.  Thus the signed-energy high-moment method reaches exactly the
project speed as an **upper concentration bound**.

## 4. Why the absolute deficit does not follow

Write

```math
Q_S=Q(A[S]),
\qquad
\delta_S^{abs}(z)=Q_S-|C_S|.
```

Its first two moments are exactly

```math
\boxed{
\begin{aligned}
\mathbb E\delta_S^{abs}
&=\mathbb E Q_S-\mathbb E|C_S|,\\
\mathbb E(\delta_S^{abs})^2
&=\mathbb E Q_S^2-2\mathbb E[Q_S|C_S|]
+\mathbb E C_S^2.
\end{aligned}
}
\tag{R41.12}
```

Only the last term is supplied by (R41.4).  Global minimality gives the coarse
pointwise interval

```math
q_m\le Q_S\le q_n.
\tag{R41.13}
```

The upper inequality follows by averaging full extensions of a child spin;
some extension has absolute energy at least its principal energy.  Neither
(R41.13), the cut cap, the spectral bound, nor row minimality controls
`E Q_S`, `E[Q_S|C_S|]`, or their higher mixed cumulants.  Using `E Q_S` at
the needed value is exactly the circular uniform-mean route (10.970).

If `sigma=sign(E)` and `Delta=q_n-|E|`, then (R41.11) says that outside an
affordable exceptional set

```math
|C_S|\ge p_2(q_n-\Delta)-u.
\tag{R41.14}
```

Combining this with a principal-norm lower tail would give

```math
Q_S\le p^{3/2}q_n-p_2\Delta+t-u
\quad\Longrightarrow\quad
\delta_S^{abs}(z)\le B_{n,m}+t.
\tag{R41.15}
```

But even the existence of one selector satisfying the left side of (R41.15)
would already be a power-saving restriction statement.  Treating that as an
input merely relocates the original problem.  Signed concentration alone
therefore cannot close the route.

## 5. A genuinely nonlinear moment criterion

The weakest clean successor exposed by the calculation is a **mixed
principal-norm/signed-energy cumulant**, not another moment of `C_S` alone.
For a low-row center `z`, set

```math
H=B_{n,m}+t,
\qquad
M_S(z)=H-Q_S+\sigma C_S,
\qquad
K_z(\theta)=\log\mathbb E_{U_m}e^{\theta M_S(z)}.
\tag{R41.16}
```

Since `|C_S|>=sigma C_S`,

```math
M_S(z)\ge0
\quad\Longrightarrow\quad
\delta_S^{abs}(z)\le H.
\tag{R41.17}
```

Let `L_0=n^(3/4-c)` and take `theta=Theta(n^(-3/4))`.  The following is an
exact sufficient theorem:

```math
\boxed{
R(z)\le2n(n-1),\qquad
K_z(\theta)\ge aL_0,qquad
K_z(2\theta)-2K_z(\theta)\le bL_0
}
\tag{R41.18}
```

for fixed `a>0,b<infinity`.  Indeed, with `W=e^(theta M_S)`, apply
Paley--Zygmund at fraction `alpha=e^(-aL_0)`.  Since
`alpha E W>=1`,

```math
\boxed{
U_m\{\delta_S^{abs}(z)\le H\}
\ge U_m\{M_S\ge0\}
\ge(1-e^{-aL_0})^2e^{-bL_0}.
}
\tag{R41.19}
```

Thus (R41.18) gives a zero-radius exceptional center, stronger than
(10.967), and proves the restriction chain.

This criterion is not a disguised hard-tail assertion: it asks for one
positive exponential pressure and its one-replica gap.  It also pinpoints why
linearization is circular.  If
`G=E Q_S-p^(3/2)q_n`, then

```math
K_z'(0)=t-G-p_2\Delta.
\tag{R41.20}
```

So the first derivative is exactly the failed mean balance.  A successful
proof of (R41.18) must use nonlinear curvature at `theta~n^(-3/4)` and the
mixed correlation between high principal norms and exceptional aligned
center energy.  Hanson--Wright controls the unweighted `C_S` factor, but
global minimality currently supplies no bound for the `e^(-theta Q_S)`
weight or for the replica gap in (R41.18).

## 6. Exact finite collision inside `A_8`

Let `A_8` be the verified exact minimizer in the repository, with `q_8=20`,
and take `m=4`.  The two projective cuts

```math
z_0=(1,-1,-1,-1,-1,-1,-1,-1),
```

and

```math
z_{12}=(1,-1,-1,-1,1,1,-1,-1)
```

both have

```math
(E,R)=(-16,72),
\qquad72<2n(n-1)=112.
```

Hence (R41.3)--(R41.5) give the same signed first two moments,

```math
\mathbb E C_S=-\frac{24}{7},
\qquad
\mathbb E C_S^2=\frac{208}{7}.
```

Both also have the same marginal principal-norm histogram:
`Q_S=8` on 60 selectors and `Q_S=12` on 10.  Nevertheless their absolute
deficit histograms are

```text
z_0:   {0:19, 4:34, 8:15, 12:2},
z_12:  {0:22, 4:28, 8:16, 12:4}.
```

Since

```math
B_{8,4}=20\left(2^{-3/2}-\frac3{14}\right)=2.7853\ldots<4,
```

their `t=0` favorable selector masses are respectively `19/70` and `22/70`.
Their third signed moments already differ.  This is a finite mechanism wall,
not an asymptotic falsifier of (10.967): both tails are large, and no unbounded
family is constructed.  It proves sharply that `E`, `R`, the complete
second-moment identity, and the marginal law of `Q_S` do not determine the
mixed favorable tail, even inside an exact minimizer.

## 7. Verification and recommendation

`tmp/exceptional_center_moments_r41_check.py` exhaustively verifies
(R41.3)--(R41.5) for every projective center and every proper slice size of
`A_6,A_8,A_9`, and verifies all `A_8` counts above exactly.

The exceptional-center route remains viable but the ordinary moment route is
now sharply delimited.  Do not retry the uniform mean or an unweighted
`C_S` high moment.  A meaningful successor must attack the mixed pressure and
replica gap (R41.18), or supply an actual-minimizer inverse theorem forcing
positive alignment of large `Q(A[S])` with rare signed slice deviations at
the Hanson--Wright boundary scale.
