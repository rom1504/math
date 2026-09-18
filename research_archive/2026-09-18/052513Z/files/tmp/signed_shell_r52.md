# Wave 52: level-two signed shells and maximal-selector port exchange

## Status

No asymptotic signed-shell theorem, controlled positive prefix, restriction
estimate, or convergence proof is obtained.  There are three rigorous outputs:

1. an exact diagonal/off-diagonal decomposition of every level-two shell;
2. signing-specific maximal-selector port identities and a quantitative
   positive-threshold exchange star; and
3. a precise scale obstruction showing that a cap-neighborhood/local-exchange
   amplification from one seed is exponentially too small for the level-two
   threshold.

The finite shell signs are also audited.  Their apparent positivity is almost
entirely a small-order self-loop artifact and must not be treated as
asymptotic evidence.  Exact calculations are reproduced by
`tmp/signed_shell_r52_check.py` and `tmp/signed_shell_r52_check.out`.

## 1. Exact level-two shell decomposition

Let `N=binom(n,m)`, let `F_z` be the favorable selector family of a center,
and put `r_z=|F_z|=Na_z`.  Set

```math
d_2=\binom m2\binom{n-2}{m-2},\qquad
h_2=\binom{n-2}{m-2}^{-1}.
```

The number `h_2` is the diagonal atom `K_2(S,S)`.  Separating diagonal and
off-diagonal ordered pairs in the common-core formula gives the **Verified
exact identity**

```math
\boxed{
p_2(z)=h_2+\frac{\Omega_2(z)}{r_zd_2},\qquad
\Omega_2(z)=\sum_{\substack{S,T\in F_z\\S\ne T}}
\binom{|S\cap T|}{2}.}
\tag{R52S.1}
```

Consequently, for the row shell `H_R={z:R_2(z)=R}`, its unnormalized
spectral excess is exactly

```math
\boxed{
N^2\sum_{z\in H_R}a_z^2\{p_2(z)-\lambda_2(2)\}
=\sum_{z\in H_R}
\left[r_z^2\{h_2-\lambda_2(2)\}
+\frac{r_z}{d_2}\Omega_2(z)\right].}
\tag{R52S.2}
```

Thus, once `h_2<lambda_2(2)`, a nonnegative shell is precisely an
off-diagonal pair-core collision theorem:

```math
\sum_{z\in H_R}r_z\Omega_2(z)
\ge d_2\{\lambda_2(2)-h_2\}\sum_{z\in H_R}r_z^2.
\tag{R52S.3}
```

There is no hidden row monotonicity in this identity; all signing-specific
content must enter through `Omega_2` on the prescribed shell.

The elementary upper bound `Omega_2(z)<=r_z(r_z-1)binom(m,2)` gives

```math
p_2(z)\le r_zh_2=\frac{a_z}{\rho_2},\qquad
\rho_2=\frac{(m)_2}{(n)_2}.
\tag{R52S.4}
```

This independently cross-checks the main-agent core-load sandwich.  In
particular a nonnegative level-two prefix already forces a center of degree
`Omega(n^{-2})` at fixed density.  The shell target is clean, but it is not a
weak local-collision statement: it manufactures exponentially many favorable
selectors at one center.

At fixed density `p in (1/2,1)`,

```math
h_2=\exp\{-nH(p)+O(\log n)\},\qquad
\lambda_2(2)=\frac{2(1-p)^2+o(1)}{p^2n^2}.
```

Hence diagonal self-retention becomes negligible.  For example, at
`n=12,m=7`, `h_2=1/252<2/189=lambda_2(2)`.

## 2. Signing-specific maximal-selector exchange identities

Let

```math
q_* = \max_{|S|=m}Q(A[S]),
```

fix a maximizing selector `S`, and orient a ground `y` so that its energy is
`q_*`.  Let `T=S^c`, `k=|T|`, and complete `y` to a full center `z`.  For
`i in S` and `j in T` define

```math
r_i=\sigma y_i\sum_{u\in S\setminus\{i\}}a_{iu}y_u,
\qquad
t_{ji}=\sigma z_j\sum_{u\in S\setminus\{i\}}a_{ju}y_u.
```

One-spin optimality inside `S` and maximality over all `m`-selectors prove
the **Verified port bounds**

```math
\boxed{0\le r_i\le q_*/2,\qquad |t_{ji}|\le r_i.}
\tag{R52S.5}
```

Indeed, flipping `i` gives energy `q_*-4r_i`.  Replacing `i` by `j`, with
the two signs at the incoming port, gives candidate energies
`q_*-2r_i+/-2t_{ji}`; both have absolute value at most `q_*`.

If `l_i=sigma z_i(Az)_i`, direct summation gives two further **Verified exact
identities**:

```math
\boxed{
\sum_{i\in S,j\in T}(r_i-t_{ji})
=(n-1)q_*-(m-1)\sum_{i\in S}l_i,}
\tag{R52S.6}
```

```math
\boxed{
\sum_{i\in S,j\in T}(r_i^2-t_{ji}^2)
=k\lVert A[S]y\rVert_2^2
-(m-2)\lVert A[T,S]y\rVert_2^2-mk\ge0.}
\tag{R52S.7}
```

For (R52S.7), if
`h_j=sigma z_j sum_(u in S) a_(ju)y_u`, then
`sum_i t_(ji)^2=(m-2)h_j^2+m`.

There is a useful positive-threshold collision consequence.  For each port
`(i,j)`, choose between `z` and the outside flip `z^j` so the incoming sign
maximizes the exchanged energy.  For `S'=S-i+j`, its deficit obeys

```math
\boxed{D_{z\text{ or }z^j}(S')\le2(r_i-|t_{ji}|).}
\tag{R52S.8}
```

Moreover

```math
\sum_{i,j}2(r_i-|t_{ji}|)\le2kq_*.
```

Thus, for every `H>0`, at least

```math
mk-\frac{2kq_*}{H}
```

ports give an `H`-favorable exchanged selector, distributed among only the
`k+1` centers `z,z^j`.  All those centers satisfy

```math
R_2(z^j)\le\{\sqrt{R_2(z)}+2\sqrt n\}^2.
```

Taking `H=4q_*/m` proves that some one-flip center has at least

```math
\boxed{\frac{mk}{2(k+1)}}
\tag{R52S.9}
```

distinct near-ground exchanged selectors (up to the harmless integer
rounding).  This is a genuine signing-specific theorem for the
**positive-threshold family** `D_z(S)<=H`; it is not a zero-deficit overlap
theorem for `F_z={D_z=0}`.

At zero threshold, however, (R52S.8) only certifies an exact neighboring
`q_*`-ground when `|t_{ji}|=r_i`; neither (R52S.6) nor (R52S.7) forces such
saturation.

## 3. Why the local star cannot reach the level-two threshold

The star (R52S.9) has only `Theta(n)` selectors.  Its normalized degree is
`Theta(n)/binom(n,m)=exp{-Theta(n)}`, and (R52S.4) then makes its level-two
retention exponentially smaller than `lambda_2(2)=Theta(n^{-2})`.

The same obstruction persists under sublinear iteration.  Suppose an
amplification scheme starts from one seed selector, uses at most `u` physical
outside flips, and only records selectors within `u` element exchanges of
that seed.  Every reached center then has at most

```math
B_{m,k}(u)=\sum_{s\le u}\binom ms\binom ks
=\exp\{O(u\log(n/u))\}
\tag{R52S.10}
```

available selector labels.  The row-neighborhood inequality guarantees the
controlled target cap `O(R_0+n^2)` only for

```math
u=O\left(\sqrt{\frac{R_0+n^2}{n}}\right)
=O(n^{5/8-c/2})=o(n)
```

when `R_0=O(n^{9/4-c})`.  Hence (R52S.10) is `exp{o(n)}`, while
`binom(n,m)=exp{Theta(n)}`.  Such a path-local cap-neighborhood construction
has `p_2=exp{-Theta(n)}` and cannot establish a nonnegative level-two prefix.

This is a **Falsification only of local cap/exchange amplification from one
seed**.  It does not rule out a theorem forcing linear-distance selector
closure on the same low-row center, a global supply of exponentially many
seeds, or another exact-minimizer mechanism not charged one physical flip per
exchange.

## 4. Exact-minimizer finite audit and its scope

For the exact order-ten minimizer used in Waves 50--51 at `m=6`, every active
level-two shell is positive and every prefix has `P_2-lambda_2>0`.  At cap ten
the gap is `1/75`.  But here `h_2=lambda_2=1/70`, so nonnegativity follows
from the diagonal atom alone; only strictness sees an off-diagonal collision.

Across all sixteen stored `(A_6,A_8,A_9,A_10,m)` cases, every active
level-two shell is nonnegative.  In fifteen cases this is automatic from
`h_2>=lambda_2`; the sole nonautomatic case is `A_10,m=5`, at density exactly
`1/2`, outside the target compact window in `(1/2,1)`.  Therefore this atlas
provides no asymptotic signing-specific evidence for shell positivity.

The same exact order-ten minimizer gives sharp local-exchange warnings:

- its two globally row-optimal cap-ten centers have five favorable selectors
  each, with pair intersections only `2,4,6` and no one-exchange pair;
- all forty one-coordinate outside flips of a row-ten incidence have row
  `42`, and each supplies exactly one adjacent favorable selector to the
  seed;
- the maximal-selector value at `m=6` is `q_*=22`; its 160 ground incidences
  occur only at rows `106,138`; and all 3,840 audited maximal ports have
  strict `|t_{ji}|<r_i`.

Thus the equality-saturation implementation of exact exchange is false even
on an exact minimizer, while the positive-threshold theorem remains valid
(its audited minimum assigned star is 12 at `H=44/3`).  A surviving
level-two shell proof must obtain genuinely global off-diagonal selector
collisions, not only a one-flip star or a sublinear exchange ball.
