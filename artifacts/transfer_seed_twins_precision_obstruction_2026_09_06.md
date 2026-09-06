# Near-optimal twin seeds defeat every untruncated precision-determinant bound

Date: 2026-09-06. Seed-transfer track. The construction and inequalities
below are proved. They concern ACTUAL hollow signings that remain
asymptotically near-optimal, not weighted or conference-only models.

The result does not disprove convergence. It excludes deriving vanishing
full-rank Gaussian-elimination prefactor loss from Boolean near-optimality
ALONE, even after optimizing every sign diagonal, precision, and coordinate
ordering. Selecting another near-optimal seed sequence, or discarding
low-dimensional defects in a quantitatively controlled way, remains open.

## 1. A small row set with small bilinear norm

For a real symmetric hollow sign matrix `A` of order `n`, put

```math
\mathcal B(A)=\max_{x,y\in\{\pm1\}^n}|y^TAx|,
\qquad
R_J(A)=\max_{x\in\{\pm1\}^n}\sum_{i\in J}|(Ax)_i|.
```

Polarization on the disjoint supports of `(x+y)/2` and `(x-y)/2`
gives `B(A)<=4Q(A)`. The use of fractional vectors here is exact:
the absolute maximum of a hollow multilinear quadratic polynomial on
`[-1,1]^n` equals its Boolean maximum.

For every integer `m>=3` with `4m<=n`, there exists `J`, `|J|=2m`, with

```math
R_J(A)\le 8\frac mn\mathcal B(A)+8n\sqrt m.                 (1)
```

Here is a proof with constants. Let `eta_i` be independent Bernoulli
variables of mean `p=4m/n`, and set `S={i:eta_i=1}`. Symmetrization and
the Rademacher contraction inequality for the one-Lipschitz function
`u -> |u|` give

```math
\mathbb E R_S(A)
\le p\mathcal B(A)
 +2\mathbb E_{\eta,\epsilon}
           \sup_x\sum_i\epsilon_i\eta_i|(Ax)_i|
\le p\mathcal B(A)
 +2\mathbb E_{\eta,\epsilon}\|A^T(\eta\epsilon)\|_1
\le p\mathcal B(A)+2n\sqrt{pn}.
```

The last estimate is Cauchy--Schwarz for each coordinate and uses
`sum_i A_ij^2=n-1<=n`. Markov gives `R_S<=2 E R_S` with probability at
least one half. Chebyshev gives
`Pr{|S|<2m}<=1/m<1/2`. Thus some `S` satisfies both inequalities.
Choose any `2m` elements of it; row norm is monotone in the row set.
Substituting `p=4m/n` proves (1).

For completeness, the contraction inequality used above requires no
external regularity theorem. Replace one coordinate at a time. With
arbitrary finite numbers `a_x,b_x`, the average over one Rademacher sign is
`(1/2)max_(x,y)[a_x+a_y+phi(b_x)-phi(b_y)]`. The Lipschitz property bounds
this by `(1/2)max_(x,y)[a_x+a_y+|b_x-b_y|]`. Interchanging `x,y` shows
that the latter equals the expression for `phi(b)=b`. Iterating proves
the required contraction. The usual ghost Bernoulli sample and symmetry
of its difference prove the preceding symmetrization bound.

## 2. Duplicate m pairs at a cap-negligible price

Partition the set from (1) into pairs `(r_i,s_i)`, `1<=i<=m`. Let
`pi(s_i)=r_i` and let `pi` fix every other coordinate. Define a hollow
signing `A'` by

```math
A'_{uv}=A_{\pi(u),\pi(v)}\quad\text{if }\pi(u)\ne\pi(v),
\qquad A'_{r_i,s_i}=d_i\in\{\pm1\}.
```

Thus the two members of each pair have identical links to every outside
vertex. They are genuine twin vertices in an actual complete signing.
For a Boolean vector `x`, collapse each twin pair into the representative
coordinate: `y_(r_i)=x_(r_i)+x_(s_i)`, `y_(s_i)=0`, and `y_j=x_j`
elsewhere. Then exactly

```math
P_{A'}(x)=P_A(y)+\sum_i d_i x_{r_i}x_{s_i}.
```

Write `y=x+z`. The vector `z` is supported on `J` and has entries in
`[-1,1]`. Therefore

```math
|P_A(y)-P_A(x)|
\le |x^TAz|+\frac12|z^TAz|\le\frac32R_J(A).
```

Consequently

```math
Q(A')\le Q(A)+12\frac mn\mathcal B(A)+12n\sqrt m+m
\le Q(A)+48\frac mnQ(A)+12n\sqrt m+m.                       (2)
```

In particular, if `Q(A_n)=O(n^(3/2))` and `m_n=o(n)`, this modification
costs `o(n^(3/2))`. Starting from any asymptotically minimizing sequence,
the resulting twin sequence is still asymptotically minimizing.

## 3. Every sign-diagonal completion retains many tiny singular values

The `m` orthonormal twin-difference vectors
`v_i=(e_(r_i)-e_(s_i))/sqrt(2)` satisfy `A'v_i=-d_i v_i`.
Let `D` be ANY sign diagonal and `B=A'+D`, a full symmetric sign matrix.
On the span of the `v_i`,

```math
\|Bv\|\le\|A'v\|+\|Dv\|\le2\|v\|.
```

The singular-value min--max principle gives at least `m` singular values
of `B` at most two. Set `G=B^TB/n`. Then `tr G=n` exactly, because all
`n^2` entries of `B` are signs, and at least `m` eigenvalues of `G` are
at most `4/n`. Arithmetic--geometric mean on the remaining eigenvalues
therefore gives

```math
\frac1n\log\det G
\le \delta\log(4/n)-(1-\delta)\log(1-\delta),
\qquad \delta=m/n.                                        (3)
```

If `G` is singular, the left side is `-infinity` and the statement remains
valid. This bound holds for EVERY diagonal completion, regardless of how
it was chosen from the full seed.

For `m=floor(Cn/log n)`, (3) tends to `-C`. For the more discriminating
choice `m=floor(n/sqrt(log n))`, it tends to `-infinity`, at least as
fast as `-sqrt(log n)+o(1)`. Yet (2) still preserves the normalized cap,
with error `O((log n)^(-1/4))` for bounded original normalized caps.

## 4. Every Gaussian precision choice pays the determinant defect

The Gaussian transport potential has the scalar precision representation

```math
g_t(v)=\sup_{0<\lambda\le t}
           \{c_t(\lambda)-\lambda v\},
\qquad
c_t(\lambda)=\frac14\log\left[
          \frac\lambda t\left(2-\frac\lambda t\right)\right].
```

This follows directly from `lambda=-g_t'(v)=t(1-rho)` and
`2tv=rho/(1-rho^2)`. The endpoint `v=0` has `lambda=t`.

For `M=B/sqrt(n)`, choose ANY diagonal precision
`Lambda=diag(lambda_i)` with `0<lambda_i<=t`. Set `P=M^T Lambda M`,
and let `delta_i` be its sequential Cholesky pivots in ANY coordinate
ordering. Since every entry of `M` has squared magnitude `1/n`,

```math
P_{jj}=\frac1n\sum_i\lambda_i\le t,
\qquad 0<\delta_j\le P_{jj}\le t.
```

Assume first that `B` is nonsingular. Exact Gaussian elimination has the
scalar-prefactor log difference

```math
\mathcal E=\sum_i c_t(\lambda_i)-\sum_i c_t(\delta_i).
```

Determinant cancellation, `prod delta=det(P)=det(G)prod lambda`, gives

```math
\mathcal E
=-\frac14\log\det G
 +\frac14\sum_i\log\frac{2-\lambda_i/t}{2-\delta_i/t}
\ge-\frac14\log\det G-\frac n4\log2.                       (4)
```

Thus neither unequal precisions, dependence of those precisions on a
child label, nor a different elimination order can remove the normalized
determinant payment. The inequality is pointwise in every such choice.
For a singular seed the untruncated full-rank formula already has an
infinite penalty; no finite bound is inferred by pretending a zero pivot
is positive.

Combining (2)--(4), there exist asymptotically minimizing ACTUAL sign
sequences for which this normalized prefactor payment tends to infinity
under every sign-diagonal completion and every positive precision choice.
At `m=floor(n/log n)` the lower bound is already the fixed positive
constant `(1-log2)/4+o(1)`. At `m=floor(2n/log n)` it is at least
`(2-log2)/4+o(1)>1/4`.

## 5. Exact scope and the surviving repair

This is a universal-seed obstruction. The sufficient convergence
criterion only needs ONE well-chosen near-optimal sequence, so the
original problem remains open. The result says that a loss estimate
depending only on how close `Q(A_n)/n^(3/2)` is to the original liminf
cannot control this full-rank prefactor.

The twins occupy only `o(n)` dimensions and are removable by `o(n)`
principal deletions. Consequently a plausible repair is a precision
factorization which is insensitive to such low-dimensional defects,
plus a theorem that every remaining obstruction can be removed or
controlled while preserving near-optimality. One must prove the latter;
the counterexample does not establish it. Replacing the full determinant
by a truncated determinant without redoing the kernel and information
inequalities would not be a valid repair.
