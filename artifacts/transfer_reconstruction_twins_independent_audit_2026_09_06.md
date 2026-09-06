# Independent reconstruction of the near-optimal twin precision obstruction

Date: 2026-09-06. The proof in
`transfer_seed_twins_precision_obstruction_2026_09_06.md` survives independent
reconstruction. This note supplies the normalization checks and a separate
integer/rational verifier. It does not promote numerical tests to an
asymptotic proof.

## 1. Row selection: the square-root term and both constants

Let `A` be an actual hollow symmetric sign matrix of order `n`,
`Q(A)=max_x |x^TAx|/2`, and

```math
\mathcal B(A)=\max_x\|Ax\|_1,
\qquad R_J(A)=\max_x\sum_{i\in J}|(Ax)_i|.
```

Polarization with `u=(x+y)/2`, `v=(x-y)/2` gives
`y^TAx=2(P_A(u)-P_A(v))`, where `P_A(z)=z^TAz/2`. Since `A` is hollow,
this polynomial is affine in each coordinate separately. Its absolute
maximum on `[-1,1]^n` is its Boolean maximum. Thus
`B(A)<=4Q(A)`, with the stated convention for `Q`.

For independent Bernoulli selectors of mean `p`, the ghost-sample
symmetrization bound, followed by contraction, is

```math
\mathbb E R_S(A)
\le p\mathcal B(A)
+2\mathbb E\sup_x\sum_i\epsilon_i\eta_i|(Ax)_i|
\le p\mathcal B(A)+2\mathbb E\|A^T(\eta\epsilon)\|_1
\le p\mathcal B(A)+2n\sqrt{pn}.
\tag{1}
```

There is no omitted extra contraction factor: for one Rademacher
coordinate, the averaged supremum is one half of
`max_{x,y}(a_x+a_y+phi(b_x)-phi(b_y))`. A one-Lipschitz `phi` bounds this
by the analogous expression with `|b_x-b_y|`; interchanging `x,y`
identifies that maximum with the identity-function expression. Apply this
one coordinate at a time with `phi=abs`. The final inequality in (1)
uses zero cross terms from the independent signs and
`sum_i A_ij^2=n-1`; hence the square-root scale is `n sqrt(pn)`, not
`n sqrt(p)` or `n^(3/2) p`.

Take `p=4m/n`, with `m>=3` and `4m<=n`. Markov gives
`R_S<=2 E R_S` with probability at least one half. Chebyshev gives
`P(|S|<2m)<=1/m<1/2`, so the events have a common outcome. Restricting
that outcome to any `2m` rows only decreases `R`. Therefore

```math
R_J(A)\le8(m/n)\mathcal B(A)+8n\sqrt m,
\qquad |J|=2m.
\tag{2}
```

The row set exists for every signing; no distributional or near-optimality
assumption was used here.

## 2. Actual sign edits and the cap payment

Partition `J` into pairs `(r_i,s_i)`. Copy each representative row/column
to its twin, except that the newly created twin edge is an arbitrary sign
`d_i`. All off-diagonal entries remain signs and the diagonal stays zero.
This operation is symmetric even between two different twin pairs.

For any Boolean `x`, collapse its two entries to
`y_(r_i)=x_(r_i)+x_(s_i)`, `y_(s_i)=0`, keeping every other entry. The
resulting signing `A'` obeys the exact polynomial identity

```math
P_{A'}(x)=P_A(y)+\sum_i d_i x_{r_i}x_{s_i}.
\tag{3}
```

Although `y` can have entries of magnitude two, its increment `z=y-x`
is supported on `J` and has entries of magnitude at most one. Thus
`|x^TAz|<=R_J(A)`. Also the function
`w -> sum_{i in J}|(Aw)_i|` is convex, so its maximum on `[-1,1]^n`
is attained at Boolean points. Consequently
`|z^TAz|<=R_J(A)`. This is the fractional-input step needed to justify

```math
|P_A(y)-P_A(x)|\le\frac32R_J(A).
```

Combining (2)--(3), including all `m` new twin edges, gives

```math
Q(A')\le Q(A)+48(m/n)Q(A)+12n\sqrt m+m.
\tag{4}
```

For `Q(A)=O(n^(3/2))` and `m=o(n)`, the normalized error is
`O(m/n)+12sqrt(m/n)+m/n^(3/2)=o(1)`. Starting with exact minimizers or
any asymptotically minimizing sequence therefore produces actual signings
with the same asymptotic normalized cap. This conclusion uses the original
universal `O(n^(3/2))` upper bound, not a spectral-flatness assumption.

## 3. Every diagonal completion, with no eigenvector assumption on it

The orthonormal twin differences `v_i=(e_(r_i)-e_(s_i))/sqrt(2)` obey
`A'v_i=-d_i v_i`. For every sign diagonal `D`, not necessarily constant
on any pair, `B=A'+D` satisfies

```math
\|Bv\|\le2\|v\|\qquad(v\in\operatorname{span}\{v_i\}).
```

The subspace need not be invariant under `D` or `B`. The singular-value
min--max principle alone gives at least `m` singular values of `B` at
most two. In fact `B(e_r-e_s)` remains supported on its pair, and its
squared norm is one of `0,4,8`, making this estimate particularly direct.

The Gram matrix `G=B^TB/n` has trace exactly `n`, not `n-1`: its parent
`B` is a full sign matrix. At least `m` of its eigenvalues are at most
`4/n`. Bounding the sum of all remaining eigenvalues by `n` and applying
arithmetic--geometric mean gives, uniformly over every sign diagonal,

```math
\frac1n\log\det G
\le\frac mn\log(4/n)
-\left(1-\frac mn\right)\log\left(1-\frac mn\right).
\tag{5}
```

Singularity means the left side is `-infinity`. For `m=floor(Cn/log n)`
the right side tends to `-C`; for `m=floor(n/sqrt(log n))` it equals
`-sqrt(log n)+o(1)`, whereas (4) has normalized error
`O((log n)^(-1/4))`.

## 4. Precisions and elimination order do not remove this prefactor

For a nonsingular completion, let `M=B/sqrt(n)` and choose any positive
diagonal `Lambda<=tI`. Every diagonal entry of `P=M^T Lambda M` equals
`n^-1 sum_i lambda_i<=t`, because all squared entries of `M` are `1/n`.
Every sequential Cholesky pivot therefore lies in `(0,t]`, in every
coordinate ordering. With

```math
c_t(u)=\frac14\log[(u/t)(2-u/t)],\qquad
\mathcal E=\sum_i c_t(\lambda_i)-\sum_i c_t(\delta_i),
```

the determinant identity is `prod delta=det(G)prod lambda`, not its square
or square root. Hence exactly

```math
\mathcal E=-\frac14\log\det G
+\frac14\sum_i\log\frac{2-\lambda_i/t}{2-\delta_i/t}
\ge-\frac14\log\det G-\frac n4\log2.
\tag{6}
```

This is pointwise in all positive precisions and all coordinate orderings;
allowing a label to choose them does not change the lower bound. A singular
completion has a zero pivot, for which the untruncated positive-precision
full-rank prefactor is infinite; it cannot be evaluated by silently deleting
that pivot.

Equations (4)--(6) prove the stated obstruction. In particular a vanishing
normalized Boolean-cap error is compatible with a divergent normalized
prefactor, uniformly over all diagonal completions, precisions and orders.
The claim concerns this prefactor payment; it is not a claim that every
possible negative quadratic term in every different kernel argument is
incapable of compensating it.

## 5. Reproducible exact checks and logical scope

Run

```sh
.venv/bin/python computations/transfer_reconstruction_twins_exact_audit_2026_09_06.py
```

The independent checker uses only integers and rational arithmetic. It
enumerates 256 row subsets and 6561 signed Bernoulli vectors for the
expectation/contraction bounds. Its order-six, eight and ten twin cases
check every spin up to global reversal, and every one of their 64, 256
and 1024 sign-diagonal completions. Bareiss determinants verify the
exponentiated form of (5). Seventy-two rational precision/order cases
verify every pivot range, determinant cancellation, and the exponentiated
form of (6), without numerical logarithms or eigenvalues. All passed.

The obstruction is universal-seed, not existential-seed. Convergence only
needs one suitable near-optimal sequence; this construction shows that
asymptotic scalar near-optimality alone cannot certify every sequence.
Moreover deleting just the copied vertex `s_i` from each pair makes `A'`
exactly the corresponding principal submatrix of the original seed `A`.
Thus the exhibited obstruction is explicitly removable in `m=o(n)`
vertices. It does not prove that all possible determinant defects admit
such a removal, nor that a truncated kernel argument is valid.
