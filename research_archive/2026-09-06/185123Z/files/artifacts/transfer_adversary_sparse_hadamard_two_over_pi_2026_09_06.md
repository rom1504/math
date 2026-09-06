# Sparse Hadamard compression has an actual cap floor of 2/pi

Date: 2026-09-06. Corollary of the independently audited uniform local
moment theorem and finite Gaussian covariance inequality. This is an
actual Boolean-witness lower bound for uniform random restrictions, not
an obstruction only to a certificate, and not a lower bound on the
original minima.

## Theorem

Let `H_D` be any sequence of real symmetric full sign Hadamard matrices,
`H_D^2=D I`. Select an exactly uniform `n`-principal subset, with
`n -> infinity` and `n/D -> 0`, and hollow the resulting signing `A_n`.
Then

```math
\frac{Q(A_n)}{n^{3/2}}\ge\frac2\pi-o_P(1).
```

Precisely, for every fixed `epsilon>0`, the probability that the left
side is below `2/pi-epsilon` tends to zero. The assertion is uniform
over the symmetric Hadamard parents, and requires no condition such as
`n^2/D -> 0`. It makes no statement about exceptional deterministic
selectors or the exact limiting cap.

## 1. Fixed integer polynomials approach the semicircle edge

Put `L=A_n/sqrt n`. Let

```math
V_0(x)=1,\quad V_1(x)=x,\quad
V_{j+1}(x)=xV_j(x)-V_{j-1}(x),\qquad
P_q(x)=\sum_{j=0}^q V_j(x),\qquad R_q=I+P_q(L)^2.
```

Every polynomial has integer coefficients, and `R_q` is positive
semidefinite with every diagonal entry at least one. The `V_j` are the
orthonormal polynomials of the variance-one centered semicircle law.
For a direct check, write `x=2 cos(theta)` and
`V_j(x)=sin((j+1)theta)/sin(theta)`. Its measure is
`(2/pi) sin(theta)^2 dtheta` on `[0,pi]`, so ordinary sine orthogonality
gives `E V_i(S)V_j(S)=1_(i=j)`. The recurrence gives

```math
\mathbb E P_q(S)^2=q+1,\qquad
\mathbb E S P_q(S)^2=2q.                                 (1)
```

The second identity counts the `q` adjacent pairs of indices in the
three-term recurrence, once in each orientation. Since `E S=0`, (1)
implies

```math
v_q:=\mathbb E[1+P_q(S)^2]=q+2,\qquad
w_q:=\mathbb E[S(1+P_q(S)^2)]=2q.                        (2)
```

## 2. The local moment theorem supplies every needed hypothesis

Fix `q` before taking any order limit. The theorem in
`transfer_seed_sparse_hadamard_local_moments_2026_09_06.md` gives
mean-square concentration of each diagonal `(L^k)_ii` around its
semicircle moment, uniformly over parents, for every fixed degree.
Finite linear combinations and Cauchy--Schwarz therefore give

```math
\frac{\operatorname{tr}R_q}{n}\longrightarrow q+2,
\qquad
\frac{\operatorname{tr}(LR_q)}n\longrightarrow2q,
```

and

```math
\frac1n\sum_i\left((R_q)_{ii}-\frac{\operatorname{tr}R_q}{n}\right)^2
\longrightarrow0,\qquad
\frac{\operatorname{tr}(R_q^2)}n=O_P(1).                  (3)
```

All convergence is in probability. Only degrees through `4q` and
`2q+1` occur, so they are fixed finite degrees at this step. In fact the
trace of `R_q^2` converges to its finite semicircle value. No operator
norm bound on the sparse compression, or entrywise control of all
covariance entries, is inferred or needed.

## 3. Gaussian signs give actual Boolean witnesses

Apply the finite inequality in
`transfer_seed_finite_spectral_gaussian_cap_2026_09_06.md` to `R_q`,
with lower diagonal bound `a=1`. Equations (2)--(3) make its two error
terms vanish in probability. It follows, for each fixed `q`, that

```math
\frac{Q(A_n)}{n^{3/2}}
\ge\frac{2q}{\pi(q+2)}-o_P(1).                           (4)
```

This inequality comes from the expectation of the energy of actual
Gaussian-rounded Boolean spins. The cap dominates that expectation,
so no equality between an auxiliary norm and the Boolean cap is being
asserted.

For any desired `epsilon>0`, first choose a finite integer `q` so that
`4/[pi(q+2)]<epsilon/2`. Then (4) puts the probability of falling below
`2/pi-epsilon` at zero in the order limit. This proves the theorem.
There is no growing-degree moment estimate and no exchange of
uncontrolled limits.

## Independent audit record

The local graph proof was checked at its delicate points: parity
cancellation and component retention, the extra degree-two elimination
when the parity graph is nonempty with `b=0`, quotient closed walks in
partition-Mobius inversion, the exact correction `p^c`, and the common-root
double-tree factorization. The finite Gaussian inequality was checked
including diagonal normalization, the arcsine remainder, and the hollow
factor of two. Their exact checkers reproduced 2198 structural graphs,
150 homomorphism identities, 150 Mobius identities, and the rational
cubic constants. This corollary strengthens the previously displayed
single-cubic bound; it does not require an additional probabilistic theorem.

The integer polynomial identities (1)--(2) are independently replayed by
`computations/transfer_adversary_sparse_hadamard_two_over_pi_2026_09_06.py`.

## 4. Verified extension to every bounded-cap parent sequence

The later theorem
`transfer_seed_sparse_bounded_cap_parent_theorem_2026_09_06.md`
extends the conclusion to every hollow signing sequence satisfying
`Q(A_D)=O(D^(3/2))`, with the same uniform selector and simultaneous
`n -> infinity`, `n/D -> 0`. Its proof has been independently audited
in full; the Hadamard assumption above is no longer needed for this
stronger statement. The polynomial argument itself is unchanged.

The primary graph input was checked directly in
[Mingo--Speicher, Theorem 6, pp. 6--7](https://arxiv.org/pdf/0909.4277).
It concerns unrestricted colorings of directed graphs, permits loops
and multiple edges, and uses the product of ordinary edge operator
norms. Its exponent is the leaf weight of the forest of two-edge-connected
components, with weight one for a trivial leaf. An Eulerian component
has no bridge, by cut parity; hence each such component contributes one.
Isolated vertices contribute their free color factor separately. Setting
every edge matrix equal to the full-sign completion `B=A_D+I` therefore
gives

```math
|t_F(B)|\le(C+1)^{|E(F)|}D^{-b(F)}
\quad\hbox{when }\|A_D\|_{op}\le C\sqrt D.
```

This is the exact needed specialization, not an injective-coloring
or independent-matrix assertion. The normalization by `D^(|V(F)|)`
is applied only after the unrestricted graph sum is bounded.

The remaining steps are elementary and were also rechecked. A nonempty
zero-excess parity graph must have a nonloop edge; conditioning all
other colors bounds its average by `beta(B)/D^2<= (C+1)/sqrt D`.
Thus the zero-excess term vanishes. Quotient graphs remain closed-walk
parity graphs, so every nontrivial injectivity correction retains the
factor `(n/D)^h`; the rooted second-moment argument is unchanged.
These facts yield the required local polynomial control for every
fixed normalized operator bound.

For a bounded-cap parent, the diagonal-majorant proof in Section 1 of
`resumed_bound_audit_minimal_proof_2026_09_06.md` was reconstructed
independently: its simultaneous majorant has trace at most
`K_G beta(A_D)<=4K_G Q(A_D)`, so deleting the largest diagonal weights
leaves a `(1-epsilon)D` core with operator norm `O_epsilon(sqrt D)`.
For fixed `epsilon`, the intersection of the uniform sample with this
core has size at least `(1-epsilon-o_P(1))n` and, conditional on its
size, is uniform in the core. Principal cap monotonicity follows by
averaging omitted unbiased spins. These statements preserve the
normalization factor `(1-epsilon)^(3/2)` exactly.

For a desired final error, choose `epsilon>0` small and a finite `q`
large; then take the order limit at those fixed choices. This gives
`2/pi` for arbitrary bounded-cap parents without an unregularized
operator bound. It still concerns typical random selectors only, not
the minimum over selectors or the original minimax sequence.
