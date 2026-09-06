# Near-minimizers can exclude even approximately isotropic near-ground laws

Date: 2026-09-06. Actual signing theorem. The root proposed the clique
orientation-gap construction; this proof uses the already audited
fixed-half spectral core and needs no random-compression theorem.
It is a counterexample for ALL-near-minimizer isotropy, not a theorem
about exact minimizers or nonconvergence.

Conventions: `H_A(x)=sum_(i<j) A_ij x_i x_j`,
`Q(A)=max_x |H_A(x)|`, `K_G=pi/(2 asinh(1))`.

## 1. Exact finite construction

Let `A` be any order-`N` hollow signing, with `Q(A)=M`. Switch an
absolute maximizing spin to all ones and, if necessary, reverse every
edge sign, so that `H_A(1)=M`.

The simultaneous diagonal-majorant spectral-core lemma supplies a
principal core of size at least `N/2` and operator norm at most
`8K_G M/N`. For any integer `r<=N/2`, choose `r` vertices `S` inside
that core. Writing `u=Q(A_S)`,

```math
u\le\frac r2\|A_S\|_{op}\le4K_G\frac rN M.              (1)
```

The finite majorant/core theorem is proved in Section 1 of
`resumed_bound_audit_minimal_proof_2026_09_06.md` and independently
mapped in `transfer_seed_sparse_bounded_cap_parent_theorem_2026_09_06.md`.
No vertex outside this core is deleted from the final parent.

Overwrite the `S` principal block by an all-positive clique, keeping
every other edge unchanged. Call the resulting hollow signing `A'`,
and set

```math
P=\max_x H_{A'}(x),\qquad R=\max_x[-H_{A'}(x)].
```

The clique energy is
`J_S(x)=((sum_(i in S)x_i)^2-r)/2`, so
`-r/2<=J_S(x)<=r(r-1)/2`. Since
`H_(A')=H_A-H_(A_S)+J_S`, one obtains

```math
P\ge M+\frac{r(r-1)}2-u,\qquad
R\le M+\frac r2+u,                                    (2)
```

and therefore

```math
P-R\ge\frac{r^2}{2}-r-2u
     \ge\frac{r^2}{2}-r-8K_G\frac rN M.                 (3)
```

The first bound in (2) tests the all-ones old maximizer. The second
holds pointwise for every spin; it is not a comparison of separate
random witnesses. The absolute cap also obeys

```math
M+\frac{r(r-1)}2-u\le Q(A')
\le M+\frac{r(r-1)}2+u.                                (4)
```

The upper bound is the exact block-replacement triangle inequality.

## 2. A scalable family of actual near-minimizers

At every order choose an exact minimizer `A_N`, so `M=M_N`, and take
`r=floor(N^(2/3))`. The established upper bound gives
`M_N=O(N^(3/2))`. Hence (1) gives `u=O(r sqrt(N))=o(r^2)`.
Equations (3)--(4) imply

```math
P_N-R_N\ge(1/2-o(1))r^2,
\qquad Q(A'_N)=M_N+O(N^{4/3})=M_N+o(N^{3/2}).           (5)
```

Thus these are actual asymptotically minimizing signings in normalized
cap, and their two oriented extrema are separated by order `N^(4/3)`.
In particular `P_N>R_N` eventually and `Q(A'_N)=P_N`.
The same argument works whenever
`sqrt(N)<<r<<N^(3/4)`; the displayed power is a fixed clean choice.

The construction starts with an already-existing exact minimizer. It
does not compute one, transport a smaller seed upward, or show that
the modified parent remains exactly minimizing. In fact (4) shows
its unnormalized optimality gap is itself of order `r^2`.

## 3. Uniform-slack and mean-slack isotropy are impossible

Let `tau_N=o(r^2)`. By (5), eventually
`tau_N<P_N-R_N`. Every spin in the absolute near-ground shell

```math
|H_(A'_N)(x)|\ge Q(A'_N)-\tau_N
```

then has strictly POSITIVE energy, at least `P_N-tau_N`. A law with
`E xx^T=I` would instead satisfy `E H_(A'_N)(x)=0`, because the matrix
is hollow. Therefore no isotropic law is supported on this shell.
This includes exact grounds and every `o(sqrt(N))`-slack shell.

There is also an exact mean-slack obstruction. For any hollow matrix
whose oriented extrema satisfy `P>R>=0`, the chord bound for the convex
function `|h|` on `[-R,P]` is

```math
|h|\le\frac{2PR}{P+R}+\frac{P-R}{P+R}h.
```

Every isotropic law has mean energy zero. Consequently

```math
E[P-|H|]\ge\frac{P(P-R)}{P+R}\ge\frac{P-R}{2}.          (6)
```

Applied to (5), this forces mean slack at least
`(1/4-o(1))r^2`. Thus allowing a few distant states does not rescue
an isotropic law with mean slack `o(r^2)`.

## 4. Quantitative failure of approximate isotropy on the shell

For ANY probability law supported on the positive shell above, put
`K=E xx^T`. Since `tr A'_N=0` and
`||A'_N||_F=sqrt(N(N-1))`, Frobenius Cauchy--Schwarz gives

```math
\frac{\|K-I\|_F}{\sqrt N}
\ge\frac{2(P_N-\tau_N)}{N\sqrt{N-1}}.                  (7)
```

The operator norm is at least the left side. Writing
`ell=liminf M_N/N^(3/2)`, equations (5) and (7) imply

```math
\liminf\|K-I\|_{op}
\ge\liminf\frac{\|K-I\|_F}{\sqrt N}\ge2\ell.
```

The proved lower bound makes this at least `0.8666442233281614`.
This estimate assumes no bounded operator norm for the clique-modified
parent. Its signing Frobenius norm is enough.

## 5. Exact scope

This falsifies a principle asserting small-slack isotropy for EVERY
asymptotically minimizing sequence, even if approximate isotropy in
operator norm is allowed. It does not exclude such laws for selected
near-minimizers, and does not decide whether exact minimizers at all
sufficiently large orders have an appropriate ground/slack law.

Finite exact-minimizer failures are separately preserved in
`transfer_adversary_exact_minimizer_isotropy_finite_2026_09_06.md`.
Neither those finite failures nor the present near-minimizer family
is promoted to a nonconvergence theorem or an asymptotic exact-minimizer
extension obstruction.
