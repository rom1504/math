# Independent adversarial audit of the rooted lower bound

Date: 2026-09-05. Auditor: fresh-limit literature subagent.

## Verdict and scope

**Audit passed.** I independently reconstructed the transfer estimates in
`fresh_limit_rooted_response_2026_09_05.md`, the four contraction classes in
`fresh_limit_rooted_gaussian_2026_09_05.md`, and the exact Boolean partial-update
inequality. I read the primary Gaussian-chaos theorem statements, checked
their hypotheses, and reran and inspected the exact rational certificate.
No mathematical gap or counterexample was found.

The resulting statement is an asymptotic universal lower bound:

\[
\liminf_{n\to\infty} M_n/n^{3/2}
\ge 0.339649621211865756397462117164309583060688819678602946169260
>0.3396.
\]

This does **not** establish convergence, a finite-order bound with a specified
cutoff, or a Gaussian limit for the nonrooted response field. The last decimal
display is the outward rational lower endpoint for the fixed choices
`t=7/8`, `p=8/125`, not an assertion of optimality for the original problem.

## 1. Elementary bootstrap and transport

Write `m=n-1`, `B=A/sqrt(m)`, `Q=B²`, and assume `q(A)<=C n^(3/2)` with fixed
`C`. The eigenvector argument really gives `||A||op²<=beta(A)`, where
`beta=max_signs |xᵀAy|`; polarization on the cube gives `beta<=4q(A)`.
Consequently `Qii=1`, `Tr Q=n`, `||Q||op=O_C(sqrt(n))`, and
`Tr Q²=O_C(n^(3/2))`. No spectral flatness or bounded operator norm of `B`
has been assumed.

For an even, fixed smooth function `g` with polynomial-growth derivatives,
endpoint-spin extraction gives exactly

\[
E[S_jg(G_j)S_kg(G_k)]
=\tfrac14 E[(g(U+\beta)-g(U-\beta))
             (g(V+\beta)-g(V-\beta))],
\]

where `beta=Bjk` and `(U,V)` is independent of both endpoint spins. Taylor
expansion followed by two-dimensional Lindeberg replacement yields
`m^-1 K_g(Qjk)+O_g(m^-3/2)` off the diagonal. Hybrid linear row sums have
uniformly bounded polynomial moments, so polynomial-growth test derivatives
are sufficient. The missing variance `1/m` in each coordinate can be added
as independent Gaussian noise, costing only `O_g(m^-2)` in this formula;
there is no problem when the correlation is close to one.

Because `g'` is odd, its Gaussian kernel is a nonnegative sum of odd Schur
powers of `Q`. Schur multiplication by a correlation matrix contracts
operator norm. Sandwiching the covariance matrix by the unit row `B_i`
therefore proves, uniformly in `i`,

\[
E\{B[Sg(BS)]\}_i^2=E g(Z)^2+o(1).
\]

In particular this applies to `g=h-P_D`, with any fixed even Hermite
truncation `P_D`. The degree is held fixed before dimension tends to infinity.
One never multiplies the polynomial approximation error by `||B||op`.

## 2. Finite-chaos transfer and the contraction audit

The discrete recurrence

\[
G_jW_{j,r}=W_{j,r+1}+r(1-(r-1)/m)W_{j,r-1}
\]

is exact by counting the colliding multiplied index. Comparing with the
Hermite recurrence gives lower-degree coefficients `O_r(1/m)`; explicitly,
`H2=W2` and `H4=W4-8W2/m-2/m`. The preceding transport estimate bounds the
transported lower terms uniformly in L².

For the resulting homogeneous multilinear polynomial of degree `r+1`, every
coefficient is bounded by `C_r m^(-(r+1)/2)`: a fixed set has at most `r+1`
possible roots. Every coordinate influence is thus `O_r(1/n)` and their sum
is bounded. Joint Lindeberg replacement, using fixed-degree hypercontractivity
on hybrid Gaussian/Rademacher products, has error `O_r(n^-1/2)` against a
fixed smooth test with bounded third derivatives. The same argument includes
the linear field `B_i S` if desired.

After replacement, repeated-index Gaussian Hermite terms have squared norm
`O_r(1/n)`: there are only `O_r(n^r)` such multi-indices, their factorial
weights are degree-bounded, and each has at most `r+1` root contributions.
The full Gaussian polynomial is the pure degree-`r+1` chaos

\[
F_{i,r}=\sum_j B_{ij} Z_j H_r(B_j Z).
\]

Hollowness is essential to this assertion. Its exact variance is

\[
r!\left[1+\frac r m
 (B_i Q^{\circ(r-1)}B_i^\top-1)\right]=r!+o(1).
\]

For `T=sum_j a_j e_j tensor b_j^(tensor r)` and `K=Sym(T)`, I checked each
possible contraction between two permuted copies of `T` independently:

| Contracted roots | Squared-norm upper bound |
| --- | --- |
| Root paired to root | `Tr(Q²)/m²` |
| Neither root contracted | `Tr(Q²)/m²` |
| Exactly one root paired to an opposite branch | `||Q||op/m` |
| Both roots paired to opposite branches | `||Q||op²/m²` |

For the third case the residual branch Gram is
`P=Q^(circ(r-l+1))`, so its norm is at most `||Q||op`; for the fourth,
the coefficient matrix is `Cjk=a_j a_k Bjk² Qjk^(l-2)`, and the exact squared
norm is `Tr(C P Cᵀ P)`. The exponent `r-l+1` is at least one throughout.
All displayed quantities vanish, including for coordinates with
`(Q²)ii` of order `sqrt(n)`. Symmetrization only averages finitely many
terms, each related to these forms by an isometry on the remaining tensor
coordinates. Triangle inequality therefore suffices; no cancellation has
been silently assumed.

## 3. Imported theorem scope and all limiting quantifiers

I read Theorem 1 in
[Nualart--Peccati, *Central limit theorems for sequences of multiple stochastic
integrals*](https://arxiv.org/pdf/math/0503598), Annals of Probability 33
(2005), 177--193. At fixed chaos order, asymptotic unit variance and vanishing
of every nontrivial self-contraction imply convergence to a standard normal.
Embedding all finite coordinate spaces in one `ell²(N)` supplies the common
isonormal Hilbert space in the theorem. Apply it to `K/sqrt(r!)`.

I also read setup (9), Lemma 6, and Theorem 7 in
[Nualart--Ortiz-Latorre, *Central limit theorems for multiple stochastic
integrals and Malliavin calculus*](https://arxiv.org/pdf/math/0703240).
The finite-dimensional statement allows chaos orders starting at one.
Componentwise Gaussian convergence and covariance tending to the identity
give joint standard Gaussian convergence. Here distinct chaos orders are
exactly orthogonal after normalization. In particular the linear field may be
included; no unverified mixed-contraction assumption is needed.

These statements apply to **every** sequence of dimensions, permitted
signings, and selected coordinates. Since all estimates depend only on the
fixed cap `C`, fixed degree, and fixed smooth functions, a contrary sequence
would contradict the same theorem. This proves the claimed uniformity over
coordinates and matrices, without needing a quantitative chaos CLT.

For every fixed `tau>0`, first fix a Hermite degree `D`, let `n` tend to
infinity, and then let `D` tend to infinity. The transported L² tail is
bounded in the limsup by `E(h-P_D)²`, which tends to zero. Thus

\[
\{B[S h(BS)]\}_i\Rightarrow N(0,E h(Z)^2)
\]

uniformly. Its second moments are uniformly bounded by the transport lemma,
so the linearly growing Lipschitz test `max(a,|z|)` has converging
expectations. Only uniform integrability of the first moment is needed:
`E[|Y|;|Y|>K]<=EY²/K`. The weak proof's residual moment-passage and cubic
witness are unnecessary in the stronger proof.

The hard threshold is obtained only **after** these dimension and polynomial
limits: the final liminf is at least the smooth formula for every fixed
`tau>0`, hence at least its limit as `tau` decreases to zero. No estimate
uniform in vanishing dither has been used. Likewise `t` and `p` are fixed,
not dimension-dependent optimizers.

## 4. Own-spin Jensen, energy, and exact Boolean update

With `U_j=G_j-Bji S_i`, Taylor's deterministic remainder after multiplication
by row `B_i` is `O_tau(m^-1/2)`. The coefficient of `S_i` in the nonrooted
field is `m^-1 sum_(j!=i) f'(U_j)`, which tends to `a` in L². Its Gaussian
comparison covariance is bounded by the square of the removed-coordinate
Gram; the sum is `o(m²)` by `Tr Q²=O(n^(3/2))`. Removing the coordinate
changes Gram entries by at most `1/m`.

The corresponding change in the rooted field has leading term
`S_i m^-1 sum_(j!=i) S_j h'(U_j)`. Endpoint extraction bounds its off-diagonal
summand covariances by `||h''||infinity²/m`, so the term vanishes in L².
The two remaining fields are jointly independent of `S_i`. Averaging only
this unreplaced Rademacher spin gives the required convex inequality.
There is no need to prove a joint Gaussian limit with this spin, nor a CLT
for `B f(BS)`.

The one-probe energy can also be derived directly from the present bootstrap,
without importing an older theorem. The paired energy is exactly
`e=E uᵀ Bv`. Extracting the endpoint `S_j` gives

\[
e=\frac1m\sum_{i\ne j} E[f'(G_i-Bij S_j)h(G_j)]+o(n).
\]

Two-dimensional Lindeberg comparison and even-function Gaussian covariance
give `e=ab n+o(n)`: the summed covariance correction is bounded by
`O(Tr Q²/m)=O(sqrt(n))`. A pointwise assertion that the displayed expectation
is at least `ab` at every correlation would be false near correlation one;
the averaged square-Gram estimate, not that assertion, is the valid proof.

For actual Boolean initial probes `X^sigma`, independently overwrite each
coordinate with `Y^sigma=sigma sign(B X^sigma)` with probability `p`.
Hollowness makes the expected new energy exactly its quadratic expansion.
The cross term is the field norm; bounding the new signed energy above by
`q(B)` and the signed `Y` energy below by `-q(B)` gives exactly

\[
(1+p^2)q(B)\ge (1-p)^2 e+p(1-p)\ell.
\]

All factors agree with the convention `q(B)=max|xᵀBx|/2`. Conditional-dither
Jensen only improves `ell`. The conversion from `q(B)/n` to
`q(A)/n^(3/2)` multiplies by `sqrt((n-1)/n)`.

To infer the universal liminf, it suffices to suppose a contradicting
minimizing subsequence lies below, for example, the fixed cap `C=1`.
The preceding uniform proof then applies. No exchange of a minimum with an
infinite supremum, and no external all-order upper construction, is required.

## 5. Counterexample attempts and numerical certificate

The all-positive matrix is a decisive counterexample if the low-cap
assumption is deleted: its cubic rooted field tends to `H3(Z)`, not
`N(0,2)`. Its cap is of order `n²`, and `||Q||op` is of order `n`, so it
does not satisfy the audited hypotheses. This check confirms that the
contraction argument has not accidentally proved an unrestricted universality
claim.

A more serious stress family starts with a symmetric Sylvester Hadamard
matrix of order `N`, deletes its diagonal, and attaches `k=floor(sqrt(N))`
all-positive twin vertices, with all-positive cross edges and twin edges.
It has the rigorous cap

\[
q(A)\le \tfrac12N(\sqrt N+1)+kN+\tfrac12k(k-1)=O(n^{3/2}).
\]

The twin Gram block has correlations tending to one and size `sqrt(n)`;
therefore it tests the worst permitted operator-norm scale, not just a
conference-like bounded-spectrum case. All contraction bounds still vanish.
An additional, nonproof Monte Carlo check used seed `88021`, 6144 samples
per distribution/order, and the cubic response `h=H2`. At `N=4096`, the
Rademacher twin-root second/fourth moments were `2.0526, 12.8414`, and the
Gaussian-input values were `2.0795, 13.0915`; the target is `2,12`.
The ordinary-core values were `2.0298,12.2599` and `1.9887,12.2213`.
These are diagnostics, not a substitute for the uniform contraction proof.

Finally I inspected and ran
`computations/fresh_limit_rooted_lower_certificate.py`. Its interval
operations enclose exact real operations. The integer-square-root endpoints
are outward bounds; the alternating arctangent remainder, exponential
Lagrange remainder, and integrated alternating Gaussian-series remainder
are all valid for the asserted argument ranges. Machin's identity is exact
(the tangent subtraction gives one, with the angle in the first quadrant).
The program uses no floating-point arithmetic, and reproduces the displayed
rational interval. It certifies the fixed formula, not any global numerical
optimization. No repair to the mathematical argument or certificate was
needed in this audit.
