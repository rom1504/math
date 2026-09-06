# Independent audit of zero-first masked covariance closure

Date: 2026-09-06. The proposed theorem is in
`continued_feedback_zero_first_masked_covariance_2026_09_06.md`.
This note independently reconstructs the kernel, star, collision, input-law,
approximation, and energy arguments. It records three necessary proof-order
or scope clarifications, with their repairs. No involution assumption is used.

## 1. Audited statement and verdict

Let `B=A/sqrt(n-1)` be a symmetric hollow signing with fixed
`||B||op<=L`, `Q=B²`, and independent Boolean seed `S`. For bounded odd
Gaussian-a.e.-continuous scalar `f` with `E[Nf(N)]=0`, write

```math
f=\sum_{p\ge3,\ p\ {m odd}} f_p h_p,\quad
R_f=\sum_p f_p^2 Q^{\circ p},\quad T=BR_fB.
```

Let `H` be bounded even Gaussian-a.e.-continuous and `psi` bounded odd
Lipschitz. Define the ACTUAL vectors `G=BS`, `Y=Bf(G)`, and
`C_i=H(G_i)psi(Y_i)`. If independent Gaussian vectors `U,Z` have covariance
`Q,T`, respectively, put

```math
K_{ij}=\mathbb E[H(U_i)H(U_j)]\,
       \mathbb E[\psi(Z_i)\psi(Z_j)].
```

The proof reconstructs:

```math
\frac1n\|\mathbb E[CC^{\mathsf T}]-K\|_*\to0. \tag{1}
```

Hence multiplication by `B` on both sides also gives normalized-nuclear
control for the covariance of the actual next field `BC`, and contraction
against `B` identifies its current self-energy. This is a real covariance
closure for this restricted nonlinear family. It is NOT a joint-law theorem
for `BC` and its history, or a license for thresholding `BC` again.

For endpoint gain assertions, additionally require `|psi|<=1` (or the
stronger feasibility condition `|f|+H||psi||infinity<=1`). Mere boundedness
of `psi`, with `|f|+H<=1`, does not put `f(G)+C` in the cube.

## 2. Kernel and collision estimates

Let `b_j` be row j of B, and temporarily use Gaussian seeds. The normalized
pure p-chaos kernel is

```math
K_{p,i}=\sum_j B_{ij} b_j^{\otimes p}.
```

Its root-map Gram matrix is `B Q^(circ p) B`; row Hilbert norms and its
operator norm are bounded. For a split of the p slots into r and p-r,
the flattening is

```math
M_r^{\mathsf T}\operatorname{diag}(b_i)M_{p-r},
\qquad (M_r)_{j,\cdot}=b_j^{\otimes r}.
```

Since `||M_r||op²=||Q^(circ r)||op<=L²`, its norm is at most
`L²/sqrt(n-1)`. In particular every marked-slot influence is `O(1/n)`.
The old scalar branch has degree one and the same influence bound by the
flat entries of B.

Writing `m=n-1`, forcing any two slots to the same label a gives the EXACT
identity

```math
K_{p,i}(a,a,\mathbf c)
=m^{-1}K_{p-2,i}(\mathbf c)
 -m^{-1}B_{ia} b_a^{\otimes(p-2)}(\mathbf c). \tag{2}
```

After summing a, the first term has squared Hilbert norm `O(n/m²)=O(1/n)`;
the second has squared norm at most `m^-2 sum_a B_ia²=O(n^-2)`.
The lower kernel is still bounded, including `p=3`, where it is row i of Q.
Deleting within-branch collisions therefore costs uniformly `O(n^-1/2)`
in row Hilbert norm. A finite union over slot pairs suffices.

For two branches, a cross-branch repeated slot has squared mass bounded by
`sum_a influence_a(K) influence_a(K')`, hence by a maximum influence times
a bounded total influence, or `O(1/n)`. Symmetrization is an orthogonal
projection and costs only the fixed-degree factors. This verifies all
forest collision deletions without assuming flat coefficients for the
transported kernels themselves.

## 3. Local Wick replacement and input-law transfer

Separate the finitely many input degrees p of a polynomial f. Under Gaussian
seeds, distinct p channels are orthogonal by homogeneous chaos degree.
Every proper contraction is small by Section 2; a full contraction of a
smaller branch into a larger one is proper on the latter. The old degree-one
branch has zero covariance with every transported degree p>=3 and small
proper contraction into it. The product formula therefore replaces local
variance-Hermites by their pure forests with uniformly vanishing L2 error.
Use unnormalized variance-Hermites; no lower variance bound is needed.

For an explicit sign-input transfer, first delete all repeated slots from
every branch and forest. Treat these squarefree polynomials as finitely many
multilinear summary coordinates. They have bounded moments under every
Gaussian/sign hybrid law and seed derivative Lp norms `O(n^-1/2)` for every
fixed p. Apply one-coordinate replacement to the SQUARED local replacement
error, a fixed polynomial in these summaries. Since every summary is affine
in the replaced seed, four derivatives give four small derivative factors.
Their expectation is `O(n^-2)` by fixed-degree moment bounds. Signs and
Gaussians match moments through degree three; summing fourth-order errors
over n coordinates costs `O(1/n)`. The local L2 replacement therefore holds
for signs as well.

Covariance of the squarefree forests is exactly the same for Gaussian and
Boolean inputs: every surviving input label occurs once in each copy.
Restore the UNDELETED Gaussian forests at this point, charging their row-L2
deletion errors in normalized nuclear norm. This explicit restoration is
necessary before the unrestricted tensor factorization in Section 4;
injectivity restrictions themselves must not be silently factorized.

The covariance error bound used throughout is

```math
\|\mathbb E[UV^{\mathsf T}]\|_*
\le(\mathbb E\|U\|^2)^{1/2}(\mathbb E\|V\|^2)^{1/2}. \tag{3}
```

## 4. Exhaustive partial-pairing and lone-star audit

Each local variance-Hermite term contains an even number of old degree-one
branches and an odd total number of transported branches, each of odd input
degree at least three. A pairing of two forests is a bipartite branch graph.
A component that is neither a whole edge nor a star contains a path on four
vertices. Its two inner branches are split, hence transported; contracting
the two disjoint outer edges first gives two `O(n^-1/2)` gains. Two separate
nontrivial stars give the same two gains. All these entries are `O(1/n)`,
so their matrices have Frobenius norm `O(1)` and nuclear norm `O(sqrt(n))`.

The only other partial case is one nontrivial star plus whole edges. If its
center is a transported degree-p branch on the left, let s be the number of
transported leaves and t the number of old leaves. Matching all remaining
branches forces

```math
s\text{ odd},\qquad t\text{ even},\qquad s+t\ge3. \tag{4}
```

This excludes a lone cubic-to-three-old-branches star, but does NOT exclude
all stars when several transported input degrees occur. Their correct
formula, before collision deletion and up to fixed normalization factors,
is

```math
BW,\qquad W_{aj}=Q_{aj}^t\prod_{\ell=1}^s
                  (Q^{\circ q_\ell}B)_{aj},\qquad q_\ell\ge3. \tag{5}
```

To reconstruct it, expand the center as `sum_a B_ia b_a^(tensor p)`.
Every old leaf contributes `Q_aj`; every transported degree-q leaf
contributes `sum_b B_jb Q_ab^q`. This is exactly (5); permutations only
change fixed factorial constants.

For `S_q=Q^(circ q)B`, sign-entry flatness gives

```math
|(S_q)_{aj}|\le\frac1{\sqrt{n-1}}\sum_b|Q_{ab}|^q
\le\frac{L^2}{\sqrt{n-1}},\qquad \|S_q\|_F=O(\sqrt n). \tag{6}
```

If s>=3, use two factors in maximum-entry norm and one in Frobenius norm:
`||W||F=O(n^-1/2)`. If s=1, (4) forces t>=2, so both absolute row and column
sums of W are `O(n^-1/2)`, because `sum_j |Q_ij|^t<=L²`. Thus in either
case `||BW||op=O(n^-1/2)`. Remaining whole-edge factors are Gram matrices
with bounded diagonal, so their Schur multiplier norms are bounded. They
preserve this small operator estimate; a right-centered star is its
transpose. Every partial-pairing case is now accounted for.

Whole-branch matches are precisely Gaussian Wick covariance pairings for
independent old and transported-degree Gaussian vectors. Their sum is K.
All symmetrization and normalized-Hermite factorials cancel as in the
ordinary local Gaussian Wick formula. This proves (1) for the squarefree
polynomial model.

## 5. Necessary approximation order for the bounded theorem

First extend the squarefree polynomial-model result to bounded H and
Lipschitz bounded psi. Use a finite partition of the bounded variance range,
polynomial approximation on bins away from zero, and
`|psi(y)|<=Lip(psi)|y|` on small-variance bins. Bounded deterministic
row-dependent coefficients are allowed in Sections 2–4. Fixed-root moments
are bounded in this model, so local convergence transfers the polynomial
approximations. Approximate H by fixed even Gaussian polynomials, also in
L2. Remove polynomial/bin errors only after the matrix limit.

Only NOW substitute the actual raw `Bf(BS)` for the squarefree channel.
An averaged L2 error alone does not pass through an unbounded polynomial
psi. It does pass through the now bounded-H/Lipschitz-psi response.

For completeness the required scalar raw replacement is elementary. For a
single p, let W_p be the squarefree degree-p part of `h_p(G_i)`. Its Gaussian
Hermite normalization makes `E W_p² ->1`. Because h_p has leading term
`z^p/sqrt(p!)`, its top Boolean Fourier component is exactly W_p; hence

```math
\mathbb E[h_p(G_i)W_p]=\mathbb E W_p^2,
\quad
\mathbb E|h_p(G_i)-W_p|^2
=\mathbb E h_p(G_i)^2-\mathbb E W_p^2\to0. \tag{7}
```

This is uniform in i, since every signed row sum has the same distribution
as a sum of n-1 unbiased signs divided by `sqrt(n-1)`. Polynomial moment
counting gives the convergence. Multiplication by the fixed-cap B preserves
averaged L2 convergence. This supplies the raw substitution without a
missing covariance-operator theorem.

Finally approximate bounded zero-first f by finite Hermite sums with degrees
at least three. The actual error is at most `L||f-P||2` in limiting
normalized L2. On the comparison side couple independently, for each p,
Gaussian fields of covariance `B Q^(circ p)B`. Multiplying them by the two
sets of Hermite coefficients gives a normalized mean-square discrepancy at
most `L^4||f-P||2²`. Bounded H, Lipschitz psi, and (3) transfer covariance.
This completes the bounded theorem in fixed, explicitly ordered limits.

## 6. Energy collapse and endpoint normalization

Let `mu=E H(N)`, `sigma_i²=T_ii`, and for positive sigma_i set
`a_i=E[N psi(sigma_i N)]/sigma_i`; set a_i=0 at zero. A vanishing diagonal
of the positive semidefinite T makes its whole row vanish, so that convention
is immaterial. Oddness and Lipschitz continuity give `|a_i|<=Lip(psi)`.

In K's variance-Hermite expansion, every nonconstant even-H term contains
`Q_ij^r`, r>=2; every nonlinear odd-psi term contains `T_ij^s`, s>=3.
Both matrices have bounded operator norm and `O(sqrt(n))` Frobenius norm.
Since `|B_ij|=(n-1)^-1/2` off the diagonal, these terms have normalized trace
`O(n^-1/2)`: use squared-entry sums for pure T terms, and Frobenius
Cauchy--Schwarz for mixed Q/T terms. Thus

```math
\frac{\mathbb E C^{\mathsf T}BC}{2n}
=\frac{\mu^2}{2n}\operatorname{Tr}(B D_a T D_a)+o(1). \tag{8}
```

For bounded-function passage, project the Gaussian approximation onto its
linear local chaos before applying (3). Its weighted coefficient error is
bounded by the original local L2 error. This justifies (8) even when some
sigma_i approach zero, without a termwise infinite expansion.

If H>=0, `|f|+H<=1`, `|psi|<=1`, and `y psi(y)>=0`, the endpoint identity
with `u_+=f(G)+C`, `u_-=-f(G)+C` gives

```math
\Lambda(B)\ge
\mu\frac1n\sum_i\mathbb E[\sigma_i N\psi(\sigma_iN)]
+\frac{\mu^2}{2n}\left|\operatorname{Tr}(B D_a T D_a)\right|-o(1). \tag{9}
```

Indeed the exact cross term is `f(G)^T BC=Y^T C`; the scalar old and
Gaussianized transported channel are locally independent. The old
self-energy is negligible since f has zero first chaos. Both normalized
energies are bounded in absolute value by Lambda(B), and
`max(|a+b|,|a-b|)=|a|+|b|` fixes the factor of two. Hollowness licenses
cube-to-Boolean rounding exactly.

Hard threshold psi=sign requires a separate approximation justification;
a uniform positive variance floor is sufficient. The explicit .8528
variance example only disproves a unit variance floor. It does not by
itself disprove every possible smaller uniform positive floor.

## 7. Audited scope

The repaired theorem and energy formula are verified by this independent
reconstruction. They use actual sign-entry flatness in (2) and (6), not
only generic covariance convergence. A first Gaussian coefficient in f
would introduce the exceptional `Q S` return and degree-one transported
branches, destroying the present lone-star estimate. Multicoordinate old
responses likewise require new branch identities. Neither extension follows
from this proof, and neither convergence nor a higher universal lower
coefficient is established by the theorem alone.

## 8. Subsequent actual-signing vanishing-variance audit

The later construction in
`continued_feedback_apex_vanishing_variance_2026_09_06.md` really does
disprove EVERY pointwise positive variance floor, unlike the earlier .8528
example. Its proof and program were independently read and replayed.

With a Sylvester H_m, let `R=[[1,-1],[-1,1]]`, put the hollowed
`R tensor H_m` in the bulk, and adjoin an apex joined positively to all
2m bulk coordinates. Normalize by `sqrt(2m)`. Before diagonal deletion,
the bulk acts on the antisymmetric twin subspace with norm sqrt(2), while
the apex star acts on the orthogonal apex/symmetric subspace with norm one.
The diagonal correction costs at most `1/sqrt(2m)` in operator norm.

Writing d for the diagonal of H and V=A², the bulk square entries are

```math
V_{(a,u),(b,v)}=2m r_ar_b1_{u=v}
-r_ar_bH_{uv}(d_u+d_v)+1_{a=b,u=v}+1.
```

For u=v, summing all four twin choices gives
`2(2m)³+2(3-2m)³`. For u!=v it gives `28+24d_ud_v`.
Since `sum_u d_u=0`, the full cube sum is
`72m³-80m²+2m=2m(36m²-40m+1)`. The apex row is all ones on the bulk,
and the four powers of the normalization in `B(Q^(circ3))B` give

```math
T_{00}=\frac{36m^2-40m+1}{8m^3}\to0.
```

The integer/Fraction executable
`computations/continued_feedback_variance_falsifier_2026_09_06.py --max-lift 64`
passed, including the exact value `144897/2097152` at m=64.
Only one exceptional root is required to falsify a pointwise floor.
This says nothing yet against a small-variance exceptional-set estimate
in normalized mass, which could still suffice for a hard-threshold
extension of the covariance theorem.
