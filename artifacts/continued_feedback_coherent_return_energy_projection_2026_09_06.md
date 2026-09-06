# A retained-coherent-return energy projection

Date: 2026-09-06. Status: proved at the stated scope; independent audit
is in `continued_audit_coherent_return_energy_projection_2026_09_06.md`.
The purpose is an actual-energy statement for nonzero-first scalar
responses. The coherent return `QS` is retained literally; it is not
replaced by an independent Gaussian or by the own spin.

## 1. Theorem

Let `B=A/sqrt(n-1)` be an arbitrary symmetric hollow signing with fixed
`||B||op<=L`, put `Q=B²`, and let `S` have independent sign coordinates.
Let `f` be a fixed bounded odd Gaussian-a.e.-continuous scalar function.
Write

```math
G=BS,\quad V=QS,\quad b=\mathbb E[Nf(N)],\quad
r=f-bh_1=\sum_{p\ge3\ {m odd}}f_p h_p,
\quad Z=Br(G),
\quad T=B\left[\sum_{p\ge3\ {m odd}}f_p^2Q^{\circ p}\right]B,
\quad\sigma_i^2=T_{ii}.
```

Let `H` be fixed bounded even and Gaussian-a.e.-continuous, and let
`psi` be fixed odd `C²`, with `psi,psi',psi''` bounded. The actual
mixed response is

```math
C_i=H(G_i)\psi(bV_i+Z_i).
```

Define, using a scalar Gaussian independent of the ACTUAL coherent
variables,

```math
c_i^0=H(G_i)\mathbb E_N\psi(bV_i+\sigma_iN),\qquad
a_i=\mathbb E_{S,N}[H(G_i)\psi'(bV_i+\sigma_iN)],
\qquad D_a=\operatorname{diag}(a_i).
```

The `a_i` are deterministic. The proposed conclusion is

```math
\frac{\mathbb E[C^{\mathsf T}BC]}{2n}
=\frac{\mathbb E[(c^0)^{\mathsf T}Bc^0]}{2n}
 +\frac{\mathbb E[(c^0)^{\mathsf T}B D_a Z]}{n}
 +\frac{\operatorname{Tr}(B D_a T D_a)}{2n}+o(1).                (1)
```

Thus only for the purpose of quadratic energy, the nonlinear noisy
response can be replaced by

```math
c^0+D_aZ.
```

This is NOT an `L²` approximation of the response, and its covariance
need not converge in normalized nuclear norm to that of this simpler
vector. The omitted noise Hermite levels can have positive variance;
their energy trace vanishes. The right side retains the literal
coherent response and its exact cross with `Z`.

The same argument gives the one-root tested comparison

```math
\frac1n\mathbb E\sum_i H(G_i)(bV_i+Z_i)\psi(bV_i+Z_i)
=\frac1n\sum_i\mathbb E_{S,N}
 H(G_i)(bV_i+\sigma_iN)\psi(bV_i+\sigma_iN)+o(1).                (2)
```

No variance floor is assumed. The apex signing in the companion note
shows why such a floor cannot be silently imported.

## 2. Consequence for feasible endpoints

Suppose also `H>=0`, `|f|+H<=1`, `||psi||infinity<=1`, and
`y psi(y)>=0`. Define `j_n` as the right side of (2), and let `e_n^0`,
`k_n`, `t_n` be the three displayed terms on the right of (1). The
known scalar old-self-energy computation gives

```math
\frac{\mathbb E[f(G)^{\mathsf T}Bf(G)]}{2n}
=\frac{b^2}{2}\frac{\operatorname{Tr}(B^3)}n+o(1).
```

The exact two-endpoint identity therefore yields

```math
\Lambda(B)\ge j_n+
 \left|\frac{b^2}{2}\frac{\operatorname{Tr}(B^3)}n
          +e_n^0+k_n+t_n\right|-o(1).                         (3)
```

This retains the common actual feedback energy. It does not assert that
the absolute term is uniformly positive, or that its terms have the
same sign. The new reduction is from a threshold of a transported
nonlinear channel to one static coherent response plus a LINEAR tested
nonlinear channel. It does not close a second arbitrary threshold.

## 3. Coherent Walsh kernels have bounded global root maps

For the polynomial proof, allow finitely many coherent linear fields
`M^{(j)}S`, where every matrix `M^{(j)}` has uniformly bounded operator
norm. Here only `M=B,Q` is needed. All row and column Euclidean norms
are then bounded.

Every fixed-degree coherent polynomial has an exact Boolean Fourier
(Walsh) decomposition. Its kernel in each nonzero input degree is a
finite linear combination of tensors whose marked slots are distinct,
with row factors of the form

```math
N_{ia}=\prod_{j=1}^t M^{(j)}_{ia},\qquad t\ge1,
```

and bounded deterministic row scalars. To see this directly, partition
the input positions of each monomial by their repeated seed label.
Odd-multiplicity blocks remain marked; even blocks disappear after
`S_a²=1`. Inclusion-exclusion removes distinctness restrictions involving
the disappeared blocks, producing bounded row sums of products and
possibly merging blocks. There are finitely many partitions at fixed
degree. Keep the distinctness projection on the remaining marked slots.

For `t=1`, `N` has bounded operator norm by hypothesis. For `t>=2`,
both its absolute row sums and its absolute column sums are bounded:
retain two factors in Cauchy--Schwarz, and bound the remaining entries
by their operator norms. Thus `N` again has bounded operator norm.
The disappeared even blocks have at least two factors, so the same
bound controls every row scalar.

The global root map of a product of row tensors has a Gram matrix
equal to a Schur product of bounded-op Gram matrices with bounded
diagonal. It therefore has bounded operator norm. Marked-slot
distinctness is an orthogonal projection in their common tensor space
and cannot increase that norm. Finite sums preserve the bound.

Consequently EVERY fixed-degree nonconstant coherent Walsh component
has a bounded global root map. Its cross-covariance with another such
component, or with another product row map, is a bounded Schur multiplier
and has Frobenius norm `O(sqrt(n))`.

Internal coherent-slot collisions are NOT declared small. They are
exactly the Walsh reduction just described; high entries of `Q`,
including `Q_ii=1`, are fully retained.

## 4. Local separation of the nonlinear channel from coherent inputs

For a fixed odd input degree `p>=3`, the transported Gaussian kernel is

```math
K_{p,i}=\sum_k B_{ik}b_k^{\otimes p}.
```

Its global root map is bounded, and every proper fixed-root flattening
has norm `O(n^{-1/2})`, by the exact argument in the zero-first note.
Repeated marked slots cost `O(n^{-1/2})` in row Hilbert norm.

Proper contractions into a coherent Walsh kernel are small by this
flattening bound and the coherent row Hilbert bound. There is one
additional case: ALL `p` slots of the nonlinear kernel contract into
a coherent tensor. In a product-row term from Section 3, its coefficient
is a finite sum of expressions containing

```math
\sum_k B_{ik}\prod_{\ell=1}^p (B N_\ell^{\mathsf T})_{ki}.
```

Each matrix `B N_l^T` has bounded operator norm. Since `p>=3`, retain
two column factors in Cauchy--Schwarz, bound all others entrywise, and
use the flat factor `|B_ik|=1/sqrt(n-1)`. The result is `O(n^{-1/2})`.
Uncontracted coherent row factors have bounded Hilbert norm. Exact
marked-slot exclusions can be handled explicitly. An exclusion between
two contracted slots is already enforced by the squarefree nonlinear
kernel. Forbidding a contracted label from equaling a residual coherent
label changes the contraction by `O(n^{-1/2})` in Hilbert norm: fix
that slot, use maximum nonlinear slot influence `O(1/n)`, and sum the
bounded squares of the residual coherent row factors. Restrictions
only among residual coherent slots remain as their exact orthogonal
Walsh projection. Alternatively the finite partition expansion gives
matrices `N_l` of the same bounded-op class.

This proves the full-contraction bound as well, without Gaussianizing
the coherent linear fields. It also proves local moment independence
of the Gaussianized nonlinear channels and the actual coherent fields:
after local noise Hermite subtraction, every remaining mixed diagram
has a small proper contraction or the just-bounded full contraction.

For products of coherent and nonlinear kernels, delete all marked
collisions TOUCHING a nonlinear branch. Their square Hilbert cost is
`O(1/n)`, using the nonlinear maximum influence and bounded total
coherent influence. Keep internal coherent Walsh collisions exactly.
As in the zero-first proof, work on squarefree nonlinear channels
throughout the polynomial stage; transfer raw `Br(G)` only after the
outer response is bounded and Lipschitz.

## 5. Expand only in local NOISE Hermites

For polynomial `H,psi` and a finite Hermite sum `r`, expand the local
response into variance-Hermites of the finitely many nonlinear channels.
Their covariance at a root is diagonal in original input degree.
The coefficients are exact coherent polynomials in `(G_i,V_i)`.

The zero-noise component is precisely `c_i^0`. The coefficient of one
copy of the total channel `Z_i` is

```math
A_i(G_i,V_i)=H(G_i)\mathbb E_N\psi'(bV_i+\sigma_iN),
```

and its constant Walsh component is `a_i=E A_i`. Every other Walsh
component of this coefficient has positive input degree and a bounded
global root map. All higher noise terms have at least two transported
branches. Local noise/coherent contractions are removed at the cost
proved in Section 4.

Gaussian and sign covariance of the resulting squarefree forests
agree exactly. Restore unrestricted nonlinear tensors at their small
row `L²` cost, but retain the coherent Walsh projection.

## 6. Complete covariance-pairing classification for the ENERGY trace

Expand coherent Walsh tensors into their product-row factors. Every
coherent marked factor is a degree-one leaf, whose coefficient matrix
belongs to the bounded-op class of Section 3. Transported branches have
odd degree at least three and are the only split vertices of a branch
pairing graph.

The coherent marked-slot exclusions can be kept entirely inside the
coherent-to-coherent bridge factors. An exclusion touching a label
paired to a nonlinear branch is, on the opposite side, a collision
touching that nonlinear branch; restoring it has already been charged
by the row `L²` deletion estimate. Thus it does not obstruct the star
factorizations below. Exclusions involving only coherent bridge slots
remain exact and preserve their bounded root-map estimates.

There are three partial-pairing possibilities.

1. A non-star connected component contains a four-vertex path. Its
   split inner endpoints are nonlinear branches. Contract its two
   disjoint outer edges first to obtain two `O(n^{-1/2})` gains.
2. Two nontrivial stars give the same two gains.
3. Exactly one nontrivial star, with all other components whole matches.

The first two cases have entrywise `O(1/n)` bounds, Frobenius `O(1)`,
and negligible normalized nuclear norm. Consider the last case. Suppose
its center is a nonlinear branch, with `s` transported leaves and `t`
coherent degree-one leaves. Its total leaf count is odd and at least
three. After summing leaf kernels, its matrix has the form

```math
B W,\qquad
W_{kj}=\prod_{u=1}^{s}(Q^{\circ q_u}B)_{kj}
       \prod_{v=1}^{t}(B N_v^{\mathsf T})_{kj}.                 (4)
```

Every `S_q=Q^{circ q}B`, `q>=3`, is entrywise `O(n^{-1/2})`, has
bounded operator norm, and has Frobenius norm `O(sqrt(n))`.
Every `B N_v^T` has bounded operator norm and bounded row and column
Euclidean norms.

If `s>=3`, retain one factor in Frobenius norm and two in their
entrywise-small norms: `||W||F=O(n^{-1/2})`. If `s=2`, odd leaf count
forces `t>=1`; retain a coherent factor in Frobenius norm and the two
small transported factors. If `s=1`, then `t>=2`; the product of two
coherent factors has bounded absolute row and column sums, so the small
transported factor gives `||W||op=O(n^{-1/2})`. Thus ALL stars with
at least one transported leaf are small in operator norm.

The only remaining case is `s=0`: all leaves are coherent. Its matrix
need not have small operator norm, but the flat leading `B` and column
Cauchy--Schwarz give an ENTRYWISE `O(n^{-1/2})` bound. If any whole
component remains elsewhere in the pairing, its bounded-op covariance
factor has Frobenius norm `O(sqrt(n))`; their Schur product has
Frobenius norm `O(1)` and negligible normalized nuclear norm.

Therefore the only unremoved partial contribution is a BARE star with
one nonlinear branch on one side and only a coherent polynomial on
the other side. This is exactly the cross-covariance of `D_aZ` with
`c^0`. A nonconstant coherent coefficient accompanying the nonlinear
branch would require an additional coherent bridge, which has just
been shown negligible. The bare star is retained, not estimated away.

## 7. Whole-branch terms

The wholly coherent covariance is exactly that of `c^0` and is retained.
One whole nonlinear pair with no other component contributes
`D_a T D_a`. Every other whole-pairing term has at least two genuine
covariance factors, one for a nonlinear pair and one for another
nonlinear pair or a positive-degree coherent bridge.

Each factor has bounded operator norm and bounded diagonal. For any
two such factors `Gamma_1,Gamma_2`, bounded further entrywise factors
do not change the estimate

```math
\frac1n\sum_{i\ne j}|B_{ij}|\,
 |(\Gamma_1)_{ij}(\Gamma_2)_{ij}|=O(n^{-1/2}),
```

by Frobenius Cauchy--Schwarz. Thus every multi-component whole term
has vanishing normalized energy. Combining Sections 6--7 proves (1)
for fixed polynomial responses on squarefree nonlinear channels.

## 8. Ordered bounded-function extension

The coherent fields `(G_i,V_i)` are uniformly subgaussian because they
are bounded-row-norm linear combinations of independent signs. Their
distributions may be non-Gaussian and may depend on the root; no
Gaussian replacement is used for them. The local result of Section 4
identifies all fixed moments of the nonlinear channel with those of
an independent Gaussian of variance `sigma_i²`, uniformly in the root.

First approximate `H` in its universal one-dimensional Gaussian
marginal `L²`; boundedness of `psi` and `psi'` controls that error
without any continuity assumption on the joint law of `(G,V)`.
Approximate `psi` and its first derivative on fixed compact regions
by finite trigonometric polynomials, and then approximate those by
ordinary Taylor polynomials. Uniform subgaussian
moments of the coherent fields and the Gaussianized noise control the
Taylor tails. For the actual squarefree noise, take its matrix limit
for each fixed Taylor polynomial first; its limiting fixed moments
are the same Gaussian moments. Tails of bounded responses are charged
before this polynomial approximation. Variance bins on a bounded
interval permit deterministic root-dependent polynomial coefficients.
At variance zero, use the literal function and derivative at `bV_i`;
no inverse standard deviation is introduced.

This extends (1)--(2) first to bounded `H` and the stated smooth bounded
`psi` on squarefree polynomial nonlinear channels. Their averaged
`L²` distance from raw `Br(G)` then transfers the result by the fixed
operator cap and the Lipschitz bounds. Approximate a general bounded
`f` by fixed finite odd Hermite sums while keeping its first coefficient
`b` exact. The residual approximation error is controlled in averaged
`L²` after multiplication by `B`.

Couple the comparison Gaussians using independent channels with
covariances `B Q^{circ p}B`, and change only their scalar Hermite
coefficients. This controls the changes in `c^0` and in `a` without
matrix square-root continuity. For the weighted raw-channel cross,
use the normalized nuclear covariance theorem for `Z`: its row variances
are close in averaged absolute value to the uniformly bounded `T_ii`.
Consequently averaged squared changes in bounded deterministic `a_i`
also control `D_a Z` in averaged `L²`. This completes the stated
ordered extension, subject to the current independent audit.

## 9. Checks, falsifiers, and limitations

When `b=0`, `c^0=0` by oddness and `a_i=(EH) E psi'(sigma_iN)`.
Equation (1) reduces exactly to the already audited zero-first theorem.

On an involution comparator `Q=I`, the coherent return is the actual
own spin. Then `c^0` is an own-spin times an even old function, and
all three common-energy terms vanish in the established scope. This
does not license substituting that result at arbitrary `Q`.

For the actual Steiner family `Q=I+gamma B`, local coherent return is
`S+gamma G`. The predicted specialization for
`c(s,g,z)=H(g)psi(b(s+gamma g)+z)`, `z~N(0,tau²)`, is

```math
e(C)=\frac\gamma2(a_G^2+\tau^2a_Z^2)+a_Ga_S+a_Za_r+o(1),
```

where `a_G=E partial_g c`, `a_S=E Sc`, `a_Z=E partial_z c`, and
`a_r=E r(G)c`. This specialization itself still requires its coherent
energy/cross derivation; it is recorded as an independent falsification
target, not used in the proof of (1).

The program
`computations/continued_feedback_coherent_energy_projection_2026_09_06.py`
tests (1) directly on actual signings, using `f=sin`, `psi=tanh`, and
optionally the feasible mask `H=1−|sin|`. It reports paired Monte Carlo
errors and explicitly excludes pilot-coefficient uncertainty and
Gaussian-quadrature error from those errors. Numerical agreement is
not a proof step.

No universal lower constant, minimizer spectral flatness, second-query
law, convergence, or nonconvergence conclusion is asserted here.
