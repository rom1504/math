# Energy projection for the first genuine marked history

Date: 2026-09-06. Status: proved at the stated scope, independently
reconstructed by the director and in
`continued_audit_marked_energy_projection_2026_09_06.md`.
This is not a claim about arbitrary-depth feedback.

## 1. Frozen theorem

Let `B=A/sqrt(n−1)` be an arbitrary symmetric hollow signing with fixed
`||B||op<=L`. Put `Q=B²`, take independent Boolean seeds `S`, and set

```math
G=BS,\qquad D_i=S_i h_2(G_i),\qquad Y=BD.
```

Let `f(g,y)` be fixed, bounded, odd under `(g,y)→(−g,−y)`, and
Gaussian-a.e.-continuous. Its Gaussian first projection and residual are

```math
b_0=\mathbb E[N_0f(N_0,N_1)],\quad
b_1=\mathbb E[N_1f(N_0,N_1)],\quad
r=f-b_0h_1(g)-b_1h_1(y),
```

where `N0,N1` are independent standard normals. Define the ACTUAL fields

```math
V=b_0QS+b_1QD,\qquad Z=Br(G,Y),\qquad Bf(G,Y)=V+Z.
```

Write the normalized bivariate Hermite expansion

```math
r=\sum_{p+q\ge3\ {\rm odd}}r_{pq}h_p(g)h_q(y),\qquad
T=B\left[\sum_{p,q}r_{pq}^2Q^{\circ(p+q)}\right]B,
\qquad\sigma_i^2=T_{ii}.
```

Let `H(g,y)` be fixed bounded even and Gaussian-a.e.-continuous. Let
`psi` be fixed odd `C²`, with `psi,psi',psi''` bounded. Set

```math
C_i=H(G_i,Y_i)\psi(V_i+Z_i),
\qquad
c_i^0=H(G_i,Y_i)\mathbb E_N\psi(V_i+\sigma_iN),
\qquad
a_i=\mathbb E_{S,N}[H(G_i,Y_i)\psi'(V_i+\sigma_iN)].
```

The `a_i` are deterministic; every coherent return inside `c0` and `a`
is the LITERAL Boolean-seeded return. The claimed energy identity is

```math
\boxed{
\frac{\mathbb E C^{\mathsf T}BC}{2n}
 =\frac{\mathbb E(c^0)^{\mathsf T}Bc^0}{2n}
 +\frac{\mathbb E(c^0)^{\mathsf T}B D_a Z}{n}
 +\frac{\operatorname{Tr}(B D_a T D_a)}{2n}+o(1).}            (1)
```

This is an ENERGY identity only. It does not assert that `C` is close
in `L²` or covariance trace norm to `c0+D_a Z`. The omitted local
noise Hermite levels can have positive variance. No `Q=I`, spectral
flatness, or Gaussianization of `QD` is used.

The already proved local comparison also gives

```math
\frac1n\sum_i\mathbb E H(G_i,Y_i)(V_i+Z_i)\psi(V_i+Z_i)
 =\frac1n\sum_i\mathbb E_{S,N}
 H(G_i,Y_i)(V_i+\sigma_iN)\psi(V_i+\sigma_iN)+o(1).           (2)
```

## 2. Proof modules and the degree bookkeeping

The proof uses the following independently reconstructed modules.

- `continued_feedback_marked_local_noise_separation_2026_09_06.md`:
  local Gaussian separation of `Z` from the ACTUAL coherent primitives
  `W=(S,G,Y,QS,QD)`, exact source cuts, the degree-three exception,
  and uniform exponential tails of these coherent primitives.
- `continued_feedback_all_global_walsh_cuts_2026_09_06.md`:
  all global cuts of positive Walsh kernels in fixed centered polynomial
  computations are at most a fixed polylogarithm.
- `continued_director_two_factor_flat_transport_2026_09_06.md`:
  flat transport of a product of TWO or more globally controlled rooted
  tensors has every proper local cut `O(n^{-1/2} polylog n)`.
- Section 2 of `continued_feedback_marked_energy_star_working_2026_09_06.md`:
  exact old nonlinear source forests have bounded absolute covariance
  row and column sums; their covariance followed by `B` is entrywise
  `O(n^{-1/2})`.

In this note an individual noise branch is a fixed original-Walsh-degree
component of the transported residual forest. It has degree at least
three. A coherent primitive in `W` has degree at most three. The only
noise branch of degree three is transported `h3(G)`. Its complete cross
with an old `Y` is `O(n^{-1/2})` in operator norm, by

```math
\mathbb E[R_{3,j}D_k]
 =\sqrt3 B_{jk}\left[Q_{jk}^2-
                           \frac{n-2}{(n-1)^2}\right],      (3)
```

followed by the bounded root transports. These facts hold for actual
squarefree Boolean kernels. Different original degrees are orthogonal,
but SAME original degrees from different local Hermite monomials are
not declared orthogonal without a covariance calculation.

For each fixed polynomial stage, write `epsilon_n=n^{-1/2}`. Constants
below can be multiplied by fixed powers of `log(n+1)`; every displayed
vanishing power of `n` dominates these factors.

## 3. Polynomial expansion in the noise coordinate only

First fix polynomial `f,H,psi`. Work with the exact squarefree residual
forest main. Temporarily use its actual local variance `v_i` when
normal-ordering the noise; it is uniformly bounded. Replacing `v_i`
by `T_ii` is postponed until the complete energy identity is proved.
The same argument allows an arbitrary fixed globally odd polynomial
response in `(W,Z)`, not only the initially displayed factored form:
only its noise-Hermite expansion and coefficient parities are used.
This permits the finite-catalog approximation in Section 8.

Expand the response in local variance-Hermites of the total noise:

```math
C=c^0+A(W)Z+R,\qquad
R=\sum_{m\ge2}A_m(W)\operatorname{He}_m(Z;v_i),\qquad
a_i=\mathbb E A_i(W_i).
```

The sum is finite at this stage. The coefficient `c0` is odd in the
coherent variables, `A` is even, and `A_m` has parity `m+1`. Decompose

```math
U=(A-a)Z,\qquad C=c^0+D_a Z+U+R.                            (4)
```

Every positive Walsh component of `A−a` has EVEN degree at least two.
After local normal ordering and exact squarefree source projection,
each term of `U` has original degree at least five and has at least
two rooted factors: a positive coherent Walsh kernel and one noise
branch. Every term of `R` has at least two noise branches. Its total
original degree is greater than three (in fact at least seven here).

The product replacement requires a genuine check. Full contractions
of a noise branch into an aggregated high-degree coherent coefficient
are not automatically proper. Expand that coefficient into its primitive
coherent factors of degree at most three and use the Boolean hypergraph
first-merge estimate from the local-separation proof. Every noise-touching
contraction gains a vanishing factor, including the degree-three full
pair through (3). In a squared product error, one such gain already
gives a vanishing row `L²` error; an `O(sqrt(epsilon_n))` norm rate is
enough. Internal coherent-to-coherent contractions are retained exactly.
Noise-to-noise complete contractions are those subtracted by the
variance-Hermites; the remaining partial or higher-multiplicity noise
collisions have the same vanishing bound.

All product-main replacements therefore cost `o(1)` in normalized
`L²`. Bounded-op root transport preserves that cost. Throughout the
diagram calculation one can work with exact squarefree main kernels;
source distinctness preserves every cut bound and is not silently
restored as a small error when it touches only coherent slots.

## 4. Transported nonlinear products become locally Gaussian test fields

Apply the two-factor flat-transport theorem to each exact product-main
term of `U` or `R`. All its global factors have polylogarithmic cuts.
Consequently each original-degree component of `BU` and `BR` has every
proper local cut `O(epsilon_n polylog n)` and original degree greater
than three.

At a fixed root, any finite family consisting of such transported
components and the original noise components is asymptotically jointly
Gaussian, independently of the actual coherent primitives `W`. Their
FULL mutual covariances are retained: they need not vanish. To check
this moment statement, the same hypergraph first-merge proof applies.
Proper partial matches give a small cut, while isolated complete pairs
are exactly the Gaussian covariance terms. No transported component
of `BU,BR` can form a complete pair with a primitive in `W`, whose
degree is at most three; the original degree-three noise exception is
handled by (3).

Rows of `BU,BR` may have POLYLOGARITHMIC variance, not a uniform constant
variance. This causes no gap in the present polynomial, linearly tested
uses: every fixed moment error is bounded by `epsilon_n` times a fixed
polylogarithm, using their row Hilbert bounds and the small-cut gains.
Their averaged variances are uniformly bounded by the operator cap and
the bounded averaged source second moments. Alternatively truncate
their variances at a fixed `K`, take `n→infinity`, and remove the
discarded linearly tested terms with `O(K^{-1/2})` Cauchy--Schwarz cost.

It follows immediately that

```math
\frac{\mathbb E(c^0)^{\mathsf T}BU}{n}\to0,\qquad
\frac{\mathbb E(c^0)^{\mathsf T}BR}{n}\to0.                  (5)
```

Also

```math
\frac{\mathbb E U^{\mathsf T}BU}{n}\to0,
\qquad \frac{\mathbb E U^{\mathsf T}BR}{n}\to0,\qquad
\frac{\mathbb E R^{\mathsf T}BR}{n}\to0.                    (6)
```

For the first two, the only possible leading local moment is
`E(A−a)` times a Gaussian covariance, and `E(A−a)=0` exactly. For
the last, a Gaussian linear variable has zero expectation against
`He_m(Z;v_i)` for every `m>=2`, even if it is correlated with `Z`.
These are local moment identities after transporting the left factor;
no cross-root covariance closure is being assumed.

The two remaining crosses are `U` or `R` against `D_a Z`. The next
sections handle them at the original source, including all degree aliases.

## 5. Higher-noise terms against the mean linear noise

Write `Z=B P`, where `P` is the old residual source forest main, and
put `M=B D_a B`. Its operator norm is bounded for each fixed polynomial
stage. Thus

```math
\mathbb E R^{\mathsf T}B D_a Z
 =\sum_{i,j}M_{ij}\mathbb E[R_iP_j].                         (7)
```

Each term of `R_i` contains at least TWO noise branches, and `P_j`
is a source forest of old primitive `G,Y` branches, each of degree at
most three. Local source product projection has already removed
noise/coherent collisions at root `i`. Its covariance pairing is
therefore bipartite between the left noise/coherent components and the
right old primitive branches.

There is a vertex-disjoint way to match every left noise branch to a
right old branch that it meets. For any set of `k` left noises, at least
`3k` slots must meet right primitives, each of which has at most three
slots. Thus its neighbor set has size at least `k`; Hall's elementary
matching condition applies. Contract along the matched pairs FIRST.
Each pair gains `O(epsilon_n)`: either the noise cut is proper, or it
is a complete degree-three match with an old `Y`, covered by (3).
The pairs are vertex-disjoint, so at least two gains multiply. All
remaining contractions are Hilbert-contractive.

Consequently `|E R_iP_j|<=epsilon_n² polylog n` for the exact main.
Since every absolute row sum of `M` is at most `sqrt(n)||M||op`,
(7), divided by `n`, tends to zero. This proves

```math
\mathbb E R^{\mathsf T}B D_aZ/n\longrightarrow0.             (8)
```

This argument includes unequal individual noise degrees, their exact
within-degree covariances, and any original-degree equality with the
right source. It is not restricted to a single residual chaos level.

## 6. The centered one-noise coefficient, including degree aliases

Again pull the right noise to its old source:

```math
\mathbb E U^{\mathsf T}B D_aZ
 =\sum_{i,j}M_{ij}\mathbb E[(A_i-a_i)Z_iP_j].                (9)
```

Fix one positive EVEN Walsh component of `A−a`, treat its exact kernel
as ONE coherent vertex, and fix a noise branch. If the total original
degree differs from the right source degree, the covariance is zero.
Otherwise apply the following complete dichotomy to the right old
primitive branches.

If some old branch splits between the noise and the coherent coefficient,
it must be an old `Y` (a `G` has only one slot). Contract this split
`Y` with the coherent coefficient first. Its proper cut gains
`epsilon_n`. The noise must meet some OTHER right old primitive as well:
it has at least three slots but receives fewer than three from the split
`Y`. Contract the noise with that other primitive; its cut is proper
because it still meets the first `Y`, and gives another `epsilon_n`.
These two first merges use disjoint pairs of vertices. The entire entry
is therefore `epsilon_n² polylog n`, which is negligible against `M`
exactly as in Section 5.

If NO old branch splits, the right branches partition wholly into a
noise group and a coefficient group. The noise group has odd local
branch count, because its original degree is odd. If that count is at
least three, its covariance with `Z` is an entrywise-small matrix `S`
by the exact source summability theorem. If the count is one, it is
the complete old `Y`/transported `h3(G)` exception, which also has
entrywise size `O(epsilon_n)`.

The coefficient group has a nonempty EVEN number of old branches,
since the coefficient's Walsh degree is positive and even. Its
cross-covariance with the coherent coefficient has a polylogarithmic
global root-map factorization, hence a matrix `Gamma` with
`||Gamma||F<=sqrt(n) polylog n`. This pairing contribution factors as

```math
S\circ\Gamma,
\qquad \|S\circ\Gamma\|_F\le\operatorname{polylog}(n).       (10)
```

The factorization retains each group's exact source projection. Any
additional distinctness restriction BETWEEN the two groups pulls back
to a local coincidence between the noise and coherent coefficient at
root `i`; restoring just those restrictions has already been charged
as the normalized `L²` product error in Section 3. No internal
coherent projection is discarded.

Finally `||M||F<=sqrt(n)||M||op`; Frobenius Cauchy--Schwarz in (9)
and (10) gives `o(n)`. Together with the split case, this proves

```math
\mathbb E U^{\mathsf T}B D_a Z/n\longrightarrow0.            (11)
```

This is the missing mixed-degree term. In the single original-noise-
degree case it also vanishes by degree orthogonality, but (9)--(11)
explicitly handle every alias `q_other=q_noise+q_coefficient`.

## 7. Assemble the polynomial energy

Expanding (4), equations (5)--(6), (8), and (11) leave only

```math
e_n(C)=e_n(c^0)+\mathbb E(c^0)^{\mathsf T}B D_a Z/n
                          +e_n(D_a Z)+o(1).
```

The last term is
`Tr(B D_a Cov(Z) D_a)/(2n)`. The previously proved normalized nuclear
residual covariance comparison replaces `Cov(Z)` by `T`, because `a`
is a bounded deterministic row sequence at this fixed polynomial stage.
The source normal-ordering variances differ from `T_ii` in averaged
absolute value. Gaussian coupling and polynomial moment bounds first,
and the bounded-function continuity below afterward, replace them in
`c0,a` without asserting an operator bound on the raw covariance error.
This proves (1) for fixed polynomial data in the stated ordered sense.

## 8. Bounded closure: a finite polynomial catalog, not a common Gaussian law

The coherent primitive vectors `W_i` have uniformly bounded exponential
moments near zero. In particular their row laws form a weakly compact
family after closure, with uniform moments of all orders. The independent
Gaussian noise coordinate has bounded variance. Thus the corresponding
coherent-plus-independent-noise product laws have the same compactness
and exponential-moment property. None is assumed to be a common Gaussian
law in its coherent coordinates.

Polynomials are dense in `L²` for each such law. One direct proof is to
take a function orthogonal to all polynomials: its product with the law
is a finite signed measure with a nearby exponential moment by
Cauchy--Schwarz. Its Fourier transform is analytic near zero with all
derivatives zero, hence vanishes identically, so the function is zero.

For a fixed bounded continuous target response, choose an approximating
polynomial separately at each law. The squared approximation error is
continuous in that law, by weak convergence and uniform polynomial
moments. A finite cover of the compact law family therefore produces
a FINITE CATALOG of polynomials that approximates all row laws to any
fixed requested `L²` tolerance. Assign a polynomial from this finite
catalog deterministically at each root. Degrees and coefficients remain
uniformly bounded, so all root-map, global-cut, and diagram estimates
above still apply. Symmetrization respects the required parities.

Gaussian conditioning is an `L²` contraction. Thus approximation of
the full response automatically approximates its zero-noise coefficient
`c0` and its first-noise coefficient in the weighted norm
`E sigma_i² |A−A_approx|²`. The audited local separation transfers this
weighted norm to the actual product with `Z`, after the fixed polynomial
matrix limit. Averaging in the coherent law similarly controls the
deterministic coefficient `a`; its weighted raw-channel norm follows
from the normalized nuclear covariance comparison, with bounded
coefficient approximants fixed BEFORE the matrix limit.

Approximate `H` in the established old bivariate Gaussian marginal
first; boundedness of the feedback and its first derivative controls
this error without any continuity assumption on the joint law of
`(G,Y,QD)`. Apply the finite catalog to the remaining continuous bounded
response. The energy functional is continuous in normalized `L²` under
the fixed operator cap. This extends the polynomial identity to the
stated bounded smooth feedback for a fixed polynomial residual main.
Its averaged `L²` distance from the raw residual then transfers the
result using the Lipschitz bound on `psi`.

Finally approximate bounded regular `f` by fixed odd bivariate Hermite
polynomials, keeping its first coefficients `b0,b1` exact. The old
Gaussian marginal gives averaged `L²` convergence; bounded transport
preserves it. Couple the Gaussian residual comparisons by independent
channels with covariances `B Q^{circ k}B`, changing only their Hermite
coefficients. This controls their variance and response changes without
matrix square-root continuity. In fact the ideal covariance error is
uniformly operator-controlled by the coefficient tail: the Schur bound
`||Q^{circ k}||op<=L²` and coefficient Cauchy--Schwarz give
`||T_r−T_P||op<=C_L ||r−P||_2(||r||_2+||P||_2)`. In the channelwise
coupling the difference variance at EACH root is at most
`L^4 ||r−P||_2²`. Since the coherent first coefficients are kept exact,
the bounds on `psi',psi''` therefore make the changes in `c0,a`
uniformly small. This directly controls deterministic-coefficient
weighted raw-channel crosses using their bounded averaged variance;
no merely averaged coefficient error is multiplied by an uncontrolled
raw row variance. All approximation
parameters are fixed before `n→infinity` and removed afterward.

## 9. Actual feasible-endpoint consequence

Suppose `H>=0`, `|f|+H<=1`, `||psi||infinity<=1`, and `t psi(t)>=0`.
Let `j_n` be the comparison on the right of (2), and `e0_n,k_n,t_n`
the three terms on the right of (1). The old-frame self-energy is

```math
e_n(f(G,Y))=\frac{b_0^2+b_1^2}{2}\frac{\operatorname{Tr}(B^3)}n+o(1).
```

The exact pair of feasible endpoints `f+C` and `−f+C` therefore gives

```math
\Lambda(B)\ge j_n+
 \left|\frac{b_0^2+b_1^2}{2}\frac{\operatorname{Tr}(B^3)}n
                    +e0_n+k_n+t_n\right|-o(1).              (12)
```

This retains actual new feedback energy on a true dependent marked
history. It does not assert a positive uniform common term or improve
the banked universal constant by itself.

When `tau²=||r||_2²>0`, the hard threshold can be obtained in the ordered cutoff form from
`continued_feedback_threshold_without_variance_floor_2026_09_06.md`:
the same odd-Schur residual kernel obeys the exact inverse-variance
bound. Remove the low-variance roots before taking the unbounded
threshold derivative coefficient, keep the raw `Z` cross cutoff until
after `n→infinity`, and then remove the cutoff. The Gaussian trace
can be written cutoff-free. This is not an unbounded diagonal
multiplier applied silently to a trace-norm covariance error.
The degenerate case `tau=0` is included in the smooth theorem, but
this variance argument does not justify a discontinuous threshold
replacement for an asymptotically vanishing raw residual there.

## 10. Scope and falsifier record

The noise mixture and its original-degree aliases are fully retained
in the proof. The coherent return is `b0 QS+b1 QD`, not a fresh Gaussian
and not `b0 S+b1 D`. The proof is specific to the primitive original
degrees `1,3` of the first marked history: the old-source split and
Hall-matching arguments use that degree bound and exception (3).

The numerical program
`computations/continued_feedback_marked_old_tree_projection_2026_09_06.py`
tests (1) and a centered-coefficient cross on actual Steiner and positive
twin-Steiner signings. Those tests motivated and stress-tested the
statement; none is a proof step. No arbitrary-depth closure, spectral
flatness of minimizing signings, improved universal constant, or
convergence/nonconvergence conclusion follows automatically.
