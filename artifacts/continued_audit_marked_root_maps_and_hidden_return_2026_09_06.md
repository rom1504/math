# Audit of marked coherent root maps and a hidden higher return

Date: 2026-09-06. Independent reconstruction of the feedback agent's
`continued_feedback_marked_coherent_walsh_root_maps_2026_09_06.md` and
`continued_feedback_hidden_third_return_2026_09_06.md`.

## 1. Exact marked covariance and gradients

For G=BS, D_i=S_i h_2(G_i), Y=BD, V=QD, the source's exact covariance
formula Cov(D)=(1-3/m)I+(2/m)Q was independently derived by this audit.
Distinct roots must occur together in the same marked triple, leaving
one common free label; the covariance is 2B_ij²Q_ij. At equal roots the
variance is (E G_i^4-1)/2=1-1/m. This also verifies the source's V
energy Tr(B^5)/(2n)+O_L(1/n): its global Boolean degree remains three,
so the earlier first-Walsh energy theorem for bounded-op LINEAR coherent
fields is inapplicable.

With Delta_a F=(F(S)-F(S^a))/(2S_a), direct expansion gives exactly

```math
J_D=\operatorname{diag}(h_2(G))
+\sqrt2\operatorname{diag}(SG)B
-\sqrt2\operatorname{diag}(S)(B^{\circ2})\operatorname{diag}(S).
```

The diagonal derivative is h2(G_i), since B_ii=0. For a!=i, the last
two terms cancel precisely the S_a dependence of G_i. Thus each entry
is independent of its differentiated seed, as required. The row and
column sums of B^{circ2} are one, giving the claimed logarithmic
gradient-operator moments from the subgaussian maximum of G.

## 2. Polynomial differences and Poincare

For multivariate polynomial P and fields W, the identity
`W(S^a)=W(S)-2S_a Delta_a W` makes the source's finite Taylor formula
exact, including the sign (-1)^(|alpha|+1), factor 2^(|alpha|-1), and
the factorial. Every derivative matrix term is a row-diagonal polynomial
coefficient times a Hadamard product of base gradients and a column sign
diagonal. The Hadamard operator bound is a tensor compression and holds
for rectangular matrices too, so independent seed colors introduce no gap.

The base fields have fixed degree at most three and bounded row second
moments. Hypercontractivity gives Lp maximum O(log(n+1)^(3/2)); the base
gradient maximum has moments O(log(n+1)). A derivative term of order ell
therefore has operator moments bounded by
`log(n+1)^(3(d-ell)/2+ell)`. The largest exponent is (3d-1)/2. Squaring
gives the source's E||J_C||op²=O(log(n+1)^(3d-1)). Hölder uses only
fixed higher moments, so the bound is uniform in matrix order.

Discrete Poincare has exactly the stated normalization:
Var(v^T C)<=sum_a E[Delta_a(v^T C)]². Applying this to every v gives
Cov(C)<=E[J_C J_C^T] in positive semidefinite order. Different Boolean
Walsh degrees are exactly orthogonal across every pair of roots. Thus
every positive-degree covariance is below Cov(C), yielding the claimed
polylogarithmic global root maps. Degree-zero row means are separated.
No internal coherent collision is dropped in this proof.

The source's fixed-depth polynomial-circuit corollary is also valid with
its explicit centering hypothesis before each matrix transport. Induct
on fixed Boolean degree, polylogarithmic gradient moments, and row
moments. Poincare controls transported coordinate variances and centering
removes transported means. Without that hypothesis an actual bounded-op
apex signing sends the constant vector to an apex coordinate of size
sqrt(n); multiplying that coordinate by a fresh own spin then produces
a positive-Walsh root map of size sqrt(n). Thus centering is not cosmetic.

The reproducible script
`computations/continued_audit_marked_boolean_gradients_2026_09_06.py`
enumerates all seeds at n=6 and n=8, checking J_D, the exact cube-gradient
formula for W=G+Y+V, the marked covariance, and matrix Poincare. These
are floating checks of exact algebra; the bounds are proved above.

## 3. The hidden B-cubed return construction

The separate four-block construction also survives reconstruction.
Normalized symmetric Sylvester W and a bent diagonal sign D give flat
orthogonal WDW, with fixed blocks A12=W, A23=DW, A34=WDW. Their two-step
products are flat, but A12 A23 A34=I exactly. Scaling each block by 1/2
therefore inserts I/8 into block (1,4) of B0³.

Every unspecified off-diagonal block has independent two-sided sign
gauges and every diagonal block has an independent symmetric gauge PWP.
Block orthogonality makes each diagonal block of B0² exactly I, and the
four-block operator bound gives ||B0||op<=2. Hollowing and renormalizing
changes the operator by O(m^-1/2), so fixed matrix powers retain the
claimed asymptotic entries.

The short-path random estimates have no missing deterministic large term.
For a path using a unique random off-diagonal block, condition on other
blocks and choose an interior gauge. One neighboring factor is a flat
row/column and the other is orthogonal, so the coefficient-square sum is
O(1/m). If the random block occurs in the first or last position, choose
the gauge on its interior side; the remaining two-block product need
only be orthogonal, not flat. Repeated off-diagonal edges are backtracks
and reduce to a flat block.

For an interior random diagonal block PWP, the centered quadratic-chaos
coefficient-square sum is O(1/m). The repeated-label bias equals an entry
of U diag(W) V and is O(m^-1/2), since diag(W) has that operator norm.
The same endpoint reasoning applies when that block is first or last.
Adjacent repeated diagonal blocks square to the identity. These cases
cover the finite length-two and length-three paths. Bounded-degree
hypercontractivity at moment order proportional to log(m), followed by
a union bound over O(m²) entries, makes all nonexceptional entries o(1).
Every triangle on three distinct block vertices uses a random edge,
which also verifies the vanishing diagonal of B0³ and cubic trace.

Hence the source produces actual hollow bounded-op signings with uniformly
vanishing off-diagonal Q entries and vanishing normalized cubic trace,
yet a positive linear number of off-diagonal B³ entries tend to 1/8.
The integer fixed-path identity was replayed through n=256. This disproves
the Q-only higher-return delocalization shortcut; it does not construct
minimizers or refute a theorem imposing all required higher-return bounds.

## 4. Independently checked marked local-noise separation

The complete local proof in
`continued_feedback_marked_local_noise_separation_2026_09_06.md`
passes after its explicit Gaussian-a.e.-continuity condition on general
bounded f was added. A merely measurable bounded response would not be
justified by Gaussian marginal convergence: it could be supported on the
countable union of the actual lattice inputs and vanish Gaussian-a.e.

The transported mixed G/Y tensor cuts follow from the tensor-power audit,
including the exact distinctness pinching. The moment contraction proof
now retains rather than discards all coherent Boolean hyperedges. First
remove distinctness between equality blocks by finite inclusion-exclusion;
coarser blocks that meet a primitive diagonal-free kernel twice vanish.
All remaining blocks have even cardinality. Merge two primitive tensors
along all common labels, retain labels meeting future vertices, and sum
completed labels. The inequality

```math
\sum_x\|A_xB_x^T\|_F^2
\le\sum_x\|A_x\|_F^2\|B_x\|_F^2
\le\|A\|_F^2\|B\|_F^2
```

makes every subsequent merge Hilbert-contractive. A first merge of a
small-cut noise kernel with a lower-degree coherent primitive gains
O(n^-1/2). When shared labels are retained, slicing their coordinates
uses the small influence; completed labels use a proper flattening.
Equal degrees have only the isolated complete-pair exception. It is
handled by the exact formula

```math
\mathbb E[h_3(G_j)D_k]
=\sqrt3B_{jk}\left(Q_{jk}^2-\frac{n-2}{(n-1)^2}\right),
```

whose absolute row and column sums are O(n^-1/2). The raw h3's lower
Boolean component is orthogonal to D, so it does not alter that identity.
This resolves the earlier entangled-coherent-network objection at the
polynomial level. Same original degrees of different residual monomials
are correctly kept with their actual covariance, not falsely declared
orthogonal.

The bounded-test passage has an independently verified special tail
argument. For any row weights w with bounded Euclidean norm, random
three-color decoupling gives

```math
\sum_jw_jD_j
=\frac{27}{\sqrt2}\mathbb E_{\rm color}
\sum_{j\in C_0}w_jS_j
 \left(\sum_{a\in C_1}B_{ja}S_a\right)
 \left(\sum_{b\in C_2}B_{jb}S_b\right).
```

The coefficient is exact: the ordered tetrahedral expansion has factor
1/sqrt(2), and the prescribed colors have probability 1/27. Conditional
on center-color spins, the other two colors form a bilinear Rademacher
chaos with matrix Frobenius norm at most L||w||_2. Its Lp norm is at most
(p-1) times that norm. Thus ||sum w_jD_j||_p<=C_L p||w||_2, not the
generic cubic p^(3/2) bound. It yields nearby exponential moments for
Y,QD and hence a moment-determinate coherent subsequential law. Gaussian
noise marginals plus Cauchy--Schwarz give joint exponential moments of
any limiting law before assuming independence.

Uniform local row moment bounds and averaged source-error L² bounds
justify raw replacement; averaged variance agreement then replaces the
exact squarefree variance by T_ii. These arguments establish the stated
averaged bounded local comparison while retaining literal QS,QD.

There is also a separate self-contained alternative:
`continued_audit_boolean_stable_noise_stein_2026_09_06.md` proves the
stable finite-cube characteristic-function comparison from small mixed
derivative contractions. It permits high coherent influences and does
not require coherent moment determinacy. The director independently
reconstructed that proof. Neither local route supplies cross-root energy
without the additional mixed-star classification.

## 5. Weak regressions and the recovered response gain

The director's `continued_director_marked_return_regression_2026_09_06.md`
also passes. In fixed old Hermite tests S^delta h_p(G)h_q(Y), the leading
source degree is delta+p+3q. Therefore QS has only its S and G matches,
with exact coefficients 1 and gamma_i=(B³)_ii. The only degree-three
tests for QD are D,Y,h3(G). The exact input covariance gives coefficients
1+O(1/n), gamma_i+O(1/n), and the h3 formula above bounds the last by
O(n^-1/2). This proves the displayed weak tested regressions on every
fixed polynomial test. Old Gaussian-marginal L² approximation and the
uniform second moments of QS,QD extend them to bounded continuous tests.

The average of the absolute tested error is important: it allows bounded
root-dependent coefficients in finite test catalogs. Since gamma_i and
sigma_i stay in compact intervals, finite parameter bins justify the
conditional-mean interpretation in each limiting local experiment with
these parameters retained. No finite-n conditional-L² convergence follows.

For a nonnegative feasible old mask H, the verified local noise separation
and weak regressions permit conditional Jensen for the absolute gain.
Its coherent mean is

```math
S[b_0+b_1h_2(G)]+\gamma_i[b_0G+b_1Y].
```

The independent comparison noise adds sigma_i N. Conditional on G,Y,
the S term plus Gaussian noise is symmetric, so adding the deterministic
gamma shift cannot lower its expected absolute value. After removing
that drift, the common function
`F(t)=E H(G,Y) E_N|b0+b1 h2(G)+tN|` is convex and nondecreasing for
t>=0. Jensen and the proved average sigma_i>=tau yield the director's
gain (4). This recovers the earlier scoped two-coordinate creation gain;
it does not improve the universal coefficient or identify feedback energy.
