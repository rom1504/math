# Independent audit: energy projection with the coherent return retained

Date: 2026-09-06. Source theorem:
`continued_feedback_coherent_return_energy_projection_2026_09_06.md`.
This audit reconstructs its exact Walsh kernels, collision bookkeeping,
pairing classification, and bounded approximation. The normalized-nuclear
covariance counterexample in the Steiner audit is kept as a scope check.

## 1. Verified statement and its narrower conclusion

Let `B=A/sqrt(n-1)` be a bounded-op hollow signing, `Q=B²`, and S fresh
Boolean seed. For bounded odd Gaussian-a.e.-continuous scalar f, write
`b=E[Nf(N)]`, `r=f-bz`, `G=BS`, `V=QS`, `Z=Br(G)`, and

```math
T=B\left[\sum_{p\ge3\ {m odd}}f_p^2Q^{\circ p}\right]B,
\qquad\sigma_i^2=T_{ii}.
```

For bounded even Gaussian-a.e.-continuous H and odd C² psi with psi, psi',
psi'' bounded, define

```math
C_i=H(G_i)\psi(bV_i+Z_i),
\quad c_i^0=H(G_i)\mathbb E_N\psi(bV_i+\sigma_iN),
\quad a_i=\mathbb E_{S,N}[H(G_i)\psi'(bV_i+\sigma_iN)].
```

The proposed conclusion reconstructs:

```math
\frac{\mathbb E C^{\mathsf T}BC}{2n}
=\frac{\mathbb E(c^0)^{\mathsf T}Bc^0}{2n}
 +\frac{\mathbb E(c^0)^{\mathsf T}B D_aZ}{n}
 +\frac{\operatorname{Tr}(B D_aTD_a)}{2n}+o(1). \tag{1}
```

The one-root tested comparison in the source also follows. The theorem
identifies energy, not response L2 or normalized-nuclear covariance.
The residual noise Hermite levels generally retain nonzero variance.
The coherent field QS, its probability law, and the cross term in (1)
are all kept literally. This is why the actual Steiner nonzero-first
covariance falsifier does not defeat (1).

## 2. Exact coherent Walsh root maps

For a monomial in finitely many bounded-op coherent linear fields
`M_j S`, group its input positions by repeated seed labels. An odd block
leaves one marked spin; an even block disappears because S_a²=1.
Inclusion-exclusion over the disappeared labels expresses each fixed Walsh
degree as a finite sum of products of marked row vectors

```math
N_{ia}=\prod_{j=1}^t (M_j)_{ia},\qquad t\ge1,
```

times row scalars obtained by summing products of at least two entries.
Keep the distinct-label projection on the remaining marked slots exactly.
Only finitely many equality partitions occur at fixed polynomial degree.

For t=1 the matrix N is bounded-op by assumption. For t>=2, its absolute
row and column sums are bounded by Cauchy--Schwarz on two factors and
bounded maximum entries on the rest. Every disappeared even block is
bounded by the same argument. Product-row tensor Gram matrices are Schur
products of bounded-op Gram matrices with bounded diagonal; their root
maps are therefore bounded. The common distinct-label projection cannot
increase the root-map operator norm.

This proves the bound for EVERY POSITIVE Walsh degree. The degree-zero
row mean is separated explicitly: a constant vector has uncentered
covariance J and is not a bounded-op root map. This distinction is needed
for the coefficient a_i in (1).

No internal coherent collision is called small. In particular the diagonal
entries Q_ii=1, and any finite large off-diagonal coherent influences,
remain inside the exact Boolean Walsh algebra.

## 3. Full nonlinear contractions and the delicate exclusions

For one fixed odd p>=3, the transported nonlinear kernel is
`K_(p,i)=sum_k B_ik b_k^(tensor p)`. Its proper flattenings are
`O(n^-1/2)`, its row Hilbert norms are bounded, and its maximum influence
is `O(1/n)`, by the independently audited zero-first theorem.

A proper contraction into a coherent Walsh kernel is small by the
flattening and coherent Hilbert bounds. For all p slots contracted into
an unrestricted product of coherent rows, the scalar factor is

```math
\sum_k B_{ik}\prod_{\ell=1}^p(BN_\ell^{\mathsf T})_{ki}. \tag{2}
```

Every `BN_l^T` is bounded-op. Pull out `max_k|B_ik|=(n-1)^-1/2`, retain
two column vectors in Cauchy--Schwarz, and bound the other factors
entrywise. Thus (2) is `O(n^-1/2)`. Any uncontracted coherent row tensor
has bounded Hilbert norm.

Here is an explicit check that coherent distinctness does not break this
factorization. Nonlinear marked slots are already distinct. In a full
contraction, the remaining restriction is that a contracted nonlinear
label a_s not equal a residual coherent label b_t. Fixing a_s=b_t gives a
slice of the nonlinear kernel of squared norm at most its maximum
influence `O(1/n)`. Contract the other slots in Hilbert norm and sum the
remaining bounded row-factor squares. The error kernel has squared norm
`O(1/n)`. A finite union over such equalities suffices; restrictions among
residual coherent labels remain exact. This also covers partially
contracted kernels using the proper-flattening bound.

Thus products of noise Hermites and coherent polynomials can be replaced
by squarefree forests after deleting only collisions TOUCHING noise.
One never deletes an internal coherent collision at a claimed small cost.
Local nonzero noise Hermites paired with coherent Walsh polynomials have
a small proper contraction or a full star of type (2). Their mixed moments
vanish, proving Gaussian-noise factorization against the ACTUAL coherent
law without Lindeberg-replacing QS.

## 4. Cross-root exclusions and exhaustive star classification

In a covariance pairing, a coherent distinctness exclusion involving a
slot paired to a nonlinear branch pulls back to a collision touching that
branch in the opposite forest. If the other slot is also nonlinear-paired,
it is a noise/noise collision there; otherwise it is a noise/coherent
collision. Both have the local small Hilbert cost proved above. The only
restrictions left are between coherent-to-coherent bridge slots, retained
inside their exact bounded root-map covariance. Hence unrestricted star
factorization is legitimate after the charged error.

Only transported branches can split: coherent product-row factors have
one marked slot each. A non-star component contains a four-vertex path
and gives two disjoint proper-flattening gains. Two stars also give two
gains. Such matrices have O(1) Frobenius norm and negligible normalized
nuclear norm.

A single star with s transported leaves and t coherent leaves has matrix

```math
BW,\qquad W_{kj}=\prod_{u=1}^s(Q^{\circ q_u}B)_{kj}
                     \prod_{v=1}^t(BN_v^{\mathsf T})_{kj}. \tag{3}
```

Its total leaf count is odd and at least three. Each transported-leaf
factor is entrywise O(n^-1/2) and O(sqrt(n)) in Frobenius norm. Each
coherent factor is bounded-op, with bounded row and column Euclidean
norms. The cases are exhaustive:

- s>=3: two small entry factors and a third Frobenius factor give
  `||W||F=O(n^-1/2)`.
- s=2: t>=1, and the third factor can be coherent; the same bound holds.
- s=1: t>=2, so two coherent factors have bounded absolute row and
  column product sums; `||W||op=O(n^-1/2)`.
- s=0: the flat leading B and column Cauchy--Schwarz give only an
  entrywise O(n^-1/2) bound. Any remaining whole component supplies
  a bounded-op covariance bridge of Frobenius norm O(sqrt(n)), making
  their Schur product O(1) in Frobenius norm.

Thus the sole potentially nonnegligible partial contribution is a BARE
star with exactly one nonlinear branch against a purely coherent
polynomial. It is precisely the cross term between D_a Z and c^0 in (1).
A nonconstant coherent coefficient beside that nonlinear branch would
require another bridge and has already been shown negligible.

## 5. Whole-pairing terms and factors of two

The wholly coherent covariance is exactly that of c^0 and is retained.
One whole nonlinear pair with no additional component gives D_a T D_a.
Every other whole-pairing term contains at least two bounded-op covariance
factors, one of which is nonlinear. Frobenius Cauchy--Schwarz and the
flat entry magnitude give

```math
\frac1n\sum_{i\ne j}|B_{ij}||\Gamma_{1,ij}\Gamma_{2,ij}|
=O(n^{-1/2}).
```

The two orientations of the retained bare-star covariance contribute
`E(c0)^T B D_a Z/n` to normalized HALF-energy; the two self terms each
retain denominator 2n. This checks every factor in (1).

## 6. Bounded closure without replacing coherent fourth moments

For fixed polynomial residual r, coherent fields G,QS are uniformly
subgaussian; their laws may depend on i and n. The local moment argument
above compares the noise to an independent Gaussian with its actual
deterministic variance profile, not the coherent fields to Gaussians.
Approximate bounded smooth responses and their first noise derivatives
on compact sets by trigonometric polynomials, then by fixed Taylor
polynomials. Uniform subgaussian moments control Taylor tails in the
comparison model; take the matrix limit first for each polynomial in
the actual model. This avoids assuming a finite-n subgaussian tail for
an arbitrary high-degree noise chaos.

H can be approximated separately in Gaussian L2 since each flat row G_i
has the same standardized Rademacher-sum law and converges to normal.
The boundedness of psi and psi' controls this error. Variance bins
permit fixed row-dependent polynomial coefficients without inverse
variance factors. At variance zero use the function and derivative
literally.

Only after H and the outer response are bounded and Lipschitz is the raw
channel Br(G) substituted for the squarefree channel in averaged L2.
For general bounded f, use fixed finite Hermite sums preserving b exactly.
The residual r may be unbounded through its subtracted linear term; it
has at most linear growth, and the scalar L2 approximation and moment
uniform integrability still apply.

The comparison Gaussians for two residual approximants are coupled via
independent degree-p channels. In particular bounded psi'' gives

```math
|a_i-a_i'|\le C|\sigma_i-\sigma_i'|,
```

and the averaged squared variance-profile difference is controlled by
the residual coefficient L2 difference. The normalized-nuclear raw
noise covariance estimate gives averaged absolute agreement of raw
variances and T_ii. Since the coefficients a_i are uniformly bounded,
this converts coefficient errors into averaged L2 errors in D_a Z.
Thus all three right-hand terms of (1) pass through the approximation.

For the one-root tested expression, `y psi(y)` has at most linear
growth and a difference bound controlled by `|y-y'|(1+|y'|)`.
Normalized L2 Cauchy--Schwarz supplies the needed uniform integrability
and raw-channel transfer. No variance floor or small-noise division
is required.

## 7. Consequence and exact surviving gap

With H>=0, `|f|+H<=1`, `|psi|<=1`, and `y psi(y)>=0`, both `+/-f(G)+C`
are feasible. The old self-energy remains `b² Tr(B³)/(2n)+o(1)`.
The exact endpoint identity therefore gives the source's common-energy
lower bound with (1) retained. Its sign is not assumed positive.

The theorem is verified at the scalar, fixed-response, bounded-operator
scope stated above. It removes nonlinear noise Hermite terms from the
ENERGY while retaining all coherent Boolean information. It does not
identify c0's energy in general, remove its exact cross with Z, identify
the joint law of BC, or prove another nonlinear iteration, a new universal
coefficient, convergence, or nonconvergence.

## 8. Independently audited finite-color extension

The statement in `continued_feedback_colored_seed_energy_projection_2026_09_06.md`
also passes reconstruction. Fix d before order and regard a seed coordinate
as a pair (color,vertex). Each old linear row is an embedded copy of a row
of B in R^(dn). Product-Hermite multiindices with different color counts
are exactly orthogonal, even when their total degrees coincide. A matched
multiindex of total degree p has covariance Q_ij^p, with the usual product
factorials canceled.

Every proper colored flattening splits into finitely many two-group maps
with middle diagonal b_i; the group Grams are Schur powers of Q or zero
when colors do not match. Repeated marked slots are possible only in the
same color, where the scalar repeated-slot formula applies unchanged.
Same vertex labels in different colors are NOT collisions: their input
seeds are independent coordinates.

A matched transported leaf of degree q still contributes Q^(circ q)B.
A coherent leaf is restricted to its matching color block before forming
BN^T. Restriction cannot increase operator norms; all small-star estimates
therefore persist, with constants depending on fixed d. Exact Walsh
bookkeeping takes place on the colored seed coordinates. The bounded
extension uses the product-Gaussian marginal of the d old row sums and
keeps all first coefficients exactly.

The old self-energy has coefficient ||b||² Tr(B³)/(2n): distinct colors
give zero first-chaos covariance, and all higher total degrees vanish in
the flat-B trace. This verifies the colored energy and endpoint formulas.
The extension does not cover shared-seed higher-degree canonical trees,
whose coherent returns are not bounded-op linear maps of the independent
seed coordinates. That remains a separate mathematical obligation.

## 9. Independently audited intentionally Gaussian-seeded certificate

The director's separate algorithm in
`continued_director_gaussian_seed_energy_certificate_2026_09_06.md`
starts with an ACTUAL standard Gaussian vector xi. Thus G=Bxi and V=Qxi
are exactly jointly Gaussian, rather than replacements for Boolean laws.
The same noise-star proof applies, now without Boolean input replacement
or internal Walsh collision bookkeeping.

For these seeds, Gaussian Hermite covariance gives Cov(r(G))=R and
Cov(Z)=BRB=T EXACTLY, including the infinite L2 Hermite sum. The linear
and residual input-chaos cross covariance is also exactly zero. The
proper nonlinear flattenings and bare-star estimates are unchanged.

Assuming the director's bounded first two derivatives of H and psi,
ambient Gaussian integration by parts gives

```math
\mathbb E[\xi_k c_i^0]=u_iB_{ik}+b a_iQ_{ik},\qquad
u_i=\mathbb E[H'(G_i)\psi(bV_i+\sigma_i\eta)].
```

Hence the first GLOBAL Gaussian projection is Mxi with
`M=D_u B+b D_a Q`. No inverse of Cov(G_i,V_i) is used. Singularity of that
local covariance creates no gap.

For a polynomial coherent response, every higher global Gaussian chaos
has degree at least three. Its covariance is a finite sum of products
of at least three entries of Q,B³,B⁴, with bounded deterministic row
coefficients. The corresponding matrix powers have bounded operator
norm and O(sqrt(n)) Frobenius norm. Two factors in Cauchy--Schwarz and
the flat B entries make every higher-chaos normalized energy
O(n^-1/2). Bounded passage follows from L2 approximation and contraction
of the global first-chaos projection. Thus the coherent energy is
`Tr(BMM^T)/(2n)+o(1)`.

Finally `Z=Br(G)` gives the EXACT remaining cross

```math
\mathbb E(c^0)^{\mathsf T}B D_aZ
=\langle BD_aB,K\rangle_F,\qquad K_{ij}=\mathbb E[c_i^0r(G_j)].
```

The covariance matrix of (G_i,V_i,G_j) is exactly

```math
\begin{pmatrix}
1&(B^3)_{ii}&Q_{ij}\\
(B^3)_{ii}&(B^4)_{ii}&(B^3)_{ij}\\
Q_{ij}&(B^3)_{ij}&1
\end{pmatrix}.
```

Its positivity follows from the literal linear representation. Adding
one independent eta produces the stated at-most-four-dimensional
Gaussian integral for each K_ij. The resulting certificate, with all
three terms and the residual cross retained, is verified. It is a new
randomized algorithm and must not be relabeled as the Boolean-seeded
law. Gaussian seeding does not itself eliminate the nontrivial bare-star
cross or supply a universal numerical improvement.

## 10. Independently audited actual Boolean sine certificate

The director's `continued_director_boolean_sine_feedback_certificate_2026_09_06.md`
passes an independent reconstruction. Put f(g)=kappa sin(g), b=kappa
exp(-1/2), H=h, psi(y)=sin(ty). The residual kernel is exactly the
Gaussian-Hermite comparison `R=b²(sinh^circ(Q)-Q)`, with T=BRB. This
does not assert that R is the exact finite Boolean covariance.

The coherent comparison is `c0_i=d_i sin(tb(QS)_i)`, where
`d_i=h exp(-t²T_ii/2)`. The product characteristic function of the
ACTUAL Boolean seed gives its first-Walsh coefficient and its mean
noise derivative exactly. The identity sin(u)sin(v)=(cos(u-v)-cos(u+v))/2
gives K_ij=E[c0_i r(G_j)] as the two cosine products minus b(MB)_ij.
The matrix orientation follows directly from

```math
\mathbb E(c^0)^T B D_a B r(G)
=\sum_{ij}(BD_aB)_{ij}\mathbb E[c_i^0r(G_j)]
=\langle BD_aB,K\rangle_F.
```

Consequently the displayed three-term energy formula uses precisely the
retained theorem plus the exact first-Boolean-Walsh energy reduction.
Differentiating the characteristic product with respect to the outer
frequency t, while keeping T fixed, gives the stated tested gain. This
also follows from one-dimensional Gaussian integration by parts in the
independent comparison noise. Both calculations give the same factor
`b sum_k Q_ik L_ik + t T_ii P_i`.

For 0<tb<pi/2 all cosine factors are positive and each
Q_ik sin(tbQ_ik) is nonnegative; Q_ii=1 makes the gain strictly positive.
Outside this regime its sign cannot be assumed. The source correctly
retains the absolute value in the general statement. Sine does not obey
the optional pointwise sign-preserving condition, but that condition is
not needed after the exact gain is evaluated. Feasibility follows from
kappa+h<=1. The endpoint identity has cross normalization 1/n and common
self-energy normalization 1/(2n), as displayed.

The independent script
`computations/continued_audit_boolean_sine_identities_2026_09_06.py`
enumerates all 64 seeds for the order-six Steiner signing and all 1024
seeds for a fixed random order-ten signing. At tb=0.4, pi/2, and 2.4 it
checks the first-Walsh matrix, residual cross kernel, mean derivative,
cross-energy orientation, and tested gain (the last by 80-node Gaussian
quadrature). All discrepancies are below 2e-15 in the recorded replay.
The tests at zero cosine and at negative gain verify the division-free
product formulas and the need for the absolute value. These are floating
checks of exact comparison identities, not finite-order feedback error
certificates. The original code's actual-feedback simulations are also
correctly labeled asymptotic diagnostics.

The matrix products and finite cosine products cost O(n³). Excluding one
factor by prefix/suffix products remains valid when cosine factors vanish;
division by a full product would not. No full Boolean-seed distribution
is needed for evaluating this explicit matrix functional. None of these
algebraic reductions by itself proves a new universal lower constant or
cross-order convergence.
