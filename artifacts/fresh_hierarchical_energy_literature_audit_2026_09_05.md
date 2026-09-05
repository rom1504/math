# Independent audit of hierarchical tree energy

Date: 2026-09-05. Auditor: literature track. This audits the complete files
`fresh_limit_hierarchical_tree_energy_2026_09_05.md` and
`fresh_ownfree_sobolev_transport_audit_2026_09_05.md`, using the independently
proved odd-tree moment theorem in `fresh_tree_chaos_literature_audit_2026_09_05.md`.
No external state-evolution theorem is imported.

## Verdict

The polynomial energy identity, its fixed-family smooth extension, the
Gaussian isometry, and the stable-mask recursion pass this audit. Their
scope is a universal lower-bound certificate for the original minima, not
convergence of those minima and not an actual infinite-depth matrix
iteration. The existing paired-certificate ceiling remains applicable.

## 1. Polynomial identity: critical labels and the bridge

Expand fixed monomials in the fields at roots `i,j`, include the bridge
`A_ij`, and include the explicit spin `S_j`. If the tree copies contain `D`
free positions, the normalization after averaging over roots is
asymptotic to `n^(-(D+3)/2)`. Every nonzero spin expectation has even label
multiplicities. Thus there are at most `(D+1)/2` spin labels and one
additional label for `i`.

The critical case has exactly `(D+3)/2` distinct labels: every spin label
occurs twice, and `i` occurs in no spin position. Lower cases lose at least
one power of `n`. Cases with odd `D+1` vanish identically. The bridge is
nonzero only when `i != j`.

For any critical quotient, reduce all edge multiplicities modulo two. If
a surviving edge has endpoints `u,v`, both labels are summed variables in
this energy calculation. Fixing all other labels leaves `A_uv f(u)g(v)`
with bounded unary factors. Distinct-label exclusions are handled by zero
unary weights, and the exclusion `u=v` is already supplied by `A_uu=0`.
Consequently the contribution is bounded by `beta(A)/n^2`, with a constant
depending only on the fixed monomials. No parity argument for the
remaining graph is needed here.

If no edge survives, the connected quotient has `V=(D+3)/2` vertices and
exactly `2(V-1)` edge occurrences. Connectivity and even multiplicities
force an underlying tree with every edge used exactly twice.

The bridge partner must be a top edge of an `F` tree: its endpoint `i`
cannot be internal in any copy. The first child of this selected copy is
`j`, and its spin is paired with the explicit `S_j`. In particular **no
other free position has label `j`**.

Here is the separation detail needed for the subsequent Wick count. Delete
the bridge from the quotient tree. Each unselected `F` copy stays wholly
on the `i` side, since crossing the bridge would require another internal
occurrence of `j`. Each `H` copy stays wholly on the `j` side, since crossing
would require the forbidden internal label `i`. The selected `F` tree's
child branches remain on the `j` side. Whole-copy propagation therefore
pairs the remaining `F` copies among themselves. At `j`, every selected
child branch pairs with an `H` copy, never with another child branch of the
selected copy, because that copy is injectively embedded. Remaining `H`
copies pair freely.

Selecting the top `F` copy gives `partial_T F`. For each child type `tau`,
the rule forbidding internal pairings of its `m_tau` selected branches is
exactly the Wick polynomial `He_m_tau`. Since

`a(T) = product_tau m_tau! a(tau)^m_tau`,

the normalized coefficient is `1/sqrt(m_tau!)`. This proves precisely

`lim n^-1 E F(X)^T B[S H(X)]
 = sum_T E partial_T F(Z) E h_T(Z) H(Z)`.

For example, selecting the three-free-vertex rooted star against two edge
fields yields `2/sqrt(2)=sqrt(2)`, agreeing with
`E[(Z^2-1) Z^2]/sqrt(2)` and checking the potentially delicate factorial.

## 2. Own-free Sobolev transport

For a fixed `d`-vertex injective tree field, each Walsh coefficient is at
most `C n^(-d/2)`. A fixed spin belongs to `O(n^(d-1))` supports. Hence its
influence has squared `L2` norm `O(1/n)`, and fixed-degree
hypercontractivity supplies the stated `Lp` bounds. Own-spin independence
is exact, not asymptotic.

For `j != k`, own-freeness gives the exact identity

`E S_j S_k g(X_j) g(X_k)
 = E D_k g(X_j) D_j g(X_k)`.

Taylor expansion at the actual endpoint field has remainder `O_L2(1/n)`
for every fixed smooth polynomial-growth function. No independence of
that derivative and the field is asserted or needed. Holder plus the
`L4` influence bound yields the Gaussian `L2 + W^(1,4)` transport estimate.
Its leading constant depends on the fixed field family, not on the
chosen approximation remainder. Constants in the vanishing error may
depend on that remainder, so the order of limits in the artifact is
essential and sufficient.

The provided polynomial-density proof is valid: first smooth cutoff,
then truncated Fourier integral/Riemann sums in Gaussian `W^(1,4)`, then
Taylor approximation of each plane wave dominated by `exp(C|x|)`.
Real parts and parity projection preserve the needed approximation.
On the Gaussian right-hand side, `E partial_T F = E Z_T F`, so continuity
in the first argument needs only Gaussian `L2`, as claimed.

## 3. Gaussian isometry and finite realization

The correspondence between a tree and its even-sized multiset of child
types is bijective. Thus the `h_T` are exactly the normalized finite
Hermite monomials of even total degree, including the constant for the
single-edge tree. They form an orthonormal basis of the globally even
subspace of the countable Gaussian product space. This proves the stated
isometry `U` and its Gaussian-linear image, including for arbitrary `L2`
even masks by completion.

For a smooth finite-coordinate mask `0 <= H <= 1`, enlarge the chosen
coordinate set to be downward closed and use the finite projection
`V_J` of `UH`. The response `F = psi_epsilon(V_J)(1-H)` is bounded smooth
and satisfies `|F| + H <= 1`. The energy identity gives
`E V_J psi_epsilon(V_J)(1-H)` exactly in the dimension limit. Projection
convergence and then sign smoothing give `E |UH|(1-H)`.

For general masks, conditional expectation and finite-dimensional
Gaussian smoothing preserve range and global parity. The continuity
bound follows from

`|J(H)-J(K)| <= ||U(H-K)||2 ||1-H||2
               + ||UK||2 ||H-K||2 <= 2 ||H-K||2`.

Every approximation is fixed before the matrix dimension limit. Nothing
here requires realizing an infinite field family on a finite matrix.

## 4. Stable recursion and its limits

If `v=E f(Z)^2 > 0`, `0 <= f <= 1`, and `D=E f'(Z)^2 < v`, then
`H0=sqrt(v)` is an admissible mask. Inductively `V_t=U H_t` is a centered
Gaussian linear form of variance `v`, and `H_(t+1)=f(V_t/sqrt(v))` has
squared norm `v`. Every finite list of the `V_t` is jointly Gaussian
because its coefficients in the fixed Gaussian coordinate family are
deterministic.

The initial consecutive correlation is `q1=E f/sqrt(v) in [0,1]`.
Thereafter `q_(t+1)=K(q_t)`, with
`K(q)=E f(Z)f(Z')/v`. Squared Hermite coefficients give
`K(1)=1` and `0 <= K'(q) <= D/v < 1` for `0 <= q <= 1`.
Thus `q_t -> 1` geometrically, and for `t >= 1`

`||H_t-f(V_t/sqrt(v))||2^2 = 2v [1-K(q_t)] -> 0`.

Applying Cauchy--Schwarz to the certificate at `H_t` proves the artifact's
formula (12). It does not require that the masks themselves converge or
that an infinite recursion exists on a finite signing. Finite `t`, finite
coordinate approximation, and dimension limit occur in that order.

Finally, this is still the paired, orientation-cancelled energy. The
Bessel bound on the coefficients of `UH` gives the same scalar ceiling
already established for arbitrary finite Gaussian feature families. No
claim of reaching `1/2` or proving convergence is justified by this module
alone.

## 5. Independent audit of the unbounded-response fixed point and certificate

The subsequently added Sections 5–6 also pass. An even intermediate
`g in W^(1,2)(gamma)` with `Eg^2=1`, `Eg>0`, and `Eg'^2<1` is allowed:
it is used only on the Gaussian probability space. Starting from the
single-edge Gaussian, all iterates remain Gaussian linear forms. Their
consecutive correlations obey the squared-Hermite-coefficient kernel
recursion. Since `K'(q)<=Eg'^2<1`, consecutive `L2` distances are
geometrically summable, proving convergence in the closed first Gaussian
chaos to a standard Gaussian `V`. Hermite expansion proves
`g(V_t)->g(V)` in `L2`, so `V=U[g(V)]`. There is no unjustified invocation
of matrix state evolution for the unbounded polynomial.

For `H=1{|V|<=alpha}`, isometry gives `Var(UH)=p` and
`Cov(V,UH)=w=E g(Z)1{|Z|<=alpha}`. Thus write `UH=wV+sZ`, where
`s^2=p-w^2`. Independent integration of the folded-normal expectation
gives the artifact's formula (17): the Gaussian-tail contribution from
the conditional exponential term is
`4s^2/sqrt(p) phi(0) tailPhi(alpha sqrt(p)/s)`, and integration by parts
of the conditional linear term contributes
`4w^2/sqrt(p) phi(0) tailPhi(alpha sqrt(p)/s)`. Their sum is exactly the
stated `4sqrt(p)` coefficient. This checks the normalization directly.

I read both the entire new certificate script and its imported interval
primitive, then independently ran

`.venv/bin/python -B computations/fresh_limit_hierarchical_fixed_point_certificate.py`.

The run succeeded and reproduced

`Eg'^2 in [0.997762540769427182369956738328581413730480359774925265750549,
           0.997762540769427182369956738328581413730480359774925265751941]`,

with `Eg>0.8493298981`, and

`J in [0.426090354752424224090748824426603635740460994016689775818265,
       0.426090354752424224090748824426603635740460994016689775819420]`.

The recurrence computes the probabilists' Hermite polynomials. Its three
rational tail sums are precisely the norm, derivative norm, and covariance
of the degree-200 resolvent polynomial; the common factor `4phi(alpha)^2`
and degree-zero mass terms are correct. Normalization is exact in the
definition of g, so no finite truncation tail is missing.

The Fraction interval primitives round outward, division checks exclusion
of zero, and integer square roots bound both endpoints correctly. Machin's
identity and explicit alternating arctangent remainders enclose pi. The
extended Gaussian integral sums its finite Taylor series exactly at each
rational endpoint. Its omitted integrated terms decrease from that point
onward: for `x<=8`, `u=x^2/2<=32`, while the first omitted index is 257.
The absolute first omitted term therefore bounds the alternating tail.
Endpoint monotonicity of Phi then encloses the entire argument interval.
The actual arguments are below 2.5. No floating-point optimization or
floating-point feasibility check enters this certificate.

The positive margin above `213/500` is consequently certified. The exact
one-parameter optimization in Section 5 is not needed for this numerical
conclusion. Its weighted Cauchy--Schwarz dual and covariance-monotonicity
calculation are also consistent with the stated scope: they optimize only
the stable fixed-point covariance subclass at a fixed threshold.
