# From diffuse row spectra to a global stability loss in the actual weave

2026-09-07. **PASS, including the director's convex-product row lemma and
random-selector sparse-spectrum counting lemma.** The independent
reconstruction is in
`principle_invent_2026_09_07_stability_selector_independent_audit.md`.
This is a strict positive symbolic improvement of an actual all-order
upper certificate, not a changed displayed decimal or convergence.

## 1. Deterministic row-to-column thermal-variance transfer

Let u=(u_i(j)) be a real m by m array with

    sum_j u_i(j)^2=m for every i.                          (1)

Fix t>0, 0<delta,epsilon<1 and V>=1. A row i is diffuse if

    sum_(j: |u_i(j)|<=V) u_i(j)^2 >=epsilon m.             (2)

Suppose at least delta m rows are diffuse, and write d=delta epsilon.
Choose fixed constants

    U^2>=4 V^2/d,      C>=4 V^2/d,
    eta=d/(8V^2),      kappa=(d/8) sech^2(2t U V).          (3)

For all sufficiently large m, at least eta m columns i simultaneously
satisfy

    q_i:=sum_j u_j(i)^2 <=C m,
    sum_(j!=i, |u_j(i)|<=V)
       u_j(i)^2 sech^2(2t u_i(j)u_j(i)) >=kappa m.         (4)

Proof. The total mass on directed entries of magnitude at most V is at
least d m^2, by (2). Delete those entries whose reverse entry has magnitude
greater than U. There are at most m^2/U^2 reverse-heavy entries, since
the total squared mass is m^2, and each deleted light entry has weight
at most V^2. This loses at most V^2 m^2/U^2<=d m^2/4.

There are at most m/C columns with q_i>Cm. Each column has at most V^2m
light mass, so deleting them loses at most V^2m^2/C<=d m^2/4. Finally,
deleting all diagonal entries loses at most V^2m. The remaining light
mass is at least 3dm^2/8 for m>=8V^2/d. A column carries at most V^2m
such mass. Columns carrying less than dm/8 account for at most dm^2/8;
therefore at least dm/(4V^2), and in particular eta m, columns carry
at least dm/8. On these remaining entries both magnitudes are bounded
by V and U in the required directions, so the sech-squared factor is
at least sech^2(2tUV). This proves (4).

This statement is uniform over every row permutation, every choice of
bases, every Boolean configuration and both energy orientations. No
random-matching approximation or independence of the entries of u is
used. The diffuse-row classification itself is invariant under arbitrary
coordinate permutations within rows.

## 2. Consequence of the convex row stability lemma

The exact stability event and tilted signs are defined in
`principle_invent_2026_09_07_stability_aware_weave.md`.
Assume the director's row lemma in the following precise form: for fixed
p=k/m in a compact subset of (0,1], C,V,kappa>0, there is c_row>0 such
that (4) implies

    P_(t,u,sigma)(A_i)<=exp(-c_row k),                     (5)

uniformly over all selected Hadamard rows and all physical words.
Heavy incoming signs may be conditioned arbitrarily; the physical
diagonal is retained in (5). This is an actual finite-dimensional
orthant statement, not Gaussian replacement at exponential accuracy.

By the exact tilted-edge Finner bound, every configuration with at
least delta m diffuse rows gains the uniform factor

    product_i P(A_i)^(1/2)
       <=exp[-c_row eta k m/2]
       =exp[-Delta_good m^2],
    Delta_good=c_row eta p/2>0.                           (6)

Because (6) holds before any basis/permutation average, it can be pulled
out and the remaining integrand bounded by the original graph-CS/orbit
moment theorem. This is important: no decoupling of the p_i themselves
under the basis average is required.

## 3. Why random physical selectors dispose of the complementary region

Choose each retained physical selector T_i independently and uniformly
among the k-subsets of [m], independent of its basis. This is allowed in
the original actual-sign construction. Conditional on any orthogonal
basis, averaging T and summing its 2^k row words is constant-type ternary
counting divided by binom(m,k):

    (1/binom(m,k)) sum_(y: #zeros=m-k,
                       y_a in {0,+1/sqrt(p),-1/sqrt(p)}) .

A non-diffuse row has all but epsilon m energy in coordinates of
magnitude greater than V. There are at most m/V^2 such coordinates.
The director's sparse-subspace net count gives, uniformly over the
orthogonal basis, at most exp[b(epsilon,V)m+o(m)] such ternary words,
where b can be made arbitrarily small by first taking epsilon small
and then V large. One explicit elementary choice, with
rho=1/V^2 and theta=16p epsilon<1/2, is

    b(epsilon,V)=h2(rho)+rho log(1+2/sqrt(epsilon))
                       +h2(theta)+theta log 2,            (7)

up to harmless rounding/o(m). It follows by an epsilon-scale net on
the unit ball in each sparse coordinate subspace, followed by a
ternary Hamming-ball bound. The exact constants in (7) are not needed
for the combination, only b tending to zero in the indicated order.

The orbit factor L_t is at most one. The maximum deleted-coordinate
factor in the original row theorem is at most sqrt(2m)L_t, so its
polynomial loss is harmless. Thus the restricted one-row partition
function for non-diffuse rows obeys

    limsup_m (1/m) log E Z_spiky <=-h2(p)+b(epsilon,V).     (8)

The unrestricted one-row theorem at fixed depth r gives

    limsup_m (1/m) log E Z_all <=a_r,
    a_r=p log 2+B^r Phi_t(nu_p),
    limsup_(r->infinity) a_r<=a:=p log 2+E_t(nu_p).        (9)

All row indicators here are spectral-permutation invariant. They can
be retained while performing the original conditional row-permutation
graph-CS bound. The independent selectors and bases then retain the
product of independent row partition functions. Summing the at most
2^m choices of which rows are spiky costs only exp(O(m)), not exp(O(m^2)).

If fewer than delta m rows are diffuse, at least (1-delta)m rows are
spiky. Combining (8)--(9), this region has exponent at most

    a_r-(1-delta)[a_r+h2(p)-b(epsilon,V)]                 (10)

per m^2, provided the bracket is positive.

## 4. Positive phase gap is the exact closing condition

Suppose

    g:=H(nu_p)+E_t(nu_p)
       =h2(p)+p log 2+E_t(nu_p)>0.                       (11)

Choose epsilon and V so b<g/4. For example fix delta=1/2. Given any
small zeta>0, choose a sufficiently large but fixed depth so the common
unrestricted row upper bound is a+zeta. Use THIS common upper bound
in both regions, even if the true finite-depth exponent is smaller.
The bad-region loss is at least (1-delta)(g-b), whereas the good-region
loss in (6) is also a fixed positive constant. Combining
the two regions gives a strict improvement over the all-vector exponent:

    limsup_m (1/m^2) log E sum_(stable x)
          exp((t/k)sigma x^T C x)
       <=t+a-Delta                                      (12)

for some Delta>0, after choosing a sufficiently large but FIXED depth r.
The harmless physical-diagonal prefactor costs O(m), not O(m^2): indeed
the folded-kernel identity bounds it by exp(tm^2+tm). The conclusion
uses no growing-depth limit theorem.

The stable-vector first moment therefore proves an actual all-order cap
bound strictly below

    [t+p log 2+E_t(nu_p)]/(2t sqrt(p)).                   (13)

More conservatively, an explicit certified upper endpoint e_bar for E_t
may replace it throughout. For improvement over THAT certified expression,
only h2(p)+p log2+e_bar>0 is needed: use its value plus zeta as the common
unrestricted row upper in both sectors, regardless of whether the true
E_t is smaller. Finite depth is chosen after the strict exponent margin,
just as in the original all-order theorem.

At p=24/25 and t=97/20, the Gaussian lower bound E_t>=g_t(1) gives

    H(nu_p)+g_t(1)=0.0088920962596074...>0.

This is a robust positive margin, not the tiny stability constant. For the
certified endpoint e_bar=-5151/6250, the required gap is even larger and
has been independently checked with outward intervals: it exceeds .009,
while the director's choice epsilon=1/100000,V=100 gives b<.004 (using
the slightly larger p-independent Hamming radius 16epsilon). The resulting
new cap improvement is extraordinarily small because c_row deteriorates
rapidly through U,V,kappa. No displayed decimal improvement is claimed.

The same H2/H12 terminal orders, arbitrary fixed depth, and principal
restriction fill all large orders. Random selectors do not alter the
dimension N=mk or the actual full-sign property of C.
