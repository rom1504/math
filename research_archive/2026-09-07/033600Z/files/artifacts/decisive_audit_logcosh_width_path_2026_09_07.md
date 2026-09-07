# Independent audit: optimized width's exact log-cosh path

The theorem in `decisive_bridge_width_logcosh_interpolation_2026_09_07.md`
passes independent algebraic and endpoint audit.  It is an exact finite
identity, not a convergence theorem: its weighted flip cost has no proved
sign.  The recurrence obtained if that cost is controlled is almost
**superadditive**, not subadditive.

## 1. Flip cost and cavity optimality

Write `F=.5 log(Z_+ Z_-)`, `d_e=.5 A_e(C^+_e-C^-_e)`, and
`P_e=C^+_e C^-_e`.  Directly in the *full* Gibbs laws, flipping edge `e`
gives

    exp(2 Delta_e)
      =(cosh(2lambda)-A C^+ sinh(2lambda))
       (cosh(2lambda)+A C^- sinh(2lambda))
      =1+sinh(2lambda)^2 [1-P_e-2 d_e coth(2lambda)].

Thus the stated `a_e` is nonnegative at every actual global minimizer.
No Gaussian replacement or small-coupling expansion enters.

Alternatively let `c_+,c_-` be the correlations with this edge deleted,
`t=tanh(lambda)`, `b=|c_+-c_-|`, and `p=c_+c_-`.  The product insertion
factor is `cosh(lambda)^2(1+A t(c_+-c_-)-t^2 p)`.  Since this factor is
strictly increasing in `A(c_+-c_-)`, positive-weight edge optimality is
exactly `A(c_+-c_-)=-b`.  The denominator `q=1-t b-t^2 p` is positive,
and subtracting the two possible factors gives

    a_e=b(1-t^2)^2/(2t q).

This also proves the endpoint temperature estimate used in the source.
Indeed

    d_e=t-(1-t^2)(b+2t p)/(2q) <= t.

If `p>=0` the last sign is immediate.  If `p<0`, write the two opposite
absolute cavity correlations as `u,v` in `[0,1]`; then
`b=u+v>=2uv=-2p>=-2t p`.  Hence for homogeneous optimized order `d`
pressure,

    0 <= Psi_d'(beta) <= (d-1) sqrt(d) tanh(beta/sqrt(d))/2
                       <= beta(d-1)/2

at differentiability points.  The lower bound follows because every
fixed-signing paired pressure is even and convex, and therefore
nondecreasing on positive temperatures; their finite minimum remains
nondecreasing.  Convexity of that minimum is not needed or asserted.

## 2. Unequal blocks, diagonal terms, and the PSD sign

Let `z` take values `sqrt(n/m)` and `-sqrt(m/n)` on the two blocks.
Then `sum z=0`.  Put `D=zz^T`, retaining its actual diagonal values
`n/m` and `m/n`; it is **not** a hollow zero-row-sum matrix.
Nevertheless `P_ii=1`, so

    sum_(i<j) D_ij(1-P_ij)
      =.5 sum_(i,j) D_ij(1-P_ij)
      =-.5 z^T P z.

Both branch correlation matrices are PSD, irrespective of their mutual
relationship.  Their Schur product `P` is PSD as well, so this quantity
is nonpositive for every unequal or equal split.  With
`tau=log cosh(2 beta/sqrt(N))` and
`log cosh(2lambda_e(u))=tau(1+u D_e)`, the chain rule is precisely

    F_*'(u)=-tau z^T P z/8 -tau sum_e D_e a_e/4.

There is no diagonal contribution missing and no Taylor remainder.

## 3. Active minima and the zero-cross endpoint

On every closed interior subinterval each fixed-signing pressure is
smooth, and there are finitely many signings.  At almost every point
the finite minimum is differentiable, and every active smooth branch
has the same derivative there.  Thus one may use the flip inequalities
of an actual active minimizer in the derivative identity.

On the whole interval each `lambda_e` is absolutely continuous.  A
vanishing cross coefficient has derivative `O((1-u)^(-1/2))`, which
is integrable.  Each fixed-signing `F` is coordinatewise 1-Lipschitz,
and the same common modulus controls their finite minimum.  This proves
absolute continuity through `u=1` and licenses integration of the
almost-everywhere identity.  In particular, apparent division by a
vanishing cross coefficient is not an assertion of positive-weight
optimality at zero.  The two terms are individually integrable: `P`
is bounded, while the flip identity and `|Delta_e|<=2lambda_e` give
`a_e=O(1/lambda_e)` as `lambda_e` tends to zero.

For fixed positive block fractions, expansion of `log cosh` gives
`sqrt(m)lambda_1(1)=beta+O(1/N)` and the analogous second block.
The preceding derivative bound gives an `O(1)` total endpoint error.
At zero cross weight the paired partition factors *separately in each
branch*, and minimization over the two internal signings separates.
Therefore

    F_*(1)=Psi_m(beta)+Psi_n(beta)+O_(beta,m/N)(1).

If the integrated weighted slack is bounded below by `-K N^(1+alpha)`
with `alpha<1`, the nonpositive PSD term yields

    Psi_N(beta) >= Psi_m(beta)+Psi_n(beta)-O(N^alpha)-O(1).

This is the claimed sufficient almost-superadditive direction.  The
known pointwise nonnegativity of `a_e` does not prove that weighted
bound because `D_e=-1` across the split.  Moreover convergence of
optimized paired width pressure would still require a separate bridge
to the original two-sided absolute cap.

## 4. Reproducible finite audit and relocation scope

`computations/decisive_audit_logcosh_width_path_2026_09_07.py` enumerates
all switching-gauged signings at orders 3 through 6.  It tests actual
weighted global minimizers at unequal and equal splits, three
temperatures and four interior interpolation times.  It verifies both
slack formulas, every edge-flip cost, the PSD identity and finite
differences of the optimized envelope.

The same replay independently checks all 8192 projective spins of the
archived order-14 conference completion, yielding cap 21 and squared
row norm 13.  Thus both elementary counterexamples in
`decisive_director_variance_relocation_scope_2026_09_07.md` also pass:
bounded amplitude alone permits the concentrated principal block with
cap below .350; balanced row variance alone permits repeated scaled
order-14 blocks with cap below .417.  Neither example satisfies both
conditions simultaneously, and no stronger relocation impossibility
is claimed.

## 5. Endpoint variance formula and exact finite nonmonotonicity

The source's subsequently appended endpoint reduction also checks.  For
fixed equal-size child Gibbs systems, independent child spin reversal
makes the cross statistic centered.  Its paired log-pressure gain is
`epsilon^2 V(B)/4+O(epsilon^4)`, with
`V(B)=sum_s Tr(B^T C_1^s B C_2^s)`.  Since
`epsilon^2=tau(1-u)/2+O((1-u)^2)`, the cross contribution to the left
derivative at `u=1` is `-tau V(B)/8`.  Taking a finite minimum just
before the endpoint selects the **largest** endpoint derivative; when
the child radial derivatives agree this is the smallest `V(B)` among
all active child pairs, including relative global polarities.  The
generic PSD strengthening in the source is correctly marked false and
is not the actual response-energy inequality.

There is a simple exact finite example showing that even the complete
optimized path derivative need not be nonpositive.  At order 4 split
into two pairs, write `x>=y` for the two internal and four cross
magnitudes, and `t=tanh x`, `v=tanh y`.  Enumeration of the eight
switching classes gives exactly

    F_*(x,y)=4log2+2log cosh x+4log cosh y
               +log(1+v^4-2t^2 v^2).

Indeed the three possible squared paired high-temperature loop factors
are

    (1+v^4-2t^2v^2)^2,
    (1-v^4)^2,
    (1+v^4+2t^2v^2)^2-16t^2v^4.

The first is smallest: its positive square root differs from `1-v^4`
by `2v^2(v^2-t^2)<=0`, and the third minus the first is
`8t^2v^2(1-v^2)^2>=0`.  These cases follow also by classifying the
four triangle signs: they have product one; two negative triangles
give one positive and two negative four-cycle products, and four
equal triangle signs give the last displayed factor.

Along the equal-split log-cosh path its left derivative at `u=1` is

    tau (3 tanh(x_end)^2-1)/2.

For the exact temperature `beta=log3`,
`tau=log(5/3)` and `tanh(x_end)^2=8/17`, giving the strictly positive
derivative `7log(5/3)/34`.  Here the weighted slack tends to `-32/17`;
its unfavorable contribution exceeds the favorable PSD contribution.
This rules out a pointwise finite monotonicity claim even for actual
global minimizers.  It does **not** rule out a sublinear integrated
error or any asymptotic convergence theorem.
