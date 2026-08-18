# Second independent audit: speed-`r` conference basin candidates

**Frozen source:**
`extremal_information/drafts/conference_speed_r_basin_candidate_audit.md`

**SHA-256:**
`688bcfd5303932728f6b8fc6ec1675dbb1e3c00ceece034fad2b6512e033d13a`

**Verdict: PASS WITH ONE NOTATIONAL REPAIR.**  Theorems CB.1, CB.3,
CB.5 and Corollaries CB.2, CB.4 are valid, including their constants and
their stated limitations.  In CB.24--CB.25 the variance profile must be
held fixed while the translation variable is varied.  With that repair the
convex-even comparison is rigorous and no theorem statement changes.

The report does **not** rule out every family of relative probability
`exp(-Cr)`.  It rules out the listed affine, product-biased, type-shell, and
regular adaptive-gauge mechanisms as whole favorable basins.  An additional
pressure-correlated `exp(-Theta(r))` fraction of a type shell, or the
operator-irregular `exp(-Theta(r))` remainder of a switching cross-section,
can still have total probability `exp(-O(r))`.  The frozen source says this
explicitly; none of its proved displays should be read as excluding those
subfamilies.

## 1. Normalization and edit modulus

For

```math
G_{epsilon,B}(x,y)
=H_A(x)+epsilon H_A(y)+x^TBy,
\qquad t={\beta\over\sqrt{2r}},
```

one has, for every state,

```math
|G_{epsilon,B}(x,y)-G_{epsilon,B'}(x,y)|
\le\|B-B'\|_1.
```

The inequality `cosh(u+v)<=e^{|v|}cosh(u)`, applied in both directions,
therefore proves

```math
|f(B)-f(B')|\le t\|B-B'\|_1.
```

A sign-coordinate flip has `l_1` cost two, hence oscillation `2t`.  Thus
CB.7--CB.8 have no missing factor from the two symmetric blocks in the
parent matrix.  If `s_r=o(r^(3/2))`, then

```math
{2ts_r\over r}
={\sqrt2\beta s_r\over r^{3/2}}=o(1),
```

which is exactly the scale used in CB.1.

## 2. Uniform retractions and affine fibres

Let `U_r` be uniform on the bridge cube.  The hypotheses of CB.1 give the
literal coupling

```math
B_r=\pi_r(U_r),
\qquad B_r\sim\operatorname {Unif}(F_r),
```

and the preceding deterministic bound gives

```math
|f(B_r)-f(U_r)|/r=o(1).
```

The archived convergence in probability of `f(U_r)/r` therefore transfers
to `F_r`.  No injectivity, cardinality estimate, or concentration theorem is
being silently assumed.

For a rank-`q_r` affine parity system, choose `q_r` pivot columns of its
binary check matrix.  Every free assignment determines exactly one pivot
assignment in the chosen syndrome.  The overwrite map ignores the input
pivots, so every output has exactly `2^(q_r)` preimages and the uniform cube
pushes forward to the uniform fibre.  At most `q_r` coordinates change.
The fibre has exactly

```math
2^{r^2-q_r}
```

members.  This verifies both the retraction and entropy assertions in CB.2.

For prescribed row and column products, the product of all row equations
equals the product of all column equations and this is the only dependency.
The rank is `2r-1`, so the dimension is

```math
r^2-(2r-1)=(r-1)^2.
```

The corollary therefore applies to arbitrary supports of the checks, not
only local parities.

## 3. Lindeberg comparison and the convex-even direction

The auxiliary-sign representation makes `f` a log partition function in
the observable

```math
Z_e=\sigma x_i y_j\in\{+-1\}.
```

Differentiating three times in a real bridge coordinate gives exactly

```math
\partial_e^3f=t^3\kappa_3(Z_e).
```

If `mu=E Z_e`, then

```math
\kappa_3(Z_e)=-2\mu(1-\mu^2),
```

so the source's bound `|kappa_3|<=2` is valid (and non-sharp).  The biased
sign `B_e` and

```math
Y_e=m_e+\sqrt{1-m_e^2}V_e
```

have the same first two raw moments:

```math
E B_e=E Y_e=m_e,
\qquad E B_e^2=E Y_e^2=1.
```

Both are uniformly bounded.  Taylor expansion to second order in each
replacement coordinate therefore costs `O(t^3)`, uniformly in all already
replaced coordinates.  Summing over `r^2` entries gives

```math
r^2t^3=O_\beta(\sqrt r),
```

as in CB.23.

There is one necessary notational repair in the next step.  For the already
chosen vector `M`, freeze

```math
a_e=\sqrt{1-m_e^2},
\qquad W_M=(a_eV_e),
```

and define a function of a *new* translation variable `Z` by

```math
G_M(Z)=E_V f(Z+W_M).
```

Then `G_M` is convex.  Since `W_M` is symmetric and `f(-B)=f(B)`, it is
also even, so

```math
E f(M+W_M)=G_M(M)\ge G_M(0)=E f(W_M).
```

If instead the notation `G(M)` in CB.24 is read as recomputing
`a_e=sqrt(1-m_e^2)` while its argument varies, convexity does not follow.
The frozen-profile reading above is plainly the intended one and proves the
claimed comparison direction.

Couple `W_M` to a uniform bridge `V` with the same coordinate signs.  Since

```math
1-\sqrt{1-u^2}\le u^2,
```

CB.7 gives

```math
|f(W_M)-f(V)|
\le t\sum_e(1-a_e)
\le t\sum_em_e^2
=O_{\beta,C}(\sqrt r).
```

Combining this with the uniform mean asymptotic proves CB.19 with the
correct lower-bound direction.

Under the biased product law, every coordinate oscillation is `2t`, and

```math
\sum_ec_e^2=r^2(2t)^2=2\beta^2r.
```

McDiarmid's lower-tail inequality is therefore

```math
\Pr\{f\le Ef-u\}
\le\exp\{-u^2/(\beta^2r)\}.
```

For large `r`, CB.19 gives `Ef>=(h_beta-eta/2)r`.  Taking
`u=eta r/2` proves exactly the exponent
`eta^2r/(4beta^2)` in CB.20.

## 4. Exact type-shell count and conditioning

Write `N=N_r`, `m=c/sqrt(r)`, and `p=(1+m)/2`.  Stirling expansion at
`k=floor(pN)` gives

```math
\begin{aligned}
\log {N\choose k}
&=N\log2-{Nm^2\over2}+O(Nm^4+\log N)\\
&=N\log2-{\theta c^2\over2}r+O(\log r).
\end{aligned}
```

Here `Nm^4=O(1)`.  Adding the `r^2-N` free unbiased coordinates proves
CB.31--CB.32.

Under the associated product law, the number of template agreements is
`Bin(N,p)`, and `k` is within one of its mean.  Since `p->1/2`, the local
central-binomial estimate is

```math
\Pr\{\operatorname {Bin}(N,p)=k\}=\Theta(N^{-1/2})=\Theta(1/r).
```

Conditioned on that event, every shell member has the same probability and
the outside coordinates remain free uniform.  Dividing CB.20 by this
polynomial conditioning probability preserves an `exp(-c r)` estimate,
which proves CB.4.  A fixed number of disjoint block types costs only a
fixed power of `r` and gives the entropy loss in CB.34.

This result says that a uniformly chosen shell member is not a fixed lower
deviation with probability `1-exp(-Omega(r))`.  It does not say that the
shell contains no `exp(-Theta(r))` exceptional fraction.

## 5. Factorization and one-sided witnesses

Pair `y` with `-y`.  The child term

```math
C=H_A(x)+\epsilon H_A(y)
```

is unchanged and the bridge term `D=x^TBy` changes sign.  The identity

```math
{\cosh(t(C+D))+\cosh(t(C-D))\over2}
=\cosh(tC)\cosh(tD)
```

proves CB.36 exactly.

For `B=uv^T`, the state `(x,y)=(u,v)` has `x^TBy=r^2`.  Its one term in
the normalized average, together with `cosh(tC)>=1`, gives

```math
f(uv^T)\ge tr^2-(2r+1)\log2
={\beta\over\sqrt2}r^{3/2}-2r\log2-O(1).
```

Thus CB.37 has the correct leading coefficient and normalization.  This
only proves that the coarse halfspace or row-majority constraint contains a
high member; the source correctly avoids inferring that every member is
high.

## 6. Adaptive row/column gauges

If `D_sBD_u=B`, then `s_i u_j=1` for every pair `(i,j)`.  Hence `s` and
`u` are both constant with the same sign.  After quotienting by the common
global sign, the row/column action is free and has size `2^(2r-1)`.
Row-only switching is free with size `2^r`.  The cross-section cardinalities
in CB.38--CB.39 are consequently exact.

Fix `delta>0` small enough that

```math
{\beta(3+\delta)\over\sqrt2}<\kappa<1/2.
```

On the common event

```math
\|B\|_{op}\le(2+\delta)\sqrt r,
```

every row/column switched bridge has the same norm, and the conference
triangle bound gives

```math
\left\|t
\begin{pmatrix}A&D_sBD_u\\D_uB^TD_s&\epsilon A\end{pmatrix}
\right\|_{op}
\le{\beta\over\sqrt2}
 \left(\sqrt{1-1/r}+2+\delta\right)<\kappa
```

for all large `r`.  Thus all gauges lie simultaneously in the
operator-regular class.  The complement of this common event has probability
`exp(-c_0r)` by the rectangular Rademacher norm tail.

For each fixed gauge, `D_sBD_u` is again a uniform sign bridge and preserves
the bridge norm.  Equivalently, conjugation gives CB.42.  Hence the
orientation-specific regular-sector theorem bounds

```math
\Pr\{\|B\|_{op}\le(2+\delta)\sqrt r,
       f(D_sBD_u)\le(h_\beta-\eta)r\}
\le e^{-c_\eta r^2}
```

with constants independent of `(s,u)`.  Union over at most `2^(2r-1)`
gauges changes this to

```math
\exp\{-c_\eta r^2+O(r)\}\le e^{-c_1r^2}.
```

Adding the norm-tail complement proves CB.40.  If a canonical selector
maps every orbit to one cross-section, all orbits have the same size and
the uniform cube pushes forward to the uniform law on that cross-section.
The asserted vanishing lower-deviation fraction follows.

The `e^{-c_0r}` term is essential.  CB.5 proves a quadratic lower tail only
on the common regular event; it neither proves that the adaptive output is
regular with superexponential probability nor excludes a speed-`r`
favorable subset of the irregular outputs.

## 7. Archive comparison and final scope

1. A prior task-local audit at a different frozen hash records CB.1--CB.4
   and the same frozen-profile clarification.  The current source adds the
   adaptive-gauge theorem CB.5.  Thus CB.1--CB.4 are not new relative to
   that immediate draft history, although they are not duplicated as a
   general theorem in the canonical ledger.
2. Canonical Theorem 37.6 supplies the operator-regular quadratic lower
   tail and Hamming-collar localization used by CB.5.  CB.5 is a new
   switching-orbit corollary of that result, not a new concentration
   mechanism.
3. The switching covariance and orbit counts overlap algebraically with
   `switching_gauge_quotient_for_optimized_bridges.md`.  That source studies
   optimized bridge fibres; it does not prove the adaptive conference
   pressure union bound.
4. The generic entropy-transport theorem permits an `O(r)`-entropy law to
   change pressure by order `r`.  CB.3 is genuinely sharper on independent
   weak biases because convex evenness fixes the direction and removes the
   possible leading gain.
5. The mesoscopic overwrite and Hamming-collar results already identify
   `Theta(r^(3/2))` as the critical edit scale.  CB.1 gives the complementary
   pushforward statement below that scale; it does not replace the collar
   theorem.

Accordingly the final research judgment CB.43 is valid only with the
source's stated quantifiers: a family that is itself a certified
fixed-lower-pressure basin of probability `exp(-O(r))` cannot be one of the
tested affine/product/regular-gauge families.  It may still be an
additional nonlinear subset inside one of them, and the present results do
not determine whether such a subset exists.
