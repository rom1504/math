# A fixed-degree spiked response tensor for the actual-child pressure

Status: **rigorous cross-row closure obstruction, with a decisive scope
qualification**.  Exact nonnegative degree-two row densities of fixed
`L^2` norm expose a Boolean response landscape on all left-child row signs.
The associated row tensor has no small singular multiplier in its fully
active channel, and the entropy charge is independent of the queried sign
word.  Thus replacing positive-part carriers by exact polynomial densities
does not make fixed row degree into bounded global row order.

The theorem applies to the exact pressure of every pair of actual optimizing
children and uses no surrogate.  It does **not** prove that the exposed
landscape has linear range for a sequence of actual minimizers.  Therefore it
is a real algebraic/precision ceiling for fixed-row-degree closure, but not a
scalable fixed-density collision and not a resolution of the current SML.

## 1. A diffuse degree-two sign interface

Fix a row width `n`, a direction `y in {+-1}^n`, and put

```math
 z_y(b)={1\over\sqrt n}\langle y,b\rangle,
 \qquad e_y(b)={1+z_y(b)^2\over2}.                 \tag{ST.1}
```

For `v in {+-1}`, define

```math
 q_{v,y}(b)=e_y(b)+v z_y(b)
 ={(1+v z_y(b))^2\over2}.                          \tag{ST.2}
```

**Lemma ST.1 (uniform polynomial carrier).**  The function `q_(v,y)` is a
probability density with respect to the fair row law `U_n`, is nonnegative,
has Walsh degree at most two, and satisfies

```math
 \boxed{
 \|q_{v,y}\|_2^2={5\over2}-{1\over2n},
 \qquad
 E_{q_{v,y}}B={v y\over\sqrt n},
 \qquad
 \operatorname {Cov}_{q_{v,y}}(B)
 =\left(1-{1\over n}\right)I_n.}                  \tag{ST.3}
```

Its entropy `D(q_(v,y)U_n || U_n)` is independent of `v` and `y`.

*Proof.*  Under `U_n`, `z_y` is symmetric and

```math
 E z_y=0,\qquad E z_y^2=1,\qquad
 E z_y^4=3-{2\over n}.                             \tag{ST.4}
```

Equations (ST.2)--(ST.4) give nonnegativity, unit mean, and

```math
 E q_{v,y}^2={1+6E z_y^2+E z_y^4\over4}
 ={5\over2}-{1\over2n}.
```

Parity gives `E[B_jz_y^2]=0`, while
`E[B_jz_y]=y_j/sqrt(n)`, proving the mean formula.  For `j ne k`,

```math
 E_{q_{v,y}}B_jB_k={y_jy_k\over n},               \tag{ST.4a}
```

whereas the diagonal entries are one.  Subtracting the outer product of
the mean gives the covariance formula in (ST.3).  Row sign flips and
coordinate gauges act transitively on `(v,y)` and preserve `U_n`; hence all
the densities are translates and have the same entropy. `square`

Unlike a point row, (ST.2) is uniformly diffuse:

```math
 q_{v,y}(E)\le
 \sqrt{{5\over2}}\sqrt{U_n(E)}.                    \tag{ST.5}
```

Nevertheless its mean row is the coherent vector `v y/sqrt(n)`.  Across
`m` rows this gives a rank-one mean bridge at exactly the physical scale.
But the centered row covariance is `(1-1/n)I_n`: the carrier retains a
full-dimensional random residual.  Consequently “rank-one mean” is not an
exact response reduction for the quenched pressure.

## 2. The exact `2^m` response tensor

Let

```math
 L:\{+-1\}^{m\times n}\longrightarrow\mathbb R   \tag{ST.6}
```

be the exact bridge pressure of two actual contracted-temperature minimizing
children, in either orientation.  The argument below only uses its exact
point values and therefore holds uniformly over those choices.  For
`v=(v_1,...,v_m) in {+-1}^m`, let

```math
 P_{v,y}=\bigotimes_{i=1}^m q_{v_i,y}U_n,
 \qquad R_{L,y}(v)=E_{P_{v,y}}L.                   \tag{ST.7}
```

**Theorem ST.2 (unattenuated fixed-row-degree response tensor).**  The Walsh
coefficients of the Boolean response landscape `R_(L,y)` are exactly

```math
 \boxed{
 \widehat R_{L,y}(S)
 =E_{U_B}\left[
   L(B)\prod_{i\in S}z_y(B_i)
       \prod_{i\notin S}e_y(B_i)\right],
 \qquad S\subseteq[m].}                            \tag{ST.8}
```

The `2^m` test functions on the right are mutually orthogonal.  Their squared
norms are

```math
 \left({3\over2}-{1\over2n}\right)^{m-|S|}.        \tag{ST.9}
```

In particular, the fully active row channel has norm exactly one:

```math
 \left\|\prod_{i=1}^m z_y(B_i)\right\|_2=1.        \tag{ST.10}
```

Thus the passage from the `2^m` feasible carrier values to the fully active
cross-row coefficient has no exponentially small row singular multiplier.

*Proof.*  The product density in (ST.7) is

```math
 {dP_{v,y}\over dU_B}(B)
 =\prod_i\{e_y(B_i)+v_i z_y(B_i)\}.                \tag{ST.11}
```

Multiply by `v_S` and average over uniform `v`.  In row `i`, the average
selects `z_y(B_i)` if `i in S` and `e_y(B_i)` otherwise, proving (ST.8).
The function `e_y` is even under `b -> -b`, while `z_y` is odd, so they are
orthogonal.  Independence of the fair rows then proves mutual orthogonality.
Finally,

```math
 E z_y^2=1,\qquad
 E e_y^2={1+2E z_y^2+E z_y^4\over4}
 ={3\over2}-{1\over2n},                            \tag{ST.12}
```

which proves (ST.9)--(ST.10). `square`

This is stronger, for the present precision question, than exact
table-completeness caused by a tiny positive-part leakage.  Every density in
(ST.7) is itself a degree-two polynomial; its norm is fixed; and the odd row
mode in (ST.2) has unit `L^2` norm.  The exponential tensor comes solely from
composition across rows.

## 3. Correct physical scale and constant entropy charge

Put `u=beta/sqrt(N)`, `N=m+n`.  Changing one bridge bit changes the actual
pressure by at most `2u`.  If all rows except `i` have any genuine product
laws and

```math
 F_i(b)=E[L(B)\mid B_i=b],                          \tag{ST.13}
```

then fair-cube Poincare gives

```math
 \|F_i-E_UF_i\|_2\le u\sqrt n.                    \tag{ST.14}
```

Since `q_(+,y)-q_(-,y)=2z_y`, (ST.14) proves

```math
 \boxed{
 |R_{L,y}(v)-R_{L,y}(v^{(i)})|
 \le2u\sqrt n
 =2\beta\sqrt{n/N}.}                              \tag{ST.15}
```

Equivalently, Parseval on the query cube gives

```math
 \boxed{
 \sum_{S\subseteq[m]}|S|\widehat R_{L,y}(S)^2
 \le m u^2n={\beta^2mn\over N}.}                  \tag{ST.16}
```

Thus this exposed response is a Boolean landscape with order-one coordinate
oscillations and a potentially order-`N` diameter, while its typical
`L^2(v)` fluctuation is only order `sqrt(N)`.  The distinction is exactly a
rare-extreme question: (ST.16) does not exclude a linear minimum excursion
among `2^m` queries.

The entropy term does not simplify that query.  If

```math
 d_n=D(q_{+,y}U_n\Vert U_n),                       \tag{ST.17}
```

then every product in (ST.7) has entropy `m d_n`.  Consequently the full
row-product objective restricted to this explicit subfamily is

```math
 \boxed{
 \mathcal F(P_{v,y})=R_{L,y}(v)+{m d_n\over\lambda},
 \qquad
 \min_v\mathcal F(P_{v,y})
 =\min_vR_{L,y}(v)+{m d_n\over\lambda}.}           \tag{ST.18}
```

Hence no separately paid scalar channels, point rows, or hidden entropy
variation create the tensor.  A linear excursion of `R_(L,y)` would be one
coherent fixed-degree branch-(iii) certificate at the correct leading scale.

## 4. Consequence and exact limitation

The theorem rules out the following inference:

> a fixed Walsh degree in every row implies that polynomially many
> bounded-total-row-order pressure observables evaluate the carrier optimum.

Even degree two exposes all row orders, and the fully active tensor mode is
not attenuated by its presenting density.  A valid polynomial child closure
must therefore prove one of the following genuinely optimizer-specific
statements:

1. the high-row-order coefficients in (ST.8) synchronize to a smaller child
   order parameter;
2. their contribution to `min_v R_(L,y)(v)` is `o(N)` uniformly for actual
   minimizing children; or
3. a tractable recursive algebra evaluates that minimum without listing its
   `2^m` responses.

No current flip, deletion, sector--Gram, or one-row best-response identity
implies any of these statements.  The identities constrain the child Gibbs
laws or selected low-order tangents, whereas (ST.8) is a coherent positive-
density row response of the exact two-child pressure.

The theorem is **not** a robust actual-minimizer no-go at fixed density.
Equations (ST.15)--(ST.16) allow, but do not force, a linear extreme.  Proving
that some sequence of actual minimizing children has

```math
 \max_vR_{L,y}(v)-\min_vR_{L,y}(v)\ge cN           \tag{ST.19}
```

would be the required scalable obstruction; proving the left side is `o(N)`
uniformly would instead be precisely a useful synchronization theorem.  At
present neither conclusion follows.  Accordingly this audit sharpens the
macroscopic finite-degree child-closure SML but does not reset it.
