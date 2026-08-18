# Switching-flat child optimality versus the cross-row response tensor

Status: **rigorous actual-child identity audit and exact finite
falsifier**.  The degree-two query surface in ST.2 has an exact
interpretation which is useful for judging optimizer identities: its Boolean
row word and carrier direction are switching gauges of the left and right
children.  Contracted-temperature minimality is completely flat along that
two-sided orbit.  The full internal
edge-replacement/contraction atlas is identical at every point of the orbit,
and the deletion--reinsertion atlas is only gauge-permuted.  Hence those
identities provide no directional stationarity condition on the query word.

There is one universal cross-row identity: the response is antipodally even,
so all odd row-order Walsh coefficients vanish.  An order-two pair of actual
minimizing children has every presently known internal optimality defect equal
to zero but has a strictly positive even coefficient of maximal row order.
Thus neither switching symmetry, convex endpoint stationarity, the complete
flip/contraction box, nor the one-vertex Bellman identities suppresses the
first admissible cross-row mode.

The example is finite.  It does not prove a scalable linear response range,
and therefore does not decide `L_balanced-product-phase`.  It rigorously
excludes the known optimizer identities as a source of a coherent query
direction without an additional gauge-covariant multirow observable.

## 1. Exact two-sided switching-orbit representation

Let `A,D` be the two actual contracted-temperature minimizing children of
orders `m,n`, fix a relative orientation, and let

```math
 L_{A,D}(B)
```

denote their exact bridge log pressure at the physical bridge amplitude.
Additive constants in `L` will never matter.  For `v in {+-1}^m`, write

```math
 A^v=\operatorname {diag}(v)A\operatorname {diag}(v).
                                                               \tag{SO.1}
```

Fix `y in {+-1}^n` and use the exact degree-two row densities of ST.1,

```math
 q_{v_i,y}(b)={\{1+v_i\langle y,b\rangle/\sqrt n\}^2\over2},
 \qquad
 P_{v,y}=\bigotimes_{i=1}^m q_{v_i,y}U_n.           \tag{SO.2}
```

The corresponding response is

```math
 R_{A,D;y}(v)=\mathbb E_{P_{v,y}}L_{A,D}(B).        \tag{SO.3}
```

Let `mathbf1` denote the all-positive word of length `n`, and let
`P_(+,mathbf1)` denote the product whose every row has density
`q_(+,mathbf1)`.

**Theorem SO.1 (the carrier is a two-sided child switching orbit).**  One
has the exact identity

```math
 \boxed{
 R_{A,D;y}(v)
 =\mathbb E_{P_{+,\mathbf1}}L_{A^v,D^y}(B),}
 \qquad
 D^y=\operatorname {diag}(y)D\operatorname {diag}(y).         \tag{SO.4}
```

Every `A^v` and `D^y` is again an exact contracted-temperature minimizing
child, with exactly the same sector normalizations as the original child.
In addition,

```math
 \boxed{R_{A,D;y}(-v)=R_{A,D;y}(v),}               \tag{SO.5}
```

and consequently

```math
 \boxed{\widehat R_{A,D;y}(S)=0\quad\text{for odd }|S|.}       \tag{SO.6}
```

*Proof.*  Let `C` have law `P_(+,mathbf1)` and put

```math
 B=\operatorname {diag}(v)C\operatorname {diag}(y).             \tag{SO.7}
```

For row `i`, the inverse change of variables is
`C_i=v_i(B_i\odot y)`, and therefore

```math
 z_{\mathbf1}(C_i)=v_i z_y(B_i),
 \qquad q_{+,\mathbf1}(C_i)=q_{v_i,y}(B_i).         \tag{SO.7a}
```

The fair row law is preserved by this bijection, so `B` has exactly the
law `P_(v,y)`, with no normalization factor.  Changing spin variables
`x -> diag(v)x` and `y_spin -> diag(y)y_spin` in the exact parent partition
function gives

```math
 L_{A,D}(\operatorname {diag}(v)C\operatorname {diag}(y))
 =L_{A^v,D^y}(C),                                  \tag{SO.7b}
```

which proves (SO.4).  Applying the same changes of variables to either
one-sided child partition leaves the sector sign unchanged, so every sector
partition, the augmented pressure, and the contracted-temperature
minimality property are preserved exactly.

The exact bridge pressure is even under `B -> -B`: change `y_spin` to
`-y_spin` in the parent spin sum.  Also `P_{-v,y}` is the image of
`P_{v,y}` under `B -> -B`.  This proves (SO.5), and Fourier expansion on the
query cube proves (SO.6). `square`

Equation (SO.6) is an honest cross-row identity, but it removes only half of
the exponentially many query channels.  It is caused by the global spin
symmetry of the exact parent pressure, not by child optimality.

It is useful to expose the complete group action.  Define

```math
 \mathscr R_{A,D}(v,y)
 =\mathbb E_{P_{+,\mathbf1}}L_{A^v,D^y}(B).         \tag{SO.7c}
```

For any switching words `a,b`, switching composition gives

```math
 \boxed{
 \mathscr R_{A^a,D^b}(v,y)
 =\mathscr R_{A,D}(av,by).}                        \tag{SO.7d}
```

Moreover `mathscr R` is unchanged by independently replacing `v` or `y`
by its global negative.  Hence it descends to the projective switching
groups

```math
 \Gamma_m=\{+-1\}^m/\{+-\mathbf1\},
 \qquad
 \Gamma_n=\{+-1\}^n/\{+-\mathbf1\}.              \tag{SO.7e}
```

The full exposed orbit surface still has as many as
`2^(m+n-2)` entries.  Equation (SO.4) is therefore an exact orbit
parametrization, not by itself a finite-state closure.

## 2. What the complete known optimizer atlas does along this orbit

Let `S` be an arbitrary set of internal child edges.  Switching commutes
with edge replacement:

```math
 (A^v)^S=(A^S)^v.                                  \tag{SO.8}
```

Since child pressure is switching invariant, the *whole* replacement table
is pointwise identical:

```math
 {\overline Z_{(A^v)^S}(t)\over\overline Z_{A^v}(t)}
 ={\overline Z_{A^S}(t)\over\overline Z_A(t)}
 \qquad(S\subseteq E(A)).                          \tag{SO.9}
```

The same calculation with arbitrary real edge amplitudes proves that the
complete inhomogeneous contraction box OI.3, including all convex endpoint
derivatives and edge-flip subgradients, is identical at every query word.
Thus convexity supplies no missing stationarity equation in the `v`
direction: `v` moves between gauge-equivalent vertices, whereas the convex
chords between those vertices lie in the interior of the coupling cube and
need not be flat.

The one-vertex Bellman atlas is gauge-covariant rather than pointwise
different.  Deleting vertex `i` from `A^v` switches the deletion by
`v_{-i}`; its proposed incident row is relabelled by

```math
 b\longmapsto v_i(v_{-i}\odot b).                  \tag{SO.10}
```

Under this bijection the entire function
`b -> R_C^aug(b;t,t)`, all its minimizing rows, and every reinsertion defect
are preserved.  Hence retaining the Bellman table modulo its natural gauge
gives exactly the same data for all `v`.  Retaining it with a fixed external
gauge instead retains an exponentially indexed response table and is not a
low-information statistic.

Finally, scalar envelope information such as the optimized child pressure,
its temperature subgradients, adjacent-order extension deficits, and sector
weights is switching invariant.  A parent-order envelope can of course
select among the values in (SO.4), but doing so imports a target-order bridge
optimization.  It is not an internal child-minimality identity.

We may summarize the exact information split as

```text
 internal optimizer atlas at A^v
   = one switching-invariant/covariant object;

 degree-two carrier response at v
   = one externally gauge-fixed observation of A^v.
```

Therefore the known internal identities cannot output a coherent favorable
`v` unless one adds a genuinely gauge-covariant observable coupling the
child to the external direction `y`.  The sector--Gram tangent is such an
observable at row order two; ST.2 shows that its uncontrolled continuation
contains all even row orders.

This blindness has a precise selector formulation.  Suppose `S(A,D)` is
any statistic invariant under switching both children, and a decoder using
only `S(A,D)` is claimed to output one projective query pair
`sigma(A,D) in Gamma_m times Gamma_n`.  A physically well-defined selector
must be equivariant:

```math
 \sigma(A^a,D^b)=(a,b)\sigma(A,D).                 \tag{SO.10a}
```

But switching invariance of `S` makes the decoder's left side equal to
`sigma(A,D)` for every `(a,b)`.  The regular action of
`Gamma_m times Gamma_n` on itself has no fixed point when either projective
group is nontrivial.  Thus:

**Corollary SO.1a (no invariant coherent-direction selector).**  No
switching-invariant child statistic can, by itself, output a gauge-covariant
single carrier direction on a nontrivial child switching orbit.  Such a
statistic may still decide a switching-invariant scalar such as the *range*
of the orbit response; the corollary does not rule that out.  It rules out
using invariant optimizer data alone as the coherent direction required by
branch (iii).

The qualification is important.  The switching classes of the two child
matrices are finite descriptions of the input, and in that tautological
sense determine the response `mathscr R`.  But evaluating (SO.7c) from them still means
solving the exact parent pressure on the orbit.  Alternatively, storing a
gauge-covariant scalar for every orbit point stores the response table
itself.  SO.1 supplies neither an intermediate orbit statistic nor a closed
update law.  A genuine reduction would have to prove that a strictly
smaller covariant object controls the orbit minimum or range.

## 3. Zero-defect actual children with a nonzero maximal coefficient

The preceding distinction is realized by the smallest actual child pair.
Let both children be the positive order-two edge.  They are exact
contracted-temperature minimizers for every raw temperature `t>0`.  Use
relative orientation `+`, bridge amplitude `u>0`, and put

```math
 \rho=\tanh t,\qquad r=\tanh u,qquad
 c=\rho^2r^2,\qquad d=r^4.                         \tag{SO.11}
```

For a bridge `B in {+-1}^{2\times2}`, let

```math
 \alpha=B_{11}B_{22},\qquad \zeta=B_{12}B_{21}.
```

The exact forward likelihood from OI.20 is

```math
 p(B)=1+c(\alpha+\zeta)+d\alpha\zeta.              \tag{SO.12}
```

Take the carrier direction `y=(1,1)`.  Directly from (SO.2), under the two
row laws with query word `(1,1)` the probabilities of
`(alpha,zeta)=(++,+-,-+,--)` are

```math
 {1\over16}(9,3,3,1),                              \tag{SO.13}
```

whereas under query word `(1,-1)` they are

```math
 {1\over16}(1,3,3,9).                              \tag{SO.14}
```

**Theorem SO.2 (zero optimizer defect does not kill the first even row
mode).**  For the preceding pair,

```math
 \boxed{
 \widehat R_{A,D;y}(\{1,2\})
 ={1\over4}\log
 {1+2\rho^2r^2+r^4\over1-2\rho^2r^2+r^4}>0.}      \tag{SO.15}
```

At the same time every edge-replacement comparison, every inhomogeneous
contraction comparison, every actual-row deletion/reinsertion comparison,
and the complete neutral extension table is an equality; both sector biases
also vanish.

*Proof.*  The three likelihood values in (SO.12) are

```math
 p_{++}=1+2c+d,\qquad p_{+-}=p_{-+}=1-d,\qquad
 p_{--}=1-2c+d.                                    \tag{SO.16}
```

The last value is positive since
`1-2\rho^2r^2+r^4=(r^2-\rho^2)^2+1-\rho^4>0`.
Equations (SO.13)--(SO.16) give

```math
 R(1,1)-R(1,-1)
 ={1\over2}\log{p_{++}\over p_{--}}.              \tag{SO.17}
```

By (SO.5), the response has only its constant and `{1,2}` Walsh
coefficients, so (SO.17) proves (SO.15).  Strict positivity follows from
`rho,r>0`.  The zero-defect assertions are OI.16--OI.18: the normalized
order-two contraction pressure is identically one, every sign of its only
edge is switching equivalent, deletion leaves the unique order-one child,
and every order-three signing has augmented pressure `cosh(t)^3`. `square`

This is directly stronger for the ST.2 question than merely observing a
nonzero connected fourth cumulant: the surviving object is the actual
fixed-degree query coefficient itself, at the largest allowed row order for
the example.  For fixed `y=(1,1)`, (SO.15) also makes the mixed-sign word the
unique favorable element of the two-point projective query cube.  Switching
the left child by that mixed word swaps the favorable element, while every
switching-invariant optimizer datum remains unchanged.  Thus the selector
obstruction in Corollary SO.1a is witnessed inside the actual optimizing
class, not only abstractly.

## 4. Consequence for `L_balanced-product-phase`

The audit yields exactly one new universal restriction on the ST.2 response:
odd row orders vanish.  It also proves that the first admissible even mode
survives when all existing internal optimizer defects vanish.  Therefore:

1. internal flip/contraction slack cannot control the even query tail;
2. neutral deletion/reinsertion slack cannot control it;
3. convex endpoint or scalar envelope stationarity adds no new equation;
4. switching symmetry identifies the query with a flat minimizer orbit but
   does not select a favorable point on that orbit.

A useful next identity would have to be genuinely **multirow and
gauge-covariant**, with state smaller than the full switching-orbit response.
For example, it would need either to synchronize the even coefficients of
(ST.8), prove their contribution to `min_v R(v)` is `o(N)`, or give a closed
recursion for that minimum.  None of the optimizer identities audited here
does so.

Accordingly this is a rigorous route-specific no-go, not a RESET.  Its
strongest conclusion is the group-theoretic one: switching-invariant
optimizer data cannot supply the coherent gauge required by branch (iii),
while the first admissible even response mode survives zero optimizer
defect.  It does not exclude an invariant low-information statistic deciding
branch (i) or (ii), nor an as-yet-unknown covariant statistic deciding branch
(iii).  It therefore does not change the SML: bounded-row-degree cross-row
synchronization remains missing, and a scalable actual-minimizer bound or
collision is still needed.
