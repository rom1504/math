# Optimizer identities do not coerce the actual-child cluster mass

Status: **rigorous actual-child identity audit and exact zero-defect
falsifier**.  This note asks only whether the presently available exact
edge-, flip-, and vertex-minimality identities can force a power-saving
connected-cluster bound for the balanced actual-child law.  They cannot do
so through their optimality defects.  Two order-two exact minimizers make
every such defect vanish, including the complete inhomogeneous flip box and
the complete neutral one-vertex extension table, while their balanced
two-row bridge has a nonzero fourth connected cumulant and strictly positive
physical canonical interaction.

The example is finite.  It rules out a coercive deduction from the existing
optimizer defects, not a new asymptotic theorem using genuinely multirow
optimizer structure.  In particular, it is a route-specific strike, not an
asymptotic counterexample to `J_N=o(N)`.

## 1. Inventory of the exact optimizer identities

Let `A` be a signing of order `d`, let `E=binom([d],2)`, and put

```math
 \overline Z_A(t)=E_{x,\tau}e^{t\tau H_A(x)},
 \qquad
 d\nu_A={e^{t\tau H_A}\over\overline Z_A(t)}dU.
 \tag{OI.1}
```

The exact identities which actually use pressure minimality are the
following.

1. **Arbitrary edge replacement.**  If `A` is an exact minimizer, then for
   every `S subset E`,

   ```math
   {\overline Z_{A^S}(t)\over\overline Z_A(t)}
   =E_{\nu_A}\exp\left{-2t\tau
     \sum_{e\in S}a_ex_e\right}\ge1.               \tag{OI.2}
   ```

   This is AC.32.  Averaging (OI.2) with independent, edge-dependent flip
   probabilities gives FC.8,

   ```math
   {\overline Z_A(s)\over\prod_e\cosh s_e}
   \ge {\overline Z_A(t\mathbf1)\over(\cosh t)^{|E|}}
   \quad(s\in[-t,t]^E).                              \tag{OI.3}
   ```

   Conversely, the corner values of (OI.3) recover every comparison in
   (OI.2).  Thus the full inhomogeneous box is equivalent to the complete
   sign-flip pressure table; it is not a compression.  Its one-edge tangent
   is only

   ```math
   a_eE_{\nu_A}(\tau x_e)\le\tanh t.                 \tag{OI.4}
   ```

   Cut flips, including a complete star, are switching equivalences and
   give equality in (OI.2).  Star-subset inequalities are exactly the
   corresponding subfamily of (OI.2), not an independent constraint.

2. **Vertex deletion and reinsertion.**  If `C` has order `d-1` and `b` is
   a proposed incident row, the exact extension identity is

   ```math
   \overline Z_{C\oplus b}(t)
   =\overline Z_C(t)R_C^{\rm aug}(b;t,t).            \tag{OI.5}
   ```

   Every row deleted from an exact order-`d` minimizer minimizes the whole
   function `b -> R_C^aug(b;t,t)` for its own deletion.  Equivalently,

   ```math
   F_d(t)=\min_C\{\log\overline Z_C(t)
                    +\min_b\log R_C^{\rm aug}(b;t,t)\}.       \tag{OI.6}
   ```

   This is DER.2/DER.12 (equivalently EE.3).  It is exact, but its state is
   the complete one-row response table.  It does not compare two rows
   inserted jointly.

3. **A minimizing base and an arbitrary new vertex.**  If `D` itself is an
   exact order-`d` minimizer, its neutral normalized extension response
   obeys

   ```math
   z_D^0(b;t,t)\ge e^{-\delta_d(t)},
   \qquad
   \delta_d(t)=d\log\cosh t-{F_{d+1}(t)-F_d(t)\}.    \tag{OI.7}
   ```

   The adjacent deficits telescope:

   ```math
   \sum_{k=1}^{d-1}\delta_k(t)
   ={d\choose2}\log\cosh t-F_d(t).                  \tag{OI.8}
   ```

   Bias balancing converts (OI.7) into the uniform canonical-row density
   bound of Theorem 37.32.  Neither (OI.7) nor (OI.8) contains a joint
   cross-row cumulant.

There is one immediate multirow comparison, but it invokes the target-order
optimum.  For exact children `A,D` of orders `m,n`, put

```math
 L_\epsilon(B)=\log\overline Z_N(A,\epsilon D,B;t),
 \qquad
 p_\epsilon(B)
 ={e^{L_\epsilon(B)-L_\epsilon(0)}\over(\cosh t)^{mn}}.
 \tag{OI.9}
```

Uniform bridge averaging gives `E_U p_epsilon=1`.  Since every completed
signing is an order-`N` competitor,

```math
 \boxed{
 p_\epsilon(B)\ge e^{-G_\epsilon}\quad\hbox{for every }B,
 \qquad
 G_\epsilon=L_\epsilon(0)+mn\log\cosh t-F_N(t)\ge0.}
 \tag{OI.10}
```

Consequently

```math
 {1\over\lambda}\log E_U p_\epsilon^{-\lambda}
 \le G_\epsilon,
 \qquad
 D(U\Vert q_{\lambda,\epsilon})\le\lambda G_\epsilon.        \tag{OI.11}
```

The second inequality uses `E_U\log p_epsilon<=0`.  This is a genuine exact
optimizer floor, but `G_epsilon` contains `F_N(t)`.  Proving a power saving
for it is already a target-order almost-subadditivity estimate, so (OI.10)
cannot be used to manufacture the missing recurrence without circularity.
It also controls the uniform inverse escort, not the canonical
`D(r||q)` when the inverse row factors are nonuniform.

Equations (OI.2)--(OI.8) are therefore all internal-edge or one-vertex
statements.  The connected variables in Theorem 37.33,
`Q_(ij)=sX_iY_j`, first couple two different child vertices and two bridge
rows.  The next theorem shows that this mismatch is real even when every
available optimizer defect vanishes.

## 2. Zero-defect two-row synergy

Let `A=D=E_2` be the positive one-edge signing.  Take two copies, so
`m=n=2`, `N=4`, and let the common internal and bridge amplitude be

```math
 t=u={\beta\over\sqrt4}>0.
```

Both relative orientations are balanced because both sector biases vanish.

**Theorem OI.1 (exact zero-defect synergy).**  For every `t,lambda>0`:

1. both children are exact contracted-temperature pressure minimizers;
2. every edge-set flip comparison (OI.2), every inhomogeneous contraction
   comparison (OI.3), every one-edge tangent (OI.4), every actual-row
   deletion/reinsertion comparison (OI.5)--(OI.6), and the complete neutral
   extension table (OI.7) are equalities;
3. nevertheless the balanced zero-bridge child channel has

   ```math
   K_\epsilon=2\tanh^4t>0                           \tag{OI.12}
   ```

   and the fourth connected cross-row cumulant

   ```math
   \operatorname {cum}(Q_{11},Q_{12},Q_{21},Q_{22})
   =1-\tanh^4t>0;                                   \tag{OI.13}
   ```

4. its exact physical canonical interaction is

   ```math
   \boxed{
   \mathcal J(t,\lambda)
   =\log\left\{{e^{-3\lambda d/4}
                         +3e^{\lambda d/4}\over4}\right\}>0,
   \qquad
   d=\log{1+3q\over1-q},\quad q=\tanh^4t.}          \tag{OI.14}
   ```

   In particular, as `t -> 0` at fixed `lambda`,

   ```math
   \mathcal J(t,\lambda)
   ={3\over2}\lambda^2t^8+O_\lambda(t^{10}).        \tag{OI.15}
   ```

Thus no universal estimate of either the physical interaction or the
connected cluster mass by a coercive function of the presently available
edge/vertex optimality defects can hold.  This remains true if one retains
the *complete* inhomogeneous flip box and the *complete* neutral
one-vertex extension table: both objects are constant in this example.

### Proof

For an order-two signing,

```math
 \overline Z_2(t)=\cosh t,
 \qquad
 {\overline Z_2(s)\over\cosh s}=1.                 \tag{OI.16}
```

Changing the only edge sign does not change the augmented pressure, so
(OI.2)--(OI.3) are equalities.  Under its augmented Gibbs law,
`E(\tau x_1x_2)=tanh t`, proving equality in (OI.4).  Deleting either
vertex leaves the unique order-one signing and a constant normalized
reinsertion response.

Every order-three signing has

```math
 \overline Z_3(t)
 ={\cosh(3t)+3\cosh t\over4}
 =\cosh^3t.                                        \tag{OI.17}
```

Hence

```math
 F_1(t)=0,
 \quad F_2(t)=\log\cosh t,
 \quad F_3(t)=3\log\cosh t,
 \quad \delta_1(t)=\delta_2(t)=0,                  \tag{OI.18}
```

and `z_(E_2)^0(b;t,t)=1` for every proposed two-bit row `b`.  This proves
the complete zero-defect assertion.  Also
`Z_(E_2)^+(t)=Z_(E_2)^-(t)=cosh t`, so both sector biases are zero.

It remains to compute the joint bridge.  Put

```math
 \rho=\tanh t,\qquad r=\tanh u,
 \qquad
 \alpha=B_{11}B_{22},\qquad
 \zeta=B_{12}B_{21}.                                \tag{OI.19}
```

The exact high-temperature expansion on the completed `K_4` has only
three nonconstant even Eulerian subgraphs: the two four-cycles using both
internal edges and the all-bridge rectangle.  Therefore the forward bridge
likelihood is

```math
 \boxed{
 p_{\epsilon,u}(B)
 =1+\epsilon\rho^2r^2(\alpha+\zeta)+r^4\alpha\zeta.}
 \tag{OI.20}
```

The two characters `alpha,zeta` are independent fair signs under the fair
bridge law.  Averaging (OI.20) over either complete bridge row gives one,
so both exact row likelihoods are constant and the canonical inverse-row
product is exactly `r_row^(tensor2)=U_B`.

At physical amplitude `u=t`, put `q=rho^4`.  For either orientation the
likelihood has the same distribution,

```math
 p(B)=
 \begin{cases}
 1+3q,&\text{with probability }1/4,\\
 1-q,&\text{with probability }3/4.
 \end{cases}                                      \tag{OI.21}
```

Since the canonical product is uniform,

```math
 \mathcal J
 =D(U_B\Vert q_\lambda)
 =\log E_Up^{-\lambda}+\lambda E_U\log p.          \tag{OI.22}
```

Centering the two values in (OI.21) gives (OI.14).  It is strictly positive
because the centered likelihood is nonconstant and exponential Jensen is
strict.  Taylor expansion first in `d` and then in `t` gives

```math
 \mathcal J={3\lambda^2d^2\over32}+O_\lambda(d^3),
 \qquad d=4q+O(q^2),
 \qquad q=t^4+O(t^6),                               \tag{OI.23}
```

which is (OI.15).

For completeness, under the zero-bridge sector law let
`Q_(ij)=sX_iY_j`.  The sector is fair and

```math
 E(X_1X_2\mid s)=s\rho,
 \qquad
 E(Y_1Y_2\mid s)=\epsilon s\rho.                  \tag{OI.24}
```

All same-row and same-column distinct pair correlations of `Q` vanish,
whereas both opposite-corner correlations equal `epsilon rho^2`.  This
gives (OI.12).  Moreover the product of all four distinct `Q` variables is
identically one.  Of the three Wick pairings, only the opposite-corner
pairing is nonzero, and it contributes `rho^4`; hence (OI.13).  In the
ordered-tuple normalization of Theorem 37.33, the `4!` permutations of this
tuple alone give the unconditional fourth-order lower bound

```math
 \mathfrak C_4(u)\ge u^4(1-\rho^4)>0.              \tag{OI.25}
```

Whenever the full absolute cluster series converges, the same is a lower
bound on `mathfrak C_(>=4)(u)`.  This finishes the proof. `square`

### Normalization cross-check

The three normalizations in (OI.12)--(OI.15) agree with the canonical
definitions in EO.8, SP.2, and CC.6.

- EO.8 sums `i<k` and then all **ordered** column pairs `(j,l)`.  Here there
  is one row pair.  The two off-diagonal column pairs `(1,2)` and `(2,1)`
  each have `Gamma^2=rho^4`, while `(1,1)` and `(2,2)` vanish.  Hence
  `K=2rho^4`, with no factor of two missing.
- SP.2 sums ordered edge tuples and divides the order-four term by `4!`.
  The tuple of the four distinct bridge edges has `4!` permutations, all
  with absolute cumulant `1-rho^4`.  Their total contribution is therefore
  exactly `u^4(1-rho^4)`, as in (OI.25).
- CC.6 defines the canonical interaction as the centered negative log MGF
  of `G=p/prod_i p_i` under the canonical row product.  In this example
  every `p_i` is one, so `G=p` and the canonical row product is `U_B`.
  Thus CC.6 is precisely (OI.22), not an uncentered inverse work.  The
  high likelihood `1+3q` occurs on four of the sixteen bridges and the low
  likelihood `1-q` on the other twelve, yielding (OI.14).  Finally
  `d=4q+O(q^2)`, `q=t^4+O(t^6)`, and the centered two-point variance is
  `3d^2/16`; the cumulant prefactor `1/2` gives
  `3lambda^2d^2/32=(3/2)lambda^2t^8+O(t^10)`.

## 3. What the falsifier excludes

Theorem OI.1 is stronger than observing that a generic bridge need not be
row-product.  On the two actual minimizing children it simultaneously has

```text
complete edge-flip slack                 = 0,
complete inhomogeneous-contraction slack = 0,
all neutral Bellman/reinsertion slack    = 0,
adjacent extension deficits              = 0,
sector biases                            = 0,
balanced row-factor complexity           = 0,
fourth cross-row connected cumulant      > 0,
physical canonical interaction J         > 0.
```

Accordingly, none of the following can be a valid route to a power saving
without an additional genuinely multirow theorem:

- summing one-edge or star-flip slacks;
- iterating the neutral one-vertex lower envelope;
- charging the cluster tail to adjacent Bellman deficits;
- treating equality in the full contraction box as a weak-dependence
  certificate.

This also explains why the sign of (OI.2) is unusable for cluster
superconcentration.  It is a lower-MGF inequality; large fluctuations make
it easier to satisfy.  At equality, cross-row synergy can still first appear
when two vertices are inserted together.

The stored exact-minimizer falsifiers point in the same direction at richer
orders.  DER.26 shows that neutral deletion optimality need not survive the
other child's sector bias.  The two order-eight exact minimizers in EO.4
have identical complete radial pressure data but different sector--Gram
tangents.  Theorem OI.1 isolates the still sharper issue relevant here:
even *zero* edge/vertex optimality defect does not suppress the first joint
cluster.

## 4. Frontier disposition

No power-saving connected-cluster bound follows from the exact optimizer
identities audited here.  The only immediate multirow optimizer floor,
(OI.10), imports the target-order optimum and is therefore not a Level-6
mechanism.  Any continuation of preferred route 1 would need a new
multi-vertex stability identity whose information footprint is provably
smaller than the bridge landscape; edge/vertex slack cannot be its control
parameter.

This is a **finite, route-specific strike** for coercing the cluster mass
from the existing optimizer defects.  It is not an overall campaign strike:
it does not exclude a new asymptotic multirow identity, a directional
product certificate, or an orientation-uniform promotion theorem.  By
itself it does not change `L_balanced-product-phase` and does not make Level
6 credible.  It justifies freezing further edge/vertex-slack refinements;
the disposition of the broader adversarial-statistical-mechanics branch
depends on whether another route in the same campaign supplies a genuine
reset.
