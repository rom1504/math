# Boolean spectral deficit for cyclic two-fiber conferences

Status: proved universal arithmetic improvement over the conference spectral
bound; exact Boolean audits at orders 6 and 18; an exact order-38 witness; and
an exact certificate that the order-38 cyclic example is non-Paley.  The
improvement is only lower order and does not yet supply a scalable landing
theorem.

## 1. Exact cyclic Fourier decomposition

Let `s=2k^2+1`, `q=2s-1=4k^2+1`, and let

```math
S=\begin{pmatrix}A&C\\C^{\mathsf T}&-A\end{pmatrix},
\qquad S^2=qI_{2s},                                  \tag{BF1}
```

where `A` is symmetric circulant and `C` is circulant.  If `alpha_j` and
`gamma_j` are their cyclic Fourier multipliers, then every frequency has the
Hermitian symbol

```math
\widehat S(j)=
\begin{pmatrix}
\alpha_j&\gamma_j\\
\overline{\gamma_j}&-\alpha_j
\end{pmatrix},
\qquad \alpha_j^2+|\gamma_j|^2=q.                   \tag{BF2}
```

For a Boolean pair `w=(x,y)` and unitary Fourier transforms `X,Y`, its energy
is exactly

```math
H_S(w)={1\over2}\sum_j
\left[\alpha_j(|X_j|^2-|Y_j|^2)
+2\operatorname{Re}(\overline{X_j}\gamma_jY_j)\right]. \tag{BF3}
```

Write `lambda=sqrt(q)` and `P_+ = (I+S/lambda)/2`,
`P_-=(I-S/lambda)/2`.  For `epsilon=sign(H_S(w))`,

```math
\boxed{
\|P_{-\epsilon}w\|_2^2
=s-{|H_S(w)|\over\lambda}.}                         \tag{BF4}
```

Thus concentration in one sign of the conference spectrum is not a vague
Fourier condition: it is exactly the distance in (BF4).  Frequency by
frequency, exact positive concentration would require

```math
(\lambda-\alpha_j)X_j=\gamma_jY_j                  \tag{BF5}
```

and its conjugate companion wherever the Fourier pair is supported.  No
Boolean vector can satisfy this exactly, because `Sw` is integral while
`lambda` is irrational.  The quantitative question is how closely it can be
satisfied simultaneously across the Boolean Fourier support.

## 2. Universal arithmetic deficit identity

There is an exact answer sharper than ordinary Cauchy--Schwarz, although not
sharp enough by a constant proportion.  For any symmetric conference signing
of order `N=4k^2+2`--not only the cyclic examples--define the signed local
fields

```math
t_i=w_i(Sw)_i,
\qquad u_i=\operatorname{sign}(H_S(w))t_i.           \tag{BF6}
```

Every `u_i` is odd, because each row has `N-1` signs.  Conference orthogonality
gives

```math
\sum_i u_i=2|H_S(w)|,
\qquad \sum_i u_i^2=Nq.                             \tag{BF7}
```

For an odd integer `u`, there is no integer strictly between the consecutive
odd numbers `2k-1` and `2k+1`, hence

```math
(u-(2k-1))(u-(2k+1))\ge0.                           \tag{BF8}
```

Summing (BF8) and using (BF7) proves the exact identity

```math
\boxed{
2ks-|H_S(w)|
={1\over8k}\sum_{i=1}^{2s}
(u_i-(2k-1))(u_i-(2k+1)).}                          \tag{BF9}
```

In particular,

```math
\operatorname{cap}(S)\le2ks.                       \tag{BF10}
```

There are an odd number of edges, so every Boolean energy is odd, whereas
`2ks` is even.  Therefore

```math
\boxed{\operatorname{cap}(S)\le2ks-1.}             \tag{BF11}
```

This improves the spectral bound `s sqrt(4k^2+1)` by

```math
s\bigl(\sqrt{4k^2+1}-2k\bigr)+1
={s\over\sqrt{4k^2+1}+2k}+1=\Theta(k).              \tag{BF12}
```

It also completely characterizes near equality.  The summands in (BF9) are
nonnegative multiples of eight.  A coordinate pays no deficit precisely when
its oriented local field is `2k-1` or `2k+1`; all other local fields pay at
least eight.  If `D=2ks-|H|`, the total integer penalty is exactly `8kD`.

Combining (BF4) and (BF11) gives

```math
\|P_{-\epsilon}w\|_2^2
\ge s-{2ks-1\over\sqrt{4k^2+1}}
={s\over\sqrt{4k^2+1}(\sqrt{4k^2+1}+2k)}
+{1\over\sqrt{4k^2+1}}.                             \tag{BF13}
```

The right side tends only to `1/4`.  It is a constant squared-distance
obstruction, not a positive fraction of the Boolean norm `2s`.  Consequently
the arithmetic theorem explains exactly why a uniform normalized improvement
over `1/2` cannot follow from local-field quantization alone.

## 3. Exact Boolean audits of the cyclic ASDS certificates

The embedded cyclic certificates give the following exact results.

For `k=1`, order six, exhaustive enumeration gives cap five.  An extremizer
has oriented local fields

```text
1 (five times), 5 (once),
```

and deficit `2ks-|H|=1`.

For `k=2`, order 18, exhaustive enumeration of all `2^17` projective spins
gives

```math
\operatorname{cap}(S)=33=2ks-3.                    \tag{BF14}
```

An extremizer has six oriented local fields equal to one and twelve equal to
five.  Its penalty sum in (BF9) is exactly 48.

For `k=3`, order 38, the explicit little-endian positive-bit spin

```text
1e1d4bf05d
```

has exact energy 109:

```math
109=2ks-5.                                          \tag{BF15}
```

Its oriented local-field multiset is

```text
1 (once), 5 (30 times), 9 (six times), 13 (once),
```

whose total penalty is `120=8*3*5`.  Equation (BF15) is a certified lower
bound, not an exact cap claim.  The deficits `1,3,5` at `k=1,2,3` suggest an
`O(k)` rather than project-scale spectral deficit, but three finite cases do
not establish an infinite family.

## 4. Paley equivalence and the genuinely new order-38 example

Exact signed-permutation certificates show that the `k=1` and `k=2`
two-circulant matrices are presentations of the prime Paley conferences for
`q=5` and `q=17`.  Their near-spectral Boolean behavior therefore supplies no
new route beyond the established prime-Paley saturation mechanism.

The `k=3` certificate is different.  Gauge a chosen root so all incident
edges are positive, delete it, and form the graph of negative edges.  For the
cyclic two-fiber order-38 matrix this graph has

```text
615 four-cliques and 65 five-cliques.
```

For the Paley `q=37` conference graph it has

```text
555 four-cliques and no five-cliques.
```

The Paley switching automorphism group is vertex-transitive, so the root
choice does not affect its isomorphism class.  Clique counts are graph
isomorphism invariants.  This proves that the cyclic `k=3` certificate is not
switching/permutation equivalent to Paley.

## 5. Research judgment

The arithmetic identity (BF9) is a theorem-level Boolean improvement and an
exact characterization of spectral-sign concentration.  It does not give a
constant normalized deficit: the order-38 non-Paley witness already lies
within `O(1)` squared distance of one spectral sign.

### 5.1 Exact self-indexed ASDS mapping and product obstruction

Put `alpha=a+delta_0`, so that `alpha` and `c` are both sign sequences.
The complementary identity for `a,c` becomes

```math
R_\alpha(h)+R_c(h)=2a_h\in\{\pm2\}\qquad(h\ne0),   \tag{BF16}
```

If `P={j:a_j=-1}` and `Q={j:c_j=-1}`, expanding (BF16) gives precisely

```math
N_P(h)+N_Q(h)+\mathbf1_P(h)=k^2-k\qquad(h\ne0),    \tag{BF17}
```

with `P=-P`, `|P|=k(k-1)`, and `|Q|=k^2`.  Thus this is an almost
supplementary difference set with the additional **self-indexing** term
`1_P(h)` and symmetry of `P`.

General ASDS and the amicable subclass corresponding to optimal quaternary
sequences are developed by Armario--Flannery
([primary source](https://arxiv.org/abs/1911.08828)).  Their general existence
statements do not assert (BF17), `P=-P`, and the prescribed weights
simultaneously.  In fact the `k=2,3` certificates fail the extra cross-
correlation symmetry used in the OQS equivalence, so importing the OQS
families here would be an invalid mapping.

There is an exact obstruction to the most obvious product of the actual
binary pair.  For two certificates, direct tensoring gives

```math
R_{\alpha\otimes\alpha'}(h,h')
+R_{c\otimes c'}(h,h')
=R_\alpha(h)R_{\alpha'}(h')+R_c(h)R_{c'}(h').       \tag{BF18}
```

At an axial shift `(h,0)` with `h\ne0`, the right side is
`s'(R_alpha(h)+R_c(h))=2s'a_h`, whose magnitude is `2s'`, rather than the
required two.  Therefore direct tensor/CRT multiplication does not preserve
the self-indexed ASDS equation.  The two missing mixed autocorrelation terms
show exactly why a four-channel interleaving or phase-cancellation operation
would be required.

The smallest Boolean four-channel repair also fails, for a more structural
reason.  Given sign pairs `(alpha,c)` and `(beta,d)`, complex multiplication
followed by a 45-degree phase rotation gives the Boolean arrays

```math
\begin{aligned}
e&={\alpha\otimes\beta-c\otimes d
        -\alpha\otimes d-c\otimes\beta\over2},\\
f&={\alpha\otimes\beta-c\otimes d
        +\alpha\otimes d+c\otimes\beta\over2}.
\end{aligned}                                       \tag{BF19}
```

At every coordinate each numerator in (BF19) is `+/-2`, so both outputs are
Boolean.  Define

```math
S_{\alpha,c}(h)=R_\alpha(h)+R_c(h),
\qquad
K_{\alpha,c}(h)=R_{c,\alpha}(h)-R_{\alpha,c}(h).    \tag{BF20}
```

Expansion of all four tensor channels, with the mixed terms retained, gives
the exact correlation law

```math
\boxed{
R_e(h,h')+R_f(h,h')
={S_{\alpha,c}(h)S_{\beta,d}(h')
 -K_{\alpha,c}(h)K_{\beta,d}(h')\over2}.}           \tag{BF21}
```

At zero shift `K(0)=0` and `S_{beta,d}(0)=2s'`.  Hence every nonzero axial
shift satisfies

```math
R_e(h,0)+R_f(h,0)=s'S_{\alpha,c}(h)=2s'a_h.         \tag{BF22}
```

The mixed four-channel term cannot cancel an axial shift because its second
factor is exactly zero there.  Thus the violation grows linearly with the
other factor size.  The same obstruction appears after identifying a cyclic
CRT product: every nontrivial direct product has shifts supported in one
factor.  The reproducer verifies (BF21) at every shift for the `k=1` and
`k=2` certificates, obtaining axial magnitude 18 instead of two.

Equations (BF19)--(BF22) rule out the natural Gray/complex multiplication
closure.  In fact the axial argument covers the full pointwise two-channel
architecture.  Write `p(u)=(alpha_u,c_u)^T`.  Any orthogonal map carrying the
two-dimensional Boolean square to itself is a signed permutation matrix.  If
an interleaving uses an arbitrary signed permutation `T_v` at each coordinate
of the second factor,

```math
q(u,v)=T_vp(u),                                      \tag{BF23}
```

then at an axial shift

```math
\sum_{u,v}q(u,v)^Tq(u+h,v)
=\sum_v\sum_up(u)^TT_v^TT_vp(u+h)
=s'S_{\alpha,c}(h).                                 \tag{BF24}
```

Thus **every fiberwise signed-permutation/phase interleaving has the same
linear axial obstruction**, independently of how the second certificate
chooses `T_v`.  This includes (BF19).

Equations (BF23)--(BF24) do not rule out every conceivable construction.  A
surviving law must make the channel transformation depend nontrivially on the
first coordinate as well, mix group coordinates so that endpoints of an
axial shift see different transformations, or add further channels whose
cross terms remain nonzero on the axial subgroups.  Such a law is no longer a
two-state pointwise composition of the two input certificates.

The next parameter after `k=3` is `k=4`, `s=33`, conference order 66.  A
general order-66 symmetric conference matrix is known via Gritsenko's
conference graph, constructed using a more elaborate array of circulant
blocks ([primary source](https://arxiv.org/abs/2102.05432)).  This does not
supply the two-circulant form (BF1), the quotient row sums `(8,-8,1)`, or the
self-indexed ASDS.  Thus general conference existence is strictly weaker than
the family needed here.

### 5.2 Current route classification

The surviving question is narrow: determine whether the non-Paley
self-indexed ASDS at `k=3` belongs to an infinite construction closed under a
difference-set, interleaving, or composition operation, and whether (BF15)'s
`O(k)` deficit persists along it.  Such a family would be a scalable
obstruction to obtaining Boolean slack from two-fiber conference structure.

No such closure follows from the conference identity, general ASDS
existence, or direct tensoring.  Without an infinite family or a
cross-order operation, the order-38 example remains finite evidence rather
than primary convergence progress.

## Reproduction

```bash
.venv/bin/python computations/two_fiber_cyclic_conference.py \
  --output computations/results/two_fiber_cyclic_conference.json

.venv/bin/python computations/audit_two_fiber_boolean_spectral_deficit.py \
  --output computations/results/two_fiber_boolean_spectral_deficit.json
```

The second program independently verifies the cap and witness claims, the
local-field deficit identity, the explicit signed-permutation certificates,
and the Paley/non-Paley clique invariants using exact integer arithmetic.
