# Sparse selectors cannot manufacture a sublinear cross-order defect

Status: **proved method-class obstruction for actual optimizing children**.
The results below do not prove or disprove the desired recurrence.  They
show quantitatively that two broad bridge-selection mechanisms cannot create
an `o(N)` defect unless the uniform exact bridge identity already has one.

Fix `N=m+n`, `d=mn`, `t=beta/sqrt(N)`, and exact pressure-minimizing
children `A,D`.  For orientation `epsilon in {+-1}` and sign bridge `B`, put

```math
L_\epsilon(B)=\log\overline Z_N
 \left(\begin{pmatrix}A&B\\B^{\mathsf T}&\epsilon D\end{pmatrix},t\right),
\qquad T=P_m(\beta)+P_n(\beta).
```

Let `U` be uniform on the `2^(d+1)` pairs `(epsilon,B)`.  The normalized
positive joint output law is

```math
\Pi(\epsilon,B)={e^{L_\epsilon(B)}
 \over 2^{d+1}(\cosh t)^d\overline Z_A(t)\overline Z_D(t)},
\qquad D=D_{\rm KL}(U\Vert\Pi).
```

The centered identity from
[`cross_order_centered_channel_identity.md`](cross_order_centered_channel_identity.md)
holds after the **joint uniform** orientation/bridge average:

```math
\boxed{
\mathbb E_{(\epsilon,B)\sim U}L_\epsilon(B)-T
=C^0_{m,n}(\beta)+\Gamma_A+\Gamma_D-D.}                \tag{1}
```

It does not generally hold orientation by orientation: averaging a fixed
orientation over bridges leaves a term
`epsilon E sinh(tH_A) E sinh(tH_D)`.  Nor is the output orientation
marginal generally uniform.  All selectors below therefore keep the
orientation exactly fair and apply their transport separately within each
orientation.  An orientation bias can have linear leverage and is outside
these no-go statements.

## 1. Fixed parity selectors have only sparse Hamming leverage

For either fixed orientation and arbitrary real bridges, differentiation of
log partition gives

```math
|L_\epsilon(B)-L_\epsilon(B')|
\le t\|B-B'\|_1.                                      \tag{2}
```

In particular one flipped sign coordinate costs at most `2t`.
Consequently, for arbitrary bridge laws `q_epsilon` and the Hamming
Wasserstein distance,

```math
\boxed{
|\mathbb E_{q_\epsilon}L_\epsilon-
  \mathbb E_{U_B}L_\epsilon|
\le2tW_1^{\rm Ham}(q_\epsilon,U_B).}                  \tag{2a}
```

After averaging a fair orientation, (1) gives the direct cross-order arrow

```math
\boxed{
{1\over2}\sum_\epsilon\mathbb E_{q_\epsilon}L_\epsilon-T
\le\omega_N
\Longrightarrow
E_{m,n}(\beta)
\le\omega_N
 +{\beta\over\sqrt N}\sum_\epsilon
 W_1^{\rm Ham}(q_\epsilon,U_B).}                     \tag{2b}
```

Thus every selector at average Hamming transport `o(N^(3/2))` is
quantitatively redundant.  Conversely, improving a uniform certificate by
`cN` requires average transport at least
`(c/(2beta)-o(1))N^(3/2)`.  This conclusion is independent of how the law is
encoded; the parity and product results below make it checkable for two
important representations.

Marton's transport--entropy inequality for the product fair cube gives

```math
W_1^{\rm Ham}(q,U_B)
\le\sqrt{{d\over2}D(q\Vert U_B)}.                     \tag{2c}
```

Thus a fair-orientation selector with average relative entropy `o(N)` is
also redundant in (2b), irrespective of its internal dependence.  In
particular, any law which improves a genuinely linear uniform defect must
pay `Omega(N)` KL.  Linear KL is the sharp unresolved scale: Section 5
shows that an arbitrary correlated law can have linear leverage there.

Let `Q:{+-1}^d -> {+-1}^k` be a rank-`k` parity map and suppose a bridge
law `q` has density `dq/dU_B=f(Q(B))`; equivalently, it is uniform within
each fiber of one fixed rank-`k` parity system, with an arbitrary law on the
syndromes.  Choose systematic coordinates for `Q`.  Couple `B~q` to a
uniform bridge `B'` by keeping the `d-k` free coordinates and replacing the
`k` pivot coordinates by independent fair signs.  Then `B'~U_B` and

```math
\mathbb E d_H(B,B')={k\over2}.
```

Combining this coupling with (2) proves

```math
\boxed{
|\mathbb E_qL_\epsilon-\mathbb E_{U_B}L_\epsilon|
\le tk={\beta k\over\sqrt N}.}                         \tag{3}
```

The same statement holds for a law whose density is constant on the affine
leaves of an adaptive depth-`k` parity decision tree, even when the queried
spans change.  Sample a leaf with its `q`-mass, sample `B'~U_B`, and project
`B'` to that leaf by correcting a systematic pivot set.  The projection is
uniform on the leaf and changes `rank(leaf)/2<=k/2` coordinates on average.
It does not cover additional nonuniform structure inside a leaf.

Consequently, suppose a selector first samples a fair orientation and then
uses a fixed-parity bridge law of average rank `k` in that orientation.  If
it proves

```math
{1\over2}\sum_\epsilon\mathbb E_{q_\epsilon}L_\epsilon-T
\le C N^{1-\delta},
```

then (1)--(3), averaged over the fair orientation, already imply

```math
\boxed{
E_{m,n}(\beta)
\le C N^{1-\delta}+{\beta k\over\sqrt N}.}           \tag{4}
```

In particular, `k=O(N)` makes the selector quantitatively redundant: (4)
is `O(N^(1-min(delta,1/2)))`, and the uniform centered identity itself has
the same power saving.  Conversely, if the joint uniform defect obeys

```math
C^0+\Gamma_A+\Gamma_D-D\ge cN,
```

then any target-reaching fair-orientation pair of hard affine fibers must
have average rank

```math
\boxed{
{k_++k_-\over2}\ge(c/\beta-o(1))N^{3/2}.}             \tag{5}
```

Thus a speed-`N` affine triangle selector cannot cancel a linear centered
defect.  In particular at least one orientation uses a fiber of probability
at most `exp{-(c log(2)/beta-o(1))N^(3/2)}`.  If each orientation must make
the same linear improvement, the bound holds for both fibers.  No bound on
the probability of their union follows from average rank alone.

## 2. The complete triangle state is the full bridge

Fix roots `1` in both children.  The anchored parent-triangle parities

```math
u_{ij}=A_{1i}B_{1j}B_{ij}\quad(2\le i\le m,1\le j\le n),
```

and

```math
v_j=\epsilon D_{1j}B_{11}B_{1j}\quad(2\le j\le n)
```

are `mn-1` independent parity coordinates.  Given `B_11`, they reconstruct

```math
B_{1j}=\epsilon D_{1j}v_jB_{11},
\qquad B_{ij}=A_{1i}u_{ij}B_{1j}.                      \tag{6}
```

Hence a complete triangle table determines `B` up to global negation, and
`L_epsilon(B)=L_epsilon(-B)`.  A prescribed table is a two-point fiber of
probability `2^(-(mn-1))`.  Likewise the rectangle holonomies
`B_ij B_i1 B_1j B_11`, `i,j>1`, have rank `(m-1)(n-1)`.  A complete
triangle or rectangle dynamic state is therefore the full bridge
optimization modulo gauge, not a compressed route to (4).

## 3. Product tilts with linear entropy have only square-root leverage

For each fair orientation, let `q_epsilon` be an independent product law
on the bridge coordinates and write
`mu_(epsilon,e)=E_(q_epsilon) B_e`.  If their average entropy obeys

```math
{1\over2}\sum_\epsilon D(q_\epsilon\Vert U_B)\le CN,
```

then Pinsker's binary-coordinate bound gives, after averaging orientations,

```math
{1\over2}\sum_\epsilon\sum_e\mu_{\epsilon,e}^2
\le2CN.                                                \tag{7}
```

The following comparison holds separately in either orientation and is
uniform in both children (write `mu_e=mu_(epsilon,e)` locally):

```math
\boxed{
\mathbb E_{q_\epsilon}L_\epsilon
\ge\mathbb E_{U_B}L_\epsilon
-t\sum_e\mu_e^2-3dt^3.}                               \tag{8}
```

To prove it, let `V_e` be independent fair signs and set
`Y_e=mu_e+sqrt(1-mu_e^2)V_e`.  The coordinates of `Y` match the first two
moments of `q`.  The third derivative of `L_epsilon` in one bridge
coordinate is `t^3` times the third cumulant of a sign, hence has absolute
value at most `2t^3`.  A coordinatewise Taylor replacement about the common
mean therefore costs at most `3dt^3`.  Next,

```math
z\longmapsto\mathbb E_VL_\epsilon
 (z+(\sqrt{1-\mu_e^2}V_e)_e)
```

is convex and even, so it is minimized at `z=0`.  Finally (2) and
`1-sqrt(1-u^2)<=u^2` compare the variance-reduced fair vector with `V`,
giving (8).

Equations (7)--(8), averaged over the fair orientation and using
`d<=N^2/4`, yield

```math
\boxed{
{1\over2}\sum_\epsilon\mathbb E_{q_\epsilon}L_\epsilon
\ge {1\over2}\sum_\epsilon\mathbb E_{U_B}L_\epsilon
-\left(2C\beta+{3\beta^3\over4}\right)\sqrt N.}       \tag{9}
```

Therefore, if such a fair-orientation product selector proves
`(1/2)sum_epsilon E_(q_epsilon)L_epsilon-T<=C_1N^(1-delta)`, then the
exact uniform bridge identity gives

```math
\boxed{
E_{m,n}(\beta)
\le C_1N^{1-\delta}
+\left(2C\beta+{3\beta^3\over4}\right)\sqrt N.}       \tag{10}
```

An entropy-`O(N)` product/template tilt cannot manufacture a leading
improvement which was absent under uniform bridges.  Strong biases on
`O(N)` coordinates and weak product biases spread across all `mn`
coordinates are both covered.  Mixtures are covered only when the displayed
mean-square bias and componentwise replacement bounds hold after
conditioning; unrestricted mixtures can contain essential dependence and
are not claimed.

## 4. Fractional moments are a basin certificate, not a shortcut

Regard `L=L_epsilon(B)` as a function on the joint fair
orientation/bridge space and, for `lambda>0`, define

```math
\mathcal R_\lambda
=-{1\over\lambda}\log\mathbb E_Ue^{-\lambda L}.
```

The soft minimum lies above `min L`, so the parent construction gives the
direct implication

```math
\boxed{E_{m,n}(\beta)\le\mathcal R_\lambda-T.}        \tag{11}
```

For the tilted law

```math
{dq_\lambda\over dU}=
 {e^{-\lambda L}\over\mathbb E_Ue^{-\lambda L}},
```

the Gibbs variational identity is

```math
\mathcal R_\lambda
=\mathbb E_{q_\lambda}L
{1\over\lambda}D(q_\lambda\Vert U).                  \tag{12}
```

Since `L>=0`, success at error `epsilon`, namely
`R_lambda<=T+epsilon`, necessarily entails

```math
D(q_\lambda\Vert U)
\le\lambda(T+\epsilon)=O_\beta(\lambda N).           \tag{13}
```

It also forces a quantitatively rare near-target basin.  For every
`s>epsilon`, splitting the exponential moment on
`{L<=T+s}` gives

```math
\boxed{
U\{L\le T+s\}
\ge {e^{-\lambda\epsilon}-e^{-\lambda s}
       \over e^{\lambda T}-e^{-\lambda s}}.}          \tag{14}
```

In particular, fixed-`lambda` success requires an `exp{-O_beta(N)}`
near-target basin.  It does not make an exponentially rare bridge
disappear; it certifies the precise exponential rarity that must be paid.

If `q_lambda` is coordinate-product conditionally on each fair
orientation, (9), (12), and (13) give

```math
\boxed{
\mathcal R_\lambda\le T+\epsilon,quad
q_\lambda\text{ product},\quad q_\lambda(\epsilon)=1/2
\Longrightarrow
E_{m,n}(\beta)
\le\epsilon
 +O_\beta((\lambda+1)\sqrt N).}                       \tag{15}
```

Thus every `lambda=o(sqrt(N))` product fractional tilt is quantitatively
redundant: transporting it back to the exact uniform identity already
gives the same sublinear conclusion.  The unrestricted correlated and
orientation-biased fractional tilt is not excluded.

## 5. Linear entropy can have linear leverage in a log-partition model

The entropy hypothesis in Section 3 cannot be upgraded to arbitrary
correlated laws using convexity, evenness, a log-partition representation,
and Hamming Lipschitzness alone.  Partition `r^2` signs into `r` disjoint
groups `I_g` of size `r` and put

```math
f(B)=\sum_{g=1}^r\log\cosh\!\left(
 {\beta\over\sqrt{2r}}\sum_{e\in I_g}B_e\right).      \tag{16}
```

This is a nonnegative convex even log-partition function, and one sign
flip changes it by at most `2 beta/sqrt(2r)`.  Under the uniform law,
the central limit theorem and uniform integrability give

```math
{1\over r}\mathbb E_Uf
\longrightarrow
c_\beta:=\mathbb E\log\cosh(\beta Z/\sqrt2)>0.        \tag{17}
```

Fix `a>0` and condition every group independently on
`|sum_(e in I_g)B_e|<=a sqrt(r)`.  The resulting correlated law `q_a`
has

```math
D(q_a\Vert U)=\Theta_a(r)=\Theta_a(N),                \tag{18}
```

while, for every sufficiently small fixed `a`,

```math
\mathbb E_Uf-\mathbb E_{q_a}f
=(c_{\beta,a}+o(1))r=\Theta_{\beta,a}(N)              \tag{19}
```

with `c_(beta,a)>0`.  Equations (18)--(19) disprove any generic transport
claim of the form

```math
D(q\Vert U)=O(N),\quad \mathbb E_qL\le T+o(N)
\quad\Longrightarrow\quad \mathbb E_UL\le T+o(N)     \tag{20}
```

for this structural class.  Such a claim, if true for the actual parent
landscape, would have implied `E_(m,n)=o(N)` through (1); the example shows
that complete-graph cut constraints, not generic entropy geometry, must do
the work.

## 6. Quantitative boundary

Equations (3), (5), (9), and (15) use the **actual selected pressure-minimizing
children** and all orders.  They force any fair-orientation rare-bridge proof
which improves a linear uniform defect to use genuinely correlated nonlinear
dependence: it must have Hamming transport on the `N^(3/2)` scale or depend
on more than `O(N)` independent parity directions.  They do not exclude
such a law, an orientation-correlated selector, or an orientation bias, and
do not improve the currently proved `Theta_beta(N)` comparable-split
defect.
