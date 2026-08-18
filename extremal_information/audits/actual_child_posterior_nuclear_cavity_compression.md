# Posterior nuclear compression of the actual bridge cavities

**Status.**  Rigorous querywise compression theorem and operational scope
audit.  For every rank-one bridge channel, including every actual
contracted-temperature minimizing-child law and every inverse-disorder
escort, the complete deleted-edge cavity matrix is pointwise close to a
nuclear-norm-bounded posterior mean.  Truncating that posterior mean to rank
`R` gives physical squared cavity error

```math
O_\beta(N/R)+O_\beta(1).                          \tag{PN.1}
```

Thus `R=N^alpha` gives a genuine power saving in the correct escort-weighted
norm.  Positive cavity overlap also forces a macroscopic singular-value
witness.

This is **querywise geometric compression**, not an operational posterior
quotient.  The truncation is selected from the full posterior separately at
each bridge word.  A generic net at the accuracy in (PN.1) has logarithmic
size `Theta(RN)` up to a logarithmic upper-bound factor, and posterior
replica sampling stores `Theta(RN)` sign bits.  Therefore `R -> infinity`,
which is necessary for `o(N)` error, has superlinear state cost.  Actual
optimality would have to prove additional singular-value synchronization or
a lower-information rule generating the truncation.  Nothing below supplies
that rule or converts the frames to product gain/target reach, so this result
alone is not a RESET.  It sharpens the coordinate-lifting part of the current
Smallest Missing Lemma only when frame synchronization is bundled with that
separate directional conversion.

## 1. Channel and posterior notation

Let

```math
d=mn,qquad N=m+n,qquad t={\beta\over\sqrt N},qquad
\rho=\tanh t.                                     \tag{PN.2}
```

Let `mu` be any law on rank-one sign matrices

```math
Q=XY^{\mathsf T}\in\{\pm1\}^{m\times n}.         \tag{PN.3}
```

An optional independent global sign can be absorbed into either factor.
The normalized binary-channel likelihood is

```math
P(B)=E_\mu\prod_{e=1}^{d}(1+\rho B_eQ_e).         \tag{PN.4}
```

For a complete bridge word `B`, let `mu_B` be the Bayes posterior and put

```math
M(B)=E_{\mu_B}Q.                                  \tag{PN.5}
```

For an edge `e`, the exact deleted-edge cavity response is

```math
r_e(B_{-e})
={E_\mu[Q_e\prod_{f\ne e}(1+\rho B_fQ_f)]
  \over
  E_\mu[\prod_{f\ne e}(1+\rho B_fQ_f)]}.         \tag{PN.6}
```

We regard `r(B)=(r_e(B_(-e)))_e` as an `m by n` matrix evaluated at the
complete word `B`; its `e`th coordinate simply does not use `B_e`.

All conclusions are pointwise in `B`.  They therefore remain true after
averaging against the fair bridge law, the actual inverse escort

```math
{dq_\lambda\over dU}={P^{-\lambda}\over E_UP^{-\lambda}},     \tag{PN.7}
```

or any other disorder law.  Actual child minimality is not assumed in the
proof; it enters only by certifying that (PN.3)--(PN.4) is the exact child
channel.

## 2. Posterior means lie in a nuclear ball

### Theorem PN.1 (posterior nuclear norm and cavity proximity)

For every complete bridge word `B`,

```math
\boxed{\|M(B)\|_*\le\sqrt{mn}=\sqrt d.}           \tag{PN.8}
```

Moreover, coordinatewise,

```math
\boxed{
M_e(B)={r_e(B_{-e})+\rho B_e
        \over1+\rho B_er_e(B_{-e})},
\qquad |M_e(B)-r_e(B_{-e})|\le2\rho.}             \tag{PN.9}
```

Consequently

```math
\boxed{\|M(B)-r(B)\|_F^2\le4\rho^2d.}            \tag{PN.10}
```

*Proof.*  Every rank-one sign matrix has one nonzero singular value
`||X||_2||Y||_2=sqrt(mn)`.  Convexity of the nuclear norm and (PN.5) give

```math
\|M(B)\|_*
\le E_{\mu_B}\|Q\|_*=\sqrt d.
```

Bayes insertion of the deleted bit gives the first identity in (PN.9).  If
`z=B_er_e`, then

```math
|M_e-r_e|={\rho(1-z^2)\over1+\rho z}.
```

For `z>=0` this is at most `rho`; for `z=-a<=0`,

```math
{\rho(1-a^2)\over1-\rho a}
=\rho(1+a){1-a\over1-\rho a}\le2\rho.
```

Summing the coordinate bound proves (PN.10).  `square`

The use of the **full** posterior rather than `d` incompatible deleted-edge
posteriors is harmless at physical scale: (PN.10) contributes only `O(1)`
after multiplication by `t^2` at comparable splits.

## 3. Pointwise low-rank cavity approximation

### Theorem PN.2 (best-rank-`R` physical error)

Let `M_R(B)` be the rank-`R` singular-value truncation of `M(B)`.  For every
integer `R>=1` and every bridge word,

```math
\boxed{
\|r(B)-M_R(B)\|_F^2
\le {2d\over R+1}+8\rho^2d.}                     \tag{PN.11}
```

A slightly sharper norm form is

```math
\boxed{
\|r(B)-M_R(B)\|_F
\le\sqrt{d/(R+1)}+2\rho\sqrt d.}                 \tag{PN.12}
```

Hence, uniformly over the actual inverse escort and all actual optimizing
children,

```math
\boxed{
t^2E_{q_\lambda}\|r-M_R\|_F^2
\le {2t^2d\over R+1}+8t^2\rho^2d.}               \tag{PN.13}
```

At comparable splits and fixed `beta`, this is

```math
\boxed{
t^2E_{q_\lambda}\|r-M_R\|_F^2
=O_\beta(N/R)+O_\beta(1).}                       \tag{PN.14}
```

In particular, `R=N^alpha` gives error `O(N^(1-alpha))+O(1)` for every
fixed `alpha in(0,1)`.

*Proof.*  Let the singular values of `M` be
`sigma_1>=sigma_2>=...`.  Since

```math
(R+1)\sigma_{R+1}\le\sum_{j\le R+1}\sigma_j\le\|M\|_*,
```

one has

```math
\begin{aligned}
\|M-M_R\|_F^2
&=\sum_{j>R}\sigma_j^2\\
&\le\sigma_{R+1}\sum_{j>R}\sigma_j\\
&\le{\|M\|_*^2\over R+1}\le{d\over R+1}.        \tag{PN.15}
\end{aligned}
```

Combine (PN.15) with (PN.10), first by the triangle inequality to obtain
(PN.12), and then by `(a+b)^2<=2a^2+2b^2` to obtain (PN.11).  Averaging is
legitimate because the estimate is pointwise.  Finally,

```math
t^2d=\Theta_\beta(N),
\qquad t^2\rho^2d=O_\beta(1),                    \tag{PN.16}
```

which proves (PN.14).  `square`

This theorem is in the correct physical and inverse-escort norm; it does
not use the ordinary fair-law `L^2` tail from the separate Walsh-degree
audit.

## 4. Posterior replicas give the same approximation rate

There is a useful equivalent randomized presentation.  Conditional on `B`,
draw `R` independent rank-one posterior replicas

```math
Q^1,\ldots,Q^R\sim\mu_B,
\qquad \overline Q_R={1\over R}\sum_{a=1}^RQ^a.   \tag{PN.17}
```

### Corollary PN.3 (posterior-replica Maurey approximation)

For every `B`,

```math
\boxed{
E[\|\overline Q_R-M\|_F^2\mid B]
={d-\|M\|_F^2\over R}\le {d\over R}.}           \tag{PN.18}
```

Consequently

```math
\boxed{
t^2E_{q_\lambda}E[\|r-\overline Q_R\|_F^2\mid B]
\le {2t^2d\over R}+8t^2\rho^2d
=O_\beta(N/R)+O_\beta(1).}                      \tag{PN.19}
```

*Proof.*  Expand the conditional variance of the empirical mean and use
`||Q^a||_F^2=d`.  The second assertion follows exactly as in Theorem PN.2.
`square`

Thus a growing posterior-replica object really does approximate all
cavities with a power rate.  But drawing those replicas requires sampling
the posterior `mu_B`; it is not supplied by a child overlap law under zero
external field.  Storing the signed rank-one replicas naively costs
`R(N-1)` sign bits (`N-2` projective coordinate bits plus one orientation
bit per replica).  This is already `omega(N)` as soon as `R -> infinity`.

## 5. Positive overlap forces a spectral witness

The raw negative-overlap obstruction from Theorems 37.52 and 37.56 has a
direct geometric consequence.

### Corollary PN.4 (overlap-to-leading-singular-value conversion)

Let `nu` be any disorder law.  If

```math
{1\over d}E_\nu\|r(B)\|_F^2\ge\eta,              \tag{PN.20}
```

then

```math
\boxed{
E_\nu\|M(B)\|_F^2
\ge\left({\eta\over2}-4\rho^2\right)d}          \tag{PN.21}
```

and

```math
\boxed{
E_\nu\sigma_1(M(B))
\ge\left({\eta\over2}-4\rho^2\right)\sqrt d.}   \tag{PN.22}
```

For fixed positive `eta` and physical `rho=O(N^(-1/2))`, this is a
macroscopic rank-one spectral witness on a positive-`nu`-mass set.

*Proof.*  The inequality

```math
\|r\|_F^2\le2\|M\|_F^2+2\|r-M\|_F^2
```

together with (PN.10) proves (PN.21).  Pointwise,

```math
\|M\|_F^2=\sum_j\sigma_j^2
\le\sigma_1\sum_j\sigma_j
\le\sigma_1\sqrt d,                              \tag{PN.23}
```

so averaging proves (PN.22).  Since `sigma_1<=sqrt d`, a lower bound of
order `sqrt d` in expectation implies such a lower bound on a set of fixed
positive mass, after decreasing the constant.  `square`

This converts positive cavity overlap into an explicit low-rank direction,
but not into a favorable product direction or target reach.  The singular
vectors are selected after seeing the full external bridge word and can
vary over a large response image.

## 6. Metric entropy and the operational ceiling

Let

```math
\mathcal B_*^{m,n}(\sqrt d)
=\{M\in\mathbb R^{m\times n}:\|M\|_*\le\sqrt d\}.              \tag{PN.24}
```

Theorem PN.1 places every posterior mean in this ball.  Standard
singular-value truncation and volumetric nets give the following scale.

### Proposition PN.5 (ambient nuclear-ball response entropy)

For `1<=R<=c min(m,n)` and comparable `m,n`, there are absolute constants
`c_i,C_i>0` such that

```math
\boxed{
c_1RN
\le
\log \mathcal N\left(c_2\sqrt{d/R},
             \mathcal B_*^{m,n}(\sqrt d),\|\cdot\|_F\right)
\le C_1RN\log(C_2R).}                            \tag{PN.25}
```

Here the lower bound is an **ambient** nuclear-ball bound.  It is not a
claim that actual-child posterior means fill this ball.

*Proof sketch with the quantitative ingredients.*

For the upper bound, truncate every matrix to rank `R`; (PN.15) loses at
most `sqrt(d/(R+1))`.  A rank-`R` matrix of Frobenius norm at most `sqrt d`
is specified by two Stiefel frames and `R` singular values.  Euclidean nets
for these three factors at relative precision `Theta(R^(-1/2))` have total
logarithmic size

```math
O(Rm\log(CR)+Rn\log(CR)+R\log(CR)),               \tag{PN.26}
```

and perturb the product by `O(sqrt(d/R))` in Frobenius norm.

For the lower bound, take rank-`R` matrices with `R` equal singular values
`sqrt d/R`.  They have nuclear norm `sqrt d` and Frobenius norm
`sqrt(d/R)`.  A constant-angle packing of the rank-`R` Grassmann manifold
(or the corresponding elementary Stiefel volumetric packing) contains
`exp{Omega(R(m-R))}` left frames whose matrices are separated by a constant
multiple of `sqrt(d/R)`.  Varying the right frame, or simply using the
larger comparable side, gives `exp{Omega(RN)}` separated matrices.  This
proves (PN.25).  `square`

At the physical target, rank `R` gives error `N/R`.  Therefore:

```math
\text{fixed relative accuracy: }R=O(1),\quad \log\mathcal N=\Theta(N),
```

whereas

```math
\text{vanishing relative error: }R\to\infty,
\quad \log\mathcal N=\omega(N)                  \tag{PN.27}
```

for the unrestricted carrier.  This matches the raw posterior-replica cost
up to logarithmic factors.

Equation (PN.25) does **not** prove an actual-child lower bound.  The only
way actual optimality can improve the conclusion is by showing that the
realized posterior means occupy a much smaller subset: for example, common
singular frames, a finite phase orbit, a summable singular tail stronger
than the nuclear bound, or a lower-information rule which generates the
frames from child data.

## 7. Why this is not the desired quotient

The results separate three statements which should not be conflated.

1. **Existence:** a rank-`R` approximation with physical error `O(N/R)` is
   unconditional and pointwise.  This is the new positive theorem.
2. **Generation:** `M_R(B)` is obtained by first computing the full Bayes
   posterior mean at `B`.  Posterior sampling has the same dependence.  No
   zero-field child statistic currently generates either object.
3. **Reuse:** even if every `M_R(B)` is individually small, its singular
   frames can range over a carrier with response entropy `Theta(RN)` up to
   logs.  A reusable quotient must control this response image, not only the
   rank of each member.

Actual pressure minimality gives exact zero-field flip inequalities, and
macroscopic conditional spread rules out hard common coordinate blocks.
Neither statement presently controls the external-field singular frames of
`M(B)`.  Generic perturbation/Ghirlanda--Guerra regularization controls
replica Gram data, while the singular vectors in (PN.22) retain precisely
the labelled coordinate embedding that Gram data forget.

Accordingly, PN.2 is a strict compression of one **already computed** cavity
matrix, but not a strict compression of the map

```math
B\longmapsto r(B).                                \tag{PN.28}
```

It neither controls the reverse-product direction nor supplies the
target-reaching recurrence.  The sharpened operational question is:

> Do actual optimizing children force the posterior singular frames in
> (PN.22) into an `exp{o(N)}` response image, or make their discarded
> singular mass summable by a child-only rule?

Without such a theorem, adaptive posterior nuclear compression is not a
RESET.
