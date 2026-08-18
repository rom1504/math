# Adversarial audit: full-range row-magnitude no-gain theorem

**Frozen source:**
`extremal_information/drafts/conference_row_magnitude_full_beta_no_gain.md`

**SHA-256:**
`eea295989cf782ac2e28bbd2904a35e318b243f9d92fc7b706e6c9e00e6abbd4`

**Verdict:** **PASS.**  The one-sided theorem is correct with the stated
scope.  I independently checked the row-layer coupling, projected
uniform-subset covariance, matrix-Bernstein scale and exceptional
probabilities, every pressure normalization, the use of real projected
bridges, convex restoration of the rank-one component, and switching
uniformity.  Attempts to turn the population spike or an extreme allowed
layer into a counterexample fail for the precise convex-supporting-line
reason used in the source.

There is one scope clarification worth retaining in any canonical summary:
the result rules out a lower *typical pressure rate under the uniform fibre
law*.  It does not say that the fibre contains no exceptionally favorable
individual bridge or exponentially smaller subfamily.  This does not alter
MF.1--MF.5.

## 1. Gauge reduction and normalization

Let `D_v=diag(v_r)`.  Right multiplication sends the conditioned bridge to
`B D_v`, whose rows satisfy an all-ones magnitude condition.  The spin
change `y -> D_v y` sends the second child to
`D_v A_r D_v`, again a symmetric conference signing.  More strongly, for
an iid bridge `W`,

```math
f_{A_r,\epsilon D_vA_rD_v}(W)
=f_{A_r,\epsilon A_r}(WD_v),
```

and `WD_v` has exactly the law of `W`.  Thus the uniform conference pressure
center is unchanged, not merely asymptotically uniform over an unspecified
class of children.  Frobenius and operator norms are also unchanged.  Both
orientations survive this conjugation.

With `t=beta/sqrt(2r)`, the full interaction corresponding to MF.2 is

```math
X(C)=t\begin{pmatrix}A_r&C\\C^T&\epsilon A_r\end{pmatrix},
```

because one half of `z^T X(C)z` is exactly
`t(H_A(x)+epsilon H_A(y)+x^TCy)`.  Hence there is no missing factor two in
MF.17.  The block triangle inequality gives

```math
\|X(C)\|_{op}
\le {\beta\over\sqrt2}
 \left(\sqrt{1-1/r}+{\|C\|_{op}\over\sqrt r}\right).
```

Consequently a projected bridge edge `2+o(1)` produces the coefficient
`beta(3+o(1))/sqrt(2)`.  A fixed `delta>0` and `kappa<1/2` satisfying MF.18
exist exactly when `beta<sqrt(2)/6`, so the advertised full conference
high-temperature interval is normalized correctly.

## 2. Layer coupling and the exact projected covariance

For fixed plus counts `k,k'`, nesting the smaller uniform plus set inside
the larger preserves both uniform layer marginals.  The symmetric difference
has size `d=|k-k'|`.  In either nesting direction its changed-coordinate set
`T` is a uniform `d`-subset; the sign of the row difference is irrelevant
after taking an outer product.  Thus

```math
D_i^\circ=\pm2\left(1_T-{d\over r}1\right),
\qquad \|D_i^\circ\|_2^2=4d(1-d/r).
```

For `u=1_T-(d/r)1`, exchangeability and `u perpendicular 1` force
`E uu^T=c(I-P)`.  Taking traces gives

```math
c(r-1)=d(1-d/r)={d(r-d)\over r},
```

and therefore

```math
\mathbb E[(D_i^\circ)^TD_i^\circ\mid d]
={4d(r-d)\over r(r-1)}(I-P),
```

exactly MF.14.

The count-tail estimate is also uniform.  Both the ordinary count and the
conditioned count leave the `r^(3/4)` central window with probability at
most `exp(-c sqrt(r))`, with only a factor `p_0^{-1}` for the latter.  A
union bound over `r` rows is absorbed into the same stretched-exponential
form.  This proves MF.12.

## 3. Matrix Bernstein really gives the claimed subcritical cost

Condition on all `d_i` in the good count event and set

```math
X_i=(D_i^\circ)^TD_i^\circ.
```

The `X_i` are independent positive semidefinite matrices and

```math
\|X_i\|\le8r^{3/4},
\qquad
\left\|\sum_i\mathbb E X_i\right\|
\le {4\over r-1}\sum_i d_i\le(8+o(1))r^{3/4}.
```

Moreover `X_i^2=||D_i^circ||_2^2 X_i`, so the matrix variance proxy for
`X_i-E X_i` is `O(r^(3/2))`.  Matrix Bernstein at

```math
t=C r^{3/4}\log r
```

has exponent

```math
{t^2\over O(r^{3/2})+O(r^{3/4})t}=Omega(C\log r).
```

After the dimension factor `r`, increasing `C` makes the failure probability
`O(r^(-10))`.  It follows that

```math
\|D^\circ\|_{op}^2
=\left\|\sum_iX_i\right\|
\le C r^{3/4}\log r,
```

and hence `||D^circ||op=O(r^(3/8)sqrt(log r))=o(sqrt r)`.
This checks the point at which a superficially plausible but false
`r^(3/4)` operator estimate could otherwise have entered.

Projection is a Frobenius contraction, while every changed sign contributes
four to the unprojected squared Frobenius norm.  Since
`E d_i=O_(p_0)(sqrt r)`, Jensen gives

```math
\mathbb E\|D^\circ\|_F
\le(4r\mathbb E d_i)^{1/2}=O_(p_0)(r^{3/4}).
```

Thus both assertions of MF.2 have the needed uniformity.

## 4. Projected real bridges and `L^1` exceptional events

The comparison theorem applies to real symmetric interactions; it does not
require the off-diagonal block itself to remain a sign matrix.  On the event

```math
\|W\|_{op}\le(2+\delta/2)\sqrt r,
\qquad
\|D^\circ\|_{op}\le(\delta/2)\sqrt r,
```

both `X(W^circ)` and `X(B^circ)` lie in the same closed operator ball of
radius `kappa<1/2`.  Convexity of that ball contains their interpolation
segment.  The archived nuclear stability estimate and
`||D^circ||_*<=sqrt(r)||D^circ||_F` give exactly

```math
|f(B^\circ)-f(W^\circ)|
\le {K_\kappa\beta\over\sqrt2}\|D^\circ\|_F.
```

For the iid projection,

```math
\|WP\|_*={\|W1\|_2\over\sqrt r},
\qquad \mathbb E\|W1\|_2\le
(\mathbb E\|W1\|_2^2)^{1/2}=r.
```

The symmetric off-diagonal dilation doubles this nuclear norm and the
stability theorem contributes its factor `1/2`; therefore the good-event
expected difference between `f(W)` and `f(W^circ)` is `O(1)`, as claimed.

The relevant exceptional probability is

```math
O(r^{-10})+e^{-c\sqrt r}+e^{-c_\delta r}.
```

Every projected bridge has entries bounded in absolute value by two, so its
pressure is `O_beta(r^(3/2))`; sign bridges satisfy the same crude bound.
Multiplying by the displayed probability is `o(r)`.  This justifies the
`L^1` statements in MF.19--MF.22, not only convergence in probability.

## 5. Convex restoration of the population direction

For

```math
Y=t\begin{pmatrix}0&BP\\PB^T&0\end{pmatrix},
\qquad t={\beta\over\sqrt{2r}},
```

the rank-one identity and symmetric dilation give

```math
\|BP\|_*={\|B1\|_2\over\sqrt r},
\qquad
\|Y\|_*=2t{\|B1\|_2\over\sqrt r}
={\sqrt2\beta\over r}\|B1\|_2.
```

At the regular base point, the archived covariance estimate therefore
implies

```math
g'(0)\ge-{K_\kappa\over2}\|Y\|_*
=-{K_\kappa\beta\over\sqrt2 r}\|B1\|_2.
```

The pressure is a log-sum-exp of affine functions of the interaction, so
`g` is globally convex even if `X(B)` is far outside the high-temperature
ball.  Its supporting line at zero is consequently valid at `s=1`; no norm
control on the endpoint is being smuggled into MF.25.

Under the conditioned product law,

```math
\mathbb E\|B1\|_2
\le\left(r\,\mathbb E[S^2\mid E_r]\right)^{1/2}
\le {r\over\sqrt{p_0}},
```

so the expected one-sided loss is `O_(beta,kappa,p_0)(1)`.  As a robustness
check, even the deterministic bound `||B1||_2<=r^(3/2)` would make the
one-sided loss only `O(sqrt r)=o(r)`; the theorem does not rely on delicate
tail cancellation in this final step.

On failure of the regular base event,
`(h_beta-f(B)/r)_+<=h_beta` because every cosh pressure is nonnegative.
Combining this with MF.22 and the preceding supporting-line loss proves
MF.4.  Markov's inequality then gives MF.5.

## 6. Counterexample attempts

I tested the three evident adversarial regimes against the proof.

1. **A fixed Gaussian-tail magnitude band.**  This gives a constant-mass
   event with a population covariance eigenvalue strictly larger than one.
   It can make the restored rank-one component thermodynamically active, but
   convexity permits only an upward nonlinear departure after its bounded
   initial slope.  It cannot lower the pressure by `Theta(r)`.
2. **A union containing an extreme layer.**  Rare rows with sum of order
   `r` can make individual endpoint interactions very large.  They do not
   affect the projected operator proof, and the pointwise `O(sqrt r)`
   supporting-line loss still prevents a linear downward shift.  Such rows
   can only threaten a two-sided equality statement, which MF.1 explicitly
   does not assert.
3. **A gauge chosen adversarially with the conference child.**  The exact
   spin change above converts it to right multiplication of an iid bridge,
   preserving its law and all norm estimates.  There is no orientation- or
   gauge-dependent pressure center to exploit.

None produces a counterexample to the one-sided conclusion.  A genuine
escape must use information absent from this theorem's hypotheses, such as
cross-row dependence, a growing family of distinguished directions, or a
row event whose mass tends to zero.

## 7. Required corrections

None.  In a canonical synopsis, say “the uniform row-product fibre law has
no lower leading pressure rate” rather than the potentially stronger-sounding
“the fibre contains no favorable bridge.”
