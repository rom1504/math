# Independent audit: regular conference tilts and the sublevel-gradient obstruction

**Verdict: PASS.**  Both frozen task-local sources are mathematically valid
with their stated scope.  I found no normalization, constant, conditioning,
convex-extension, conference-pressure, planted-edit, or Hamming-counting
error that requires a repair.

This audit applies only to the following byte-for-byte sources:

```text
extremal_information/drafts/conference_regular_conditioned_all_tilts.md
sha256 65e93b956549af59cbf5e41585e691747cfa67d844b9b29dbe0c084ceb0c886a

extremal_information/drafts/conference_sublevel_gradient_audit.md
sha256 27c452cd7446020cf165a97e15dec616947f594aafda084af6f1d1ee83903859
```

The audit is independent of the sources and does not modify either one.
There is one optional strengthening, not a defect: the regular-conditioned
argument can be recentered after conditioning to prove the pressure-rate
conclusion for every positive `lambda_r=o(r)`, rather than only the stated
`o(sqrt(r))` range.  Section 8 records the short proof and its scope.

## 1. Operator event and conference scaling: PASS

Let

```math
S_{\epsilon,B}=\begin{pmatrix}A_r&B\\B^T&\epsilon A_r\end{pmatrix},
\qquad t={\beta\over\sqrt{2r}},
```

where `A_r^2=(r-1)I`.  Splitting `S` into its block-diagonal and
off-diagonal parts gives exactly

```math
\|tS_{\epsilon,B}\|_{op}
\le {\beta\over\sqrt2}
 \left(\sqrt{1-1/r}+{\|B\|_{op}\over\sqrt r}\right).
```

Thus the standard event `||B||_op<=(2+delta)sqrt(r)` lies in `K` for all
large `r` whenever

```math
{\beta(3+\delta)\over\sqrt2}<\kappa<\frac12.
```

Such a positive `delta` and `kappa` exist precisely in the source's strict
range `beta<sqrt(2)/6`.  The rectangular Rademacher norm tail then gives
`P(K^c)<=2e^{-cr}` uniformly in the two orientations.  No factor of two is
missing from either the conference diagonal cost or the bridge norm.

The set

```math
K_{\epsilon,r}=\{B:\|tS_{\epsilon,B}\|_{op}\le\kappa\}
```

is a convex real set because it is a sublevel of the norm of an affine map.

## 2. Gradient, boundary, and convex extension: PASS

For a real bridge increment `E`, Theorem 1.3 of
`artifacts/high_temperature_frobenius_pressure_stability.md` gives

```math
|df_B[E]|
\le {K_\kappa\over2}
 \left\|t\begin{pmatrix}0&E\\E^T&0\end{pmatrix}\right\|_*
=K_\kappa t\|E\|_*
\le {K_\kappa\beta\over\sqrt2}\|E\|_F.
```

The equality uses

```math
\left\|\begin{pmatrix}0&E\\E^T&0\end{pmatrix}\right\|_*=2\|E\|_*,
```

and the last step uses `||E||_*<=sqrt(r)||E||_F`.  Therefore the Frobenius
gradient constant is dimension-free.

There is no boundary gap.  The archived covariance theorem assumes the
closed bound `||tS||_op<=kappa` with the strict numerical parameter
`kappa<1/2`; it consequently applies at a point satisfying equality in the
definition of `K`.  Equivalently, conditioning the auxiliary orientation
`sigma` gives the two zero-field Ising laws with interactions `+tS` and
`-tS`.  Each has covariance operator norm at most `K_kappa`, so its
`r`-by-`r` cross block has Frobenius norm at most `K_kappa sqrt(r)`.
Their signed mixture has the same bound, and multiplication by `t` yields
the displayed bridge-gradient constant.

Since `f` is differentiable and convex, the supremum of its supporting
affine functions based at points of `K` is convex, has the same Lipschitz
constant, and agrees with `f` on all of `K`, including its boundary.  An
anchor in `K`, together with the Euclidean diameter `2r` of the bridge
cube, gives `|g|=O(r)` there.  The elementary complete-sign cap gives
`0<=f=O_beta(r^(3/2))`.  Multiplying these bounds by `P(K^c)=e^{-Omega(r)}`
indeed proves

```math
|\mathbb Eg-\mathbb Ef|=o(1),
```

which is stronger than the `o(r)` accuracy eventually needed.

## 3. Convex concentration and the speed-`r^2` tail: PASS

Talagrand convex-Lipschitz concentration on the Rademacher cube applies to
the convex extension `g`.  Its dimension-free Frobenius Lipschitz constant
gives, with constants independent of `r`,

```math
\log\mathbb E e^{-s(g-\mathbb Eg)}\le Cs^2
```

and a lower-tail bound `exp(-cu^2)`.  The audited conference-pressure input
gives, separately for both orientations,

```math
\mathbb Eg_{\epsilon,r}=h_\beta r+o(r),
\qquad
h_\beta=2\psi(\beta/\sqrt2)+{\beta^2\over4}.
```

On `K`, `f=g`; substituting a deviation `u=eta r/2` gives
`exp(-c_eta r^2)`.  Division by `P(K)=1-e^{-Omega(r)}` does not change this
speed.  Although RC.13 is written for the joint orientation law, the proof
also gives the same conclusion for either fixed orientation, which is the
version used later in the Hamming-collar theorem.

## 4. Fixed, growing, and tiny tilts in the frozen statement: PASS

For the full extension, Jensen and concentration give

```math
1\le\mathbb E e^{-\lambda(g-\mathbb Eg)}\le e^{C\lambda^2}.
```

Cauchy--Schwarz bounds the omitted complement by

```math
\mathbb E[1_{K^c}e^{-\lambda(g-\mathbb Eg)}]
\le \exp\{-cr/2+2C\lambda^2+O(1)\}.
```

For `lambda_r=o(sqrt(r))` this is `e^{-Omega(r)}`.  If
`lambda_r>=e^{-c_1r}`, choosing `c_1` smaller than that exponent makes both
the omitted-moment error and the conditioning normalization harmless even
after division by `lambda_r`.  The two orientation means have the same
`h_beta r+o(r)` rate; taking their two-term soft minimum remains between
their minimum and average up to the same `O(lambda_r)` centered-moment
cost, so no hidden `log(2)/lambda_r` loss is necessary.

For `0<lambda_r<e^{-c_1r}`, the source correctly changes arguments.  On
`K`, order `N=2r` and `||tS||_op<=kappa` imply

```math
0\le f\le {1\over2}\kappa\|z\|_2^2=\kappa r.
```

Hoeffding's bounded-range lemma under the already conditioned law gives

```math
0\le \mathbb E_Kf+{1\over\lambda_r}
 \log\mathbb E_Ke^{-\lambda_rf}
\le {\lambda_r\kappa^2r^2\over8}=o(r).
```

The conditioned mean is `h_beta r+o(r)` because the discarded set has
probability `e^{-Omega(r)}` and the unconditioned pressure is
`O(r^(3/2))`.  This closes the exponentially small-tilt gap without
dividing an unspecified error by `lambda_r`.

Finally, the regular event is an actual bounded-cap sign class.  With
`N=2r`,

```math
Q(S)\le {N\over2}\|S\|_{op}
\le r{\kappa\sqrt{2r}\over\beta}
={\kappa\over2\beta}N^{3/2}.
```

The coefficient in RC.22--RC.23 is therefore exact.

## 5. Exact differential and entropy identities: PASS

Under the auxiliary Gibbs law on `(sigma,x,y)`, direct differentiation
gives

```math
\nabla_Bf=tM_B,
\qquad
\nabla^2f_{ij,k\ell}
=t^2\operatorname{Cov}(\sigma x_i y_j,\sigma x_k y_\ell).
```

Because every coordinate square is one, tracing the Hessian gives
`t^2(r^2-||M_B||_F^2)`.  Expanding `||M_B||_F^2` using two replicas gives
SG.8, with no sign or normalization loss.  The relative-entropy and bridge
Euler identities SG.9--SG.10 follow directly from the Gibbs density.

For every real `C`, conditioning on `(sigma,x)` and using
`log cosh u<=u^2/2` proves

```math
\log\mathbb E_\nu e^{\sigma x^TCy}
\le {1\over2}\|C^Tx\|_2^2
\le {r\over2}\|C\|_{op}^2.
```

Entropy duality with `C=M_B/r` gives

```math
D(\mu_B\|\nu)
\ge {\|M_B\|_F^2\over r}
 -{\|M_B\|_{op}^2\over2r}
\ge {\|M_B\|_F^2\over2r}.
```

Since `t^2=beta^2/(2r)`, this is exactly
`||nabla f||_F^2<=beta^2D(mu_B||nu)`.  The trivial entropy ceiling is
`(2r+1)log 2`, so this route alone yields only the stated
`O_beta(sqrt(r))` global gradient bound.

## 6. Universal double and planted mesoscopic edit: PASS

For `epsilon=-1` and `B_r^0=A_r+I`, the bridge is entrywise signed and

```math
(S_r^0)^2
=\operatorname{diag}((2r-1)I+2A_r,(2r-1)I+2A_r).
```

Diagonalizing `A_r` makes the check transparent.  For each eigenvalue
`a=+-sqrt(r-1)`, the corresponding `2`-by-`2` block of `S_r^0` has
eigenvalues

```math
+-\sqrt{2r-1+2a}.
```

Hence the empirical law of `tS_r^0` converges to equal atoms at
`+-beta`.  Every fixed matrix power has blocks in the algebra spanned by
`I` and `A_r`, so the archived fixed-power delocalization hypothesis holds;
also `||tS_r^0||_op->beta<1/2`.  The archived pressure theorem therefore
does give

```math
f_-(B_r^0)/r\longrightarrow2\psi(\beta)=\tau_\beta.
```

The elementary norm estimate
`||S_r^0||_op<=2sqrt(r-1)+1` also follows by splitting into the conference
diagonal block and the `A_r+I` bridge block.

Now let `k=floor(c r^(3/4))` and overwrite any `k`-by-`k` block by `+1`.
Writing `D=B^1-B^0`, the conference norm gives

```math
|1_I^TB^0[I,J]1_J|
\le k(\sqrt{r-1}+1)=o(k^2).
```

Since each nonzero entry of `D` is `2`, it follows exactly that

```math
d:=1_I^TD1_J=k^2-o(k^2),
\qquad \|D\|_F^2=2d,
\qquad \sup_{x,y}|x^TDy|=d.
```

The last equality uses nonnegativity of `D` and is attained by all-one
signs.  It gives the upper pressure increment `td`.

For the reverse bound, pinning `x_I=1` and `y_J=sigma 1` makes the new
auxiliary bridge energy identically `td`.  Every pinned state has exactly
`2^(2k)` preimages.  If `z'` is its pinned image, then at most `2k`
coordinates change, and

```math
|H_{S^0}(z')-H_{S^0}(z)|
={1\over2}|(z'-z)^TS^0(z'+z)|
\le4\|S^0\|_{op}\sqrt{kr}.
```

Thus the two losses are exactly those in SG.29:
`4t||S^0||_op sqrt(kr)` and `2k log 2`.  At
`k=Theta(r^(3/4))` they are respectively
`O_beta(r^(7/8))` and `O(r^(3/4))`, both `o(r)`.  Meanwhile

```math
{td\over r}\longrightarrow{\beta c^2\over\sqrt2}=\delta
```

for `c^2=sqrt(2)delta/beta`.  This proves the claimed linear increment and

```math
\|B^1-B^0\|_F=(\sqrt2+o(1))c r^{3/4}.
```

Convexity at `B^1` gives
`nabla f(B^1):D>=f(B^1)-f(B^0)`.  Cauchy--Schwarz then gives exactly the
coefficient `delta/(sqrt(2)c)` in the `Omega(r^(1/4))` gradient lower
bound.  Choosing `tau_beta+delta<a<h_beta` puts both bridges in the stated
strict sublevel.  Any agreeing extension must pay at least the same
two-point Lipschitz quotient, independently of convexity.

## 7. Hamming collar: PASS

One bridge-bit flip changes `x^TBy` by `2`, hence changes `f` by at most
`2t`.  With

```math
s_r=\left\lfloor{\eta r\over4t}\right\rfloor,
```

a low bridge within `s_r` flips of a regular sign bridge produces a regular
center with pressure at most `(h_beta-eta/2)r`.  The orientation-specific
version of the regular-sector theorem bounds the number of centers by

```math
2^{r^2}e^{-c_0r^2}.
```

Here `s_r=Theta(r^(3/2))=o(r^2)`, and the standard binomial estimate gives

```math
\sum_{j\le s_r}{r^2\choose j}
\le (er^2/s_r)^{s_r}
=\exp\{O(r^{3/2}\log r)\}=e^{o(r^2)}.
```

The union bound therefore proves SG.43.  In the concluding sentence,
"probability `e^{-O(r)}`" must be read in its usual lower-mass sense
`P(L)>=e^{-Cr}` along the sequence; then subtracting the collar's
`e^{-Omega(r^2)}` mass leaves a subset of the same speed-`r` order outside
the collar.

## 8. Optional strengthening: the regular tilt range extends to `o(r)`

This is not needed to validate either frozen statement.  It follows from
the same ingredients after centering under the conditioned law instead of
estimating the omitted part of the full centered moment.

For a fixed orientation let

```math
p_\epsilon=\Pr(K_{\epsilon,r}),
\qquad
m_\epsilon=\mathbb E[f_\epsilon\mid K_{\epsilon,r}].
```

The exceptional-set bounds used above show
`m_epsilon=h_beta r+o(r)` and
`|m_epsilon-Eg_epsilon|=o(1)`.  Because `f=g` on `K`, RC.8 gives, for every
`lambda>0`,

```math
\begin{aligned}
1
&\le\mathbb E[e^{-\lambda(f-m_\epsilon)}\mid K_{\epsilon,r}]\\
&\le p_\epsilon^{-1}
 \exp\{\lambda|m_\epsilon-\mathbb Eg_\epsilon|+C\lambda^2\}.
\end{aligned}
```

The lower bound is Jensen.  Thus the fixed-orientation conditioned soft
minimum lies between

```math
m_\epsilon-C\lambda-o(1)-{1\over\lambda}\log(1/p_\epsilon)
\quad\hbox{and}\quad m_\epsilon.
```

For `lambda_r>=e^{-c_1r}`, the last term is `o(r)`; for smaller tilts the
bounded-range argument RC.21 applies.  Averaging the two orientations adds
only the `o(r)` spread between `m_+` and `m_-`: bounding the two-term
exponential average by its minimum on one side and Jensen by its weighted
mean on the other avoids a `log(2)/lambda` loss.  Consequently the same
proof yields

```math
{\mathcal R^K_{\lambda_r,r}\over r}\longrightarrow h_\beta
\qquad\text{for every positive }\lambda_r=o(r).
```

This larger range remains far below the `Theta(r)` disorder tilt at which
an isolated `e^{-Theta(r^2)}` universal-double bridge can compete.  It does
not control the unconditioned or cut-cap-conditioned laws.

## 9. Archive novelty and scope

The scope classifications in the sources are accurate.

1. The archived two-temperature theorem proves an unconditioned exact rate
   only for a small fixed-tilt interval because its irregular complement
   can compete.  RC.1 conditions on the stricter operator-regular event and
   newly proves a speed-`r^2` lower tail and all-fixed-tilt typical rate
   there.
2. The broader cut-cap event from Theorem 37.4 is not replaced by the
   operator event.  The unresolved mass is precisely the difference
   between those classes, or the full operator-irregular sector.
3. SG.1 is a new scalable obstruction to deriving dimension-free response
   from scalar low pressure: genuine sign bridges in one strict sublevel
   have a macroscopic pressure separation at mesoscopic Frobenius distance.
   It does not rule out a weaker `O(r^(1/4))` gradient estimate or a
   non-gradient lower-tail method.
4. SG.2 is a new localization consequence of RC.1.  It proves that any
   speed-`r` favorable basin must be Hamming-deep in the irregular sector;
   it does not prove that such a basin exists.
5. Neither source supplies the full bridge lower LDP, contracted-temperature
   pressure-minimizer compensation, a same-temperature recurrence, or any
   convergence theorem for the original signing problem.

The resulting research boundary is therefore genuine and narrow: regular
conference outputs are controlled throughout every sublinear disorder
tilt, while the only possible finite-tilt or target-reaching phase must be
carried by pressure-lowering bridges that are deeply operator-irregular.
