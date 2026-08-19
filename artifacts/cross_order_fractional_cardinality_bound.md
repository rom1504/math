# Fractional-cardinality recurrence for actual optimizing children

Status: **proved all-order cross-order inequality**.  It strictly improves
the annealed bridge coefficient at sufficiently large fixed temperature and
uses the exact radial payment of the actual child minimizers.  Its
comparable-split exponent remains one.  A sharp certificate-level lower
floor shows why this support-cardinality/Chernoff architecture cannot become
sublinear at large fixed temperature, even with the best possible scalar
radial payment.

Fix `N=m+n`, `L=mn`, `t=beta/sqrt(N)`, and exact minimizers `A,D` defining
`P_m(beta),P_n(beta)` at
`s_m=beta/sqrt(m),s_n=beta/sqrt(n)`.  Put

```math
\Delta_A=\phi_A(s_m)-\phi_A(t),
\qquad
\Delta_D=\phi_D(s_n)-\phi_D(t).                       \tag{1}
```

## 1. Rank-one-word fractional moment

Under the contracted child measures from the exact channel identity, define

```math
Z_\epsilon(B)=2\mathbb E_{\mu_A\mu_D}
 \mathbf1_{\{\tau\sigma=\epsilon\}}
 e^{t\tau x^{\mathsf T}By}
=\sum_{Q\in\mathcal Q_{m,n}}
 c_Q^\epsilon e^{t\langle B,Q\rangle},               \tag{2}
```

where `Q=tau x y^T` and `|Q_(m,n)|=2^(N-1)`.  All coefficients are
positive, and

```math
\sum_Qc_Q^\epsilon=1+\epsilon u_Au_D,
\qquad
u_A={\mathbb E\sinh(tH_A)\over\mathbb E\cosh(tH_A)}. \tag{3}
```

For every `0<q<=1`, subadditivity of `z^q`, Jensen, and bridge independence
give

```math
\begin{aligned}
q\,\mathbb E_B\log Z_\epsilon(B)
&\le\log\mathbb E_B Z_\epsilon(B)^q\\
&\le L\ell(qt)+\log\sum_Q(c_Q^\epsilon)^q.            \tag{4}
\end{aligned}
```

The support-cardinality inequality

```math
\sum_Q(c_Q^\epsilon)^q
\le2^{(N-1)(1-q)}
 \left(\sum_Qc_Q^\epsilon\right)^q                  \tag{5}
```

and averaging the two orientation logarithms yield

```math
{1\over2}\sum_\epsilon\mathbb E_B\log Z_\epsilon(B)
\le {L\ell(qt)\over q}
 +{(N-1)(1-q)\log2\over q}
 +{1\over2}\log(1-u_A^2u_D^2).
```

Dropping the favorable final term and using the exact bridge identity
`G=F_B-Delta_A-Delta_D` proves

```math
\boxed{
E_{m,n}(\beta)\le G_{m,n}
\le {L\ell(qt)\over q}
 +{(N-1)(1-q)\log2\over q}
 -\Delta_A-\Delta_D.}                                \tag{6}
```

This holds for every split, every order, every fixed `beta>0`, and every
selected pair of exact child minimizers.  It is a genuine cross-order
inequality, not a surrogate-channel statement.

Since convexity and `phi_A(0)=0` imply

```math
\phi_A(t)\le\sqrt{m/N}\,\phi_A(s_m),
```

one obtains the fully scalar actual-pressure corollary

```math
\boxed{
\begin{aligned}
E_{m,n}(\beta)\le{}&{L\ell(qt)\over q}
 +{(N-1)(1-q)\log2\over q}\\
&-(1-\sqrt{m/N})P_m(\beta)
 -(1-\sqrt{n/N})P_n(\beta).
\end{aligned}}                                       \tag{7}
```

Every estimate of the right side is therefore already an estimate of the
requested cross-order defect.

## 2. A second exact radial--Renyi form

The coefficients in (2) can be written, for any representative of `Q`, as

```math
c_Q^\epsilon
=2^{1-N}e^{-\phi_A(t)-\phi_D(t)}
 \cosh\{t(H_A(x)+\epsilon H_D(y))\}.                  \tag{8}
```

Using

```math
(\cosh z)^q\le2^{1-q}\cosh(qz)                       \tag{9}
```

and then the arithmetic/geometric inequality on the two orientations gives
the complementary bound

```math
\boxed{
\begin{aligned}
E_{m,n}(\beta)\le{}&{L\ell(qt)\over q}
 +{N(1-q)\log2\over q}\\
&+{\phi_A(qt)+\phi_D(qt)\over q}
 -\phi_A(s_m)-\phi_D(s_n).
\end{aligned}}                                       \tag{10}
```

The entropy coefficient here is `N`, not `N-1`: the two common-orientation
weights in a rank-one-word fibre need not be equal.  Bound (6) is usually
the cleaner statement; (10) records exactly how the Renyi radial data enter.

## 3. Optimized coefficient and a sharp architecture floor

For `h=N-1` in (6), define

```math
B_h(q)={L\ell(qt)\over q}+{h(1-q)\log2\over q}.
```

Its derivative has the sign of

```math
L[qt\tanh(qt)-\ell(qt)]-h\log2.                      \tag{11a}
```

Thus the finite-order optimizer is `q=1` when the expression at `q=1`
is nonpositive, and otherwise it is the unique root of (11a) in `(0,1)`.

At the equal split, the positive part of (6), divided by `N`, has the limit

```math
f_\beta(q)={q\beta^2\over8}+{(1-q)\log2\over q}.      \tag{11}
```

For `beta>=sqrt(8log2)`, its optimizer is

```math
q_*={\sqrt{8\log2}\over\beta},
\qquad
f_\beta(q_*)=\beta\sqrt{\log2\over2}-\log2.          \tag{12}
```

This is strictly smaller than the annealed coefficient `beta^2/8`, so (6)
is a real coefficient improvement for actual children.  It does not improve
the exponent.

If one uses only the convexity payment displayed in (7), the known upper
bound on the child pressures shows that the optimized **certificate** is at
least

```math
[0.442258401851\ldots\,\beta-\log2-o(1)]N             \tag{12a}
```

in the interior-`q` regime.  Thus even the endpoint-pressure version cannot
be sublinear.

There is a certificate-level obstruction which cannot be repaired by a
stronger scalar estimate of `Delta_A,Delta_D`.  Since `phi_A(t)>=0` and the
rigorous all-order construction gives
`P_m(beta)+P_n(beta)<=(beta/2+o(1))N` at an equal split,

```math
\Delta_A+\Delta_D
\le(\beta/2+o(1))N.                                  \tag{13}
```

Consequently even granting the largest radial subtraction consistent with
the complete scalar pressure data, the optimized right side of (6) is at
least

```math
\boxed{
\left[
 \left(\sqrt{\log2\over2}-{1\over2}\right)\beta
 -\log2-o(1)\right]N
=\left[(0.0887050112577\ldots)\beta-
 \log2-o(1)\right]N.}                                \tag{14}
```

It is positive for `beta>7.81407015...`.  Thus the entire
support-cardinality/subgaussian-Chernoff plus scalar-radial architecture has
an unavoidable linear certificate floor in that range.  This does not lower
bound the true defect and does not exclude a score detector using exact
large deviations or a genuinely correlated rare bridge.  It proves that
merely tuning the fractional exponent or sharpening a scalar child pressure
estimate cannot make (6) sublinear.

## Quantitative verdict

Equation (6) is the strongest new unconditional cross-order inequality in
this campaign so far.  It moves a high-temperature coefficient but leaves
the exponent `1`.  Equation (14) forces the next architecture to use more
than latent support size, generic Chernoff control, and a scalar radial
payment; it must exploit the geometry of which rank-one words receive which
child weights, or optimize an exponentially rare correlated bridge directly.
