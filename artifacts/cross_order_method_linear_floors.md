# Linear floors for four cross-order methods

Status: **proved method-class obstructions**.  These results do not show that
the desired cross-order inequality is false.  They show that four broad
ways of proving it have an unavoidable linear defect.  In every statement
the children may be actual optimizers.

Write

```math
H_A(x)=\sum_{i<j}a_{ij}x_ix_j,
\qquad
Q(A)=\max_x|H_A(x)|,
\qquad
M_n=\min_AQ(A),
\qquad
b_n=M_n^{2/3}.
```

## 1. Separately paid child and bridge envelopes

Let

```math
S=\begin{pmatrix}A&C\\C^{\mathsf T}&D\end{pmatrix},
```

where `A,D` are hollow symmetric sign matrices and `C` is a rectangular
sign matrix.  For fixed `x,y`, replacing `y` by `-y` preserves both internal
quadratic forms and reverses the bridge.  Therefore

```math
\boxed{
Q(S)=\max_{x,y}
\left\{|H_A(x)+H_D(y)|+|x^{\mathsf T}Cy|\right\}.}       \tag{1.1}
```

In particular, the certificate which pays the three channels separately is

```math
Q(S)\le Q(A)+Q(D)+\|C\|_{\infty\to1}.                   \tag{1.2}
```

This certificate can never have a sublinear `b`-scale defect on comparable
blocks.  Indeed, for every `n` by `n` sign matrix `C`, averaging over a
uniform sign vector `y` gives

```math
\begin{aligned}
\|C\|_{\infty\to1}
 &=\max_y\sum_i|(Cy)_i|\\
 &\ge n\,\mathbb E|\varepsilon_1+\cdots+\varepsilon_n|\\
 &\ge {n^{3/2}\over\sqrt2},                             \tag{1.3}
\end{aligned}
```

and the central limit theorem sharpens the last constant asymptotically to
`sqrt(2/pi)`.  Let

```math
\mu_n=\min_{C\in\{\pm1\}^{n\times n}}
       \|C\|_{\infty\to1}.
```

Using exact optimizing children in (1.2), the best separately-paid
`b`-scale statement that this method can certify at an equal split is

```math
b_{2n}\le (2M_n+\mu_n)^{2/3}.
```

Consequently its certified additive defect is at least

```math
E^{\rm sep}_{n,n}
=(2M_n+\mu_n)^{2/3}-2M_n^{2/3}.                         \tag{1.4}
```

The right side decreases as `M_n/n^(3/2)` increases.  The rigorous
asymptotic upper bound `M_n<=(1/2+o(1))n^(3/2)` and the CLT refinement
following (1.3) therefore imply

```math
\boxed{
E^{\rm sep}_{n,n}\ge
\left[
 \left(1+\sqrt{2/\pi}\right)^{2/3}-2^{1/3}+o(1)
\right]n
=(0.218646\ldots+o(1))n.}                              \tag{1.5}
```

This is `0.109323...` times the parent order.  Hence cancellation between
the bridge and the two child energies must occur **before** the absolute
values in (1.1); a scalar bridge norm cannot prove `E_N=o(N)`.

## 2. Independent-edge softmax annealing

Put

```math
e_n={M_n\over\sqrt n},
\qquad
p_n(\beta)={1\over\beta}\min_A
 \log\mathbb E_x\cosh\!\left({\beta H_A(x)\over\sqrt n}\right).
```

The exact entropy sandwich is

```math
e_n-{n\log2\over\beta}\le p_n(\beta)\le e_n.           \tag{2.1}
```

For `N=m+n`, `r=m/N`, averaging one relative child orientation and an
independent sign bridge gives

```math
\boxed{
p_N(\beta)\le
 \sqrt r\,p_m(\beta\sqrt r)
 +\sqrt{1-r}\,p_n(\beta\sqrt{1-r})
 +{mn\over\beta}\log\cosh{\beta\over\sqrt N}.}         \tag{2.2}
```

This already uses actual optimizing children for the two soft pressures.
If one now relaxes every exact soft child term using only the upper half of
(2.1), then at the balanced split converting this certificate back to `e`
incurs, per parent vertex, at least

```math
\inf_{\beta>0}\left\{
 {N\over4\beta}\log\cosh{\beta\over\sqrt N}
 +{\log2\over\beta}
\right}
\longrightarrow\sqrt{\log2\over2}=0.588705\ldots .     \tag{2.3}
```

The contraction coefficients in (2.2) can save at most

```math
{1\over2}\left(1-{1\over\sqrt2}\right)N
=0.146447\ldots N                                      \tag{2.4}
```

because `M_n/n^(3/2)<=1/2+o(1)`.  Thus this
**entropy-sandwich-relaxed certificate** leaves the linear floor

```math
\boxed{E^{\rm ann}_{n,n}
\ge(0.442258\ldots+o(1))N.}                            \tag{2.5}
```

This is a floor for that relaxed certificate, not a lower bound on the true
cross-order difference and not a no-go for an argument retaining additional
information in the exact soft child pressures.

The temperature conflict cannot be removed by an additive scalar
counterterm.  At raw edge temperature `t`, define

```math
\Psi_n(t)={1\over t}\min_A\log\mathbb E_x\cosh(tH_A(x)).
```

Independent bridge averaging yields

```math
\Psi_{m+n}(t)\le\Psi_m(t)+\Psi_n(t)
 +{mn\over t}\log\cosh t.                              \tag{2.6}
```

Every scalar `D_n(t)` which cancels the displayed bridge term for every
integer split has the form

```math
D_n(t)=\binom n2{\log\cosh t\over t}+n c(t).             \tag{2.7}
```

Indeed, subtracting the particular quadratic solution leaves a Cauchy
additive function on the positive integers.  Hence

```math
D_{2n}(t)-2D_n(t)=n^2{\log\cosh t\over t}.               \tag{2.8}
```

Making (2.8) `o(n^(3/2))` requires `t sqrt(n)->0`, whereas the entropy
error `n log(2)/t=o(n^(3/2))` requires `t sqrt(n)->infinity`.  No schedule
and no additive scalar recentering closes (2.6).

## 3. Covariance-only switching-orbit interpolation

For a signing `A` and a uniform switch `s in {+-1}^n`, consider the random
process

```math
X_A^s(x)=H_A(s\odot x).
```

Its full mean and covariance kernel are independent of `A`:

```math
\mathbb E_sX_A^s(x)=0,
\qquad
\mathbb E_sX_A^s(x)X_A^s(y)
=\sum_{i<j}x_ix_jy_iy_j
={(x\cdot y)^2-n\over2}.                               \tag{3.1}
```

Nevertheless these data miss a linear amount of fixed-temperature
pressure even inside the project-scale-cap class.  Fix `beta>0`, let `A_n`
minimize

```math
\log\overline Z_A(\beta/\sqrt n),
\qquad
\overline Z_A(t)=\mathbb E_x\cosh(tH_A(x)),
```

and put `k=floor(C n^(3/4))`.  Form `A'_n` by setting all edges of one
principal `k`-set to `+1` and leaving all other edges unchanged.

Uniformly averaging the original edge signs gives

```math
\log\overline Z_{A_n}(\beta/\sqrt n)
\le\binom n2\log\cosh(\beta/\sqrt n)
\le {\beta^2n\over4}.                                  \tag{3.2}
```

The entropy lower bound on one extremal spin then shows

```math
Q(A_n)\le
\left({\log2\over\beta}+{\beta\over4}+o(1)\right)n^{3/2}.
                                                                    \tag{3.3}
```

Editing the block changes every energy by at most `2 binom(k,2)`, so
`Q(A'_n)=O_{\beta,C}(n^(3/2))` as well.  On the other hand, condition the
first `k` spins to be all `+1`.  Averaging the remaining spins makes every
term outside the planted block vanish.  Convexity of `cosh` gives

```math
\overline Z_{A'_n}(\beta/\sqrt n)
\ge2^{-k}\cosh\!\left({\beta\binom k2\over\sqrt n}\right),
```

and therefore

```math
\log\overline Z_{A'_n}(\beta/\sqrt n)
\ge\left({\beta C^2\over2}+o(1)\right)n.               \tag{3.4}
```

Combining (3.2)--(3.4), any `C` with `C^2>beta/2` gives

```math
\boxed{
\log\overline Z_{A'_n}(\beta/\sqrt n)
-\log\overline Z_{A_n}(\beta/\sqrt n)
\ge\left({\beta C^2\over2}-{\beta^2\over4}+o(1)\right)n
=\Omega_{\beta,C}(n).}                                 \tag{3.5}
```

Switching does not change either partition function, while (3.1) says the
two switching-orbit processes have identical complete covariance kernels.
Consequently any common-reference Guerra/replacement argument which tries
to recover endpoint pressure from only this mean/covariance state has
worst-case error `Omega(N)` on the `O(N^(3/2))`-cap class and cannot imply
`E_N=o(N)`.

The scope is important.  Equation (3.5) does **not** rule out an
interpolation that retains the exact endpoint pressures and uses an
optimizer-specific higher-order identity in its derivative.  It rules out
the tempting universal covariance closure, despite the fact that its
entire covariance kernel—not merely a few moments—is retained.

## 4. Balanced switched-twin blow-ups

Fix an actual order-`m` optimizer `A`, an integer `r>=2`, and replace every
vertex `i` by a fibre of `r` switched twins.  Thus, for arbitrary signs
`sigma_(i,a)`, every cross-fibre edge is

```math
d_{(i,a),(j,b)}=a_{ij}\sigma_{i,a}\sigma_{j,b},       \tag{4.1}
```

while all within-fibre signs may be chosen arbitrarily or adaptively.  Let
`D` be any resulting signing of order `N=rm`.  For a parent spin `z`, put

```math
s_i=\sum_{a=1}^r\sigma_{i,a}z_{i,a}\in[-r,r].         \tag{4.2}
```

The complete cross-fibre energy is exactly `H_A(s)`.  Since a multilinear
quadratic form attains its absolute maximum on the box `[-r,r]^m` at a
vertex,

```math
|H_A(s)|\le r^2M_m.                                   \tag{4.3}
```

The fibre interiors contribute at most `m binom(r,2)`.  Conversely, at a
ground state `x` of `A`, choosing
`z_(i,a)=sigma_(i,a)x_i` makes `s_i=rx_i`; the interiors can change the
result by at most the same amount.  Therefore, for **every** choice of
switches and fibre fillings,

```math
\boxed{
|Q(D)-r^2M_m|\le m{r\choose2}.}                       \tag{4.4}
```

For fixed `r`, (4.4) gives the exact `b=M^(2/3)`-scale asymptotic

```math
Q(D)^{2/3}=r^{4/3}b_m+O_r(\sqrt m).                  \tag{4.5}
```

The desired `r`-copy almost-subadditive target is `r b_m+o(m)`.  Using the
rigorous lower constant `M_m>=(c_*-o(1))m^(3/2)`, every balanced
switched-twin certificate consequently has linear defect

```math
\boxed{
Q(D)^{2/3}-r b_m
\ge\left[c_*^{2/3}(r^{1/3}-1)+o(1)\right]N.}          \tag{4.6}
```

For `r=2` the coefficient is `0.1257455010...` per parent vertex.  Random
switches, pseudorandom fibre fills, and adaptive within-fibre optimization
cannot change this conclusion: (4.4) is pointwise.  The underlying
lexicographic blow-up identity was known in the archive; (4.6) records its
previously unstated quantitative cross-order consequence for actual
optimizers.  Any viable repeated-block architecture must alter a positive
density of cross-fibre edges rather than merely switch cloned children.

## Cross-order consequence

The four statements point to the following quantitative boundaries:

```text
separate absolute channels      => E_N >= c N for the certificate;
entropy-sandwich scalar softmax  => E_N >= c N for that relaxed certificate;
universal covariance closure    => endpoint error >= c N.
balanced switched-twin blow-up   => b-defect >= c_r N.
```

Thus an `o(N)` cross-order proof cannot separately pay the absolute
channels, cannot close the endpoint from only a universal second-order
disorder state, and cannot use the entropy-sandwich relaxation of the scalar
soft pressure.  This is a forced architecture change, not a characterization
of the children.
