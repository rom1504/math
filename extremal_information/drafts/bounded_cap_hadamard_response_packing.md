# A cap-`1/2` Hadamard family with `Theta(sqrt(n))` response information

Status: rigorous theorem, with a self-contained construction and exact finite
verifier.  This is a response-complexity theorem on the subsequence
`n=2^(2m)`.  It is **not** a convergence theorem for the original signing
problem.

## 1. Result

For a hollow symmetric signing `A` put

```math
H_A(x)=\frac12x^TAx=\sum_{i<j}a_{ij}x_ix_j,
\qquad
Q(A)=\max_x|H_A(x)|.                                      \tag{BH.1}
```

For a fixed sign bridge `B` define

```math
(P_BH_A)(y)=\max_x\{H_A(x)+x^TBy\}.                       \tag{BH.2}
```

The projective response metric, which permits one child-dependent additive
calibration, is

```math
d_{\rm proj}(A,A')
=\frac12\operatorname{osc}_y
  \bigl(P_BH_A(y)-P_BH_{A'}(y)\bigr).                     \tag{BH.3}
```

### Theorem BH.1 (bounded-cap Hadamard response packing)

There is an absolute constant `g_0=1/8` such that, for every sufficiently
large

```math
n=2^{2m},
```

there are a common bridge `B_n in {+-1}^{n times n}` with
`||B_n||_(2->2)=sqrt(n)` and at least

```math
N_n\ge \exp(c\sqrt n)                                     \tag{BH.5}
```

pairwise distinct hollow symmetric signings `A_1,...,A_N` satisfying

```math
Q(A_c)=\frac12n^{3/2}                                     \tag{BH.6}
```

and

```math
d_{\rm proj}(A_c,A_d)\ge g_0 n^{3/2}
\qquad(c\ne d).                                           \tag{BH.7}
```

The ordinary sup-response distance obeys the same lower bound.  Therefore,
for every fixed `epsilon<g_0/2`, any summary which answers all future
continuations through this fixed bridge to error `epsilon n^(3/2)` needs

```math
\exp(\Omega(\sqrt n))\text{ states},
\qquad\text{or }\Omega(\sqrt n)\text{ bits}.              \tag{BH.8}
```

The construction itself is indexed by `sqrt(n)` Boolean bits, so this
particular family has response description complexity
`Theta(sqrt(n))` bits at the declared scale.  The lower bound is for the
entire cap-`1/2` class as well, whereas the matching upper description is
only asserted for this explicit family.

## 2. Walsh and Maiorana--McFarland coordinates

Put `q=2^m` and `n=q^2`.  Index the coordinates by
`(u,v) in F_2^m times F_2^m`.  Let

```math
W_{(a,b),(u,v)}=(-1)^{a\cdot u+b\cdot v}.                  \tag{BH.9}
```

Then `W=W^T`, `W^2=nI`, and `||W||_(2->2)=q=sqrt(n)`.
For every Boolean function `g:F_2^m->F_2`, define the sign vector

```math
s_g(u,v)=(-1)^{u\cdot v+g(v)}.                             \tag{BH.10}
```

Summing first over `u` gives the exact transform

```math
(Ws_g)(a,b)=q(-1)^{a\cdot b+g(a)}.                         \tag{BH.11}
```

Thus `s_g` is bent, and

```math
y_g:=q^{-1}Ws_g\in\{-1,1\}^n,
\qquad
Wy_g=q s_g.                                                \tag{BH.12}
```

No external bent-function theorem is needed: (BH.11) proves every property
used below.

Take `b=s_0` and regularize the Walsh matrix by

```math
\mathcal H=D_bWD_b.                                       \tag{BH.13}
```

Equation (BH.11) with `g=0` says `Wb=qb`, hence

```math
\mathcal H\mathbf1=q\mathbf1,
\qquad
\mathcal H^2=nI.                                          \tag{BH.14}
```

Moreover

```math
\operatorname{tr}\mathcal H=\operatorname{tr}W
=\bigl(\operatorname{tr}
\begin{pmatrix}1&1\\1&-1\end{pmatrix}\bigr)^{2m}=0.      \tag{BH.15}
```

Remove the diagonal:

```math
A=\mathcal H-\operatorname{diag}(\mathcal H).             \tag{BH.16}
```

This is a hollow symmetric `+-1` signing.  Since `x_i^2=1`, (BH.15) gives

```math
H_A(x)=\frac12x^T\mathcal Hx.                             \tag{BH.17}
```

The spectral bound and the all-one eigenvector in (BH.14) now prove

```math
Q(A)=\frac12qn=\frac12n^{3/2}.                            \tag{BH.18}
```

For every `g`, switch the child by

```math
A_g=D_{s_g}AD_{s_g}.                                      \tag{BH.19}
```

Switching preserves hollowness, sign coefficients, and the exact cap
(BH.18).

## 3. An exponentially large low-bias code

Write

```math
S(g,h)=\sum_{v\in F_2^m}(-1)^{g(v)+h(v)}.                 \tag{BH.20}
```

There is a collection `G` of at least `exp(cq)` Boolean functions such that

```math
|S(g,h)|\le q/2\qquad(g\ne h).                            \tag{BH.21}
```

Indeed, sample

```math
N=\lfloor\exp(q/32)\rfloor
```

independent uniform sign tables.  For each pair, Hoeffding gives

```math
\Pr\{|S(g,h)|>q/2\}\le2e^{-q/8}.                         \tag{BH.22}
```

The union bound over fewer than `N^2` pairs is less than one for all
sufficiently large `q`.  The absolute-value condition also excludes equal
or complementary tables.  In particular, the switchings in (BH.19) are
pairwise distinct: if `D_sAD_s=D_tAD_t`, then the nonzero off-diagonal
entries force `s_i t_i` to be constant in `i`.

## 4. The exact pair Rayleigh coordinate

Fix distinct `g,h in G` and put

```math
w=s_g\odot s_h=(-1)^{g(v)+h(v)}.                           \tag{BH.23}
```

The vector `b odot w` is `s_{g+h}`.  Taking the inner product of (BH.11)
with `s_{g+h}` gives

```math
\begin{aligned}
w^T\mathcal Hw
&=(b\odot w)^TW(b\odot w)\\
&=q\left(\sum_v(-1)^{g(v)+h(v)}\right)^2
=qS(g,h)^2.                                                \tag{BH.24}
\end{aligned}
```

Consequently

```math
0\le\rho:=\frac{w^T\mathcal Hw}{qn}
=\frac{S(g,h)^2}{q^2}\le\frac14.                         \tag{BH.25}
```

This rooted pair coordinate is the only part of the `sqrt(n)`-bit table
used by the response upper bound.

## 5. Resolvent trust-region audit

Let

```math
P_\pm=\frac12(I\pm\mathcal H/q).                          \tag{BH.26}
```

These are the orthogonal projections onto the `+-q` eigenspaces.  Equation
(BH.25) is equivalently

```math
||P_+w||_2^2=\frac n2(1+\rho),
\qquad
||P_-w||_2^2=\frac n2(1-\rho).                            \tag{BH.27}
```

There is a particularly short spherical trust-region bound.  For every
`u` with `||u||_2^2=n`, put `K=2qI-\mathcal H`.  Since the eigenvalues of
`mathcal H` are `+-q`, the matrix `K` is positive definite.  Completing the
square gives

```math
\begin{aligned}
\frac12u^T\mathcal Hu+(qw)^Tu
&=qn-\frac12u^TKu+(qw)^Tu\\
&\le qn+\frac12(qw)^TK^{-1}(qw).                          \tag{BH.28}
\end{aligned}
```

The involution identity `mathcal H^2=q^2I` makes the inverse explicit:

```math
K^{-1}=\frac{2qI+\mathcal H}{3q^2}.                       \tag{BH.29}
```

Consequently

```math
\begin{aligned}
\frac12u^T\mathcal Hu+q w^Tu
&\le qn+\frac{2qn+w^T\mathcal Hw}{6}\\
&=nq\left(1+\frac{2+\rho}{6}\right)
\le\frac{11}{8}nq.                                       \tag{BH.30}
\end{aligned}
```

The completion uses only the sphere containing the Boolean cube.  Thus
(BH.30) proves a strict `1/8` gap from `3/2`; no Boolean relaxation or
unproved optimizer assertion remains.

## 6. Response separation

Use the common bridge

```math
B=W.                                                       \tag{BH.32}
```

At query `y_h`, substitute `u=D_{s_g}x`.  Equations
(BH.12), (BH.17), and (BH.19) give the exact identity

```math
(P_WH_{A_g})(y_h)
=\max_{u\in\{-1,1\}^n}
 \left\{\frac12u^T\mathcal Hu+q(s_g\odot s_h)^Tu\right\}.
                                                                    \tag{BH.33}
```

On the diagonal `g=h`, the vector `u=1` attains the spherical upper bound:

```math
(P_WH_{A_g})(y_g)=\frac32nq=\frac32n^{3/2}.               \tag{BH.34}
```

Off the diagonal, (BH.30) gives

```math
(P_WH_{A_g})(y_h)
\le\frac{11}{8}n^{3/2}.
                                                                    \tag{BH.35}
```

Thus, at `y_g`, the response of child `g` exceeds that of child `h` by at
least `g_0n^(3/2)`.  At `y_h`, the reverse difference has the same magnitude.
The oscillation of the response difference is at least `2g_0n^(3/2)`, which
proves both the projective and sup-norm claims in (BH.7).

The usual coordinate-pinning continuation proves that sup-norm response
distance is exactly all-future contextual distance.  The projective packing
then gives (BH.8) even if every decoded response may be recalibrated by its
own additive constant.

## 7. Why generic random fields do not prove this theorem

The construction was motivated by regular Hadamard/conference trust-region
bounds, but the ordinary random-field attempt has a sharp leading-constant
ceiling.  Suppose a field has asymptotic Euclidean norm `n`, its mass is
equally split between the two eigenspaces, and its planted-pole gain is the
generic Gaussian absolute-mean value

```math
\mu=\sqrt{2/\pi}.                                         \tag{BH.36}
```

The off-pole spherical extra gain is

```math
\max_{0\le r\le1}
\left\{-r^2+\frac{\sqrt{1-r^2}+r}{\sqrt2}\right\}
=\frac{3\sqrt3-2}{4}
=0.7990381056\ldots,                                      \tag{BH.37}
```

attained at `r=sin(pi/12)`.  But

```math
\sqrt{2/\pi}=0.7978845608\ldots.                          \tag{BH.38}
```

So a proof using only typical `l_1` gain, typical `l_2` norm, and balanced
spectral projection misses separation by

```math
0.0011535448\ldots\ n^{3/2}.                              \tag{BH.39}
```

This is a ceiling for that coarse trust-region certificate, not a theorem
that random bridges cannot pack bounded-cap children.  Rare fields with
larger `l_1/l_2` ratio could still evade it.  The exact flat fields in
(BH.12) do evade it; their pole gain is `1` instead of `sqrt(2/pi)`.

## 8. Exact verification

Run

```bash
./.venv/bin/python \
  extremal_information/experiments/verify_bounded_cap_hadamard_response.py
```

The verifier uses integer arithmetic for all Walsh, cap, transform,
Rayleigh, and Boolean-response identities at `n=4,16`, checks a larger set
of transform/Rayleigh identities at `n=64`, and separately audits the two
constants in (BH.37)--(BH.39).

## 9. Interpretation and scope

This theorem closes one previously open qualitative gap.

1. The planted-pole packing for arbitrary sign quadratics relied on children
   of quadratic internal scale.  Here every child has the exact natural cap
   `n^(3/2)/2`.
2. Bounded cap therefore does **not** imply bounded, logarithmic, or
   `o(sqrt(n))` all-future response information.
3. On this structured family, a `sqrt(n)`-bit rooted table is both sufficient
   and necessary up to constants.  This is a genuine sub-landscape quotient:
   it is far smaller than the `Theta(n^2)` signing coefficients and the
   `2^n`-entry response landscape.
4. The theorem is only on the subsequence `n=2^(2m)`, only gives
   `Omega(sqrt(n))` bits, and says nothing about convergence of
   `M_n/n^(3/2)`.  It also does not prove that exact or near minimizers below
   the `1/2` cap have the same response complexity.

The sharp next discriminator is whether the `sqrt(n)` lower bound can be
raised to `Omega(n)` within an asymptotic cap-`1/2` class, or whether every
such class admits an `O(sqrt(n))`- or `o(n)`-bit quotient after imposing an
appropriate synchronization invariant.  The generic-field ceiling says
that either direction needs information beyond one planted pole and one
two-eigenspace projection statistic.
