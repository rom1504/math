# Adaptive product clusters: field dispersion, a variance counterexample, and spectral counting

Date: 2026-09-06. Three finite lemmas are proved below. The first gives a
strictly stronger partition lower bound when the aligned local fields have
nonvanishing average absolute deviation. The second shows that a variance
bound alone cannot force such an improvement. The third turns the compressed
weave's small middle eigenspace into a counting restriction, not an exclusion
of a nearly constant-field maximizer. No new numerical signing bound is
claimed in this note.

## 1. Exact adaptive product law

Let `A` be a real symmetric hollow matrix of order `N`, with
`||A||op<=C sqrt N`. Write `H_A(x)=x^TAx/2` and `Q=max_x|H_A(x)|`.
Choose an objective sign `sigma` and a spin `x` with `sigma H_A(x)=Q`.
Gauge-transform to

```math
B=\sigma\operatorname{diag}(x)A\operatorname{diag}(x),\qquad
\ell_i=\frac{(B1)_i}{\sqrt N},\qquad
\bar\ell=\frac1N\sum_i\ell_i=\frac{2Q}{N^{3/2}}.
```

Single-coordinate maximality gives `ell_i>=0`, and
`N^{-1}sum ell_i^2<=C^2`. For arbitrary independent flip probabilities
`delta_i in [0,1]`, let `Y_i=x_i xi_i` with `P(xi_i=-1)=delta_i`.
Then, exactly,

```math
\mathcal H(Y)=\sum_i h(\delta_i),\qquad
\mathbb E\sigma H_A(Y)
=Q-2\sqrt N\sum_i\ell_i\delta_i+2\delta^TB\delta. \tag{1}
```

Thus Jensen's finite Gibbs inequality at `b=beta/sqrt N` gives

```math
\log Z_\sigma(b)\ge bQ+
\sum_i\{h(\delta_i)-2\beta\ell_i\delta_i\}
+\frac{2\beta}{\sqrt N}\delta^TB\delta.        \tag{2}
```

Replacing the last term by `-2 beta C sum delta_i^2` gives the separable
lower bound `bQ+sum_i max_{d in [0,1/2]}[h(d)-2 beta ell_i d-2 beta C d^2]`.
The following centered version retains the exact homogeneous baseline.

## 2. A quantitative extra credit from average absolute field dispersion

Fix `0<delta<=1/3` and `q=1-2delta`. Define

```math
d_1=\frac1N\sum_i|\ell_i-\bar\ell|,\qquad
K=\frac1{\delta(1-\delta/2)}+2\beta C.
```

Then

```math
\log[Z_+(b)+Z_-(b)]
\ge Nh(\delta)+bq^2Q+N\mathcal G(d_1),       \tag{3}
```

where the explicit nonnegative credit is

```math
\mathcal G(d)=
\min\left\{\frac{(\beta qd)^2}{K},
                 \frac{\beta q\delta d}{4}\right\}.       \tag{4}
```

In particular, a lower bound `d_1>=d_*>0`, at the relevant extremizers,
would improve the uniform-cluster extraction by a fixed amount at fixed
`beta,C,delta`. This is a sufficient, directly checkable field condition;
it is not proved for the selected weave signings.

A truncated-variance condition is enough: for any fixed `R>0`,
`d_1>=N^{-1}sum_i (ell_i-bar ell)^2 1_{|ell_i-bar ell|<=R}/R`.
Thus positive variance together with uniform integrability of the squared
fields implies a positive `d_1`; an operator-norm bound alone gives no such
uniform integrability.

Proof. Put `s_i=sign(ell_i-bar ell)`, with `sign(0)=0`, and
`u_i=-epsilon(s_i-bar s)`. Then

```math
\sum_i u_i=0,\quad |u_i|\le2\epsilon,\quad
\sum_i u_i^2\le N\epsilon^2,\quad
\sum_i\ell_i u_i=-N\epsilon d_1.             \tag{5}
```

For `0<=epsilon<=delta/4`, the probabilities `delta_i=delta+u_i`
lie in `[delta/2,3delta/2] subset (0,1/2]`. On this interval
`|h''(z)|=1/[z(1-z)]<=2/[delta(1-delta/2)]`. Taylor's theorem and the
zero sum in (5) therefore give

```math
\sum_i h(\delta+u_i)\ge
Nh(\delta)-\frac{\sum_i u_i^2}{\delta(1-\delta/2)}.
```

Expanding (1) around the uniform probability gives

```math
\mathbb E\sigma H_A(Y)
=q^2Q-2q\sqrt N\sum_i\ell_i u_i+2u^TBu.
```

Use `u^TBu>=-C sqrt N sum_i u_i^2` and (5). The extra term in the
Gibbs lower bound is at least
`N[2 beta q epsilon d_1-K epsilon^2]`. Choose
`epsilon=min(delta/4,beta q d_1/K)`. If the second argument is the
minimum, the gain is `(beta q d_1)^2/K`. Otherwise it is at least
`beta q delta d_1/4`. This proves (3)--(4). The proof uses only one of
the two objective signs, so no symmetry assumption or missing factor two
is hidden in the two-sided partition.

## 3. Positive field variance does not force a product-pressure improvement

There is a sequence of hollow sign matrices of bounded normalized operator
norm, with a global absolute maximizer whose aligned field variance tends
to one, while its best adaptive-product pressure has the same limiting
value as its best uniform-product pressure.

Let `H_4=J_4-2I_4`. For odd positive integers `r`, put
`n=4^r` and `H=H_4^{otimes r}`. Then `H` is symmetric, has diagonal
`-1`, satisfies `H^2=nI`, and has row sum `sqrt n`. Its hollowing
`A_0=H+I` has largest absolute eigenvalue `sqrt n+1`, attained at `1`.
Adjoin a single universal positive hub:

```math
A=\begin{pmatrix}0&1^T\\1&A_0\end{pmatrix},\qquad N=n+1.
```

For every spin, the core energy has absolute value at most
`n(sqrt n+1)/2`, and the hub term at most `n`. Equality in the sum
is attained at the all-positive spin. Hence

```math
Q(A)=\frac{n(\sqrt n+3)}2,\qquad
\|A\|_{op}\le2\sqrt n+1=O(\sqrt N).          \tag{6}
```

At that maximizer the hub field is `n/sqrt N`, while each core field
is `(sqrt n+2)/sqrt N`. Consequently

```math
\bar\ell\longrightarrow1,\qquad
\frac1N\sum_i\ell_i^2\longrightarrow2,\qquad
\operatorname{Var}(\ell)\longrightarrow1,
\qquad d_1\longrightarrow0.                \tag{7}
```

This is not merely a failure of the particular estimate (4). Define the
best product variational value at fixed `beta>0` by

```math
P_N(\beta)=\sup_{m\in[-1,1]^N}
\left\{\sum_i h\!\left(\frac{1-m_i}{2}\right)
       +\frac{\beta}{2\sqrt N}m^TAm\right\}.
```

The function

```math
\varphi(s)=h\!\left(\frac{1-\sqrt s}{2}\right)
=\log2-\sum_{j\ge1}\frac{s^j}{2j(2j-1)},\quad 0\le s\le1,
```

is concave. For the core, the spectral inequality and Jensen give

```math
\sum_{i=1}^n\varphi(m_i^2)+
\frac{\beta}{2\sqrt N}m^TA_0m
\le n\sup_{s\in[0,1]}
\left\{\varphi(s)+\frac{\beta(\sqrt n+1)}{2\sqrt N}s\right\}.
```

Equality is attained by a constant vector of means. The hub adds at most
`log2+beta n/sqrt N=o(N)` to this upper bound; giving the hub zero mean
provides a matching lower bound up to `o(N)`. Thus

```math
\lim_{N\to\infty}\frac{P_N(\beta)}N
=\sup_{0\le\delta\le1/2}
\left\{h(\delta)+\frac\beta2(1-2\delta)^2\right\}. \tag{8}
```

The right side is exactly the best homogeneous cluster value using (6).
Therefore positive variance, even at a genuine absolute maximizer of a
sign matrix with `||A||op=O(sqrt N)`, does not force additional adaptive
product credit. A variance-to-dispersion inference requires a tail condition
or additional weave-specific structure. This example is not asserted to
be a compressed weave.

## 4. What the compressed spectral identity does give

Let `K` be a principal `N`-coordinate restriction of a symmetric matrix
`W` with `W^2=m^2I_M`. Block multiplication gives
`K^2+LL^T=m^2I_N`. Thus `B=K/sqrt N` has at most `d=M-N` eigenvalues
strictly between `-s` and `s`, where `s=m/sqrt N`. Let `E` be that middle
eigenspace. Fix `|mu|<s` and `g=s-|mu|>0`.

If `||(B-mu I)x||_2<=epsilon sqrt N` for a Boolean vector `x`, then

```math
\operatorname{dist}(x,E)\le\frac\epsilon g\sqrt N. \tag{9}
```

This follows by diagonalizing `B` on `E` and its orthogonal complement.
For exact equality `(B-mu I)x=0`, all such Boolean vectors lie in `E`.
A `d`-dimensional subspace contains at most `2^d` Boolean vectors: choose
`d` coordinate functionals that give an injective map on the subspace.
If `mu` varies, there are at most `N` eigenvalues, hence at most `N2^d`
exact constant-field vectors with interior eigenvalue. This still permits
one or exponentially many candidate maximizers.

There is also an explicit entropy bound from an `L1` residual. Suppose

```math
\frac1N\|(B-\mu I)x\|_1\le\eta.
```

Choose `tau>0`, `theta=eta/tau<=1/2`, and a net radius `alpha>0`.
Put `e=sqrt(tau eta)/g`, and assume `(e+alpha)^2<=1/2`. Then the number
of these Boolean vectors is at most

```math
\exp N\left[
h(\theta)+\left(\frac dN+\theta\right)\log(1+2/\alpha)
             +h\big((e+\alpha)^2\big)\right].             \tag{10}
```

Proof. Write `w=(B-mu I)x` and let `J={i:|w_i|>tau}`. Then
`|J|<=theta N` and `||w_{J^c}||_2^2<=tau eta N`. On `E^perp`,
the inverse `R=(B-mu I)^{-1}` has norm at most `1/g`. Hence `x`
is within `e sqrt N` of the subspace

```math
U_J=E+R P_{E^\perp}\operatorname{span}\{e_j:j\in J\},
\qquad\dim U_J\le d+|J|.
```

The radius-`sqrt N` ball in `U_J` has an `alpha sqrt N` net of size at
most `(1+2/alpha)^{d+|J|}`: disjoint balls of half the net radius fit
in the ball with that half-radius added. Use the orthogonal projection of
`x` onto `U_J`, whose norm is at most `sqrt N`. Round a nearby net center
coordinatewise to a Boolean vector. Each mismatched coordinate costs at
least one squared unit of distance, so the Hamming distance is at most
`(e+alpha)^2 N`. The corresponding Hamming ball has at most
`exp[N h((e+alpha)^2)]` points. Finally sum over the at most
`exp[N h(theta)]` possible sets `J`. These binomial entropy estimates
follow directly from the binomial generating function; floors only
improve them.

For a bounded range of `mu`, use a fixed grid: replacing `mu` by a grid
point within `zeta` increases the normalized `L1` residual by at most
`zeta`, and costs only the number of grid points. The grid is fixed before
the order limit. Thus (10) applies also when `mu` is the spin's own mean
field. Hollowing `K` changes every normalized field by at most `1/sqrt N`
and its centered average absolute deviation by at most `2/sqrt N`.

For the actual retention `p=31/32`, one has `d/N<=1/31+o(1)` and
`s=1/sqrt p+o(1)`. Whenever the mean field stays a fixed distance inside
`(-s,s)`, (10) gives an entropy strictly below `log2` for sufficiently
small fixed `eta`. For instance, with `alpha=1/8` its limiting exponent
as `eta` decreases to zero is at most `log(17)/31+h(1/64)<log2`.

## 5. What is still missing

The spectral facts alone cannot exclude a constant-field Boolean vector.
Already the sign Hadamard `W=J_4-2I_4`, restricted to three coordinates,
has `K=J_3-2I_3`, `K^2=4I_3-11^T`, and `K1=1`. Its middle space is
one-dimensional and contains the all-positive spin; the hollowing is the
positive triangle and that spin is an absolute maximizer. This finite
example is not an asymptotic counterexample at retention `31/32`, but
it directly defeats the inference that small middle-space dimension
alone forces local-field spread at every maximizer.

To turn (3) into a new all-order coefficient requires either a uniform
lower bound on `d_1` at the relevant extremizers, or a separate probabilistic
argument excluding the small-field-dispersion candidates counted in (10).
In the actual weave, `beta=2t/sqrt p+o(1)` and
`C<=1/sqrt p+o(1)`. If such a fixed `d_*>0` were proved for the selected
realizations, an annealed pressure bound `A_* m^2+o(m^2)` would yield

```math
c\le\frac{A_*-p h(\delta)-p\mathcal G(d_*)}
              {2t\sqrt p(1-2\delta)^2},
```

with the limiting values of `beta,C` used in (4).
Counting their total number after seeing the matrix is not by itself an
annealed energy estimate: the candidate set is matrix-dependent. Neither
the marginal Gaussian row-profile theorem nor the rank identity supplies
that missing correlation control.

Reproduction: `computations/transfer_reconstruction_adaptive_cluster_exact_checks_2026_09_06.py`
checks the finite adaptive inequality in outward rational arithmetic,
the hub identities and the three-coordinate compressed-Hadamard example.
