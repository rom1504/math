# Independent audit: bounded-cap thin-tail pressure consequence

**Verdict: PASS.**

The finite-temperature sandwich, the integer active-dimension schedule, and
the weighted-star countermodel are all correct in the repository's
normalized-expectation and half-quadratic conventions.  In particular, the
countermodel simultaneously has the complete-edge variance, fixed-order
entropy squeeze, adjacent-order regularity, fixed-parameter monotonicity,
exact edge-centered subadditivity, uniform fixed-rate one- and two-sided
thin tails, and nonconvergent ground and every fixed-temperature SK-diagonal
pressure.

The source audited and frozen for this verdict is

```text
extremal_information/drafts/bounded_cap_thin_tail_pressure_consequence.md
sha256 25be86fb5c98514bce4578cda91f05f4c71f3a474da6d3877b3b135bfb6b3a32
```

No source repair is required.  Two harmless scope qualifications should be
kept in mind: the constructed forms are sparse with real weights of order
`sqrt(n)`, not complete signings, and (TP.13) is a scalar inequality
matching the complete-signing composition law rather than a physical
restriction/gluing realization of these particular landscapes.

## 1. Pressure normalization and Proposition TP.1: PASS

The source uses the uniform spin expectation

```math
\psi_A(\beta)
=\frac1n\log\left[2^{-n}\sum_x
 \cosh\left(\frac{\beta H_A(x)}{\sqrt n}\right)\right].
```

Thus replacing the expectation by the unnormalized sum adds exactly
`log 2`; there is no missing factor two from the edge-sum convention
`H_A(x)=sum_(i<j)a_ijx_ix_j=(1/2)x^TAx`.

Let `Delta(x)=Q(A)-|H_A(x)|`.  Assumption (TP.2) gives

```math
\mathbb E_x e^{-\beta\Delta(x)/\sqrt n}
\le e^{-\kappa n}+e^{-\beta dn}.
```

Therefore

```math
\mathbb E_x e^{\beta|H_A(x)|/\sqrt n}
\le 2e^{n\{\beta q(A)-\min(\kappa,\beta d)\}},
```

and `cosh u<=e^|u|` gives the upper half of (TP.3), including the term
`(log 2)/n`.

For the lower half, a homogeneous quadratic is even, so a maximizing spin
`x_*` and its distinct antipode both have absolute energy `Q(A)`.  Each
contributes at least `e^(beta q(A)n)/2` by
`cosh u>=e^|u|/2`.  After division by `2^n`, their combined contribution is
`2^(-n)e^(beta q(A)n)`.  Hence

```math
\psi_A(\beta)\ge\beta q(A)-\log2,
```

exactly as stated.

For the minimized pressure, the lower inequality holds for every signing,
whereas the upper inequality may be evaluated at an exact minimizer.  The
standard random-sign upper bound places exact minimizers in the `C=1`
bounded-cap class for all sufficiently large orders.  Theorem 36.26 then
gives fixed `d_1,kappa_1`, so (TP.7) follows with
`g(beta)=min(kappa_1,beta d_1)>0`.  This is a same-order strict roof deficit,
not a comparison between orders.

## 2. Lemma TP.2 and the integer schedule: PASS

Start at `(n,s)=(2,1)`.  During a growth phase both `n` and `s` increase by
one, so `s/n` increases toward one; during a decay phase only `n` increases,
so `s/n` decreases toward zero.  Every phase is finite.  Passing from ratio
`3/5` to `4/5` and back requires a number of steps proportional to the
current order, so endpoint orders and phase lengths diverge.  The first
threshold crossing overshoots by `O(1/n)`, proving

```math
\liminf s_n/n=3/5,
\qquad\limsup s_n/n=4/5.
```

The bounds `1<=s_n<=n-1` are invariant under both updates.  At a decay
endpoint of order `n`, the preceding ratio is greater than `3/5`, whence
`s_n>3(n-1)/5>=(n-1)/2`; the initial growth starts at equality, and both
updates preserve the needed inequality until the next endpoint.  Thus
`2s_n>=n-1` at every order.

Finally, if `s_(n+1)=s_n`, (TP.9) is immediate.  If it equals `s_n+1`, then

```math
\frac{s_n+1}{E_{n+1}}\le\frac{s_n}{E_n}
\quad\Longleftrightarrow\quad
(s_n+1)(n-1)\le s_n(n+1)
\quad\Longleftrightarrow\quad
n-1\le2s_n.
```

Hence `theta_n=sqrt(s_n/E_n)` is nonincreasing exactly as required later.

## 3. Weighted-star landscape and all elementary properties: PASS

For

```math
H_n(x)=b_nx_1\sum_{j=2}^{s_n+1}x_j,
\qquad b_n=\sqrt{E_n/s_n},
```

the `s_n` products `x_1x_j` are independent uniform signs.  Consequently

```math
\|H_n\|_\infty=s_nb_n=\sqrt{s_nE_n},
\qquad
F_n(t)=s_n\log\cosh(b_nt),
```

and

```math
q_n^2=\frac{s_n}{n}\frac{n-1}{2n}.
```

The endpoint subsequences therefore give the separated limit values
`sqrt(3/10)` and `sqrt(2/5)`.  The schedule in fact traverses intermediate
ratios too; the theorem only claims the existence of these two separated
limit points, not that they are the only ones.

Both `s_n` and `E_n` are nondecreasing, so the unnormalized maximum
`sqrt(s_nE_n)` is nondecreasing.  Moreover

```math
q_n^2=\frac{s_n(n-1)}{2n^2}
```

and one update changes `s_n` by at most one.  Since `q_n>=1/4`, direct
subtraction gives `|q_(n+1)-q_n|=O(1/n)` uniformly.

The formula for `F_n` proves evenness, analyticity, convexity, and `F_n(0)=0`.
It also gives the exact complete-edge variance

```math
F_n''(0)=s_nb_n^2=E_n=\binom n2.
```

Using `|u|-log2<=log cosh u<=|u|` gives

```math
a_n|t|-s_n\log2\le F_n(t)\le a_n|t|,
```

which implies (TP.12) because `s_n<=n`.  Finally, TP.9 makes
`b_n=sqrt(E_n/s_n)` nondecreasing, while `s_n` itself is nondecreasing.
Since `log cosh(bt)` is nonnegative and nondecreasing in `b` for fixed
`t>=0`, `F_n(t)` is nondecreasing in the order.

## 4. Exact centered subadditivity: PASS

For `theta=sqrt(s/E)`, one has

```math
F_n(t)=E_n\theta_n^2\log\cosh(t/\theta_n),
\qquad
R_n(t)=E_nh_{\theta_n}(t).
```

Here `0<theta_n<=1`, because `s_n<=n-1<=E_n` for `n>=2`.  Differentiating
the first term of `h_theta` gives

```math
\partial_\theta\{\theta^2\log\cosh(t/\theta)\}
=\theta\{2\log\cosh u-u\tanh u\},
\qquad u=t/\theta.
```

For `k(u)=2log cosh(u)-u tanh(u)`, one has

```math
k'(u)=\tanh u-u\operatorname{sech}^2u,
\qquad
k''(u)=2u\operatorname{sech}^2u\tanh u\ge0
```

on `u>=0`, and `k(0)=k'(0)=0`; evenness covers negative `u`.  Thus
`h_theta` is nondecreasing in `theta`, and `h_theta<=h_1=0`.

Since the schedule makes `theta_n` nonincreasing, for `N=m+n` one has
`theta_N<=theta_m,theta_n`.  Using
`E_N=E_m+E_n+mn`, first negativity and then monotonicity yield

```math
R_N(t)
\le(E_m+E_n)h_{\theta_N}(t)
\le E_mh_{\theta_m}(t)+E_nh_{\theta_n}(t).
```

This proves (TP.13) for all orders in the declared sequence (in particular
`m,n>=2`) and every real `t`, with no defect term.

## 5. Uniform two-sided and one-sided thin tails: PASS

If `k` active products carry the minority sign, then

```math
a_n-|H_n(x)|=2b_nk.
```

The schedule bound gives

```math
q_n=\sqrt{\frac{s_n(n-1)}{2n^2}}
\ge\frac{n-1}{2n}\ge\frac14.
```

With `d=1/16`, membership in the absolute `dn^(3/2)` endpoint layer forces

```math
\frac{k}{s_n}<\frac d{2q_n}\le\frac18.
```

Each active product word has exactly `2^(n-s_n)` preimages.  There are at
most two endpoint balls, so the layer size is bounded by

```math
2^{n-s_n+1}\sum_{k\le s_n/8}\binom{s_n}{k}
\le 2^{n-s_n+1}e^{s_nh(1/8)}.
```

Put `c_*=log2-h(1/8)>0`.  Since `s_n>=(n-1)/2`, its logarithm is at most

```math
n\log2-\frac{c_*}{2}(n-1)+\log2
\le n(\log2-c_*/4)
```

for all sufficiently large `n`.  This proves (TP.14) with the stated
`kappa=c_*/4`.  Retaining only the all-positive or all-negative active-word
ball removes the leading factor two and proves each one-sided version.

## 6. Every fixed positive SK temperature oscillates: PASS

Writing `alpha_n=s_n/n` gives exactly

```math
\frac1nF_n(\beta/\sqrt n)
=\alpha_n\log\cosh\left(
 \beta\sqrt{\frac{1-1/n}{2\alpha_n}}\right).
```

Along `alpha_n->alpha`, the limit is

```math
J_\beta(\alpha)
=\alpha\log\cosh\left(\frac\beta{\sqrt{2\alpha}}\right).
```

For `beta>0`, put `u=beta/sqrt(2alpha)`.  Then

```math
J_\beta'(\alpha)
=\log\cosh u-\frac u2\tanh u
=\frac12k(u)>0.
```

Hence the `alpha=3/5` and `alpha=4/5` subsequences have different limits for
every fixed positive `beta`.  At `beta=0` the pressure is identically zero,
which is why the theorem correctly restricts the nonconvergence claim to
`beta>0`.

## 7. Archive collision and exact scope

The construction is the finite-landscape realization of the same analytic
family used in Theorem 2.1 of
`artifacts/finite_temperature_scalar_no_go.md`: there the quantities
`L_n` and `theta_n` were chosen as real analytic parameters, whereas here
`L_n=s_n` is an integer active dimension and the formula is realized by an
actual weighted Boolean star.  The centered-subadditivity calculation is
therefore archived algebra, as the source states.  The genuine increment is
that the countermodel now has an honest finite quadratic landscape and the
new fixed-rate endpoint-tail property supplied by the schedule.

This does not collide with the isolated-extremizer example MH.11, which has
neither the complete-edge variance/subadditivity package nor the present
oscillating finite-star realization.

The no-go proves exactly that the following **scalar data** remain
insufficient:

- a uniform strict soft/ground roof deficit;
- exact variance and the spin-entropy squeeze;
- scalar fixed-parameter monotonicity and adjacent regularity;
- exact edge-centered subadditivity;
- fixed-rate endpoint thinness.

It does not falsify any argument using exact `+-1` completeness, an actual
cross-order bridge coupling, joint shell/overlap information, vector-valued
responses, microcanonical composition, or another signing-specific state.
It proves no statement about whether the actual minimized pressure or
`M_n/n^(3/2)` converges.  The source's final judgment is therefore properly
negative and non-overclaimed: Theorem 36.26 improves a same-order sandwich,
but the scalar deficit alone is not the missing thermodynamic compactness
mechanism.
