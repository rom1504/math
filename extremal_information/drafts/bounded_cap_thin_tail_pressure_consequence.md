# Thin tails improve the soft/ground sandwich but not scalar thermodynamic compactness

**Status.** Task-local proved theorem and abstract finite-landscape no-go.
This note is deliberately scoped to the new bounded-cap thin-tail theorem
(Theorem 36.26).  It derives its exact finite-temperature consequence and
then tests whether that consequence repairs the archived scalar-pressure
bottleneck.  It does not edit the canonical theorem files and it proves no
cross-order result for complete signings.

Throughout,

```math
H_A(x)=\sum_{i<j}a_{ij}x_ix_j,
\qquad Q(A)=\max_x|H_A(x)|,
\qquad q(A)={Q(A)\over n^{3/2}},
```

and expectation over spins is uniform.  The normalized absolute canonical
pressure is

```math
\psi_A(\beta)
 ={1\over n}\log \mathbb E_x
 \cosh\!\left({\beta H_A(x)\over\sqrt n}\right).
                                                        \tag{TP.1}
```

This is the normalized-expectation convention of
`artifacts/finite_temperature_scalar_no_go.md`; adding the unnormalized
spin factor changes (TP.1) by `log 2`.

## 1. A strict finite-temperature roof deficit

### Proposition TP.1 (thin-tail pressure sandwich)

Suppose a quadratic Boolean landscape of order `n` obeys

```math
\#\{x:Q(A)-|H_A(x)|<dn^{3/2}\}
 \le 2^n e^{-\kappa n}                         \tag{TP.2}
```

for some `d,kappa>0`.  Then every `beta>0` satisfies

```math
\boxed{
\beta q(A)-\log2
\le \psi_A(\beta)
\le \beta q(A)-\min\{\kappa,\beta d\}
                 +{\log2\over n}.}             \tag{TP.3}
```

In particular, Theorem 36.26 gives constants `d_C,kappa_C>0`, uniform over
all complete signings with `Q(A)<=Cn^(3/2)`, for which (TP.3) holds with
`d=d_C` and `kappa=kappa_C`.

#### Proof

Put `Delta(x)=Q(A)-|H_A(x)|`.  Splitting the expectation according to
`Delta<dn^(3/2)` and its complement gives

```math
\begin{aligned}
\mathbb E e^{\beta|H_A|/\sqrt n}
 &=e^{\beta q(A)n}\mathbb E e^{-\beta\Delta/\sqrt n}\\
 &\le e^{\beta q(A)n}
       \{e^{-\kappa n}+e^{-\beta dn}\}\\
 &\le 2\exp\{n(\beta q(A)-\min(\kappa,\beta d))\}.
                                                        \tag{TP.4}
\end{aligned}
```

Since `cosh u<=e^|u|`, this proves the upper bound.  If `x_*` realizes
`|H_A(x_*)|=Q(A)`, then the two distinct spins `x_*,-x_*` have the same
quadratic energy.  Their two terms, together with
`cosh u>=e^|u|/2`, give

```math
\mathbb E\cosh(\beta H_A/\sqrt n)
 \ge 2^{-n}e^{\beta q(A)n},                     \tag{TP.5}
```

which proves the lower bound. `square`

There is an immediate statement for the minimized complete-sign pressure.
Let

```math
p_n(\beta)
 ={1\over n}\min_A\log\mathbb E_x
   \cosh\!\left({\beta H_A(x)\over\sqrt n}\right),
\qquad m_n={M_n\over n^{3/2}}.                  \tag{TP.6}
```

The lower half of (TP.3), applied to every signing, and the upper half,
applied to an exact minimizer and Theorem 36.26 with `C=1`, yield

```math
\boxed{
\beta m_n-\log2
\le p_n(\beta)
\le \beta m_n-g(\beta)+{\log2\over n},
\qquad
g(\beta)=\min\{\kappa_1,\beta d_1\}>0.}        \tag{TP.7}
```

Thus bounded-cap thinness gives a genuine, uniform separation of fixed-
temperature pressure from the zero-temperature roof.  It improves the old
upper inequality `p_n(beta)<=beta m_n`; it does not by itself compare two
different orders.

## 2. A finite Boolean-landscape no-go with the same thin tail

The improvement (TP.7) might conceivably rule out the oscillating scalar
countermodel in `artifacts/finite_temperature_scalar_no_go.md`.  It does
not.  The next construction realizes a strengthened countermodel by honest
finite Boolean landscapes and gives them a uniform two-sided thin tail.

### Lemma TP.2 (an oscillating monotone active dimension)

There is an integer sequence `s_n`, `n>=2`, such that

```math
1\le s_n\le n-1,
\qquad s_{n+1}-s_n\in\{0,1\},
\qquad 2s_n\ge n-1,                              \tag{TP.8}
```

```math
{s_{n+1}\over E_{n+1}}\le {s_n\over E_n},
\qquad E_n={n\choose2},                          \tag{TP.9}
```

and

```math
\liminf_n{s_n\over n}={3\over5},
\qquad
\limsup_n{s_n\over n}={4\over5}.               \tag{TP.10}
```

#### Construction and proof

Start with `s_2=1`.  In a growth phase use `s_(n+1)=s_n+1` until the ratio
first reaches `4/5`; in a decay phase use `s_(n+1)=s_n` until the ratio first
falls to `3/5`; then repeat.  The finite initial passage from `1/2` to
`4/5` is included in the first growth phase.  Every phase is finite: adding
one drives `s_n/n` toward one, while adding zero drives it toward zero.
The phase lengths tend to infinity, so the overshoot at each endpoint is
`O(1/n)`, proving (TP.10).  During growth the ratio is at least its value at
the preceding decay endpoint.  At a decay endpoint of order `n`, its order-
`n-1` predecessor had ratio greater than `3/5`, so
`s_n>3(n-1)/5>=(n-1)/2`.  The initial growth starts at equality.  This
proves the third assertion in (TP.8) at every order.

If `s_(n+1)=s_n`, (TP.9) is immediate.  If `s_(n+1)=s_n+1`, then (TP.9)
is equivalent to

```math
(s_n+1)(n-1)\le s_n(n+1),
```

or `n-1<=2s_n`, which is (TP.8). `square`

### Theorem TP.3 (thin-tail scalar-pressure no-go)

There are quadratic Boolean landscapes `H_n:{+-1}^n->R` with all of the
following properties.

1. Their absolute normalized maxima `q_n=||H_n||_infinity/n^(3/2)` do not
   converge, even though the unnormalized maxima are nondecreasing and
   `|q_(n+1)-q_n|=O(1/n)`.
2. Their normalized partition functions

   ```math
   F_n(t)=\log\mathbb E_x\cosh(tH_n(x))             \tag{TP.11}
   ```

   are even, analytic and convex, satisfy `F_n(0)=0`,
   `F_n''(0)=E_n`, the spin-entropy squeeze

   ```math
   ||H_n||_infty|t|-n\log2\le F_n(t)
                         \le||H_n||_infty|t|,       \tag{TP.12}
   ```

   and are nondecreasing in `n` for fixed `t>=0`.
3. Their exactly edge-centered pressures satisfy the complete-signing
   scalar composition law

   ```math
   R_{m+n}(t)\le R_m(t)+R_n(t),
   \qquad R_n(t)=F_n(t)-E_n\log\cosh t.             \tag{TP.13}
   ```
4. There are absolute `d,kappa>0` such that, for all sufficiently large
   `n`,

   ```math
   \#\{x:||H_n||_infty-|H_n(x)|<dn^{3/2}\}
      \le\exp\{(\log2-\kappa)n\}.                  \tag{TP.14}
   ```

   The analogous one-sided statement holds at each endpoint.
5. Nevertheless, for every fixed `beta>0`, the SK-diagonal pressures

   ```math
   {1\over n}F_n(\beta/\sqrt n)                    \tag{TP.15}
   ```

   do not converge.

Consequently, adding a fixed-rate endpoint thin-tail axiom to evenness,
analyticity, convexity, exact variance, restriction monotonicity, the spin
entropy squeeze, adjacent-order regularity, and exact centered scalar
subadditivity still does not imply a thermodynamic limit on the SK
diagonal.

#### Construction

Take the sequence `s_n` from Lemma TP.2 and put

```math
\theta_n=\sqrt{{s_n\over E_n}},
\qquad b_n={1\over\theta_n}=\sqrt{{E_n\over s_n}},
\qquad a_n=s_nb_n=\sqrt{s_nE_n}.                  \tag{TP.16}
```

On `x=(x_1,...,x_n)` define the weighted star quadratic form

```math
\boxed{H_n(x)=b_nx_1\sum_{j=2}^{s_n+1}x_j.}       \tag{TP.17}
```

This is a genuine quadratic Boolean landscape, but it is sparse and has
real weights; it is **not** asserted to be a hollow complete sign matrix.
That distinction is exactly the scope of a scalar-axiom no-go.

The products `x_1x_j`, `2<=j<=s_n+1`, are independent uniform signs.  Thus

```math
||H_n||_infty=a_n,
\qquad
F_n(t)=s_n\log\cosh(b_nt).                       \tag{TP.18}
```

Since `s_n/n` has limit points `3/5` and `4/5`,

```math
q_n={a_n\over n^{3/2}}
 =\sqrt{{s_n\over n}{n-1\over2n}}               \tag{TP.19}
```

has the two distinct limit points `sqrt(3/10)` and `sqrt(2/5)`.  Both
`s_n` and `E_n` are nondecreasing, so `a_n` is nondecreasing.  A direct
one-step calculation from `s_(n+1)-s_n in {0,1}` gives the advertised
`O(1/n)` modulus for `q_n`.  Also `s_n<=n` gives `q_n<1`, so this is a
uniformly bounded-cap landscape family at the same `n^(3/2)` scale.

Equation (TP.18) proves analyticity and convexity, while

```math
F_n''(0)=s_nb_n^2=E_n.                            \tag{TP.20}
```

The elementary bound `|u|-log2<=log cosh u<=|u|` gives (TP.12), because
`s_n<=n`.  Lemma TP.2 makes `s_n` nondecreasing and `E_n/s_n`
nondecreasing; hence both the multiplier and the argument in (TP.18) are
nondecreasing, proving fixed-`t` monotonicity.

#### Exact centered subadditivity

For `0<theta<=1`, set

```math
h_\theta(t)=\theta^2\log\cosh(t/\theta)-\log\cosh t.
                                                        \tag{TP.21}
```

The function `h_theta(t)` is nonpositive and nondecreasing in `theta`.
Indeed,

```math
{\partial\over\partial\theta}
 \{\theta^2\log\cosh(t/\theta)\}
=\theta\{2\log\cosh u-u\tanh u\},
\qquad u=t/\theta,                                  \tag{TP.22}
```

and the expression in braces is nonnegative.  Indeed, for
`k(u)=2 log cosh(u)-u tanh(u)`, one has
`k'(u)=tanh(u)-u sech^2(u)` and
`k''(u)=2u sech^2(u)tanh(u)>=0` for `u>=0`, while `k(0)=k'(0)=0`;
evenness covers negative `u`.  Also `h_1=0`.  Now (TP.9) says that
`theta_n` is nonincreasing, and

```math
R_n(t)=E_nh_{\theta_n}(t).                         \tag{TP.23}
```

For `N=m+n`, use `E_N=E_m+E_n+mn`, negativity of `h`, and then its
monotonicity:

```math
\begin{aligned}
R_N(t)
 &=(E_m+E_n+mn)h_{\theta_N}(t)\\
 &\le(E_m+E_n)h_{\theta_N}(t)\\
 &\le E_mh_{\theta_m}(t)+E_nh_{\theta_n}(t).
\end{aligned}                                      \tag{TP.24}
```

This is (TP.13) with no error term.

#### Uniform thin tail

Write the active products in (TP.17) as `w_1,...,w_(s_n)`.  If `k` of
them have the minority sign, then

```math
a_n-|H_n(x)|=2b_nk.                                \tag{TP.25}
```

From `2s_n>=n-1`,

```math
q_n\ge{n-1\over2n}\ge{1\over4}.                  \tag{TP.26}
```

Take `d=1/16`.  The event in (TP.14) therefore forces
`k<s_n/8`.  If `h(p)=-p log p-(1-p)log(1-p)` is binary entropy, its size is
at most

```math
2^{n-s_n+1}\sum_{k\le s_n/8}{s_n\choose k}
\le 2^{n-s_n+1}e^{s_nh(1/8)}.                    \tag{TP.27}
```

Since `s_n>=(n-1)/2`, (TP.14) follows, for example, with

```math
\kappa={\log2-h(1/8)\over4}>0                    \tag{TP.28}
```

for all sufficiently large `n`.  Omitting one of the two endpoint balls
gives the one-sided version.

#### Oscillation at every fixed temperature

Put `alpha_n=s_n/n`.  Equations (TP.16)--(TP.18) give

```math
{1\over n}F_n(\beta/\sqrt n)
=\alpha_n\log\cosh\!\left(
 \beta\sqrt{{1-1/n\over2\alpha_n}}\right).       \tag{TP.29}
```

Along a subsequence with `alpha_n->alpha`, this tends to

```math
J_\beta(\alpha)
=\alpha\log\cosh\!\left({\beta\over\sqrt{2\alpha}}\right).
                                                        \tag{TP.30}
```

For `beta>0` this function is strictly increasing.  With
`u=beta/sqrt(2alpha)`,

```math
J_\beta'(\alpha)
=\log\cosh u-{u\over2}\tanh u>0;                \tag{TP.31}
```

positivity is the strict version of the calculation in (TP.22).
Consequently the `3/5` and `4/5` subsequences have distinct pressure
limits.  This proves (TP.15) and the theorem. `square`

## 3. Archive comparison and research judgment

1. **Collision with the scalar no-go.**  The analytic formula in Theorem
   2.1 of `artifacts/finite_temperature_scalar_no_go.md` already proves
   that centered scalar subadditivity, convexity, restriction monotonicity,
   and the usual entropy squeeze do not control the shrinking-temperature
   diagonal.  It was an abstract analytic pressure and did not encode the
   newly proved fixed-rate endpoint tail.  Theorem TP.3 strengthens that
   no-go by realizing the pressure with finite quadratic Boolean landscapes
   and imposing exactly such a tail.
2. **Collision with rare-state examples.**  Proposition MH.11 in
   `drafts/microcanonical_hypograph_compactness.md` shows that bounded
   temperatures can miss one isolated extremizer.  It supplies neither the
   exact complete-edge centered composition law nor oscillating diagonal
   pressures.  TP.3 is therefore not merely that isolated-state example.
3. **What is genuinely new.**  Proposition TP.1 is the direct thermal
   consequence of Theorem 36.26: a positive pressure deficit
   `min(kappa_C,beta d_C)` below the ground roof.  Theorem TP.3 shows that
   this new scalar datum is still insufficient even when coupled to all of
   the archived scalar regularity and composition axioms.
4. **Cross-order verdict.**  There is **no new cross-order mechanism**.
   Thinness improves the soft/ground sandwich at each order, but it supplies
   no relation between the entropy-bearing shells of different orders and
   does not repair the characteristic
   `beta -> beta sqrt(m/(m+n))` contraction.  Any thermodynamic-limit use of
   Theorem 36.26 must therefore add signing-specific joint information: for
   example a cross-order shell coupling, a joint bridge law, or another
   non-scalar state.  TP.3 rules out treating the fixed-rate deficit itself
   as that missing input.
