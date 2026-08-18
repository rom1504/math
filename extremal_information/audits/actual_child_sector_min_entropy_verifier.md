# Independent audit of actual-child quadratic spread

Status: **passed**.  This note checks
[`actual_child_sector_min_entropy.md`](actual_child_sector_min_entropy.md)
and Theorem 37.61, with particular attention to the explicit entropy
constant, the rank-one preimage factor, and the projective-cap exponent.

## 1. Ground cap and Hamming-sphere estimate

The normalized augmented partition function is

```math
\overline Z_A(t)=2^{-m}\sum_x\cosh(tH_A(x)).
```

If `|H_A(x_*)|=K_A`, the two terms indexed by `x_*` and `-x_*` contribute
at least `e^(tK_A)` before the factor `2^(-m)`.  Combining this with
`log Zbar_A<=d_m log cosh(t)` gives

```math
{tK_A\over m}
\le\log2+{m-1\over2}\log\cosh t
\le\log2+{\beta^2\over4}=C_\beta,                  \tag{VME.1}
```

where `t=beta/sqrt(N)` and `m<=N`.  Thus ME.2--ME.3 are correctly
normalized.

For a uniformly random `r`-set `R`, an edge is cut by `R` with probability
`2r(m-r)/(m(m-1))`.  Flipping `R` changes that edge's quadratic character
by `-1`, so its average multiplier is

```math
1-{4r(m-r)\over m(m-1)}.                            \tag{VME.2}
```

Jensen's inequality on the Hamming sphere then gives ME.6.  At `r=1`, its
exponent is exactly `4tK_A/m<=4C_beta`; hence

```math
\|\mu_{A,s}\|_\infty\le {e^{4C_\beta}\over m}
={16e^{\beta^2}\over m}.                            \tag{VME.3}
```

The finite constant in ME.6a and (37.215a) is correct (although the bound is
of course vacuous when its right side exceeds one).

## 2. Explicit exponential constant

For `0<q<1`,

```math
h(q)=q\log(1/q)+(1-q)\log(1/(1-q))
\ge q\log(1/q)+q(1-q),                              \tag{VME.4}
```

because `-log(1-q)>=q`.  Set
`q_beta=exp(-4C_beta)=exp(-beta^2)/16`, which is below `1/2`.  Then

```math
\begin{aligned}
h(q_\beta)-4C_\beta q_\beta(1-q_\beta)
&\ge4C_\beta q_\beta+q_\beta(1-q_\beta)
       -4C_\beta q_\beta(1-q_\beta)\\
&=q_\beta(1-q_\beta)+4C_\beta q_\beta^2\\
&\ge q_\beta.
\end{aligned}                                       \tag{VME.5}
```

Therefore

```math
\eta_\beta\ge q_\beta={e^{-\beta^2}\over16},       \tag{VME.6}
```

as claimed.  Applying the type-class estimate to a fixed asymptotic radius
`r/m -> q` proves the uniform liminf ME.9.  Since
`underline eta_beta=e^(-1)q_beta<eta_beta`, the fixed margin also absorbs
all polynomial factors and the later factor two, justifying the eventually
uniform version ME.13a.

## 3. Rank-one preimages and finite pair constant

Under the exact sector mixture,

```math
\nu_\epsilon(s,x,y)
=\pi_s^{(\epsilon)}\mu_{A,s}(x)\mu_{D,\epsilon s}(y),
\qquad Q=sxy^{\mathsf T}.
```

For a fixed signed `Q` and fixed `s`, there are exactly two preimages,
`(x,y)` and `(-x,-y)`, and their quadratic Gibbs weights agree.  Hence

```math
\begin{aligned}
\Pr(Q)
&\le\sum_s2\pi_s^{(\epsilon)}
       \|\mu_{A,s}\|_\infty
       \|\mu_{D,\epsilon s}\|_\infty\\
&\le2\max_s\bigl(
       \|\mu_{A,s}\|_\infty
       \|\mu_{D,\epsilon s}\|_\infty\bigr).       \tag{VME.7}
\end{aligned}
```

The second sector is already accounted for by `sum_s pi_s=1`; it does not
add another factor two.  Substitution of (VME.3) for both children gives

```math
\|\mu_\epsilon\|_\infty
\le2{16e^{\beta^2}\over m}{16e^{\beta^2}\over n}
={512e^{2\beta^2}\over mn}.                         \tag{VME.8}
```

Thus ME.11, ME.11a, (37.219), and (37.219a) have the correct preimage and
orientation factors.  On comparable splits, the two child exponents add to
`eta_beta(m+n)=eta_beta N`; collision is at most maximum atom, so the stated
Rényi-2 bound follows.

## 4. Projective cap exponent

The condition

```math
{|\langle x,u\rangle|\over m}\ge1-\delta
```

is equivalent to projective Hamming distance at most `delta m/2`.  For fixed
`0<delta<1`, its number of projective classes is

```math
\sum_{j\le\delta m/2}{m\choose j}
=\exp\{m h(\delta/2)+o(m)\}.                         \tag{VME.9}
```

A projective sector atom has twice the corresponding full-cube atom; this
constant is swallowed by `o(m)`.  Multiplying (VME.9) by the maximum-atom
bound yields

```math
\kappa_{A,s}(\delta)
\le\exp\{-[\eta_\beta-h(\delta/2)]m+o(m)\}.          \tag{VME.10}
```

The upper comparison in PA.13c is a convex sector sum of a left and a right
cap mass.  Their exponents add, so for a comparable split

```math
\Xi_\epsilon(\delta)
\le\exp\{-[\eta_\beta-h(\delta/2)]N+o(N)\}.          \tag{VME.11}
```

There is no extra orientation entropy in the leading exponent.  The union
and binary-data-processing corollaries consequently have thresholds
`xi+h(delta/2)<eta_beta` and coefficient
`eta_beta-h(delta/2)`, respectively, exactly as stated.

## 5. Scope and verdict

The theorem proves exponential atom/collision diffuseness and excludes
fixed posterior mass on a low-rate catalogue or sufficiently narrow
projective cap at sublinear KL cost.  It does **not** upper-bound posterior
KL, rule out retuning across exponentially many tiny atoms, or prove target
reach.  Both the source and Theorem 37.61 preserve this limitation.

All requested constants and normalizations pass independent verification.
No mathematical correction is required.
