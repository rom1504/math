# Independent audit: expander refresh, bounded rank, and non-Gaussian tangents

Scope: independent audit of

- `expander_phase_refresh_complexity.md`;
- `bounded_operator_rank_barrier.md`;
- `non_gaussian_tangent_closure.md`.

I reconstructed the arguments from their hypotheses rather than accepting the
draft narratives.  I also ran

```text
python extremal_information/experiments/verify_expander_phase_refresh_complexity.py
python extremal_information/experiments/verify_bounded_operator_rank_barrier.py
python extremal_information/experiments/verify_non_gaussian_tangent_closure.py
```

All three verifiers pass.  They are finite smoke tests; the verdicts below are
based on the proofs.

## 1. Expander phase refresh

**Verdict: REPAIR (core proposition, theorem, Walsh application, and constants
all pass; only endpoint hypotheses and two wording claims need tightening).**

### ER.1

The conjugation identity is exact.  With

```math
U_{r,g,h}=V_{r+1,h}^{*}L_rV_{r,g},
```

unitarity of the signed coordinate permutations gives

```math
U_{r,g,h}^{*}T_{r+1,h}U_{r,g,h}
=V_{r,g}^{*}L_r^{*}S_{r+1}L_rV_{r,g},
```

independently of `h`.  Hence an arbitrary phase kernel really may be inserted.
The `L^2` and `L^infinity` contraction properties survive conjugation, and a
signed permutation bijects the Boolean cube, so the response is gauge
invariant.

For the displayed Walsh gauge, direct calculation gives

```math
\left\|(H-DHD)/2\right\|=2.
```

Tensoring with normalized Hadamard factors preserves this norm, and diagonal
hollowing cancels because diagonal conjugation does not change diagonal
entries.  Thus the asserted distance two and the radius-one obstruction are
correct.

### ER.2 and the rearrangement in ER.13

From uniform recovery and positivity,

```math
g\le P_jg+\delta_j.
```

Because every `P_j` preserves `pi`, its mean-zero subspace is invariant, and
the product contracts that subspace by `rho^t` even without reversibility.
For mean-zero `u`, point evaluation is represented by
`1_{x}/pi(x)-1`, whose norm is `sqrt(1/pi(x)-1)`.  This proves ER.10 (with the
slightly weaker displayed factor `1/sqrt(pi(x))`).

The ceiling rearrangement is also correct.  Put

```math
A={\log(2B\sqrt{S/\kappa}/D)\over\lambda},
\qquad \lambda=\log(1/\rho).
```

Since `0<D<=B` and `pi(x)>=kappa/S`, one has `A>0` and
`t_*<=A+1`.  ER.12 gives

```math
{D\over2\delta}\le t_*
\le1+{\log(2B/D)+\tfrac12\log(S/\kappa)\over\lambda},
```

which is exactly

```math
\log S\ge\log\kappa+{\lambda D\over\delta}
-2\lambda-2\log(2B/D).
```

No factor of two or logarithm-base error is present.  For `rho=1/2` and
`delta<=C/sqrt(N)`, division by `log 2` cancels the `lambda=log 2` in the
leading term, so the bit coefficient in ER.21 is indeed `D_*/C`.

### Walsh stationary law

The law is valid.  If `S=|X|`, put mass `99/200` at each of `1,4` and spread
the remaining `1/100` uniformly over `X\setminus{1,4}`.  Every remaining atom
has mass `1/[100(S-2)]>=1/(100S)`, in particular the atom at `3`.  The endpoint
identities and the global cap imply

```math
\int\Phi\,d\pi\le .99\cdot1+.01\cdot2=1.01,
```

so

```math
D_3\ge {89\over48\sqrt3}-1.01
=0.060503624122\ldots .
```

ER.20--ER.21 therefore follow conditionally on the stated common stationary
law, uniform scrambling, and combined transfer-plus-recovery error.  The draft
correctly does **not** claim a lower bound for presentations outside that
architecture.

### Required repairs

1. State `0<rho<1` in the logarithmic form ER.11--ER.13.  The case `rho=0`
   is valid but is a separate one-step statement; `log(1/rho)` is undefined.
2. State `epsilon_j>=0`, `omega_j>=0`, and `kappa>0`.  These are intended by
   “toll/error,” but ER.13 divides by the resulting `delta`.
3. In ER.13 say `0<delta<D/2`; if `delta=0`, ER.10 already contradicts
   `D>0` after enough steps.
4. Replace “sharp qualitative dichotomy” by “quantitative dichotomy” unless a
   matching construction is supplied.  The proof uses the non-sharp bound
   `||g-Pi g||_2<=B` and does not establish sharp constants.

Novelty is correctly scoped: the mixing estimate itself is classical; the
project-level contribution is its use as a response-recovery state-size
lower bound and its gauge/semantic separation.

## 2. Bounded-operator rank barrier

**Verdict: REPAIR (BR.1 passes exactly; BR.2 must consistently describe an
operator-norm certificate, not the actual Boolean error of every fixed
matrix).**

The Frobenius argument is exact.  If

```math
r=\#\{j:\sigma_j>\epsilon\sqrt n\},
```

then the strict threshold means every omitted value satisfies
`sigma_j<=epsilon sqrt(n)`, while every retained value is at most
`C sqrt(n)`.  Hence

```math
n^2=\sum_j\sigma_j^2
\le rnC^2+(n-r)n\epsilon^2,
```

and BR.3 follows with the correct denominator and strictness.  The assumptions
`C>=1` and `epsilon<1` make the denominator positive.  BR.4 follows from the
same identity and is exact for the all-ones matrix.  Hadamard equality in BR.3
is also correct for every `epsilon<1` when `C=1`.

Eckart--Young additionally shows that **any** rank-`r` approximation with
operator error at most `epsilon sqrt(n)` needs

```math
r\ge\#\{j:\sigma_j>\epsilon\sqrt n\}.
```

Thus the claimed linear barrier for SVD/operator-norm interfaces is valid.

The needed wording repair is in BR.2.  Theorem 18.7 supplies the certified
upper bound

```math
\sup_{x,y\in\{\pm1\}^n}|x^T(R-R_r)y|
\le n\sigma_{r+1}(R).
```

For a fixed `R`, equality need not hold because its leading omitted singular
vectors need not be Boolean.  Therefore “gives worst-case Boolean response
error” should read “gives the certified worst-case upper bound,” and the next
sentence should say that BR.3 is necessary **for this operator-norm/SVD
certificate**.  The draft's later scope disclaimer has the right intent, but
the earlier equality-like wording is too strong.

BR.1 is the classical stable-rank/Frobenius calculation specialized to dense
sign bridges.  Its novelty should be presented as closing this repository's
rank-growth benchmark, not as a new general matrix theorem.  The explicit
warning that full algebraic rank need not imply contextual incompressibility
is correct and essential.

## 3. Non-Gaussian tangent closure

**Verdict: REPAIR (NG.1, NG.2, and the off-centre exponent pass; narrow the
title and add one scaling lemma before claiming the discrete carrier
consequence).**

### Roof and tangent exponents

Strict convexity for every `p>1` gives the unique split and the coefficient

```math
c=(a^{-1/(p-1)}+b^{-1/(p-1)})^{-(p-1)},
```

so NG.1 is correct.

At output zero the change of variables
`k=n^{1-1/p}u` gives the factor `n^{1-1/p}` in NG.4.  At fixed `z!=0`, both
pieces of the minimizing split are nonzero.  Consequently

```math
p(p-1)\{a|x_*|^{p-2}+b|z-x_*|^{p-2}\}
```

is positive and finite for the full range `p>1`, including `1<p<2`.
Ordinary lattice Laplace asymptotics therefore give `Theta(sqrt n)` after
the leading exponential is removed.  The two exponents in NG.5 are correct;
there is no hidden restriction to `p>=2`.

### CLT rigidity

NG.2 is valid.  Equality within the same fixed-`p` scale family says
`X_1+X_2 =_d sX`.  Finite nonzero variance forces `s=sqrt2`; dyadic iteration
makes every normalized dyadic sum have the law of `X`; and the classical CLT
forces that law to be Gaussian.  Comparing positive densities gives `p=2`,
and scaling then gives `b=a/2`.  No infinite-divisibility assumption is being
smuggled in.

### Required repairs

1. The title “Non-Gaussian tangent types do not close under convolution” is
   broader than the theorem.  What is proved is that the **fixed
   power-exponential/generalized-Gaussian family** closes under self-
   convolution only for `p=2`.  Other non-Gaussian finite convolution
   semigroups are explicitly not excluded later, so the title should match
   that scope.
2. Before the “Consequently” paragraph after NG.2, add the local-limit bridge

   ```math
   n^{-(1-1/p)}(A_n^a*A_n^a)
   (\lfloor n^{1-1/p}u\rfloor)
   \longrightarrow
   \int_{\mathbb R}e^{-a|v|^p}e^{-a|u-v|^p}\,dv
   ```

   locally uniformly in `u` (up to the stated normalizing constants).  It is
   an elementary dominated Riemann-sum argument, but without writing it the
   jump from the continuous density theorem to failure of the discrete
   reusable tangent carrier is implicit rather than proved.
3. Describe NG.2/its proof as classical Gaussian stability via the CLT.  The
   new project contribution is the falsification of the proposed tangent
   carrier, not the rigidity fact itself.

The final disclaimer is accurate: this kills the fixed-`p`
generalized-Gaussian extension of Theorem 32.2, not all finite stratified or
non-Gaussian carriers.

## Consolidated decision

No draft fails mathematically.  The theorem cores survive independent proof
reconstruction and their constants are correct.  Canonicalization should wait
for the small repairs above, especially the SVD “upper bound” wording and the
explicit discrete-to-continuous tangent scaling statement.

## Repair resolution

Before canonicalization, ER.2 was restricted to nonnegative errors and its
logarithmic form to `0<rho<1, 0<delta<D/2`, with the zero cases separated;
the unsupported word “sharp” was removed.  BR.2 now explicitly concerns the
operator/SVD certificate.  The tangent draft has a scoped title, labels the
CLT rigidity classical, and proves the locally uniform discrete tangent
limit (NG.10).  These changes implement every required repair above.
