# Independent audit: phase refresh synchronization

Audited files:

- [`phase_refresh_synchronization.md`](phase_refresh_synchronization.md)
- [`verify_phase_refresh_synchronization.py`](../experiments/verify_phase_refresh_synchronization.py)

**Verdict: PASS, with three nonfatal presentation qualifications.**  The
operator normalization, the one-sided absolute-response pullback, Markov
iteration, Doeblin maximum principle, two-state falsifier, and quantitative
Walsh obstruction are all correct.  The theorem does not assume constancy of
the limiting profile and its operator certificate is genuinely weaker than a
comparison of the full Boolean response landscapes.

The qualifications are:

1. the abstract kernel statement should include the usual measurability of
   `(y,U) -> U^*T_(r+1,y)U` (automatic in every finite presentation used here);
2. the regular-Hadamard matching example should explicitly verify uniform
   recovery, and should be described as an illustrative certificate rather
   than evidence that the theorem is needed there, because a direct operator
   comparison already proves its phase collapse; and
3. “sharp falsifier” means sharp in the error-to-refreshed-mass scale, not that
   (PR.7) is a necessary condition for every system.

No audited theorem, draft, or verifier was changed in this audit.

## 1. Probability-space and quadratic-form normalization

On the uniform `N`-point probability space, identify functions with column
vectors and let a matrix act in the usual coordinate convention.  If

```math
T_C={C\over\sqrt N},
```

then

```math
\langle f,T_Cf\rangle_{L^2(\mu_N)}
 ={1\over N^{3/2}}f^TCf.                              \tag{A.PR.1}
```

Consequently, for a symmetric hollow matrix `C`,

```math
\sup_{\|f\|_\infty\le1}|\langle f,T_Cf\rangle|
 ={1\over N^{3/2}}\max_{x\in\{+-1\}^N}|x^TCx|
 ={2Q(C)\over N^{3/2}}.                               \tag{A.PR.2}
```

The first equality follows by maximizing one coordinate at a time: with
zero diagonal, the quadratic form is affine in each coordinate while the
others are fixed, and the absolute value of an affine function is convex on
`[-1,1]`.  The theorem itself does not require hollowness.  For a nonhollow
operator, (PR.1) is simply the relaxed cube response, and every step of the
proof still holds.

For signed equal-fibre replication, write the new coordinates as `(a,i)`,
`1<=a<=h`, and take

```math
(Uf)_(a,i)=u_af_i,
\qquad u\in\{+-1\}^h.                                \tag{A.PR.3}
```

Then `U` is an isometry for the two uniform `L^2` norms and is an exact
`L^infinity` isometry.  Its probability-space adjoint is

```math
(U^*g)_i={1\over h}\sum_a u_ag_(a,i).                 \tag{A.PR.4}
```

If `H u=sqrt(h)u`, direct substitution gives the exact normalized pullback

```math
U^*{H\otimes C\over\sqrt{hN}}U={C\over\sqrt N}.       \tag{A.PR.5}
```

Indeed `u^THu=h^(3/2)`; the `1/h` in the adjoint and the
`1/sqrt(h)` in the new normalization cancel the remaining factor.  Thus no
factor `h`, `sqrt(h)`, `N`, or two is missing in the intended Hadamard
application.

## 2. Absolute response survives the common operator average

Fix `f` with `||f||_infinity<=1`.  The normalized probability-space norm
satisfies `||f||_2<=1`, so (PR.3) implies

```math
\left|\langle f,T_(r,x)f\rangle
-\int\langle Uf,T_(r+1,y)Uf\rangle\,d\Gamma_(r,x)\right|
\le\epsilon_r.                                       \tag{A.PR.6}
```

The triangle inequality is applied only after the channels have been added
at operator level:

```math
|\langle f,T_(r,x)f\rangle|
\le\int|\langle Uf,T_(r+1,y)Uf\rangle|\,d\Gamma
  +\epsilon_r.                                        \tag{A.PR.7}
```

Since `Uf` remains in the target unit cube, the integrand is bounded by
`phi_(r+1)(y)`.  The phase marginal of one witness-independent `Gamma` is
exactly `P_r`, hence

```math
phi_r\le P_rphi_(r+1)+\epsilon_r.                     \tag{A.PR.8}
```

This verifies the important absolute-value point: the argument neither
polarizes nor separately pays positive and negative channels.  What would
be invalid is choosing `Gamma` after seeing `f` or after selecting a
maximizer.  The draft explicitly excludes that dependence.

Iterating (A.PR.8) is legitimate because a Markov kernel preserves order
and fixes constants.  More explicitly,

```math
phi_r\le P_r\cdots P_(r+l-1)phi_(r+l)
 +\epsilon_r+\cdots+\epsilon_(r+l-1).                 \tag{A.PR.9}
```

No unrecorded kernel norm or product of errors occurs.

## 3. Doeblin maximum principle and constants

Let `x_*` maximize the continuous limit `phi`, with maximum `M`.  Uniform
recovery gives

```math
M-\omega_r\le\phi_r(x_*).                             \tag{A.PR.10}
```

The minorization is exactly equivalent to

```math
K_(r,l_r)(x_*,\mathord\cdot)
=\alpha_r\nu+(1-\alpha_r)\rho_r                      \tag{A.PR.11}
```

for a probability measure `rho_r` (with the evident interpretation when
`alpha_r=1`).  Combining (A.PR.9)--(A.PR.11), using
`phi<=M`, yields the finite-level inequality

```math
\alpha_r\left(M-\int\phi\,d\nu\right)
\le\omega_r+\omega_(r+l_r)+E_r.                      \tag{A.PR.12}
```

Since (A.PR.12) holds for every `r`, taking `liminf` proves (PR.8).  If its
right side vanishes, full topological support and continuity imply
`phi=M` everywhere: otherwise a nonempty open sublevel set has positive
`nu` mass.  On a finite space,

```math
M-\int\phi\,d\nu
=\sum_x\nu_x(M-\phi(x))
\ge\min_x\nu_x\,\operatorname{osc}(\phi),             \tag{A.PR.13}
```

so the coefficient `1/mu` in (PR.9) is correct.

There is no circular use of constancy in this argument.  Uniform convergence
to some continuous `phi` is assumed, but equality of its phase values is the
conclusion.  Likewise, (PR.5) is a transition-kernel condition and (PR.3) is
an operator-norm condition, neither of which refers to the maximizing
Boolean spin.

For complete formal generality, the statement should say that the pulled-
back operator field is measurable and Bochner integrable.  All spaces and
branch families in the applications are finite, where this is automatic.

## 4. Regular-Hadamard matching-flip example

Let `N=h^r`, let `D_r` be the diagonal of `H_r`, and let `E_(r,x)` be the
symmetric perturbation obtained by flipping one perfect matching.  In a
matching-adapted ordering, `E_(r,x)` is a direct sum of signed blocks

```math
\begin{pmatrix}0&2\\2&0\end{pmatrix},
```

so `||E_(r,x)||=2`; also `||D_r||=1`.  Thus

```math
\left\|{A_(r,x)\over\sqrt N}
-{H_r\over\sqrt N}\right\|\le {3\over\sqrt N}.       \tag{A.PR.14}
```

Use (A.PR.5) for every branch.  Pullback by an `L^2` contraction cannot
increase operator norm, so for arbitrary source phase `x`, target phase
`y`, and therefore for an arbitrary average over `y`,

```math
\left\|{A_(r,x)\over\sqrt N}
-U^*{A_(r+1,y)\over\sqrt{hN}}U\right\|
\le {3\over\sqrt N}+{3\over\sqrt{hN}},                \tag{A.PR.15}
```

which is the stated bound.

The draft omits one easy verification of hypothesis 1.  The unperturbed
operator `H_r/sqrt(N)` has norm one, and `u^(tensor r)` is a Boolean unit-
cube witness attaining response one.  Hence its response equals one.
Response is one-Lipschitz in operator norm, and (A.PR.14) gives

```math
\sup_x|\phi_r(x)-1|\le {3\over\sqrt N}.               \tag{A.PR.16}
```

For a fixed finite phase set this proves the required uniform convergence.
It also exposes a limitation of this example: (A.PR.16) proves phase
collapse directly, without refresh.  The example correctly demonstrates
that (PR.3) is a finite spectral certificate even after `Theta(N)` sign
flips, but it is not evidence that refresh adds power on this particular
class.  Canonical wording should avoid presenting it as such.

## 5. Two-state falsifier

For `T_0=0`, `T_1=I`, and

```math
P_r=(1-a_r)I+a_r\nu,
\qquad \nu=(\delta_0+\delta_1)/2,                     \tag{A.PR.17}
```

the identity pullback gives averaged operators `a_rI/2` from state zero and
`(1-a_r/2)I` from state one.  The defect is exactly `a_r/2` in both cases.
Each row has full support and has minorization `P_r(x,.)>=a_r nu(.)`, but

```math
{\epsilon_r\over\alpha_r}={1\over2}.                 \tag{A.PR.18}
```

For `a_r=2^(-r)`, the future refresh mass is summable and its tail tends to
zero, so the two phases retain distinct responses forever.  This correctly
falsifies “pointwise full support plus vanishing defect.”  It establishes
the correct *order* of the error-versus-refresh condition.  It does not show
that (PR.7) is logically necessary for collapse in every other system, so
the word “sharp” should be read in that restricted sense.

The use of the diagonal one-dimensional operator is harmless: this is an
abstract falsifier for PR.1, whose definition does not impose hollowness or
sign entries, not a claimed dense-sign construction.

## 6. Walsh-prefix obstruction

At phase `t` and scale `r`, take

```math
\phi_r(t)={2Q(A_(\lfloor t4^r\rfloor))
             \over\lfloor t4^r\rfloor^{3/2}}.
```

Theorem 30.1 gives uniform convergence `phi_r -> Phi=2L`.  Its explicit
Walsh certificate gives

```math
\Phi(1)=\Phi(4)=1,
\qquad \Phi(3)\ge {89\over48\sqrt3}.                 \tag{A.PR.19}
```

To check the global upper bound used in the draft, a prefix of order
`n=floor(t4^r)` is a coordinate compression of `H_(r+1)`, whose norm is
`2sqrt(4^r)`.  Hollowing adds operator norm at most one.  Therefore

```math
\phi_r(t)\le {2\sqrt{4^r}+1\over\sqrt n},
```

and passage to the limit gives `Phi(t)<=2/sqrt(t)<=2`.  For

```math
\nu_*={99\over200}(\delta_1+\delta_4)+{1\over100}\lambda,
```

the total atomic weight is `99/100`, the measure has full support because
of its Lebesgue component, and

```math
\int\Phi\,d\nu_*
\le {99\over100}+{1\over100}\,2=1.01.                \tag{A.PR.20}
```

Thus

```math
M-\int\Phi\,d\nu_*
\ge {89\over48\sqrt3}-1.01
=0.060503624121\ldots.                                \tag{A.PR.21}
```

Applying (PR.8) proves (PR.22) with the displayed direction and constant.
This is a conditional impossibility statement about *any* proposed refresh
certificate obeying (PR.21), not an assertion that such a kernel exists.
It does not assume the profile is constant; it uses its independently proved
nonconstancy to lower-bound the necessary transfer defect.

## 7. Information content and verifier coverage

The certificate is strictly weaker than full response comparison in its
query obligation: one verifies a witness-independent average of pulled-back
linear operators in `2->2` norm, rather than aligning maximizers or comparing
all Boolean response values.  It is also weaker than rowwise/path
simulation.  However, a raw operator can still contain quadratically many
coefficients.  PR.1 alone is therefore a strict **mathematical reduction**
of the response obligation, not yet a theorem that the certificate always
has sub-landscape description complexity.  That distinction should be
preserved in any portfolio claim.

The Python verifier passes.  It correctly checks:

- the exact `a_r/2` defects in the two-state example;
- the finite-state Doeblin deficit identity and the `mu` oscillation bound;
- the Walsh constant and normalization conversion.

As advertised, it does not machine-check the functional-analytic theorem,
Hadamard pullback, uniform-recovery step, or prefix operator bound.  Those
are checked algebraically above.  No additional computational certificate
is needed for the theorem.
