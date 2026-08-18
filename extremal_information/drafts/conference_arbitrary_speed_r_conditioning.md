# Arbitrary speed-`r` conditioning is invisible at sufficiently small beta

**Status.** Task-local theorem and falsifier; no canonical edits.  This note
tests the proposed statement uniformly over *every* bridge event, including
events defined using the pressure itself.  The event theorem is true.  Its
correct law-level extension uses a max-density (Renyi-infinity) budget.
An ordinary relative-entropy budget does not suffice, and an explicit
speed-mixture counterexample is given below.

## 1. Setup

Let `A_r` be a symmetric conference signing, fix an orientation `epsilon`,
and put

```math
f_r(B)=\log\left[2^{-2r}\sum_{x,y}
\cosh\left\{{\beta\over\sqrt{2r}}
\big(H_A(x)+\epsilon H_A(y)+x^TBy\big)\right\}\right].
\tag{AC.1}
```

For the uniform bridge law `U_r`, the audited conference theorem gives

```math
{f_r(B)\over r}\longrightarrow
h_\beta=2\psi(\beta/\sqrt2)+{\beta^2\over4}
\tag{AC.2}
```

in probability and mean in the campaign high-temperature range.

We use two already audited quantitative inputs.

1. For universal constants `c_*>0`, every `s>=0` satisfies, after harmless
   small-order adjustment,

   ```math
   U_r\{\|B\|_{op}>(2+s)\sqrt r\}
   \le2e^{-c_*s^2r}.
   \tag{AC.3}
   ```

2. On any fixed strict operator-temperature set

   ```math
   \mathcal K_r(\kappa)
   =\left\{B:
    \left\|{\beta\over\sqrt{2r}}
    \begin{pmatrix}A&B\\B^T&\epsilon A\end{pmatrix}
    \right\|_{op}\le\kappa<\frac12\right\},
   \tag{AC.4}
   ```

   the convex-extension/Talagrand argument is two-sided: for each fixed
   `eta>0`,

   ```math
   U_r\{B\in\mathcal K_r(\kappa):
      |f_r(B)-h_\beta r|>\eta r\}
   \le2e^{-c_{\beta,\kappa,\eta}r^2}
   \tag{AC.5}
   ```

   for all large `r`.  The center follows from the same uniform-bridge mean
   calculation used in the regular-sector theorem.

## 2. Uniform theorem for arbitrary events

### Theorem AC.1 (every event of mass `e^{-Cr}` has the same pressure rate)

For every finite `C>=0` there is `beta_C>0` such that the following holds
for every fixed `0<beta<beta_C`, separately for both orientations.  For
*every* sequence of bridge events `F_r`, allowed to depend on `A_r`, on the
orientation, and on `f_r` itself, satisfying

```math
U_r(F_r)\ge e^{-Cr},
\tag{AC.6}
```

one has, uniformly over those events,

```math
\boxed{
{f_r(B)\over r}\longrightarrow h_\beta
\quad\text{in probability and in }L^1
\quad\text{under }U_r(\,\cdot\mid F_r).}
\tag{AC.7}
```

The quantifier order is

```text
for every C, choose beta_C;
for every fixed beta<beta_C and eta>0;
uniformly over all admissible event sequences F_r.
```

It is not claimed that one positive `beta` works simultaneously for all
`C`.

**Proof.**  Choose `s_C` so large that

```math
c_*s_C^2>C+2.
\tag{AC.8}
```

Then choose `beta_C` small enough that `beta_C<sqrt(2)/6` and

```math
{\beta_C\over\sqrt2}(3+s_C)<\kappa_C<\frac12
\tag{AC.9}
```

for some fixed `kappa_C`.  On

```math
K_r=\{\|B\|_{op}\le(2+s_C)\sqrt r\},
\tag{AC.10}
```

the block triangle inequality puts the parent in
`K_r(kappa_C)` for all large `r`.  Equations (AC.3), (AC.6) give

```math
U_r(K_r^c\mid F_r)
\le2\exp\{-(c_*s_C^2-C)r\}
\le2e^{-2r}.
\tag{AC.11}
```

For every fixed `eta>0`, (AC.5) similarly gives

```math
\begin{aligned}
&U_r\{|f_r-h_\beta r|>\eta r,\ K_r\mid F_r\}\\
&\hspace{20mm}\le
2\exp\{Cr-c_{\beta,\kappa_C,\eta}r^2\}=o(1),
\end{aligned}
\tag{AC.12}
```

uniformly in `F_r`.  This proves conditional convergence in probability.

On `K_r`, the scaled interaction has norm at most `kappa_C`, so

```math
0\le f_r(B)\le\kappa_Cr.
\tag{AC.13}
```

Globally, `f_r(B)<=C_beta r^(3/2)`.  Thus (AC.11) makes the conditional
`K_r^c` contribution to `E(f_r/r)` at most
`O_beta(sqrt r)e^{-2r}=o(1)`.  On `K_r`, the normalized variables are
uniformly bounded, and the uniform convergence in probability from
(AC.12) implies convergence in `L^1`.  This proves (AC.7). `square`

The same proof permits

```math
U_r(F_r)\ge\exp\{-Cr-o(r)\}
\tag{AC.14}
```

after increasing `s_C` by a fixed margin.

## 3. The correct law-level strengthening

### Corollary AC.2 (max-density budget)

Under the choices of AC.1, let `q_r` be any bridge law satisfying

```math
\left\|{dq_r\over dU_r}\right\|_\infty\le e^{Cr}.
\tag{AC.15}
```

Then (AC.7) holds with `q_r` in place of the conditioned law, uniformly
over all such sequences `q_r`.

**Proof.**  For every event `E`, (AC.15) gives
`q_r(E)<=e^(Cr)U_r(E)`.  Apply this directly to (AC.3) and (AC.5), then use
the same boundedness argument as in AC.1. `square`

Event conditioning is the special case
`dq/dU=1_F/U(F)<=e^(Cr)`.  Thus the operational hypothesis is not the
geometric form of the event but the absence of likelihood spikes beyond
speed `r`.

The proof also extends to a fixed Renyi order `alpha>1` with a corresponding
linear Renyi-divergence budget, by Holder's inequality and a larger choice
of `s_C`.  The `alpha=1` endpoint is qualitatively different.

## 4. Ordinary entropy `D(q||U)<=Cr` is not enough

The event theorem does **not** imply the same assertion for arbitrary laws
with only a Kullback--Leibler budget.  A law may spend a small but fixed
amount of mass on a much rarer event while keeping its average information
cost linear.

### Proposition AC.3 (KL-budget counterexample)

Fix any `C>0`, any `beta>0`, and one orientation.  There are bridge laws
`q_r` with

```math
D(q_r\|U_r)\le Cr
\tag{AC.16}
```

for all large `r`, but `f_r/r` does not converge in `q_r`-probability or
in `L^1` to `h_beta`.

**Proof.**  Fix sign vectors `u,v` and choose constants `eta>0` and `a>0`
so large that

```math
{\beta a\over\sqrt2}-2\log2>h_\beta+2\eta.
\tag{AC.17}
```

Let

```math
G_r=\{B:u^TBv\ge a r^{3/2}\}.
\tag{AC.18}
```

The left side is a sum of `r^2` independent signs.  Stirling's formula at
the nearest admissible binomial layer gives

```math
U_r(G_r)\ge
\exp\left\{-{a^2\over2}r-O_a(\log r+1)\right\}.
\tag{AC.19}
```

Pairing `y` and `-y` in the partition function gives the exact factorization

```math
\overline Z
=\mathbb E_{x,y}
 \cosh\{t(H_A(x)+\epsilon H_A(y))\}
 \cosh(t x^TBy).
\tag{AC.20}
```

The single configuration `(x,y)=(u,v)` therefore yields, for every
`B in G_r`,

```math
f_r(B)\ge {\beta a\over\sqrt2}r-2r\log2-\log2
\ge(h_\beta+\eta)r
\tag{AC.21}
```

for all large `r`.

Let `P_r=U_r(.|G_r)`.  Equation (AC.19) gives

```math
D(P_r\|U_r)=\log{1\over U_r(G_r)}\le D_0r
\tag{AC.22}
```

for a finite constant `D_0=D_0(a)` and all large `r`.  Set

```math
theta=\min\left\{{1\over2},{C\over2D_0}\right\}>0,
\qquad
q_r=(1-\theta)U_r+\theta P_r.
\tag{AC.23}
```

Convexity of relative entropy gives

```math
D(q_r\|U_r)
\le\theta D(P_r\|U_r)\le Cr.
\tag{AC.24}
```

But `q_r(G_r)>=theta`, and every point of `G_r` has normalized pressure at
least `h_beta+eta`.  Therefore convergence in probability fails.  Also,
using the uniform mean limit and (AC.21),

```math
\liminf {\mathbb E_{q_r}f_r\over r}
\ge h_\beta+\theta\eta,
\tag{AC.25}
```

so `L^1` convergence fails. `square`

This counterexample is high-pressure rather than favorable.  That is enough
to disprove a two-sided conditional-pressure theorem under KL alone.  A
KL-bounded **low-pressure** counterexample would require the still-unknown
speed-`r` favorable basin and is not asserted.

## 5. Consequence for the conference basin search

At every fixed entropy speed `C`, sufficiently small beta makes all events
of that speed thermodynamically invisible.  Therefore a speed-`r` favorable
conference basin, if it exists at a given beta, must undergo a temperature
threshold: its entropy coefficient must exceed the operator-tail exponent
available at smaller beta, or it must disappear as beta decreases.

The theorem is stronger than the preceding affine, template, gauge, and
row-magnitude no-go results in its small-beta range because it assumes no
structure whatsoever.  Those structured results remain useful at larger
beta or with sharper constants.  The exact surviving question is no longer
whether an exotic event of a *fixed* speed can beat arbitrarily small beta;
it is how the minimal favorable entropy speed grows as beta tends to zero,
and whether it ever falls below the finite-tilt reward at a fixed positive
temperature.
