# A quantitative ascent lemma for hierarchical Gaussian masks

Date: 2026-09-05. This note builds on the proved Gaussianization isometry
`𝒰` in `fresh_limit_hierarchical_tree_energy_2026_09_05.md`. Its numerical
examples are deterministic diagnostics, not interval certificates.

## 1. Exact quadratic minorant

For an even mask `0≤H≤1`, put `W=𝒰H` and
`J(H)=E|W|(1-H)`. Choose another even mask `0≤H1≤1`, and set
`δ=H1-H`, `Δ=𝒰δ`. For every `0≤ε≤1`, convexity gives

\[
 |W+\epsilon\Delta|\ge |W|+\epsilon\operatorname{sign}(W)\Delta.
\]

Multiplying by the nonnegative `1-H-εδ` and taking expectations proves
the exact minorant

\[
 \boxed{J(H+\epsilon\delta)\ge
 J(H)+\epsilon\gamma-\epsilon^2 C,}
 \qquad
 \begin{cases}
 \gamma=E\delta[\mathcal U^*F-|W|],\\
 C=E\operatorname{sign}(W)\Delta\delta,\\
 F=\operatorname{sign}(W)(1-H).
 \end{cases}                                                   \tag{1}
\]

Here `𝒰*` first projects onto Gaussian first chaos and then maps each
`Z_T` to `h_T`; thus it is a bounded adjoint operator on `L²`.
The isometry also gives `|C|≤||δ||₂²`, so

\[
 J(H+\epsilon\delta)\ge J(H)+\epsilon\gamma
                                     -\epsilon^2\|\delta\|_2^2. \tag{2}
\]

Taking `H1=1{𝒰*F>|W|}` maximizes `γ` over all masks. If `γ>0`, (2)
already gives a strict, quantitatively bounded improvement. The full
minorant (1) can be much sharper, and it needs no differentiation through
the absolute-value corner.

## 2. Closed two-field calculation at a scalar fixed point

Let `V=𝒰g(V)` be a unit Gaussian fixed point, let
`H0=1{|V|≤α}`, and write `W=𝒰H0`. Thus

\[
 p=EH_0,\quad w=E[g(V)H_0],\quad s^2=p-w^2,
 \qquad \operatorname{Cov}(V,W)=w.
\]

The first-chaos projection of `F0=sign(W)(1-H0)` has the exact form
`aV+bW`, where

\[
 a=2\phi(\alpha)[2\Phi(w\alpha/s)-1],\qquad
 b={4\phi(0)\over\sqrt p}\Phi(-\alpha\sqrt p/s).       \tag{3}
\]

Indeed Gaussian regression gives
`a=(p E VF0-wJ(H0))/s²`, `b=(J(H0)-w E VF0)/s²`;
the elementary folded-normal integrals simplify to (3). Since
`𝒰*V=g(V)` and `𝒰*W=H0`, the gradient mask is

\[
 H_1=1\{|W|<T(V)\},\qquad
 T(v)=\max\{0,a g(v)+b1_{\{|v|\le\alpha\}}\}.         \tag{4}
\]

More generally any fixed measurable threshold `T(v)≥0` gives a legitimate
test mask; no optimality of (4) for the finite step is asserted.

Every quantity in the quadratic minorant reduces to a one-dimensional
Gaussian integral. Put

\[
 p_1=EH_1,\quad r_1=E[g(V)H_1],\quad r_0=E[H_0H_1],
 \qquad
 \binom{c}{d}=
 \begin{pmatrix}1&w\\w&p\end{pmatrix}^{-1}
 \binom{r_1}{r_0}.                                  \tag{5}
\]

The conditional mean of `𝒰H1` given `(V,W)` is `cV+dW`; its remaining
Gaussian component is independent of both. Consequently

\[
 \begin{split}
 \gamma&=a(r_1-w)+b(r_0-p)-B,\\
 C&=cA+(d-1)B,\\
 A&=E[V\operatorname{sign}(W)(H_1-H_0)],\\
 B&=E[|W|(H_1-H_0)].                                \tag{6}
 \end{split}
\]

For fixed `V=v`, set `m=wv`, `u=(T(v)-m)/s`, and
`l=(-T(v)-m)/s`. Conditional on `v`, the three necessary moments are

\[
 \begin{split}
 P(H_1=1\mid v)&=\Phi(u)-\Phi(l),\\
 E[\operatorname{sign}(W)H_1\mid v]
   &=\Phi(u)+\Phi(l)-2\Phi(-m/s),\\
 E[|W|H_1\mid v]
   &=m[\Phi(u)+\Phi(l)-2\Phi(-m/s)]
     +s[2\phi(m/s)-\phi(u)-\phi(l)].                 \tag{7}
 \end{split}
\]

Thus (1), (5)--(7) furnish a concrete scalar-integration lower certificate
for a two-field mask. Direct evaluation of `J((1-ε)H0+εH1)` is also
possible using the conditional folded-normal law of
`(1-ε)W+ε𝒰H1`, but is not required to prove an improvement.

## 3. Numerical diagnostic from the certified 0.426 seed

Use the finite degree-200 `g` with `α=37/50`, resolvent `a0=97/10`
from the hierarchical exact certificate. The deterministic computation
`computations/fresh_limit_fixed_point_mask_ascent.py` gives

\[
 \begin{array}{c|r}
 a&0.5962057641940863\\
 b&0.013856407262371098\\
 p_1&0.5382529964548735\\
 r_1&0.6914656684113833\\
 r_0&0.5108679108177415\\
 \gamma&0.007557744123072973\\
 C&0.007355512066963357
 \end{array}
\]

The quadratic minorant is numerically `0.42803173883655893` at
`ε≈0.5137469733`; even `ε=1/2` gives approximately `0.42803035`.
The direct Gaussian integral is numerically `0.42808794476063183` at
`ε≈0.5177838366`. Neither decimal is currently promoted to a rigorous
numerical lower bound: its scalar quadrature has not been enclosed by
outward intervals. The exact lemma and integral formulas above are proved.

This refinement remains inside the full one-mask certificate, whose
numerical ceiling is approximately `0.449555`; it does not address the
original convergence question.

## 4. A larger rational-threshold step

A two-parameter threshold search suggested the following completely fixed
mask, for the same certified degree-200 scalar fixed point:

\[
 \widetilde H(v,w)=1\left\{|v|\le1,\quad
 |w|<\max\left(0,{5\over2}a g(v)+bH_0(v)-{4\over5}\right)\right\}.
                                                               \tag{8}
\]

The imposed support `|v|≤1` is part of the exact definition. It removes
all tail integrals involving the unbounded polynomial `g`. The old mask
moments remain their exact folded-normal closed forms.

At the full step `ε=1`, the quadratic minorant also has a simpler direct
description. Keep the sign of the old Gaussian field `W` in the final
response, using `F=sign(W)(1-\widetilde H)`. If the regression coefficients
`c,d` are defined by (5) with `H1=\widetilde H`, then

\[
 \boxed{\quad
 J(\widetilde H)\ge
 c\left(\sqrt{2/\pi}{w\over\sqrt p}
                       -E[V\operatorname{sign}(W)\widetilde H]\right)
 +d\left(\sqrt{2p/\pi}-E[|W|\widetilde H]\right).
 \quad}                                                       \tag{9}
\]

Indeed `E[𝒰\widetilde H|V,W]=cV+dW`, and
`|𝒰\widetilde H|≥sign(W)𝒰\widetilde H` pointwise. This is precisely (1)
at `ε=1` after cancellation, but it has fewer repeated interval quantities.
It only needs four of the five scalar moments; the mass `p1` may still be
computed as an independent diagnostic.

Numerically, the untruncated version of (8) gives a minorant
`0.4297871946133101`; the directly optimized final sign gives approximately
`0.42983892584`, agreeing at 96 and 192 inner Gauss--Legendre nodes.
The exact-support version (8) is the one submitted to rigorous enclosure.

## 5. Rigorous one-dimensional enclosure procedure

The exact-arithmetic script
`computations/fresh_limit_mask_ascent_certificate.py` implements (9),
with a target above `0.429`. Its run is separate from the numerical search.
All arithmetic is outward Fraction interval arithmetic; its steps are:

1. Reconstruct the finite Hermite coefficients and check `Eg'^2<1` with
   outward intervals. Compute `p,w,s,a,b` from their closed forms.
2. Split `[0,1]` exactly at `α=37/50`. All integrals for the new mask are
   twice the positive-half integrals, by joint parity.
3. At each rational midpoint, evaluate `g,g',g''` in the normalized
   Hermite basis. Uniform Taylor enclosures for a bin of radius `h` use
   `|g'''|<2300` on `[0,1]`.
4. If the whole threshold interval is negative, the bin contributes zero.
   If it straddles zero, enclose it directly using
   `P(|N(m,s²)|≤T)≤2T φ(0)/s`, together with
   `|E sign(W)H1|≤P(H1)` and `E|W|H1≤T P(H1)`.
5. On every strictly positive-threshold bin, propagate interval derivatives
   through the exact formulas (7). The midpoint error for each integrand
   is at most `sup|f''| width³/24`. Adding these finite error intervals
   encloses the complete integral, including any threshold-crossing bins.

Here is the uniform polynomial derivative bound used in step 3. With
normalized Hermites and coefficient derivative energy below one,

\[
 \|g'''\|_2^2
 =\sum_r r(r-1)(r-2)u_r^2
 \le199\cdot198\sum_r r u_r^2<39402.
\]

For `t=99/100` and `|v|≤1`, Mehler's identity gives

\[
 \sum_{r=0}^{197}{\operatorname{He}_r(v)^2\over r!}
 \le {t^{-197}\over\sqrt{1-t^2}}
                \exp\!\left({tv^2\over1+t}\right)<128.
\]

The strict bound uses the elementary inequalities `t^-197<8`,
`(1-t²)^-1/2<8`, and `exp(t/(1+t))<2`; the first two are checked
rationally by the script. Cauchy--Schwarz now gives
`|g'''(v)|²<39402*128<2300²`.

The script's extended density uses range reduction to the already audited
exponential series. Its point CDF evaluations integrate the alternating
Taylor series through degree 128; for arguments at most 8 the first omitted
tail is less than `10^-24`, checked by an exact rational inequality. Coarse
CDF ranges `[0,1]` suffice when bounding second derivatives. An exact
interval square, with lower endpoint zero on a zero-straddling input,
keeps all density evaluations inside the exponential routine's domain.

The mathematical enclosure method and all formulas have passed an
independent variational-agent audit. The exact run completed successfully,
certifying 1,842 positive-threshold bins and one threshold-crossing bin.
Its stdout is saved at
`computations/results/fresh_limit_mask_ascent_certificate.json` and encloses
the lower certificate by

\[
 [0.429786450737628041355434677764144259463671564581969931627991,
  0.429787953860539027972266545724508214067476532943989947373807].
\]

Consequently

\[
 \boxed{\displaystyle
 \liminf_{n\to\infty}{M_n\over n^{3/2}}>0.429.}
                                                               \tag{10}
\]

This assertion uses the exact interval run, not the floating-point search.
The convergence of the original normalized minima remains unresolved.

The director independently read the complete new interval program and
reconstructed the midpoint, threshold-crossing, Hermite derivative, and
regression formulas. A separate rerun completed at 21:51 UTC and reproduced
the displayed enclosure exactly. The earlier interval-square correction
is present in the committed program.
