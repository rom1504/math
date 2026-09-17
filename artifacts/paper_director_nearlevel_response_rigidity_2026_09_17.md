# Forced isotropic response in every full near-extreme code

2026-09-17. **Verified:** director derivation, independently reconstructed
by all three researchers, with two alternative proofs of the mean bound.
This strengthens the discrepancy track's forced Hamming-sphere entropy
obstruction. It applies to every bounded positive-cap full signing,
not just a surrogate or a particular construction family.

## 1. Statement and quantifiers

Let A_n be hollow full signings with Q(A_n)/n^(3/2)->c>0. Define

```math
E_n(\eta)=\{x\in\{\pm1\}^n:
 Q(A_n)-|H_{A_n}(x)|\le\eta n^{3/2}\}.
```

For EACH n let nu_n be ANY probability law on the full sign cube
with E h h^T=I_n. It can depend arbitrarily on A_n and eta. Then,
for every fixed 0<eta<c,

```math
\boxed{
\liminf_{n\to\infty}\max_{x\in E_n(\eta)}
 \frac{\mathbb E_{\nu_n}|h^Tx|}{\sqrt n}
 \ge \sqrt{\frac{2\eta}{\pi c}}.}                 \tag{1}
```

The limit in n precedes eta down to zero. In particular no isotropic
full-sign law can have whole-nearlevel squared normalized response
o(eta), regardless of its support size, tensor dimension, or adaptation
to the actual signing. This says nothing against protecting a smaller
center code and separately controlling its correlated neighborhoods.

Exact isotropy can in fact be weakened to ||E_nu h h^T||op=o(n).
The same asymptotic conclusion holds. Section6 gives a finite formulation
whose covariance-dependent error explains this extension.

## 2. A dimension-free variance lemma

For any such nu, f(z)=E_nu|h^Tz| is convex and 1-Lipschitz in
Euclidean norm: |f(z)-f(z')|<=sqrt(E[h^T(z-z')]^2)=||z-z'||_2.
For ANY product law on the sign cube,

```math
\operatorname{Var} f(Y)\le4.                       \tag{2}
```

Proof: choose a subgradient g(y) of f at every cube vertex, of norm
at most1. Convexity gives f(y)-f(y^(i))<=2y_i g_i(y), and therefore
sum_i[(f(y)-f(y^(i)))_+]^2<=4. Let Y^(i,res) resample coordinate i
independently. Efron--Stein and the symmetry of the resampled pair give

```math
\operatorname{Var}f(Y)
 \le\frac12\sum_i\mathbb E[f(Y)-f(Y^{(i,res)})]^2
 =\sum_i\mathbb E[(f(Y)-f(Y^{(i,res)}))_+]^2
 \le\mathbb E\sum_i[(f(Y)-f(Y^{(i)}))_+]^2\le4.
```

No independence of the coordinates of h is assumed. The product law
is the deliberately chosen spin perturbation Y, not the law nu.

## 3. Perturbing an actual signed ground word

Choose x0 and polarity sigma with sigma H_A(x0)=Q(A). Independently
flip each coordinate with probability rho in (0,1/2), giving Y.
The deficit D(Y)=Q(A)-sigma H_A(Y) is nonnegative for every Y, and

```math
\mathbb E D(Y)=4\rho(1-\rho)Q(A).                 \tag{3}
```

Thus if 4c rho(1-rho)<eta, a fixed positive fraction of Y lies in
E_n(eta), by Markov. We will show that almost every Y also has
f(Y)/sqrt(n)>=2sqrt(2/pi)sqrt(rho(1-rho))-o(1).

For fixed h the sum h^TY has independent summands of absolute value1,
variance exactly v_rho n, v_rho=4rho(1-rho), and mean
(1-2rho)h^Tx0. Uniformly over h and x0,

```math
\liminf_{n\to\infty}\inf_{h,x0}
 \frac{\mathbb E_Y|h^TY|}{\sqrt n}
 \ge\sqrt{\frac{2v_\rho}{\pi}}.                  \tag{4}
```

For completeness, a sequential compactness argument verifies uniformity.
If the normalized mean is unbounded along a proposed violating sequence,
Jensen contradicts a bounded upper value. Otherwise take a subsequence
whose normalized mean tends to b. The centered, normalized triangular
array has total variance v_rho and summands bounded by2/sqrt(n), so
the elementary Lindeberg central limit theorem applies. Its second
moments are bounded, giving uniform integrability of the absolute value.
The limit is E|b+sqrt(v_rho)G|, minimized at b=0 by symmetry/convexity.
This contradicts any fixed violation of (4).

An independent, fully explicit proof avoids the CLT. For any pattern of
coefficients +/-1, the characteristic function of the biased sum S has
absolute value (1-v_rho sin^2(t))^(n/2). The identity
E|S|=(2/pi)integral_0^infinity (1-Re E exp(itS))/t^2 dt gives a
pattern-independent lower bound after replacing Re by absolute value.
Substitute t=u/sqrt(n) and apply Fatou to its nonnegative integrand.
The limiting integral is the absolute mean of N(0,v_rho), proving(4).
The other independent reconstruction uses smooth-absolute-value
Lindeberg replacement with O(n^(1/3)) unnormalized error, uniformly
in the mean shift and coefficient pattern.

An explicit finite version further simplifies uniformity. Put
b_n(rho)=2 E|B-B'| for independent Binomial(floor(n/2),rho) variables.
The characteristic-function modulus above is bounded above by
(1-v_rho sin^2(t))^floor(n/2), the REAL characteristic function of
sum_(i<=floor(n/2))(xi_i-xi_i'). Hence E|h^TY|>=b_n(rho)
for every coefficient pattern, and b_n(rho)/sqrt(n)->sqrt(2v_rho/pi).
For even n equality is attained by a balanced coefficient pattern.

Average (4) over the arbitrary law nu. Equation(2) and Chebyshev imply
that f(Y)/sqrt(n) differs from its mean by o(1) in probability. Hence
its high-response event intersects the fixed-positive-probability event
from(3). We obtain the lower bound sqrt(2v_rho/pi) inside E_n(eta).
Finally increase rho to the solution of4c rho(1-rho)=eta, still taking
n to infinity FIRST. This proves(1).

## 4. Matrix-response consequence and corrected applicability

For an isotropic law nu on {+-1}^k, independent fair g in{+-1}^p,
h tensor g is an isotropic physical sign law on n=kp. For every array X,

```math
\mathbb E_{h,g}|h^TXg|
 \le\mathbb E_h\|h^TX\|_2.
```

Consequently(1) applies as a LOWER bound to the supremum of the
matrix response Phi_nu on the whole actual nearlevel code, for every
aspect ratio. There is no possible whole-code hypothesis
Phi_nu(eta)^2=o(eta). In particular the earlier sufficient profile
Phi_nu(eta)^2 log(1/Phi_nu(eta)^2)=O(eta) cannot occur.

Those conditional augmentation inequalities remain valid algebraically,
but their proposed original-problem application is VACUOUS. The correct
remaining possibility is a low-response center set with a nontrivial
anchored increment metric controlling nearby words. That possibility
is not proved for actual optimizing signings and is not implied by(1).

## 5. Scope

The proof uses only an actual maximizing word, positivity of its deficit,
the quadratic flip identity, isotropy, convexity, and an elementary CLT.
It does not assume spectral flatness or a conference construction.
It proves neither convergence nor nonconvergence. It removes a genuinely
false route hypothesis and explains why a successful response mechanism
must keep local correlations rather than charge raw nearlevel entropy.

## 6. Finite and nonisotropic form

Let nu be ANY law on full signs, L=||E h h^T||op. For any actual
signing A, tolerance t>0 and rho in(0,1/2), put
a=1-4rho(1-rho)Q(A)/t. If a>0, then

```math
\boxed{
\max_{x:Q(A)-|H_A(x)|\le t}\mathbb E_\nu|h^Tx|
 \ge b_n(\rho)-2\sqrt{L/a}.}                     \tag{5}
```

Indeed f is convex sqrt(L)-Lipschitz, so the same Efron--Stein
proof gives Var f(Y)<=4L. The near-energy event has probability at
least a; its intersection with f>=E f-u is nonempty whenever
u>2sqrt(L/a). The explicit absolute-sum lemma gives E f>=b_n.
Let u decrease to its endpoint in this finite maximum to prove(5).
Taking t=eta n^(3/2), then n to infinity under L=o(n), and finally
rho to the endpoint proves the generalization. At isotropy L=1,
the error is dimension free. No covariance reduction is presumed.
