# Independent audit: weighted ports and nonlinear filters

2026-09-16 UTC. This is an independent verification of the mechanism agent's weighted-port/core-tail theorem and the root's `bh_2026_09_16_nonlinear_filter_audit.md`, including the subsequently proposed all-signing hypercontractive extension. No global project steering or git state was changed.

## Verdict

The weighted-port inequality is valid, including its simultaneous seed/mask/core quantifiers under the stated uniform frame hypotheses. Its strict-subhalf relative constant is correctly certified. However, its core-size and variance-concentration hypotheses force `V=o(N/log^4 m)`. Consequently it does not cover a positive-density variance sector or improve a leading `N^(3/2)` cap bound. The mixed core/tail relative estimate is genuine, but the old global frame bound already makes the entire covered bulk sector additively negligible.

The original bounded-operator polynomial-filter argument is valid for arbitrary real or complex polynomials normalized on a fixed real interval. Its degree and order-of-limit statements check. The later hypercontractive extension is stronger: every signing has the same full coefficient-norm obstruction for `2<=d=o(sqrt n)`, without an operator-norm hypothesis. The weighted functional in that regime is asymptotically its linear-filter contribution after division by degree.

A further bounded-operator tail argument, independently audited below, extends the obstruction for the *current polynomial-growth BH inequalities* to **every growing sublinear filter degree**, and even to a sufficiently small fixed linear-degree window. The all-signing statement remains restricted to the smaller `o(sqrt n)` regime.

Neither result improves the quadratic signing bounds or proves convergence. A proposed global linear weighted-BH slope is constrained by H6, but a slope valid only asymptotically in degree, or with an additive intercept, is a different claim and is not excluded by that finite witness alone.

## 1. Deterministic weighted-port theorem

There are `m` fibres, each of length `q`, so `N=mq`. For real frames `T_i` with balanced columns, put `y_i=P_i x_i`, where `P_i=I-J/q`, and let

\[
h_i(j)=T_i(:,j)^Ty_i,\qquad
H_W(x)=\sum_{i<j}w_{ij}h_i(j)h_j(i),\qquad |w_{ij}|\le1.
\]

Symmetric outer signs and symmetric masks are special cases. Assume a common global and restricted bound

\[
\|T_i\|_{\rm op}^2\le M,\qquad
\|T_i(:,H)\|_{\rm op}^2\le K\le M
\]

for every fibre and every subset `H` of at most `s` columns. Set

\[
V_H=\sum_{i\in H}\|y_i\|_2^2,\quad
V_L=\sum_{i\notin H}\|y_i\|_2^2,\quad V=V_H+V_L.
\]

Then

\[
\boxed{|H_W(x)|\le\frac{MV}{2}
-\frac12\left(\sqrt{(M-K)V_H}-\sqrt{KV_L}\right)_+^2.}
\]

### Independent derivation

For arbitrary positive `lambda_i`, weighted Young gives

\[
2|H_W(x)|\le\sum_i\lambda_i^{-1}\sum_j\lambda_j|h_i(j)|^2.
\]

The sum on the right may include nonnegative diagonal terms; in the actual construction self-ports vanish. Choose `lambda_i=t>=1` on the core and 1 elsewhere. For each row,

\[
\sum_j\lambda_j|h_i(j)|^2
=\|T_i^Ty_i\|_2^2+(t-1)\|T_i(:,H)^Ty_i\|_2^2
\le[M+(t-1)K]\|y_i\|_2^2.
\]

The restriction `t>=1` matters: it makes the coefficient of the restricted square nonnegative. Therefore

\[
|H_W(x)|\le\frac12[M+(t-1)K](V_H/t+V_L)
\]

for every such `t`. Expanding this expression leaves the variable terms

\[
\tfrac12[(M-K)V_H/t+KV_Lt].
\]

If `K V_L>0`, the minimizing point is
`t=max(1,sqrt((M-K)V_H/(K V_L)))`. The claimed positive-part formula follows by substitution. Cases `V_L=0`, `K=0`, `K=M`, or `V_H=0` follow directly or by taking `t` to infinity; the infimum need not be attained for the inequality to hold.

The deficit is positive exactly when `V_H/V>K/M`, assuming `V>0`. There is no random-independence assumption involving the core. Since the restricted frame bound is simultaneous over every core of size at most `s`, the core may be chosen adaptively from the physical word or from the frames.

The mechanism agent's subsequent two-level optimality observation also checks. For arbitrary positive weights, let `a` be the largest and `b` the `(s+1)`-st largest. Layer-cake and the step restricted-norm information give row cap `M b+K(a-b)`. Raising all weights above `b` to `a`, and all others to `b`, leaves that cap unchanged while decreasing `sum_i ||y_i||^2/lambda_i`. Filling up the high set to `s` entries and choosing the largest variances further decreases the certificate. Hence two levels are optimal **for this particular layer-cake relaxation**. This is not a claim that the complete joint feasible row-square problem admits no stronger certificate.

## 2. Actual construction inputs and numerical constant

I checked the relevant existing construction statement in `artifacts/principle_synthesis_2026_09_07_balanced_bulk_sparse_active_bound.md`, particularly Sections 6.3–6.5. Its conditioned-selector law and operator-controlled repair give, uniformly in the fibres and all restricted sets of size at most `s`,

\[
M=m+o(m),\qquad K=q+o(m),\qquad q/m\to p,
\qquad s\log^4(m)/m\to0.
\]

The global bound also follows from the raw selected Hadamard norm at most `sqrt(m)`, contraction by `P_i`, and the uniform `o(sqrt m)` operator repair error. The restricted bound does not follow for arbitrary repair locations; the operator-controlled repair and the changed, conditioned selector law are essential hypotheses. This audit uses that prior construction theorem as an input rather than silently transferring its conclusions to an unconditioned law.

Let `theta=V_H/V`. The normalized new coefficient is

\[
\frac{|H_W(x)|}{\sqrt N V}
\le\frac{1-[\sqrt{(1-p)\theta}-\sqrt{p(1-\theta)}]_+^2}{2\sqrt p}+o(1).
\]

For `theta>p`, the bracket is positive and increases with `theta`. Hence the coefficient decreases. At `p=31/32` and `theta>=999/1000`, its maximum is

\[
0.497236615809144\ldots<0.498.
\]

The `o(1)` is uniform on `theta` in `[999/1000,1]`: the expression is continuous in the normalized parameters on a compact set, and the frame errors are uniform. Thus the strict eventual `0.498` conclusion is justified.

I independently reran

```sh
.venv/bin/python computations/bh_mechanism_2026_09_16_weighted_ports.py \
  --output tmp/bh_2026_09_16/independent_weighted_ports_replay.json
```

The rational phase certificate and all model checks passed. The 1,920 weighted checks and 500 one-variable optimization checks use floating-point linear algebra; they are sanity checks of the analytically proved inequality, not a proof of the asymptotic restricted frame hypothesis. The numerical optimization error was at most `4.42e-12`. The strict phase coefficient itself is certified by rational upper/lower square-root bounds. The small actual-frame tests use density `3/4`, not `31/32`, which is acceptable for checking the general deterministic inequality but should not be confused with a finite construction certificate at the target density.

## 3. Why the covered sector is additively negligible

For every Boolean word, `||P_i x_i||_2^2<=q`. Therefore

\[
\theta V=V_H\le |H|q\le sq,
\qquad
\frac VN\le\frac{s}{\theta m}.
\]

If `theta>=999/1000` and `s=o(m/log^4 m)`, then necessarily

\[
V=o(N/\log^4m).
\]

Even the old global bound consequently gives

\[
|H_W|\le MV/2=o(N^{3/2}/\log^4m).
\]

Thus the new relative coefficient does not settle a previously uncontrolled leading-order dense additive-cap contribution. Conversely, if `V>=epsilon N` and the core holds a fixed fraction `theta_0` of the variance, then `|H|>=theta_0 epsilon m`: a positive-density variance sector requires a linearly sized core, outside the available sparse restricted-frame range.

An explicit physical profile illustrates why “all fibres may be active” does not evade this observation. Let `s=floor(sqrt m)`. Make the core words balanced, so `V_H=sq`. In every other fibre flip exactly one coordinate from an otherwise constant word. Then

\[
V_L=4(m-s)(1-1/q),\qquad\theta\to1,
\qquad V/N=m^{-1/2}+O(m^{-1})\to0.
\]

All fibres are nonconstant, and eventually the core-concentration hypothesis holds, yet the entire variance remains negligible at the parent scale. This is a scope diagnostic, not a counterexample to the theorem.

The theorem genuinely pays the mixed core/tail interaction in this relative microscopic sector. A later recursive or defect-weighted application might value that relative estimate, but no such propagation theorem is supplied here. It cannot presently be promoted to a new cap constant, complete-parent estimate, or convergence recurrence.

I also checked the mechanism agent's adjacent actual-frame obstructions, Section 10 of its artifact. The rank/trace floor follows from `rank(T)<=q-1`, averaging principal Frobenius squares, and `(lambda+q)(t-lambda)>=0` for eigenvalues of a zero-trace Gram deviation with largest eigenvalue `t`. Its positive-density loss `(1-p)delta m` is correctly normalized. The allowed Sylvester tensor construction is also valid: a retained affine coset of size `K=Theta(log M)`, equipped with a nontrivial character, has Fourier support exactly `m/K`, excludes DC, and stays centered after selection. Its raw restricted squared norm is exactly `m`; an `o(sqrt m)` repair preserves `m-o(m)`. The independent-coset probability estimate and bounded-density conditioning transfer check. This is correctly limited to a full-frame-uniform obstruction; it does not establish that a typical recursive frame has the same obstruction. The adjacent temperature-mixture result was not independently reconstructed in this audit.

## 4. Original nonlinear-filter audit

Fix constants `c>0` and `L<infinity`. Assume a full signing satisfies `||A_n||op<=L sqrt n`; write `Y_n=H_A/n^(3/2)`. The existing direct decoupling/Gaussian-mgf proof gives

\[
\|Y_n\|_{2j}\le C_L(\sqrt{j/n}+j/n).
\]

For `j<=n`, this is at most `2C_L sqrt(j/n)`. Let `P(t)=sum_(j=0)^d a_jt^j` be any real or complex polynomial bounded by 1 on `[-c,c]`. The coefficient estimate

\[
|a_0|\le1,\qquad |a_j|\le(ed/(cj))^j
\]

is valid. To see this without assuming real coefficients or boundedness on a complex disk, apply the Joukowski maximum-principle construction to `P(ct)`. It gives

\[
|P(z)|\le(|z|/c+\sqrt{1+|z|^2/c^2})^d\le e^{d|z|/c}.
\]

Cauchy's formula on radius `cj/d` yields the estimate. Complex polynomials are covered exactly, not with an omitted factor.

For every fixed `s>=0`, the weighted coefficient seminorm satisfies `W_s(g)<=sqrt(en)||g||_2`, where

\[
W_s(g)^2=\sum_{r\ge1}r^{-2s}\|\widehat g^{=r}\|_{2r/(r+1)}^2.
\]

The constant term is annihilated. Triangle inequality and the preceding estimates give

\[
W_s(P(Y_n))\le\sqrt{en}\sum_{j=1}^d(z/\sqrt j)^j,
\qquad z=C_{c,L}d/\sqrt n.
\]

The series bound `sum_(j>=1)z^j/sqrt(j!)<=sqrt2 z exp(z^2)` follows by factoring out `z` and applying Cauchy–Schwarz to the sequences `2^(-ell/2)` and `(sqrt2 z)^ell/sqrt(ell!)`. Since `j!<=j^j`, it follows that

\[
W_s(P(Y_n))\le C_{c,L}d\exp(C'_{c,L}d^2/n).
\]

The same support argument with `p=4d/(2d+1)` gives

\[
\|\widehat{P(Y_n)}\|_p\le1+\sqrt2z_d e^{z_d^2},
\qquad z_d=C_{c,L}d(en)^{1/(2d)}/\sqrt n.
\]

All signs and cancellations remain inside the reduced Walsh coefficients; no positivity assumption is used.

### Limits and quantifiers

- The constants `c,L` must stay fixed. Allowing `c` to shrink or `L` to grow requires new uniform estimates.
- For every sequence `2<=d=o(sqrt n)`, `z_d->0` uniformly over `P` and the bounded-operator family. Split at `d=log n`: below it, `z_d<=C log(n)n^(-1/4)`; above it, the extra factor `(en)^(1/(2d))` is bounded.
- For growing `d=o(sqrt(n log n))`, one has `d^2/n=o(log d)`. If `d<=n^(1/4)` this is immediate; otherwise `log d` is a fixed fraction of `log n`. Thus the weighted bound is `d^(1+o(1))`, and the full coefficient norm is `d^(o(1))`.
- To obtain the latter statement directly, use `log(1+sqrt2 z exp(z^2))<=C(1+z^2)` and the bounded extra factor above `log n`; below that threshold the coefficient norm tends to 1.
- Degree 1 is excluded from the full coefficient-norm limit and cannot be silently inserted into it. The linear filter has a nonvanishing degree-two coefficient norm.
- If the actual univariate degree is `delta`, with `2delta<=n`, the top Fourier level of `P(Y_n)` is nonzero. Its coefficient on a `2delta`-set is the nonzero leading coefficient times `delta!` times an odd-integer hafnian, with the normalization factor. Hence its actual Walsh degree is `2delta`. For a polynomial of degree *at most* `d`, use its true degree rather than an artificial larger parameter when invoking a degree-sensitive conclusion.

## 5. Stronger extension: all signings, without an operator bound

The later root/power-agent extension also checks. Since `H_A` is homogeneous of degree 2, real hypercontractivity gives, for every full signing,

\[
\|Y_n\|_{2j}\le(2j-1)\|Y_n\|_2
=(2j-1)\sqrt{\binom n2}/n^{3/2}
\le\sqrt2j/\sqrt n.
\]

No spectral assumption appears. Combining this with the coefficient estimate above cancels the factor `j^j`. Put `z=sqrt2 e delta/(c sqrt n)`. For `z<1`, the nonlinear remainder obeys

\[
W_s\left(P(Y_n)-a_0-a_1Y_n\right)
\le\sqrt{en}\frac{z^2}{1-z}.
\]

Therefore, uniformly over all signings and all such normalized polynomials of degree `delta=o(sqrt n)`,

\[
W_s(P(Y_n))
=2^{-s}|P'(0)|\binom n2^{3/4}/n^{3/2}
+O_{c,s}(\delta^2/\sqrt n).
\]

This uses the reverse triangle inequality for a seminorm; it is not a false linearity assertion about `W_s`. The error may grow in absolute terms when `delta` exceeds `n^(1/4)`, but divided by degree it is `o(1)` throughout the stated `o(sqrt n)` regime.

For the full coefficient norm at exponent `4d/(2d+1)`, put `z_d=sqrt2 e d(en)^(1/(2d))/(c sqrt n)`. The same argument now gives the simpler geometric estimate

\[
\boxed{\|\widehat{P(Y_n)}\|_{4d/(2d+1)}\le1+\frac{z_d}{1-z_d}=1+o(1)}
\]

uniformly over **all** full signings and all `2<=d=o(sqrt n)`. This genuinely strengthens the original bounded-operator filter obstruction.

### Linear weighted-BH slopes

Bernstein's derivative bound on the real interval is `|P'(0)|<=delta/c`, also for complex polynomials by phase rotation. Let `w_s=2^(-s-3/4)`. Using true degree and the nonzero top-level fact,

\[
\sup_{1\le\deg P\le d\atop\|P\|_{[-c,c]}\le1}
\frac{W_s(P(Y_n))}{2\deg P}
=\frac{w_s}{2c}+O_{c,s}(d/\sqrt n+1/n).
\]

The lower bound follows from `P(t)=t/c`; the upper bound is the derivative inequality plus the remainder. This supremum includes degree 1. If one instead fixes exactly an even degree, equality in the derivative bound need not hold; odd Chebyshev degrees attain its sharp derivative scale.

Suppose hypothetically that `W_s(f)<=C m||f||infinity` holds with a single global constant `C` for every positive degree bound and every dimension. Then these filters have precisely the direct degree-two threshold `c=w_s/(2C)` at this scale. The H6 witness forces

\[
C\ge\frac{2^{-s}15^{3/4}}{10},\qquad
\frac{w_s}{2C}\le\frac5{30^{3/4}}=0.390057886553453\ldots.
\]

Thus merely proving such a universal linear-in-degree weighted inequality cannot improve the current `0.4333221116640807` lower constant through these sub-square-root-degree scalar filters. This statement does not rule out an asymptotic-in-degree slope with an additive intercept, a different coefficient functional, higher-degree filters, or an optimizer-specific theorem.

## 6. Exact full-sign false-positive family

The root's explicit family is correct. Set `R=J_4-2I_4`, `H=R^(tensor k)`, `n=4^k`, `delta=(-1)^k`, and `A=H-delta I`. Then `R^2=4I`, so `H^2=nI`; the diagonal of `H` is identically `delta`, making `A` hollow and full-sign.

The all-ones vector is a Boolean eigenvector of `R` at eigenvalue 2, and any balanced four-spin vector is an eigenvector at eigenvalue -2. For odd `k`, the all-ones tensor attains eigenvalue `sqrt n+1` of `A`. For even `k`, take one balanced factor and all other factors constant to attain eigenvalue `-sqrt n-1`. Therefore

\[
Q(A)=\frac n2(\sqrt n+1),\qquad Q(A)/n^{3/2}\to1/2.
\]

For every fixed `c<1/2`, the proposed cap is false on this subsequence, but all normalized growing-degree filters in the original `o(sqrt(n log n))` bounded-operator regime still satisfy the paper's necessary coefficient inequalities. In the `o(sqrt n)` regime, the strengthened obstruction applies to every signing, not just this family. No near-minimizer claim or all-order realization is needed for the false-positive conclusion.

## 7. Stronger bounded-operator obstruction for every sublinear degree

The root subsequently supplied a stronger argument. It is valid, with a harmless factor 2 retained in the exponential tail estimate and a fixed small linear-degree cutoff to avoid degree-collapse issues.

Fix `||A_n||op<=L sqrt n` and `c>0`. The deterministic spectral bound gives `|Y_n|<=R=L/2`. The previously audited quadratic tail estimate gives constants `gamma=gamma_(c,L)>0` such that

\[
\Pr(|Y_n|>c)\le2e^{-\gamma n}.
\]

If `c>=R`, the tail event is empty and the argument is simpler. Otherwise all constants below remain fixed with `n`. For every polynomial `P` of degree `d`, bounded by 1 on `[-c,c]`, Joukowski coefficient growth gives the uniform range bound

\[
\sup_{|t|\le R}|P(t)|\le e^{Kd},
\qquad K=\operatorname{arsinh}(R/c).
\]

Split the expectation into `|Y_n|<=c` and its complement:

\[
\|P(Y_n)\|_2^2\le1+2e^{2Kd-\gamma n}.
\]

Choose

\[
0<\eta\le\min\{1/4,\gamma/(4K)\}.
\]

Then for every `d<=eta n`,

\[
\|P(Y_n)\|_2^2\le1+2e^{-\gamma n/2},
\qquad
\|P(Y_n)\|_2\le1+e^{-\gamma n/2}.
\]

The factor 2 in the squared-norm bound should not simply be dropped. The displayed unsquared bound follows from `sqrt(1+2u)<=1+u`.

For the weighted norm this immediately gives `W_s(P(Y_n))<=C_s sqrt n`. For the full coefficient norm at `q_(2d)=4d/(2d+1)`, the support estimate, valid since `2d<=n/2`, gives

\[
\|\widehat{P(Y_n)}\|_{q_{2d}}
\le\left(\sum_{r=0}^{2d}\binom nr\right)^{1/(4d)}\|P(Y_n)\|_2
\le\sqrt{en/(2d)}\,(1+e^{-\gamma n/2}).
\]

The standard binomial partial-sum bound follows, for example, by bounding the generating function `(1+t)^n` at `t=2d/(n-2d)<=1` and using the binary entropy bound; its use here includes all subsets and hence also the smaller even Fourier support. If a larger fixed linear-degree window is desired, `D<=2^n` instead gives a bounded factor when `d>=n/4`. The conservative `eta<=1/4` version already contains every sublinear degree sequence and keeps actual Walsh degree equal to `2d` for polynomials of true degree `d`.

### Consequence for weighted polynomial growth

Fix any exponent `a>1` and choose `beta` with

\[
1/(2a)<\beta<1/2.
\]

For `d<=n^beta`, the earlier low-degree estimate gives uniformly
`W_s(P(Y_n))<=C d exp(C d^2/n)=O(d)`. For `d>n^beta`, the new `L2` estimate gives `W_s(P(Y_n))=O(sqrt n)`. Therefore, for **every growing degree sequence** `d<=eta n`,

\[
\boxed{\sup_P\frac{W_s(P(Y_n))}{d^a}\longrightarrow0.}
\]

In the low branch this is `O(d^(1-a))`; in the high branch it is at most `O(n^(1/2-a beta))`. Both vanish. The supremum is uniform over every bounded-operator signing in the fixed-`L` class and every normalized complex polynomial.

### Consequence for full coefficient polynomial growth

Fix any exponent `b>1/2`, and choose

\[
1/(2b+1)<\beta<1/2.
\]

For `2<=d<=n^beta`, the full coefficient norm is `1+o(1)` uniformly by the previously audited low-degree bound. For `d>n^beta`, it is `O(sqrt(n/d))`. Thus for every growing `d<=eta n`,

\[
\boxed{\sup_P\frac{\|\widehat{P(Y_n)}\|_{q_{2d}}}{d^b}\longrightarrow0.}
\]

The high branch is bounded by `O(n^(1/2-beta(b+1/2)))`, which vanishes by the choice of `beta`. The low branch vanishes because `d` grows.

These conclusions cover the paper's weighted/full exponents 22 and 27 and the independently refined exponents 3 and `11/2`. Consequently **no growing sublinear scalar polynomial filter** can violate those necessary inequalities on the explicit bounded-operator false-positive family in Section 6. The result is stronger than merely improving a rate inside the earlier square-root-degree regime.

The thresholds `a=1` and `b=1/2` are not covered: the required interval for `beta` then closes. Bounded filter degrees are also not covered by this ratio-to-zero statement, and `c,L` must remain fixed. Section 7 alone is not a theorem about the operator norm of actual cap minimizers. Section 9 below subsequently extends the filter obstruction to actual bounded-cap families, with a smaller degree window. Neither result obstructs model-specific coefficient information beyond the universal polynomial-growth inequalities.

## 8. What is and is not certified

Certified: the deterministic port inequality, its rational numerical phase, the negligibility implication, arbitrary-complex-polynomial coefficient extraction, the two nonlinear-filter regimes with their exact limit orders, the all-signing hypercontractive extension, the global-slope H6 obstruction, and the bounded-cap extension in Section 9.

Not certified or claimed: a positive-density variance phase, a new dense-sign construction constant, an original-problem recurrence, convergence of normalized minima, a uniform theorem at degree comparable to `sqrt n` or larger for **all** signings without a cap or operator hypothesis, or an asymptotic-only weighted-BH slope. The stronger sublinear-degree result in Section 7 has a bounded-operator hypothesis. The old restricted-frame construction is used with its conditioned selector and operator-repair hypotheses intact.

## 9. Bounded-cap families, including actual near-minimizers

This section independently verifies a later root derivation. Fix `C>0`, and consider **every** hollow symmetric full signing satisfying

\[
Q(A)=\max_{x\in\{-1,1\}^n}\left|\frac12x^TAx\right|\le Cn^{3/2}.
\]

This hypothesis includes actual minimizers, or any uniformly bounded normalized-cap family. It does not assume a spectral construction or a conjecture on the law of minimizers.

### 9.1 Elementary cap-to-spectral interpolation

Write

\[
B(A)=\max_{u,v\in\{-1,1\}^n}|u^TAv|.
\]

Multilinearity makes this also the supremum over two real unit cubes. Likewise, the hollow quadratic has the same maximum absolute value on the real unit cube as on its vertices. With `s=(u+v)/2` and `t=(u-v)/2`, polarization gives

\[
u^TAv=2\{Q_A(s)-Q_A(t)\},\qquad B(A)\le4Q(A).
\]

For complex vectors `z,w` with coordinate moduli at most one, rotate `z` by a scalar phase so that `z^TAw` is nonnegative real. Its real part is
`z_R^TAw_R-z_I^TAw_I`; both real bilinear terms have magnitude at most `B(A)`. Thus the complex bilinear unit-cube norm is at most `2B(A)`.

Here is a direct three-lines reconstruction, avoiding any ambiguity in real-versus-complex interpolation constants. For real Euclidean unit vectors `u,v`, define, on `0<=Re z<=1`,

\[
F(z)=\sum_{i,j}a_{ij}\operatorname{sgn}(u_j)
 \operatorname{sgn}(v_i)|u_j|^{2z}|v_i|^{2z}.
\]

Zero-coordinate terms are identically zero. On the line `Re z=0`, the two vectors have coordinate moduli at most one, so `|F(z)|<=2B(A)`. On `Re z=1`, each vector has absolute coordinate sum one, so `|F(z)|<=max|a_ij|=1`. Three-lines at `z=1/2` yields

\[
|v^TAu|\le\sqrt{2B(A)},\qquad
\boxed{\|A\|_{\rm op}\le\sqrt{2B(A)}\le\sqrt{8Q(A)}
 \le\sqrt{8C}\,n^{3/4}.}
\]

The strip function is a finite sum of bounded exponentials on this strip, so the usual bounded-strip hypothesis is automatic. The factor two in the complexification step is sufficient; an exact real interpolation constant is not being asserted.

### 9.2 Tail control and arbitrary normalized polynomial filters

Let `Y=Q_A(x)/n^(3/2)` for uniform independent spins, fix `c>0`, and let `P` be a real or complex polynomial of true degree `d`, normalized by `sup_[-c,c]|P|<=1`. The previously reconstructed elementary quadratic tail estimate gives an absolute `c_0>0` such that

\[
\Pr(|Y|>t)\le2\exp\left[-c_0\min\left(nt^2,
 \frac{n^{3/2}t}{\|A\|_{\rm op}}\right)\right].
\]

For every `t>=c` and `n>=1`, the new spectral estimate implies

\[
\Pr(|Y|>t)\le2e^{-a_n t},\qquad
a_n=\gamma_{C,c}n^{3/4},\quad
\gamma_{C,c}=c_0\min\{c,(8C)^{-1/2}\}>0.
\]

Indeed `nt^2>=c n^(3/4)t`. The Joukowski bound already proved above gives `|P(t)|<=exp(d|t|/c)` on the real line. With `b=2d/c<a_n`, integration of the tail, including its boundary mass, gives

\[
\begin{aligned}
\|P(Y)\|_2^2
&\le1+e^{bc}\Pr(|Y|>c)
       +\int_c^\infty b e^{bt}\Pr(|Y|>t)\,dt\\
&\le1+\frac{2a_n}{a_n-b}e^{-(a_n-b)c}.
\end{aligned}
\]

Consequently, for `d<=a_n c/4`,

\[
\|P(Y)\|_2^2\le1+4e^{-a_n c/2},\qquad
\|P(Y)\|_2\le1+2e^{-a_n c/2}.
\]

There is an even shorter proof under the cap hypothesis: `|Y|<=C` deterministically, and the Joukowski bound on `[-C,C]` is `exp(K_{c,C}d)`, with `K_{c,C}=arsinh(C/c)`. Splitting at `c` gives
`||P(Y)||_2^2<=1+2 exp(2K_{c,C}d-gamma_{C,c} c n^(3/4))`, hence the same degree scale. The tail-integration form is retained because it transparently uses only the spectral tail and has no hidden support assumption.

### 9.3 Exact uniform no-go window

There exists `eta=eta(C,c)>0` such that the preceding `L2` bounds hold uniformly for every `d<=eta n^(3/4)`. For all sufficiently large `n`, `2d<=n/2`; thus both the support bound

\[
\|\widehat{P(Y)}\|_{q_{2d}}\le\sqrt{en/(2d)}\,\|P(Y)\|_2
\]

and actual Walsh degree `2d` apply. The latter follows from the odd number of signed perfect matchings in each top `2d`-set: the top coefficient is a nonzero scalar times an odd integer and cannot cancel.

The low-degree inputs do **not** need a spectral hypothesis: Section 5's all-signing hypercontractivity gives `W_s(P(Y))=O_c(d)` for `d<=n^beta`, any fixed `beta<1/2`, and full coefficient norm `1+o(1)` uniformly for `2<=d<=n^beta`. The weighted statement follows from a linear contribution `O_c(d)` and nonlinear remainder `O_c(d^2/sqrt n)`.

Repeating the two-cutoff proof in Section 7 therefore yields, for every growing sequence `d=d_n<=eta n^(3/4)`, uniformly in all bounded-cap signings and all normalized complex polynomials of degree `d`,

\[
\boxed{\frac{W_s(P(Y))}{d^a}\longrightarrow0\quad(a>1),
\qquad
\frac{\|\widehat{P(Y)}\|_{q_{2d}}}{d^b}\longrightarrow0
 \quad(b>1/2).}
\]

Use `beta in (1/(2a),1/2)` for the first statement and `beta in (1/(2b+1),1/2)` for the second. These intervals are nonempty exactly at the strict exponent thresholds stated. On the high branch the bounds are respectively `O(sqrt n)` and `O(sqrt(n/d))`; on the low branch they are `O(d)` and `1+o(1)`.

In particular, actual near-optimal cap families pass the paper's polynomial weighted/full BH tests throughout this growing-degree window. This is a theorem about the weakness of these universal norm tests, not a theorem about all higher moments of minimizers, not a claim at the critical exponents `a=1,b=1/2`, and not a no-go for every scalar filter beyond degree `n^(3/4)` or for additional structure-sensitive coefficient inequalities.

## 10. Uniform random-restriction RMS tests obey the same barriers

This independently verifies the root's subsequent conditional-Walsh extension. Fix a coordinate set `I` of size `t>=2d`, with complement `J`, and freeze `y in {-1,1}^J` uniformly. The set `I` may depend arbitrarily on the coefficient matrix, but not on the sampled frozen spin values. Define

\[
f_y(x_I)=P(Y(x_I,y)),\qquad
\mathcal R_s(P;I)=\big(\mathbb E_y W_s(f_y)^2\big)^{1/2},
\]

and, at `q=q_(2d)`,

\[
\mathcal R_q(P;I)=
 \big(\mathbb E_y\|\widehat f_y\|_q^2\big)^{1/2}.
\]

The first is a seminorm in the underlying global function, and the second is a norm. Fiberwise coefficient-to-`L2` comparison followed by Fubini gives

\[
\mathcal R_s(P;I)\le\sqrt{et}\,\|P(Y)\|_2,
\qquad
\mathcal R_q(P;I)\le D(t,2d)^{1/(4d)}\|P(Y)\|_2,
\]

where `D(t,k)=sum_(r=0)^k binom(t,r)`. For every `1<=k<=t`, not only `k<=t/2`, the elementary estimate

\[
D(t,k)\le(et/k)^k
\]

follows by putting `u=k/t<=1` in `D(t,k)<=u^{-k}(1+u)^t<=u^{-k}e^{tu}`. Thus the second bound is at most `sqrt(et/(2d)) ||P(Y)||_2`. The even-support improvement used for an unfrozen quadratic is not needed: frozen fibers usually contain odd levels.

### The linear term remains controlled after freezing

Write the normalized quadratic on a fiber as a constant, a linear field, and its full quadratic restriction. With

\[
h_i(y)=\sum_{j\in J}a_{ij}y_j,
\]

the exact weighted identity is

\[
W_s(Y_y)^2=
\frac{\big(\sum_{i\in I}|h_i(y)|\big)^2}{n^3}
 +\frac{2^{-2s}\binom t2^{3/2}}{n^3}.
\]

Constants do not enter `W_s`. Cauchy--Schwarz and the marginal second moments give

\[
\mathbb E_y\big(\sum_{i\in I}|h_i(y)|\big)^2
\le t\sum_{i\in I}\mathbb E_y h_i(y)^2=t^2(n-t).
\]

Hence `mathcal R_s(Y;I)<=1` for all `s>=0`, uniformly over full signings and sets `I`. Correlations among distinct field coordinates do not affect this calculation.

For `P=sum a_j z^j`, the previously proved coefficient and hypercontractive estimates imply

\[
\begin{aligned}
\mathcal R_s(P-a_0-a_1Y;I)
&\le\sqrt{en}\sum_{j=2}^d|a_j|\,\|Y\|_{2j}^j\\
&\le\sqrt{en}\frac{z^2}{1-z}
=O_c(d^2/\sqrt n),
\quad z=\sqrt2 ed/(c\sqrt n)<1.
\end{aligned}
\]

Since Bernstein gives `|a_1|<=d/c`, the full low-degree RMS weighted bound is `O_c(d)` uniformly for `d<=n^beta`, every fixed `beta<1/2`.

For the full coefficient norm, the `j`th Taylor term has fiber support at most `D(t,2j)`. Its RMS norm at the common exponent `q_(2d)` is at most

\[
D(t,2j)^{1/(4d)}\|Y\|_{2j}^j
\le(en)^{j/(2d)}\|Y\|_{2j}^j.
\]

The same geometric series as in the unrestricted calculation therefore gives `mathcal R_q(P;I)<=1+o(1)` uniformly for `2<=d=o(sqrt n)`, for all signings and all eligible `I`.

### Consequences and exclusions

On each fiber, the actual Walsh degree is exactly `2d`: on a `2d`-set inside `I`, its top coefficient is `a_d d!` times the signed hafnian of the full signing on that set, divided by `n^(3d/2)`. An odd number of products of signs is a nonzero odd integer. Lower Taylor terms, frozen constants, and frozen fields cannot contribute at degree `2d`.

The high-degree `L2` bounds in Sections 7 and 9 followed by the fiber estimates above now give the identical cutoff conclusions for `mathcal R_s` and `mathcal R_q`. For every growing eligible degree sequence, their ratios to `d^a` and `d^b` tend to zero, respectively for `a>1` and `b>1/2`, throughout the fixed-operator window `d<=eta n` or bounded-cap window `d<=eta n^(3/4)`.

These statements are uniform in `I`: they allow an arbitrary coefficient-dependent choice of `I`, a supremum over `I` **after** taking the frozen-spin RMS, and any independent mixture of such partitions. They do not control the maximum over frozen assignments, nor a choice of `I` depending on those sampled assignments. The latter changes the conditional measure and invalidates the simple Fubini step. Thus uniform RMS random restrictions do not rescue these coefficient tests merely by producing conditional odd levels; distribution-sensitive or extremal conditional tests remain outside the result.

## 11. Final actual-class power refinement

The complete frozen Section 4.4 of `artifacts/bh_2026_09_16_power_analysis.md` passes independent audit. Fix `C>0`; every full signing with `Q(A)<=C n^(3/2)` satisfies the cap-to-spectrum estimate of Section 9. Combining that estimate with the previously reconstructed quadratic tail bound, and integrating its Gaussian and exponential terms, gives

\[
\frac{\|H_A\|_{2k}}{n^{3/2}}
\le C_C\left(\sqrt{k/n}+k/n^{3/4}\right)
\qquad(k\ge1).
\]

The two terms arise respectively from `sqrt(k) ||A||_F` and `k ||A||op`; the full-sign Frobenius scale is `O(n)`. This is a uniform upper bound on every bounded-cap signing, not an existential average-signing or spectral-construction argument.

Let `F_(n,k)=||hat(H_A^k)||_(q_(2k))^(1/k)`. The support and weighted comparisons give

\[
\max\{F_{n,k},W_s(H_A^k)^{1/k}\}
\le(en)^{1/(2k)}\|H_A\|_{2k},\qquad s\ge0.
\]

There is no hidden constant depending on `s`: `W_s<=W_0`, and `W_0(g)<=sqrt(en)||g||_2`. For the coefficient norm, count at most `(n+1)^(2k)` supports and use finite-dimensional `lq` versus `l2` comparison.

The uniform maximization is elementary but important. For `alpha>0,b>1,K>=2`,

\[
\max_{2\le k\le K}k^\alpha b^{1/(2k)}
\le\max\{2^\alpha b^{1/4},e^\alpha K^\alpha\}.
\]

The logarithmic derivative has sign `alpha k-(log b)/2`. If `K` is below the turning point, the maximum is at two. Otherwise the only other possible endpoint is `K`, where `b^(1/(2K))<=e^alpha`. This also covers a turning point below two. Using `alpha=1/2,1` and `b=en` gives, for every integer `K>=2`,

\[
\boxed{\frac1{n^{3/2}}\max_{2\le k\le K}
 \max\{F_{n,k}(A),W_s(H_A^k)^{1/k}\}
\le C_C\{n^{-1/4}+\sqrt{K/n}+K/n^{3/4}\}.}
\]

The lower-endpoint contributions are `O(n^(-1/4))` and `O(n^(-1/2))`; the latter is absorbed by the former. Consequently all powers `2<=k<=K_n=o(n^(3/4))` vanish after normalization by `n^(3/2)`, uniformly over every bounded-cap signing. This includes every choice of exact minimizers, since elementary constructions already bound their normalized caps. It explicitly excludes `k=1`, whose coefficient norm has a different scale.

Degree collapse does not evade the argument. For any valid positive degree bound `m<=2k`, including the actual nonzero reduced degree, the support has size `D<=(n+1)^m` and

\[
\|\widehat{H_A^k}\|_{q_m}^{1/k}
\le D^{1/(2mk)}\|H_A\|_{2k}
\le(n+1)^{1/(2k)}\|H_A\|_{2k}.
\]

A constant reduced power may be assigned the valid bound `m=1`. Thus the same theorem holds for the actual-degree coefficient variant as well.

The already audited planted-clique family shows the exponent is sharp for the broader bounded-cap class, not for minimizers. Along `n=16^t`, its tail satisfies `P(H_A/n^(3/2)>=a)>=c_a 2^(-n^(3/4))` for fixed `0<a<1/2`. For `k_n=floor(tau n^(3/4))`, fixed `tau>0`,

\[
\liminf_{n=16^t\to\infty}\frac{F_{n,k_n}(A)}{n^{3/2}}
\ge\liminf_{n=16^t\to\infty}\frac{\|H_A\|_{2k_n}}{n^{3/2}}
\ge a\,2^{-1/(2\tau)}>0.
\]

The first comparison is simply `lq>=l2`; the tail prefactor disappears after the `2k_n`th root. Thus the uniform little-o conclusion cannot be extended to all bounded-cap signings at power order `Theta(n^(3/4))`. No original minimax bound or convergence statement follows.

## 12. Exact large-core obstruction after centering a selected Hadamard frame

A final bounded algebra check requested by the mechanism agent also passes. Let `F` consist of `q` rows of a normalized `m`-column Hadamard matrix, so `FF^*=mI_q` and its DC column is the all-one vector. Let `C` be `s` active columns, excluding the DC column, and let `P` project onto the mean-zero row space.

The complementary full set has `m-s` columns, including DC. Therefore

\[
\dim\ker(F_{C^c}^*)\ge q+s-m.
\]

If `s>m-q`, this space is nonzero. It lies in `1^perp`, since DC is complementary. Every vector `u` in it satisfies

\[
F_C F_C^*u=m u,\qquad
P F_C F_C^*P u=m u.
\]

The operator norm cannot exceed `sqrt m`, so

\[
\boxed{\|P F_C\|_{\rm op}^2=m,}
\]

with top-eigenvalue multiplicity at least `q+s-m`. The strict threshold is what this dimension argument guarantees; at equality it does not force positive multiplicity. If a balanced repair differs from `P F` by `o(sqrt m)` in full operator norm, restricting that repair to `C` costs no more, and the repaired core has squared operator norm at least `m-o(m)`.

This is an exact Euclidean restricted-frame obstruction. It does not produce a Boolean spin vector attaining the Euclidean eigenvalue, a scalar energy lower bound, or a typical recursive-frame statement beyond the stated selected-Hadamard model.

The frozen mechanism Sections 10.3--10.4 and complete companion source `computations/bh_mechanism_2026_09_16_frame_saturation.py` were read. Independent replay with `.venv/bin/python` passed all 81 exact rational selected-frame/core checks and the equality-threshold counterexample `m=4,q=3,s=1`, whose centered squared norm is `8/3<4`; `py_compile` also passed. These are frame identities, not Boolean cap-attainment tests or numerical verification of the asymptotic repair theorem.

### Density consequence for the unmasked two-level Schur certificate

The root/mechanism final corollary also passes. Assume every nonself repaired column has squared norm `q`, every repair has full operator error at most `epsilon sqrt m` from its centered selected-Hadamard frame, and `epsilon->0`. Write `p_m=q/m`, `N=mq`, `V=rho N`, and `M=(1+epsilon)^2 m`. A nonempty core has common restricted-Gram upper bound `K>=q`, because at least one row sees a full nonself column. Positive two-level gain requires `theta=V_C/V>K/M`, while `V_C<=|C|q`; hence

\[
|C|/m>\rho p_m/(1+\epsilon)^2.
\]

If

\[
\rho>b_m:=\frac{(1+\epsilon)^2(1-p_m+1/m)}{p_m},
\]

then `|C|-1>m-q`. Even after deleting each row's possible self column, every core exceeds the exact rank threshold above. For `epsilon<1`, the repaired core therefore has norm at least `(1-epsilon)sqrt m`, forcing `K>=(1-epsilon)^2m`. The Schur gain is bounded by

\[
\tfrac12\big(\sqrt{(M-K)V_C}-\sqrt{K V_L}\big)_+^2
\le\tfrac12(M-K)V_C\le2\epsilon mV.
\]

If `p_m->p>0`, this is `o(N^(3/2))`, uniformly over profiles since `V<=N`. In the complementary branch `rho<=b_m`, the basic bound itself gives

\[
\frac{|H_{\rm bulk}|}{N^{3/2}}
\le\frac{M b_m}{2\sqrt N}
=\frac{(1+\epsilon)^4(1-p_m+1/m)}{2p_m^{3/2}}
\longrightarrow\frac{1-p}{2p^{3/2}}.
\]

At `p=31/32`, the limit is about `0.01639`. This low-density branch is small on the intended competitive bulk scale but is **not necessarily** `o(N^(3/2))`. The conclusion is that no fixed leading improvement of this particular common-`K`, two-level Euclidean Schur certificate occurs on profiles above the explicit small-bulk threshold. It relies on unmasked full columns and the selected-Hadamard/operator-repair structure; it does not forbid other scalar-energy or profile-sensitive mechanisms.

## 13. Closing synthesis scope audit

The complete `artifacts/bh_2026_09_16_final_synthesis.md` was read after insertion of the actual-class power theorem. Its substantive mathematical summary passes: the original interval remains a reported unchanged frontier rather than a newly reconstructed result; the analytic exponent refinement is not claimed novel; power, nonlinear-filter, restriction, ellipse, port, Bohr, and generic-moment conclusions retain their distinct hypothesis classes; and the final moment-limit implication sends the fixed-`alpha` dimension limit first and then `alpha` to infinity.

One table wording qualification was requested and confirmed corrected in the synthesis: the signed-average variance bound `eta=||c||_1^2/(D||c||_2^2)` applies **when the output level is fully flat**, not to every signed average. Identity averaging on a nonflat level falsifies the unqualified wording. The underlying detailed theorem and proof already contained the condition; the final summary now preserves it as well.

No substantive proof gap or original minimax improvement was found in the final bounded audit. Generic quadrature, planted bounded-cap examples, and Euclidean frame obstructions are not substitutes for actual optimizer-specific signed moment or scalar extremum information.
