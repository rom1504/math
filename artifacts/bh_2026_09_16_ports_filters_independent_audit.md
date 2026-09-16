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

The thresholds `a=1` and `b=1/2` are not covered: the required interval for `beta` then closes. Bounded filter degrees are also not covered by this ratio-to-zero statement, and `c,L` must remain fixed. The original scope restriction remains: this is not a theorem about the moment or operator norm of actual cap minimizers, nor an obstruction to model-specific coefficient information beyond the universal polynomial-growth inequalities.

## 8. What is and is not certified

Certified: the deterministic port inequality, its rational numerical phase, the negligibility implication, arbitrary-complex-polynomial coefficient extraction, the two nonlinear-filter regimes with their exact limit orders, the all-signing hypercontractive extension, and the global-slope H6 obstruction.

Not certified or claimed: a positive-density variance phase, a new dense-sign construction constant, an original-problem recurrence, convergence of normalized minima, a uniform theorem at degree comparable to `sqrt n` or larger for **all** signings, or an asymptotic-only weighted-BH slope. The stronger sublinear-degree result in Section 7 has a bounded-operator hypothesis. The old restricted-frame construction is used with its conditioned selector and operator-repair hypotheses intact.
