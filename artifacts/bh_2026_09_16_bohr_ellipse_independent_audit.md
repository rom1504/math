# Independent audit: Boolean Bohr localization and the uniform ellipse barrier

2026-09-16. Third bounded audit by the proof-verifier agent. The results below independently check `artifacts/bh_2026_09_16_bohr_flat_separation.md` and Section 11 of `artifacts/bh_2026_09_16_mechanism_analysis.md`. No quadratic minimax improvement or novelty claim is made.

## 1. Primary-source hypotheses and coefficient classes

The primary [Defant--Mastylo--Perez Boolean Bohr paper](https://arxiv.org/html/1707.09186v1), Definitions 1.1--1.2 and all of Section 3, was read directly. Its functions and cube are real. Lemma 3.3 chooses a sign for every prescribed coefficient; setting all coefficients on the complete degree-`m` level to one therefore gives a genuinely flat-sign upper bound. Its constant `C_0=6 sqrt(log 2)` is uniform in both dimension and degree. The real all-function radius is `2^(1/n)-1`; the union of homogeneous classes instead has radius asymptotic to `sqrt(log n/n)`.

The [complex polydisc paper](https://arxiv.org/html/1310.2834), Introduction, uses a supremum on the complex polydisc and arbitrary analytic monomials. This is a different domain. The [Boolean Sidon paper](https://arxiv.org/html/2302.00233), Section 2 and Theorem 5.1, uses real function spaces and arbitrary coefficient magnitudes. Its degree-dependent comparison with a projection constant is not an identity on the flat-sign subset.

Primary Boolean Bohr downloads were preserved as `tmp/bh_2026_09_16/1707.09186v1_independent_bohr.html` and `.pdf`. SHA-256 values are respectively `623e233cf6acb23d07ad6d9325dc1191d8563528c953128acacff84f78818554` and `b8c213e102b92fd597d1850657aa34a186c592a9cc52ed491120f4164e31ac22`.

### Exact witnesses separating the domains and coefficient classes

Consider the six-term full signing

\[
F(z)=z_1z_2-z_3z_4+(z_1+z_2)(z_3+z_4).
\]

Its real cube norm is four. Writing `z_1,z_2=e^{i(alpha+-u)}` and `z_3,z_4=e^{i(beta+-v)}`, removal of the common phase gives

\[
F=2i\sin(\alpha-\beta)+4\cos u\cos v.
\]

Hence its complex torus norm is `sqrt(20)=2 sqrt 5`, attained at `(i,i,1,1)`. Complete-support real sign coefficients do not remove the real/complex norm discrepancy.

Furthermore, any six-term signing on four variables has an even-integer value at every cube point. Parseval gives norm at least `sqrt 6`; even integrality upgrades this to four. The preceding witness proves `M_(4,2)=4`, so the flat Sidon ratio is exactly `6/4=3/2`. In contrast,

\[
x_1x_3+x_1x_4+x_2x_3-x_2x_4
=x_1(x_3+x_4)+x_2(x_3-x_4)
\]

has norm two and coefficient absolute sum four. The unrestricted real quadratic Sidon constant is therefore at least two, strictly above the flat ratio. This is a finite exact counterexample to identifying the two coefficient classes.

Changing only the scalar field also matters for nonhomogeneous Boolean Bohr radii. On the real cube, `f=1+i epsilon x_1` has norm `sqrt(1+epsilon^2)`, while its weighted absolute coefficient sum is `1+rho epsilon`. Any universal radius would satisfy `rho<=epsilon/(sqrt(1+epsilon^2)+1)` for all positive `epsilon`; therefore the all-function complex-valued radius is zero. This observation does not concern the complex polydisc theorem.

## 2. The flat homogeneous variational parameter

Set `N=binom(n,m)` and

\[
M_{n,m}=\min_{\varepsilon_S\in\{-1,1\}}
 \left\|\sum_{|S|=m}\varepsilon_Sx^S\right\|_\infty,
\qquad r^{\rm flat}_{n,m}=(M_{n,m}/N)^{1/m}.
\]

Let `r_(n,m)` be the corresponding radius for all real homogeneous degree-`m` coefficient arrays. Then `r_(n,m)<=r_flat_(n,m)`. Let `R_n` and `R_n^flat` be their minima over all integers `1<=m<=n`.

The flat upper input can be reconstructed without any complex-domain theorem. Choose independent signs on all `m`-sets. For every fixed cube vertex, the sum has tail at most `2 exp(-T^2/(2N))`. A union bound over `2^n` vertices at

\[
T=\sqrt{2N(n+2)\log 2}
\]

has failure probability at most one half. Thus one full signing has norm at most `sqrt(6 log 2) sqrt(nN)` for every `n>=1`. In particular, either this constant or the larger primary-source constant `C_0` yields the fully uniform upper bound

\[
r^{\rm flat}_{n,m}\le C_0^{1/m}(n/N)^{1/(2m)}.
\]

The argument does not suppress coefficients, change the scalar field, or assume that `m` is fixed.

## 3. Uniform quantitative consequence of polynomial BH

Assume the dimension-free real Boolean BH bound `B_m<=K m^a`, with fixed `K>=1,a>=0`. This applies to the new paper's exponent 27 and the independently reconstructed exponent `11/2`; their distinction affects constants, not the rate below.

Hölder at `q_m=2m/(m+1)` gives

\[
\sum_{|S|=m}|\widehat f(S)|
\le N^{(m-1)/(2m)}B_m\|f\|_\infty,
\quad
r_{n,m}\ge B_m^{-1/m}N^{-(m-1)/(2m^2)}.
\]

Using `N>=(n/m)^m`, this implies

\[
r_{n,m}\ge
\exp\!\left[-\frac{(a+1/2)\log m+\log K}{m}\right]
(n/N)^{1/(2m)}.
\]

Put `L=log n`, `A=a+1/2`, `B=log K`, and

\[
\Phi(t)=\frac12(\log t-1+t^{-1}).
\]

The elementary upper estimate `N<=(en/m)^m` now yields

\[
\log\frac{r_{n,m}}{\sqrt{L/n}}
\ge\Phi(m/L)-\frac{A\log m+B}{m}.                 \tag{1}
\]

The function `Phi` is nonnegative, vanishes only at one, and tends to infinity at both endpoints of the positive half-line. The next argument verifies the claimed error uniformly over **every** integer degree, including degrees comparable to `n`.

If `m<=L^2`, define `L'=L-4A log L-2B`, positive for all sufficiently large `n`. The right side of (1) equals

\[
\tfrac12\{\log(m/L)-1+(L-2A\log m-2B)/m\}
\ge\Phi(m/L')+\tfrac12\log(L'/L).
\]

The inequality uses only `log m<=2 log L`. Therefore this is at least `-O(log L/L)`. If `m>L^2`, then

\[
\Phi(m/L)\ge\tfrac12(\log L-1),\qquad
\frac{A\log m+B}{m}\le\frac{2A\log L+B}{L^2},
\]

where the second bound follows because `(A log u+B)/u` is decreasing for sufficiently large `u`. The lower bound is then positive. This proves the lower estimate for both degree-minimized radii with no unbounded numerical check or exchange of a fixed-degree limit.

For the upper estimate take `m=floor L`. The factorial expansion and
`log binom(n,m)=m log n-log(m!)+O(m^2/n)` give

\[
\log\frac{C_0^{1/m}(n/\binom nm)^{1/(2m)}}{\sqrt{L/n}}
=\Phi(m/L)+\frac{\log(2\pi m)}{4m}
 +\frac{\log C_0}{m}+O(m^{-2}+m/n)
=O(\log L/L).
\]

Consequently

\[
\boxed{R_n=\sqrt{\frac{\log n}{n}}
 \left[1+O\!\left(\frac{\log\log n}{\log n}\right)\right],
\quad
R_n^{\rm flat}=\sqrt{\frac{\log n}{n}}
 \left[1+O\!\left(\frac{\log\log n}{\log n}\right)\right].}
\]

This is a valid derived rate improvement in the Boolean setting. The older subexponential BH input already gives the leading asymptotic constant; this audit does not claim the rate corollary is new in the literature.

## 4. Exact minimizing degrees and the quadratic separation

Every asymptotically minimizing degree sequence satisfies `m/log n->1`. To see this first exclude bounded degrees: their radius divided by `sqrt(L/n)` diverges. For `m->infinity`, the penalty in (1) tends to zero; any radius asymptotic to `sqrt(L/n)` forces `Phi(m/L)->0`, hence `m/L->1`.

For every **exact** minimizer, the quantitative upper bound in Section 3 then combines with (1). In a fixed neighborhood of one, `Phi(t)` is comparable to `(t-1)^2`, and the penalty is `O(log L/L)`. Therefore

\[
\boxed{m=\log n+O(\sqrt{\log n\,\log\log n}).}
\]

This sharper window also holds for near-minimizers whose relative error is `O(log log n/log n)`. Merely specifying a `1+o(1)` relative error gives only `m/log n->1`, not this rate.

For each fixed integer `D`, all degrees `1,...,D` are eventually absent from either minimum: each has a diverging normalized radius, while degree `floor log n` has bounded normalized radius. Removing these fixed degrees therefore leaves the minima **exactly unchanged** for all sufficiently large dimensions.

The original problem corresponds exactly to `m=2`:

\[
r^{\rm flat}_{n,2}=\sqrt{M_n/\binom n2},\qquad
n^{1/4}r^{\rm flat}_{n,2}
=\sqrt{\frac{2(M_n/n^{3/2})}{1-1/n}}.
\]

Thus convergence of this fixed-degree flat radius is equivalent to the original normalized-cap convergence question. Its scale `n^(-1/4)` is much larger than `sqrt(log n/n)`. The aggregate radius asymptotic does not recover the quadratic entry through its defining variational minimum, because that entry is eventually never selected.

Multiplication by `x_1...x_n` preserves the supremum and bijects complete degree-two and complete degree-`n-2` signings. It gives `r_flat_(n,n-2)=(r_flat_(n,2))^(2/(n-2))->1`, another inactive scale. Also, an `m`th-root comparison of radii at growing degree does not imply a ratio-one comparison of the underlying Sidon constants.

## 5. Uniform joint ellipses: assumptions and Paley constraints

This separate audit checks the root/mechanism ellipse argument. Suppose, uniformly for full sign quadratics, arbitrary constants `c`, and arbitrary linear fields `h`,

\[
\|c+Q_A+h\cdot x\|_\infty^2
\ge a\|h\|_1^2+b n^3-o(n^3),\qquad a,b\ge0.       \tag{2}
\]

The uniformity is essential. A theorem only about typical independently sampled fields is not assumed. It is enough to have uniformity on every fixed natural-scale class `||h||_1<=C n^(3/2)`, with the limit taken for each fixed `C`.

The square-Paley construction in Section 2 of the mechanism artifact has `A1=0` and `A^2=nI-J`. For `h=t sqrt n 1`, its normalized norm is bounded above by `(1+t^2)/2` for `0<=t<=1`, and by `t` for `t>=1`; the row-constant spin construction asymptotically attains the first branch. Setting `c=0` in the proposed uniform theorem gives necessary constraints

\[
0\le a\le1,\qquad
b\le\begin{cases}1/4,&0\le a\le1/2,\\ a(1-a),&1/2\le a\le1.\end{cases}       \tag{3}
\]

Indeed, the second branch forces `a<=1`, and the minimum of `(1+t^2)^2/4-a t^2` over `0<=t<=1` is obtained at `t=0` if `a<=1/2`, or at `t^2=2a-1` otherwise. Only this explicit deterministic field is needed to test (2); it need not be typical under restrictions.

## 6. The restriction bootstrap cannot cross the ellipse barrier

Freeze `(1-theta)n` coordinates as independent signs, leaving `theta n` free. Each restricted linear field coordinate is a sign random walk of length `(1-theta)n`; no independence between different coordinates is needed. Its first absolute moment and Jensen yield the prospective parent certificate

\[
\frac{\|c+Q_A\|_\infty^2}{n^3}
\ge\max_{0\le\theta\le1}
 \{a s\theta^2(1-\theta)+b\theta^3\}-o(1),
\qquad s=2/\pi.                                    \tag{4}
\]

If (2) is uniform only on each fixed bounded natural-scale field class, truncate before applying it and send that cutoff to infinity after the dimension limit. For `H=||h||_1`, the bound `E H^2<=k^2 q`, where `k` is the number of free and `q` the number of frozen coordinates, gives `P(H>C k^(3/2))=O(C^-2)` and `E[H;H>C k^(3/2)]=O(n^(3/2)/C)` for fixed `theta`. Thus the truncation preserves the first-moment/Jensen lower term and the zero-field term in the required order of limits. Arbitrary constant shifts remain required.

Let `B=w^2` be the existing universal half-range squared lower constant, and require `b<=B` so that the ellipse is not already an independently stronger zero-field theorem. Define

\[
B_*=(s/3)(1-s/3)=0.16717495361482146\ldots,
\qquad \sqrt{B_*}=0.4088703383895944\ldots.
\]

For every `B in [B_*,1/4]`, the right side of (4), under (3) and `b<=B`, is at most `B`. Here is a complete optimization check. Put

\[
a_+=(1+\sqrt{1-4B})/2,\qquad a_*=1-s/3.
\]

Then `a_+<=a_*`. For `K=a s`, elementary differentiation gives

\[
\max_{0\le\theta\le1}\{K\theta^2(1-\theta)+b\theta^3\}
=\begin{cases}
b,&K\le3b,\\
4K^3/[27(K-b)^2],&K>3b.
\end{cases}
\]

If `a<=a_+`, increase `b` to `B`; then `K<=3B`, so the maximum is at most `B`. If `a_+<a<=a_*`, increase `b` to `a(1-a)<=B`; then `K<=3a(1-a)`, with the same conclusion. If `a>a_*`, use `b=a(1-a)` again; the maximum becomes

\[
\frac{4a s^3}{27(a+s-1)^2}.
\]

Its derivative has the sign of `s-1-a`, hence it decreases on this interval. Its value at `a_*` is `B_*`, completing the proof.

Therefore this entire uniform norm-only ellipse bootstrap cannot improve the present `w=0.4333221116640807`, which is already above the barrier. The argument does not exclude field-law-sensitive inequalities, nonlinear joint functionals, optimizer-specific information, or an independently improved zero-field coefficient `b>B`. It also does not identify the best possible inequality inside the feasible ellipse region.

## 7. Audit verdict

The Boolean Bohr rate, complete-support flat upper bound, integer-degree uniformity, localization window with its correct quantifier, finite domain/class separations, and generic ellipse barrier all pass independent reconstruction. These are discriminating secondary results and obstructions, not progress on the original normalized minimum or its convergence. The bounded-cap nonlinear-filter extension, including actual near-minimizers, is separately preserved in Section 9 of `artifacts/bh_2026_09_16_ports_filters_independent_audit.md`.

## 8. Supplement: permutation flattening fails for Hilbert-valued BH and influence

This bounded follow-up independently checks the power agent's proposed orbit flattening. For a scalar cube polynomial `f`, set

\[
G(x)=\big(f(\pi x)/\sqrt{n!}\big)_{\pi\in S_n},
\qquad V_r=\sum_{|S|=r}|\widehat f(S)|^2.
\]

Permutation transitivity gives

\[
\|G\|_{L_\infty(\ell_2)}\le\|f\|_\infty,
\qquad
\|\widehat G(S)\|_{\ell_2}^2=V_r/\binom nr
\quad(|S|=r).
\]

The Hilbert-valued influences are therefore exactly

\[
\operatorname{Inf}_i(G)=
\sum_{S\ni i}\|\widehat G(S)\|_2^2
=\frac1n\sum_r rV_r.
\]

These identities really do produce level-flat Hilbert coefficient norms. They do not produce a scalar level-flat polynomial, and a dimension-free Hilbert-valued version of the paper's scalar BH or influence conclusion is false even within this permutation-orbit subclass.

Indeed, choose `f=chi_(S_0)` for any fixed set of size `m>=1`, and let `D=binom(n,m)`. Then `G` has pointwise Hilbert norm exactly one, mean zero, and total variance one. Every degree-`m` coefficient has Hilbert norm `D^(-1/2)`. Consequently its coefficient norm at `q_m` is

\[
\left(\sum_{|S|=m}\|\widehat G(S)\|_2^{q_m}\right)^{1/q_m}
=D^{1/(2m)},
\]

and its weighted quantity is `m^(-s) D^(1/(2m))`. Both diverge with dimension for fixed `m`, whereas every coordinate influence is exactly `m/n`, tending to zero. This supplies a concrete falsifier for both proposed dimension-free vector extensions, including any influence lower bound that depends only on fixed degree and nonzero total variance. No generic Hilbert-space counterexample outside the orbit class is needed.

The failure does not contradict scalar hypercontractivity, which can extend to certain vector settings. The missing scalar BH/influence conclusion is independently ruled out by the displayed exact coefficients and influences.

### Exact Gram matrix for a row-balanced full quadratic signing

The power agent subsequently proposed the following sharper target-class calculation; it also passes independent reconstruction. Let `n>=5`, `D=binom(n,2)`, and let the hollow full signing have every row sum zero. Write `v_e` for the Hilbert Fourier coefficient of the orbit lift at edge `e`, and `K_(e,f)=<v_e,v_f>`. Permutation averaging gives three values:

\[
K_{e,f}=\begin{cases}
1,&e=f,\\
-1/(n-2),&|e\cap f|=1,\\
2/((n-2)(n-3)),&e\cap f=\varnothing.
\end{cases}
\]

For adjacent edges, square a zero row sum: the sum of ordered distinct products is `-(n-1)`, giving the stated average. For disjoint edges, the total edge sum is zero, so each Gram row sums to zero. There are `2(n-2)` adjacent and `binom(n-2,2)` disjoint edges, which determines the third value.

Let `B` be the unsigned vertex-edge incidence matrix. Its rank is `n`: if `z_i+z_j=0` for every pair, any triangle forces all entries to vanish. Every orbit coefficient vector is row-balanced, so `K` annihilates the row space of `B`. On `ker B`, the line-graph adjacency operator is `-2I`, because the incident sums at the endpoints both vanish. The disjoint-edge adjacency operator is `I`, because the total edge sum also vanishes. Consequently

\[
\boxed{K=\lambda P_{\ker B},\qquad
\lambda=\frac{n-1}{n-3},\qquad
\operatorname{rank}K=D-n.}
\]

In particular `lambda(D-n)=D`, as required by the trace. Let `T` send a scalarization vector `u` in the orbit Hilbert space to the coefficients `b_e=<u,v_e>`. Since `TT^*=K`, every unit `u` produces scalar variance at most `lambda`. For any target `b in ker B`, the minimum squared Hilbert norm needed to produce it is

\[
\min_{Tu=b}\|u\|_2^2=\|b\|_2^2/\lambda,
\]

attained by `u=T^*b/lambda`. In particular every row-balanced flat target, including the original signing, requires and admits

\[
\|u\|_2=\sqrt{D/\lambda}=\sqrt{D-n}.
\]

This exactly quantifies the cost of certifying a scalar cap by the generic inequality `||<u,G>||_infinity<=||u||_2 ||G||_(infinity;ell_2)`. It is **not** a lower bound on the actual cap of the resulting scalar polynomial. Indeed, reproducing the original signing is consistent with its original cap even though this Hilbert-Cauchy certificate is much larger. For nonflat targets the correct formula is `||b||_2/sqrt lambda`, not the flat-target factor.

### Final appendix audit: signed averages, exact Hilbert cap, and scalar normalization

The complete final Section 9 of `artifacts/bh_2026_09_16_power_analysis.md` was independently read and checked. The exact replay `.venv/bin/python computations/bh_2026_09_16_orbit_flattening_checks.py` and its `py_compile` both passed. This replay verifies all 120 order-five permutations, all 32 spin assignments for the orbit norm formula, the stored order-six identity, order-seven Fourier histograms, and the auxiliary-variable scalarization cap.

For a signed orbit average `Tf=sum w_pi f(pi x)` with `sum|w_pi|<=1`, each level's coefficient `l1` norm contracts. If the output level has a common magnitude on all `D_r` supports, its variance is therefore at most `||c_r||_1^2/D_r=eta_r V_r`. Cauchy--Schwarz makes `eta_r<=1`, with equality for a nonzero level exactly when the original magnitudes were already fully flat. The stored `H_7^2` histograms are `(18 at 2,3 at 6)` on level two and `(30 at 2,5 at 6)` on level four; each gives `eta=27/35`. A level-flat signed average consequently loses at least `8/35` of the square's nonconstant variance. The order-six exception is real: `H_6^2=15+2XH_6`, `X=x_1...x_6`, makes every power lie in the invariant span of `1,H_6,XH_6,X`, which occupies four separate flat levels. No generic nonflatness assertion can include this witness.

For a full row-balanced quadratic orbit, let `R` denote its **actual Hilbert supremum**, not the scalar cap. The exact formula is

\[
R^2=\frac{(n-1)^2(n+1)}{2(n-2)}.
\]

Indeed, with `s=sum x_i` and `v_(ij)=x_i x_j`, the inverse of `BB^*=(n-2)I+J` gives

\[
\|P_{\ker B}v\|^2=D-s^2-\frac n{n-2}
 +\frac{(s^2-n)^2}{2(n-2)(n-1)}.
\]

The expression is convex in `s^2`. Here `n` is odd, so the feasible interval has endpoints one and `n^2`, both attained. The latter gives zero; the former is the maximum. Multiply by `lambda=(n-1)/(n-3)` to obtain the claimed norm. Therefore `G/R` has variance `n(n-2)/(n^2-1)->1`, whereas every influence is `2/n` times that variance. This refutes a dimension-free vector influence extension even on full row-balanced quadratic orbit lifts.

For comparison, take the orbit of the **scalar-cap-normalized** power `f=(H_A/Q(A))^k`, let its variance be `V>0`, common influence be `I`, and actual nonconstant degree be `m`. Even Fourier support gives `I>=2V/n`; boundedness of `|H_A/Q(A)|` gives

\[
V\le\mathbb E|H_A/Q(A)|^{2k}\le N/Q(A)^2,
\quad N=\binom n2.
\]

Thus `I/V^2>=2Q(A)^2/(nN)`, and `m>=2`. Equality throughout is attained at the first power, where `V=N/Q(A)^2` and `m=2`. For every `a>=0` and `n>=2`,

\[
\boxed{\inf_{A,k\ge1:\,V>0}m^a\frac I{V^2}
=\frac{2^{a+1}M_n^2}{nN}.}
\]

The explicit `V>0` restriction handles any constant powers; no claim that every power is nonconstant is needed. This identity does not contradict the `R`-normalized counterexample: `R` is of order `n` for row-balanced quadratic orbits, while the scalar `Q` is at least order `n^(3/2)`. A restricted influence theorem under scalar-`Q` normalization can therefore be true, but its optimal constant is exactly a reformulation of the original minimax quantity. All these normalization and scope distinctions in the final appendix pass audit.
