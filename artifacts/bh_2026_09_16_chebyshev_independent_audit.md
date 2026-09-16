# Independent audit: Chebyshev recovery and generic moment ambiguity

2026-09-16. This bounded audit checks the complete frozen draft `artifacts/bh_2026_09_16_chebyshev_moment_recovery.md`, its replay, and the relevant archived theorem. No original minimax improvement or quadratic-energy realization is claimed.

## 1. Archive and primary-source status

The archive `artifacts/blank_slate_direct_attack_2026_08_21.md`, Sections 2 through 2.3, was read directly. Its Theorem 2.1 already proves the endpoint Chebyshev factor, including the exact endpoint-mass version, the uniform factor based on `2^(1-n)`, the parity-lattice moment LP, and the warning that the remaining limiting moment problem is quantitatively equivalent to the original convergence problem. The new draft correctly treats this upper theorem as archived, not new.

The primary [Kuczynski--Wozniakowski paper](https://epubs.siam.org/doi/10.1137/0613066) was checked at its abstract and publication record. It discusses the analogous distinction between power and Lanczos error scales for random-start eigenvalue estimation. That comparison is background only: none of its matrix-algorithm conclusions are transferred to signing energies. The proof used here is reconstructed below. The old/new Boolean BH constants are from the primary [2026 paper](https://arxiv.org/html/2609.12427v1), already fully audited separately.

## 2. Exact endpoint recovery and the nonmonotone scalar test

Let a finite real random variable `Y` have positive cap `R` and `P(|Y|=R)>=p>0`. If a probability measure supported on `[-c,c]` has the same even moments through order `2d`, then the expectations of the even polynomial `T_d(Y/c)^2` agree. For `R>c`, the competing expectation is at most one and the original expectation is at least

\[
p\cosh^2\!\left(d\operatorname{arcosh}(R/c)\right).
\]

Thus

\[
\boxed{c\ge R\operatorname{sech}
 \left(\operatorname{arcosh}(p^{-1/2})/d\right).}
\]

The case `c>=R` is immediate. Only even moments are needed because `T_d(-z)^2=T_d(z)^2`. Restricting competing measures to the correct parity lattice cannot weaken this lower bound. The abstract argument is valid for every integer `d>=1`; the archive uses the moment-LP degree range relevant to its cut-code formulation.

For an actual energy `Y=H_A/n^(3/2)`, global spin reversal pairs every maximizing spin with a distinct maximizing spin, so the absolute endpoint event has probability at least `p_n=2^(-(n-1))`. Set

\[
L_n=\operatorname{arcosh}(2^{(n-1)/2}).
\]

Both the minimum compatible central support radius and

\[
C_d=\inf\{c>0:\mathbb E T_d(Y/c)^2\le1\}
\]

lie between `R sech(L_n/d)` and `R`. Every `c>=R` is feasible for the latter definition, and the endpoint argument excludes every `c` below its lower bound. No monotonicity of `c -> E T_d(Y/c)^2` is used or generally available.

As `d/n->infinity`, `L_n=(n log 2)/2+O(1)` and `1-sech z=z^2/2+O(z^4)` yield a uniformly bounded-cap recovery error `O((n/d)^2)`. If `d=alpha n` for fixed positive `alpha`, the multiplicative recovery factor tends to `sech(log 2/(2alpha))`. The single `2d`th-moment proxy only gives the endpoint factor `exp(-(n-1)log 2/(2d))`, with error `O(n/d)`. These limit orders are correctly separated in the draft.

## 3. Linear scalar degree is complete before BH is used

Fix `c,epsilon>0`, assume `R>=c+epsilon`, and write `alpha=arcosh(1+epsilon/c)`. The endpoint estimate gives

\[
\|T_d(Y/c)\|_2\ge\sqrt{p_n}\cosh(d\alpha).
\]

Hence the false cap `R<=c` is contradicted by `L2` alone whenever

\[
d>\operatorname{arcosh}(p_n^{-1/2})/\alpha
=\frac{n\log2}{2\alpha}+O(1/\alpha).
\]

Here `d` is scalar polynomial degree, not the reduced Walsh degree. The latter is at most `m=min(2d,n)`; no equality is needed in this argument. For small `epsilon/c`, `alpha~sqrt(2epsilon/c)`, with fixed `c` understood.

For the full BH test, its coefficient norm dominates `L2`. A sufficient contradiction to the **cap-conditioned** test is

\[
d>\operatorname{arcosh}(p_n^{-1/2}\overline B_n)/\alpha,
\]

where `Bbar_n` bounds all degree constants through `n`. The older subexponential bound adds `O(sqrt(n log n)/alpha)` scalar degree; a polynomial bound adds `O(log n/alpha)`. Neither is needed for the `L2` test, which already has only the displayed `O(1/alpha)` additive error. This is not a violation of a valid BH inequality: it refutes the false bound on the composed function's supremum.

For the paper's weighted seminorm, the missing constant coefficient is handled explicitly. Since `E Y^2=(n-1)/(2n^2)`, fixed `c>0` gives `P(|Y|<=c)>=1/2` for all sufficiently large dimensions. Comparing this bulk event with the endpoint event in the identity `Var f=(1/2)E(f-f')^2` gives

\[
\operatorname{Var}(T_d(Y/c))
\ge\frac{p_n}{2}\big(\cosh(d\alpha)-1\big)^2.
\]

The endpoint may have either sign, and `d` may have either parity: its filtered magnitude is at least `cosh(d alpha)`, while bulk filtered values have magnitude at most one, so the difference estimate still holds. Since `S(f)>=m^(-5)sqrt(Var f)`, this refutes `S(f)<=K m^22` once the displayed square root exceeds `K m^27`, again with only polynomial-BH logarithmic overhead. No constant-mode contribution has been inserted into the seminorm.

## 4. Positive prescribed-node quadrature: every proof step

Let `nu` be arcsine probability measure on `[-1,1]`. Its degree-`d` reproducing kernel in the orthonormal Chebyshev basis is

\[
K_d(z,x)=1+2\sum_{j=1}^d T_j(z)T_j(x).
\]

Fix `R>1`. For every polynomial `q` of degree at most `d-1`, reproduction gives

\[
\int(R-x)K_d(R,x)q(x)\,d\nu(x)=R q(R)-(xq)(R)=0.
\]

The measure `(R-x)dnu` is strictly positive on the interval, and `K_d(R,x)` has degree exactly `d`. It therefore has `d` simple roots in `(-1,1)`. One elementary verification is to multiply it by the polynomial formed from its interior sign-changing roots: if there were fewer than `d`, the product would have constant nonzero sign on the interval, contradicting the displayed orthogonality. This also forces all roots to be interior and simple.

Use these roots together with `R` as interpolation nodes `x_j`, and define weights by integrating their degree-`d` Lagrange polynomials `L_j`. For any polynomial of degree at most `2d`, division by `(x-R)K_d(R,x)` leaves a remainder of degree at most `d` and a quotient of degree at most `d-1`. The quotient part integrates to zero by orthogonality and vanishes at every node; the remainder is integrated exactly by interpolation. Thus the quadrature is exact through degree **`2d`**, without claiming degree `2d+1`.

Exactness on `L_j^2` gives

\[
w_j=\int L_j(x)^2\,d\nu(x)>0.
\]

Exactness on `L_j q` for arbitrary `q` of degree at most `d` gives `int L_j q=w_j q(x_j)`. Reproduction identifies `L_j=w_j K_d(x_j,.)`; evaluation at `x_j` yields

\[
\boxed{w_j=1/K_d(x_j,x_j).}
\]

Constants are integrated exactly, so these positive weights sum to one. On every interior node, `|T_j(x)|<=1` implies `K_d(x,x)<=2d+1` and hence weight at least `1/(2d+1)`. The prescribed outer node has weight `1/K_d(R,R)`.

The comparison measure is ordinary `(d+1)`-node Gauss--Chebyshev quadrature, with nodes `cos((2j-1)pi/(2(d+1)))` and weights `1/(d+1)`. It matches arcsine moments through degree `2d+1`, so the two finite positive measures match through degree `2d` as claimed.

## 5. Every atom satisfies the mass constraint; asymptotic constants

For `0<p<1/(2d+1)`, `K_d(R,R)` is continuous and strictly increasing for `R>=1`, starts at `2d+1`, and diverges. Hence there is a unique `R=cosh alpha>1` with `K_d(R,R)=1/p`. The outer quadrature mass is exactly `p`, all its interior masses exceed `p`, and every comparison mass `1/(d+1)` also exceeds `p`.

The kernel satisfies

\[
\tfrac12e^{2d\alpha}\le K_d(\cosh\alpha,\cosh\alpha)
\le(2d+1)e^{2d\alpha}.
\]

The lower bound uses just `2cosh^2(d alpha)>=e^(2d alpha)/2`; the upper uses each summand separately. With `L=log(1/p)`, this gives

\[
\frac{L-\log(2d+1)}{2d}\le\alpha
\le\frac{L+\log2}{2d}.
\]

Under the two stated limit assumptions `log d=o(L)` and `L=o(d)`,

\[
\alpha\sim L/(2d),\qquad
R-1\sim L^2/(8d^2).
\]

The comparison cap is `cos(pi/(2(d+1)))=1-O(d^(-2))`. Since the assumptions imply `L->infinity`, this correction is negligible compared with `L^2/d^2`. The two measures' actual cap gap is therefore also asymptotic to `L^2/(8d^2)`. This matches the leading one-sided endpoint recovery loss. Any estimator based only on the common moments and mass floor has error at least half the gap on one of the two measures.

Taking `p=2^(-(n-1))` and polynomially bounded `d>>n` satisfies both hypotheses and gives the generic order `(n/d)^2`. The construction is deliberately an **arbitrary-positive-measure** example. Its weights are not claimed dyadic, its supports are not claimed to lie on an energy parity lattice, and no matrices realizing these measures are produced. It cannot establish a lower bound for moment recovery restricted to actual full-sign quadratic energies. Extra lattice, support-cardinality, or quadratic-chaos identities can change that narrower problem.

The root's final antipodal refinement is also valid: first construct the measures with prescribed endpoint weight `2p`, assuming `2p<1/(2d+1)`, and then symmetrize both under `x -> -x`. Every resulting atom has mass at least `p` because splitting halves a mass of at least `2p` and merging only increases it; the outer measure now has at most `2(d+1)` atoms, moments through `2d` still match, and the cap gap retains the asymptotic `L^2/(8d^2)` because `log(1/(2p))=L-log2`. This matches antipodal symmetry without establishing dyadic weights, the energy lattice, or realizability by actual quadratic signings.

## 6. Replay and verdict

The companion script was read completely, then replayed with
`.venv/bin/python computations/bh_mechanism_2026_09_16_chebyshev_recovery.py --output tmp/bh_2026_09_16/independent_chebyshev_recovery_replay.json`.
It passed all 40 exact integer/Fraction signing-filter checks and all five floating quadrature regressions; `py_compile` also passed. The replay was preserved at that path. Its largest displayed quadrature moment error was below `2.8e-14`, comfortably within the script's `1e-8` regression tolerance. These floating results are not treated as exact certificates; Sections 4--5 supply the algebraic proof.

Audit verdict: **PASS**, with the archive status, degree distinction, mass-floor quantifier, asymptotic assumptions, and generic-versus-signing separation all maintained. The archived upper theorem is not new campaign progress; the new addition is a self-contained matching example for the expressly broader moment-oracle class, with no literature-novelty claim.

## 7. Final bounded-cap sharpness appendix: planted Sylvester clique

The complete frozen `artifacts/bh_2026_09_16_planted_tail_sharpness.md` and its companion code were independently read. The following is a full proof check, including the weighted variance extension.

Take `n=2^(4t)`, `k=2^(3t)=n^(3/4)`, `t>=1`. The standard binary-indexed Sylvester matrix `S_n` satisfies `S_n^2=nI`, has trace zero, and has leading `k` principal block `S_k`. After hollowing it to a full signing `B`, its spin quadratic is still exactly `x^T S_n x/2`, because the removed diagonal contributes its trace at every spin. Thus `Q(B)<=n^(3/2)/2` exactly. The operator norm of the hollow matrix is not needed for this bound.

The leading block has total matrix sum `k` and trace zero, so its original edge-sign sum is `k/2`. Overwrite this clique's edges by `+1`. The coefficient change is nonnegative on every edge and has absolute supremum equal to its total coefficient sum:

\[
\|H_A-H_B\|_\infty=\binom k2-k/2=k(k-2)/2.
\]

Since `k^2=n^(3/2)`,

\[
\boxed{Q(A)\le n^{3/2}-k<n^{3/2}.}
\]

Conditioning all planted spins to be positive leaves independent uniform exterior spins. Every cross and exterior edge then has mean zero, so `E[H_A | planted +]=binom(k,2)`. For `Y=H_A/n^(3/2)`, the conditional mean is `1/2-1/(2k)`, while `Y<=1`. Thus for fixed `0<a<1/2`, its conditional tail probability above `a` is at least

\[
\frac{1/2-1/(2k)-a}{1-a}
\ge c_a:=\frac{1/2-a}{2(1-a)}>0
\]

whenever `k>=1/(1/2-a)`. The unconditional tail is consequently at least `c_a 2^(-k)`. This single dimension subsequence rules out a uniform fixed-threshold bound `C exp(-gamma n^beta)` with fixed positive constants and `beta>3/4` on the class `Q<=n^(3/2)`. The supported clique all-one Rayleigh vector also gives `||A||op>=k-1`, so the cap-only spectral exponent is sharp on that class.

Fix `0<c<a<1/2`, let `kappa=arcosh(a/c)`, and take `P_d(y)=T_d(y/c)`. On the high-energy event, `P_d(Y)>=exp(kappa d)/2`, giving

\[
\|P_d(Y)\|_2^2\ge(c_a/4)\exp\{2\kappa d-k\log2\}.
\]

For any fixed `delta>0`, choose

\[
d=\left\lceil\left(\frac{\log2}{2\kappa}+\delta\right)k\right\rceil.
\]

The bound grows exponentially in `k`. Because `d=Theta(n^(3/4))=o(n)`, the reduced Walsh degree is exactly `2d` for all sufficiently large dimensions: each top `2d`-set has coefficient equal to a nonzero scalar times `d!` times an odd signed hafnian. Lower Taylor terms cannot reach this degree. The full coefficient norm dominates `L2`, so it exceeds any fixed polynomial BH threshold for the hypothetical filter bound `||P_d(Y)||infinity<=1` that would follow from the false cap `||Y||infinity<=c`.

For the weighted version, Parseval gives `E Y^2=binom(n,2)/n^3`; hence `P(|Y|<=c)>=1/2` eventually. The high event `Y>=a` and low event `|Y|<=c` are disjoint, with filtered values respectively at least `exp(kappa d)/2` and of magnitude at most one. The two ordered cross-events in `Var Z=(1/2)E(Z-Z')^2` yield

\[
\operatorname{Var}(P_d(Y))\ge
\frac{c_a}{2}\,2^{-k}\left(\tfrac12e^{\kappa d}-1\right)^2.
\]

For large `d` the difference is nonnegative, as required before squaring a one-sided bound. Since `W_s(g)>=m^(-s)||g-Eg||_2`, this also grows exponentially after any fixed polynomial degree penalty. Both the full and weighted claims therefore pass.

Scope: this proves the sharp **power of dimension** for a general bounded-cap tail and low-degree non-detection statement. It does not find the optimal constant multiplying `n^(3/4)`, prove a bound for all signings from the example, or establish sharpness on actual minimizers. In fact the conditioned mean already forces this family's normalized cap to have liminf at least `1/2`.

The exact replay `.venv/bin/python computations/bh_2026_09_16_planted_tail_checks.py` and its `py_compile` passed. At `n=16,k=8` it verifies actual cap `46`, proved cap upper bound `56`, modification coefficient sum `24`, conditional mean `28`, and conditional probability `115/128` above energy `16`; the resulting unconditional lower bound is `115/32768`. These finite data illustrate, rather than replace, the dimension-uniform proof above. Final appendix audit verdict: **PASS**.
