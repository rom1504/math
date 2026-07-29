# Strategic steering

Evidence cutoff: Wave 25, §10.78 (2026-07-29). Regenerate by Wave 30, or earlier after a decisive result or route change.

## User-stated research objective

Determine whether $`\lim_{n\to\infty}M_n/n^{3/2}`$ exists. The conjectural value $`1/2`$ is not proved and is not an additional user objective.

## Present judgment

No route proves convergence. The rigorous interval remains $`0.336493364431\ldots\le\liminf\le\limsup\le1/2`$.

The leading route remains **adaptive optimized principal restriction**, but Wave 25 again strictly weakened its best implementation. The target is now one arbitrary full cut with controlled row square and stretched-exponential centered selector coverage. Exact parent grounds, a codebook, and a two-ground mixture are unnecessary.

The constant-shortfall criterion (10.617) is **not** the strongest route; it is stronger than convergence requires. Its exact signed-cycle formulation remains sharply falsifiable, but all proved generic complete-signing mechanisms give only constant reward instead of `Theta(sqrt(n))`.

## Leading route and exact sufficient lemma

Put $`q_n=2M_n`$, $`a_n=q_n/n^{3/2}`$, $`p=m/n`$, and $`p_2=(m)_2/(n)_2`$. Convergence follows if fixed $`\rho\in[1/2,1)`$, $`c>0`$, and `K` exist such that, for every large `n` and every $`m\in[\rho n,n)`$, some exact order-`n` minimizer `A` has an `m`-set `S` with

```math
Q(A[S])\le p^{3/2}q_n+Kn^{3/2-c}.
```

Different targets may use different roots and subsets. This yields $`a_m\le a_n+O_\rho(n^{-c})`$; fixed-ratio iterations have a geometric error sum, and exact target landing gives the adaptive tail (10.529), hence convergence.

The sharpest current implementation is (10.792)--(10.795). For any oriented full cut `d`, define

```math
\Delta_A(d)=q_n-\langle A,d\rangle,
\qquad
B_{n,m}=(p^{3/2}-p_2)q_n,
```

```math
\widehat\ell(S,d)
=Q(A[S])-c_A(S,d)-p_2\Delta_A(d)-B_{n,m}.
```

For fixed $`0<c<1/4`$, seek target-specific exact minimizers and cuts, with constants uniform over all sufficiently large `n` and every `m` in one fixed-ratio window, such that

```math
R_2(d)=O(n^{9/4-c}),
\qquad
U_m\{\widehat\ell(S,d)\le O(n^{3/2-c})\}
\ge\exp\{-O(n^{3/4-c})\}.
```

Conditioning `U_m` on this event gives selector KL `O(n^(3/4-c))` and a constant output. The arbitrary-cut slice mgf (10.759), with $`\lambda\asymp n^{-3/4}`$, then proves the displayed power-saving restriction lemma through (10.794).

This is strictly weaker than the Wave 24 target: it retains both `-p_2 Delta_A(d)` and the leading allowance `B_{n,m}`. For parent grounds its event is `ell<=B+t`, not `ell<=t`. Integer row-square extraction (10.788) shows any feasible captured mixture yields one such cut with only `O(log n)` additional log-loss, so two-ground mixing cannot rescue a stretched-exponential failure.

## Obstructions and falsification

Known obstructions are:

- for every fixed cut, $`\mathbb E_{U_m}\widehat\ell=\mathbb E Q(A[S])-p^{3/2}q_n`$, exactly the unknown restriction excess; scalar means and Markov are circular;
- (10.795) needs a reverse-tail/coverage lower bound for a fixed low-row-square cut, whereas the proved slice mgf gives the converse upper-tail direction under a pre-existing gap;
- row cost and coverage must be selected jointly; the universal average `E_d R_2(d)=n(n-1)` does not prevent all useful selector mass from lying on costly cuts;
- exact-ground switching has sign-indefinite restricted shores, and its continuous normal cone cannot encode discrete signing minimality without changed positive-deficit witnesses;
- the parent-ground Laplace/coarea theorem converts a soft premise to a hard column but does not prove the missing Laplace lower bound;
- every comparison must work uniformly for all requested target orders in one active ratio window with summable costs; one favorable root, subset, or landing order is insufficient.

For fixed proposed `c`, ratio window, and constants `K_R,K_t,K_u`, this implementation is falsified by infinitely many target pairs in that window for which every exact minimizer and every cut with $`R_2(d)\le K_Rn^{9/4-c}`$ has coverage below $`e^{-K_u n^{3/4-c}}`$ at threshold $`K_tn^{3/2-c}`$. Failure for all finite constants is equivalently expressed by maximal coverage $`e^{-\omega(n^{3/4-c})}`$ for every fixed `K_R,K_t`. This falsifies only that proposed window; defeating the existential implementation as a whole requires the analogous failure for every candidate `rho,c` and constants. It would still not falsify all optimized restriction. A direct power-saving falsifier is a fixed ratio and infinitely many `n` where every exact minimizer has minimum restriction excess at least $`c_0n^{3/2}`$. The adaptive comparison criterion itself fails exactly when $`\Omega_{\rm ad}(N)`$ from (10.529) does not tend to zero, equivalently when its limsup is positive (or infinite).

## Audit of constant shortfall (10.617)

The exact result needed is fixed $`\rho,\beta>0`$, finite `K`, and a path from each required exact root to every $`m\in[\rho n,n)`$ whose order-`r` deletion selects `i` with

```math
\kappa_{\beta,i}(B_r)
\ge \frac{q_n}{n^{3/2}}[r^{3/2}-(r-1)^{3/2}]-K.
```

Uniform temperature, shortfall, ratio window, and landing make (10.618) cost `O(n^(-1/2))` per geometric window. For convergence, “every minimizer” may be weakened to a target-specific exact root.

With $`\rho_\beta=\tanh(2\beta)`$ and signed even-cycle polynomial $`\mathcal P_B`$, the exact one-step theorem is the deletion ratio (10.780):

```math
\frac{\mathcal P_{B[-i]}(\rho_\beta)}{\mathcal P_B(\rho_\beta)}
\le \cosh(2\beta)^{r-1}e^{-\beta(T_r-K)}.
```

Exact minimality supplies only lower floors on character twists; coefficientwise bounds destroy the decisive cancellation, and sharp plaquette curvature yields only constant reward. A literal asymptotic falsifier is fixed `beta` and an unbounded sequence containing an exact minimizer with $`\max_i\kappa_{\beta,i}=o(\sqrt n)`$; the convergence-weakened version requires all eligible roots or every required deletion cutset to be bad.

## Ranked alternatives

1. **Parent-Gibbs full interpolation.** The sharper exact endpoint package is (10.799)--(10.800): separately control parent-cut entropy and mean adjacent finite-measure Hellinger distance. The `A_4` certificate defeats the two displayed coefficient-specific absorptions, but not a new constant-factor or asymptotic minimizer-specific comparison.
2. **Constant-shortfall signed-cycle route.** Direct and sharply falsifiable, but it needs minimizer-specific cancellation stability under deletion, not generic plaquette or coefficientwise bounds.
3. **Terminal external surplus with spatial/temporal service.** Exact Hall and replacement formulations exist; the missing input is a localized exact-minimizer surplus upper tail.
4. **Universal sharp lower bound or global thermodynamic compactness.** Either could settle the problem, but heavy-tail, scalar-pressure, fixed-replica, and finite-profile methods remain below the required resolution.

## Decision rule

Prioritize direct constructions or reverse-tail theorems for the single-cut target (10.795), including mechanisms that jointly create low `R_2` and centered selector coverage. Keep one attack on the separate parent-cut entropy/Hellinger package. Test constant shortfall only through signed-cycle deletion/cancellation or another genuinely complete-signing upper mechanism. Do not restart capped boosting, hard row regularity, two-ground rescue, scalar mean truncation, base variance, star-cover, or coefficientwise-cycle arguments as if their recorded walls were absent.
