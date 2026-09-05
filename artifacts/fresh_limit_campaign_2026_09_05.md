# Fresh limit campaign, 2026-09-05

Started at 18:29 UTC from commit 3d05fddeda74e4987608e57257ce564f0a331361.
The user authorized six hours of substantive research on convergence or
nonconvergence of M_n/n^(3/2), with freedom to replace prior architecture.
No active-time claim is inferred solely from elapsed wall time.

## Initial operational check

Local write access verified by creation of this research record. The local
Python environment is Python 3.9.2; git and rg are available. There are eight
CPU cores and approximately 66 GiB free disk. Lean and Lake are not on PATH.
Three independent researchers are assigned algebraic/combinatorial,
variational/structural, and primary-literature attacks. Their initial prompts
contain the original problem and numerical frontier, not the archived route
judgments. Existing untracked files from earlier sessions are preserved.

## First-principles notebook

Only the problem, rigorous interval, small exact values, and exact augmented
cut-code identity were read from ACTIVE_STATE before this analysis. The
original objective remains unsolved. Concrete derivations and their later
archive comparisons will be recorded below.

### Frozen director candidates before historical assessment

1. **A direct universal endpoint inequality.** Seek a universal bound
   Q(A) >= (1/2)n^(3/2)-o(n^(3/2)), possibly via the cut-Laplacian at a
   maximizing spin and its negative-energy counterpart. This is stronger
   than limit existence, but might admit a shorter proof. At a positive
   maximizer, switching makes every cut weight nonnegative; bounding its
   maximum by Q imposes a simultaneous upper cut constraint. The concrete
   first question is whether those two cut constraints and |a_ij|=1 force
   sum_i h_i >= n^(3/2)-o(n^(3/2)). No such implication is asserted.
2. **Scale transfer by nonlocal orthogonal lifting.** For an arbitrary good
   child A, seek sign realizations at larger orders whose Boolean cap is
   controlled by Q(A), rather than its operator norm. A naive A tensor H
   replaces scalar spins by coupled positive/negative eigenspace vectors.
   The exact test is whether this vector relaxation has a vanishing, rather
   than fixed, gap when its inner Hadamard dimension grows. A fixed gap kills
   this version. This is a candidate, not an imported construction theorem.
3. **Direct covering threshold.** Express the optimum as the deficit of
   the augmented cut code, and seek an all-order variational characterization
   of its worst-coset radius. Unlike ordinary radial density, any candidate
   must distinguish an empty covering hole from exponentially small positive
   mass. The concrete question is whether cut-code algebra supplies a
   nonperturbative inequality for the smallest covering multiplicity; an
   entropy limit alone does not answer it.

These are initial hypotheses, not progress claims. The team is separately
developing additional candidates. Historical comparison may now be used to
test the exact mechanisms above.

## First substantive findings (campaign continuing)

### Proved: polar-Gram tradeoff

For every hollow signing of order n>=2, with undoubled Boolean cap Q,

```math
Q(A)\ge\frac{n(n-1)}{\pi}
\arcsin\!\left(\frac{n}{\|A\|_*}\right).
```

The associated absolute vector relaxation satisfies
V(A)||A||_* >= n^2(n-1)/2. Its exact minimum n sqrt(n-1)/2 is attained
precisely by conference matrices; a relative epsilon excess forces mean
squared spectral distance from the two-point spectrum at most
2 epsilon/(1+epsilon). The polar-Gram proof was derived by the director
and independently checked by the variational and literature researchers.
It is not an assertion that Boolean minimizers minimize the relaxation.

For Q(A_n)/n^(3/2)<=c+o(1), the Boolean theorem implies
||A_n||_*/(n sqrt(n-1)) >= 1/(pi c)-o(1). This supplies a quantitative
constraint, not a better universal lower constant than the existing bound.
Full proof: [variational artifact](fresh_limit_variational_2026_09_05.md).

### Proved: a tensor-stable minimum has a limit

The algebra researcher constructed explicit symmetric regular Hadamards of
orders 4 and 144. For their tensor semigroup, the same-spin regularization

```math
\mathcal R(B)=\sup_{s=4^a144^b}
\frac{Q(H_s\otimes B)}{s^{3/2}}
```

is genuinely tensor stable. The director's full-sign-seed observation
removes the leading completion payment: allow diagonal entries +/-1 in
the finite seed and erase only the final parent's diagonal. Thus

```math
\limsup_n\frac{M_n}{n^{3/2}}
\le c_{\mathcal R}
=\inf_{B\text{ full symmetric sign}}
\frac{\mathcal R(B)}{(\operatorname{ord}B)^{3/2}},
\qquad
\frac1{\sqrt{2\pi}}\le c_{\mathcal R}\le\frac12.
```

The minimum regularized ratio itself converges to c_R. This is **not** a
lower bound on the original liminf. A finite full seed with a rigorous
regularized ratio below 1/2 would improve the original asymptotic upper
bound; none has been certified. The natural PSD-majorant upper certificate
has an exact half-scale floor for full sign seeds, so it cannot certify
such an improvement. Proof and exact tests are in the
[algebra artifact](fresh_limit_algebra_2026_09_05.md).

One quantifier correction emerged in the director audit: for a *fixed*
seed order, the finite minimum over seeds does commute with the increasing
directed supremum over outer factors. The unproved interchange concerns
unbounded seed orders, or removing the tensor-family restriction.

### Exact/numerical checks and research decision

The independent `computations/audit_fresh_nuclear_2026_09_05.py` audit
enumerates all 2,131,019 root-gauged signings through order eight.
Energy and bilinear-cap calculations are exact integers; eigensolver and
arcsine comparisons are floating-point diagnostics of the analytic proof.
It also verifies the archived order-five cut-Laplacian obstruction using
the exact negative principal determinant -1, without trusting eigenvalues.

The checked coding, Hilbertian, and synchronous-XOR sources do not give an
imported convergence theorem. In particular, independent row/column
Hadamard excess is not same-switch quadratic excess. See the
[source and adversarial audit](fresh_limit_literature_2026_09_05.md).

The original interval is unchanged. The team continues with the actual
regularized Boolean norm and its realizability constraints, while checking
whether the spectral tradeoff supports a quantitatively useful Boolean
argument. Neither auxiliary result is counted as solving convergence.
