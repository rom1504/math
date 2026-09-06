# Independent supersolution audit and orthogonal entropy theorem

Date: 2026-09-06. Status: elementary theorem below proved; random finite
tests are diagnostics only. The original strict upper bound survives this
audit of its conditional-variance supersolution. This file does not audit
the separate rational numerical certificate or all-order weave realization.

**Archive comparison, added after the independent derivation:** the g_t
specialization of Sections 2--4 is already proved in
`transfer_reconstruction_variance_kernel_factorization_2026_09_06.md` and
`transfer_seed_multivariate_schur_interface_2026_09_06.md`. The present
general-g formulation and principal-minor proof are an independent
abstraction/reconstruction, not a new signing theorem. Those artifacts
also remove the raw averaged-marginal issue below through the exact
mixture-index information cost; safe input permutations are unnecessary
once that additional inequality is used.

## 1. Selected falsification target

The proposed weak point was the use of sequential correlated labels with
non-Gaussian conditional source laws in
`continued_convergence_conditional_variance_supersolution_2026_09_06.md`.
The exact local inequality is valid. No Gaussian assumption about these
conditional laws is present: conditional linear prediction bounds their
average residual variance in the correct direction, and ordinary convexity
then permits averaging the nonlinear reward.

In fact the proof extends to every finite dimension and every orthogonal
mixing matrix. This provides a general information/variance principle,
although not a characterization of the original signing optimum.

## 2. General theorem

Let g:[0,infinity)->R be continuous, decreasing, convex, with
`s -> g(exp(s))` concave on R. For any finite real-valued source X, set

```math
T_g(X)=\sup_L\{\mathbb E g(\operatorname{Var}(X\mid L))-I(X;L)\}.
```

Finite labels suffice: writing the channel as a mixture of posterior
probability vectors reduces the objective to the concave envelope of a
continuous function on the finite source simplex, with a fixed barycentre.
Caratheodory gives a finite mixture.

**Orthogonal entropy theorem.** Let A=(A_1,...,A_d) have finite support,
let O be a real orthogonal d-by-d matrix, and set U=OA. Then

```math
\sum_{i=1}^d T_g(U_i)-\operatorname{TC}(A)
\le\sum_{i=1}^d T_g(A_i),
\qquad
\operatorname{TC}(A)=\sum_i H(A_i)-H(A).
```

The input marginals need not agree or be symmetric. In particular, if
they are all nu, then the average child envelope minus TC(A)/d is at most
T_g(nu). This is a supersolution theorem for equal-marginal orthogonal
Bellman operators.

**Scope caution.** For a general O this does NOT automatically cover a
different operator imposing only an averaged input-marginal constraint.
The original two-child Hadamard operator has an additional safe
reversal/swap symmetrization that makes its two signed marginals equal
without changing either child absolute law. An arbitrary orthogonal
operator need not possess that symmetry.

## 3. Covariance lemma in arbitrary dimension

For a positive definite covariance matrix Sigma, let d_i be its sequential
linear-regression residual variances in order 1,...,d, and let lambda_i be
its eigenvalues. Thus d_i are squared Cholesky diagonal entries and
`prod_i d_i=det Sigma=prod_i lambda_i`.

For any subset S of k coordinates, regress its coordinates only on the
earlier coordinates within S. Removing regressors can only increase each
residual variance. Their product is det(Sigma_SS). Therefore

```math
\prod_{i\in S}d_i\le\det(\Sigma_{SS})
\le\prod_{j=1}^k\lambda_j^{\downarrow}.
```

The last inequality follows from the variational eigenvalue principle for
a coordinate compression. Maximizing the left side over S proves that
log(d_i) is majorized by log(lambda_i), with equality of total sums.
Concavity of `g(exp(s))` consequently gives

```math
\sum_i g(d_i)\ge\sum_i g(\lambda_i).
```

If u_i are the diagonal entries of O Sigma O^T, then
`u_i=sum_j O'_(ij)^2 lambda_j` for another orthogonal matrix O'. Its
squared-entry matrix is doubly stochastic. Ordinary Jensen gives

```math
\sum_i g(\lambda_i)\ge\sum_i g(u_i).
```

Together these prove `sum_i g(d_i)>=sum_i g(u_i)`.
Singular covariance matrices follow by adding epsilon I and taking
epsilon down to zero. The Cholesky pivots converge to the corresponding
least-squares residual variances (using a pseudoinverse when needed);
cross-covariances lie in the range of each preceding covariance block.
Continuity of g at zero then gives the same conclusion.

## 4. Sequential channels: full proof

Choose arbitrary child channels M_i|U_i independently conditional on U,
and let M=(M_1,...,M_d). Give A_i the parent label `(M,A_1,...,A_(i-1))`.
The entropy chain rule gives exactly

```math
\sum_i I(A_i;M,A_{<i})=I(A;M)+\operatorname{TC}(A).
```

Since O is invertible, I(A;M)=I(U;M). Conditional independence of the
child channels gives

```math
I(U;M)=\sum_i I(U_i;M_i)-\operatorname{TC}(M)
\le\sum_i I(U_i;M_i).
```

Fix a value m of M and let Sigma be the covariance of A conditional on
M=m. The optimal nonlinear prediction of A_i from A_<i is at least as
good as the best linear prediction, hence

```math
\mathbb E[\operatorname{Var}(A_i\mid M,A_{<i})\mid M=m]\le d_i.
```

Convexity and monotone decrease of g imply that the averaged parent
reward at m is at least sum_i g(d_i). Section 3 bounds this by
sum_i g(Var(U_i|M=m)). Finally total variance, monotone decrease, and
Jensen give the label-refinement inequality

```math
\mathbb E g(\operatorname{Var}(U_i\mid M))
\ge\mathbb E g(\operatorname{Var}(U_i\mid M_i)).
```

Combining reward and information comparisons shows that the sum of the
parent feasible-channel values is at least the sum of the chosen child
values minus TC(A). Each parent value is at most T_g(A_i). Optimizing
the child channels proves the theorem.

The exact discarded information slack is TC(M), not an assumed
independence of M_i after marginalizing U. The M_i are generally correlated.
The conditional laws of A can be arbitrary non-Gaussian finite laws and
can have nonlinear deterministic dependencies. Neither affects the proof.

## 5. Gaussian self-transport potential

For the potential in the strict upper construction,

```math
g_t(v)=-tv(1-\rho)+\tfrac14\log(1-\rho^2),\qquad
2tv={\rho\over1-\rho^2},\quad 0\le\rho<1.
```

The envelope derivative is `g'_t(v)=-t(1-rho)`. Thus g is decreasing
and convex because rho increases with v. Also

```math
{d\over d\log v}g_t(v)=vg'_t(v)=-{\rho\over2(1+\rho)},
```

which decreases with v. Hence the logarithmic concavity hypothesis holds.
Applying the theorem with d=2 and the Hadamard rotation recovers exactly
the archived supersolution, with no equal-temperature assumption.

Both curvature requirements have real content. For example
`g(v)=exp(-v)-1` is decreasing and convex but is not log-concave globally.
Take covariance eigenvalues 1 and 100 and choose its second diagonal
entry as 10 (the other is 91, off-diagonal squared 810). The sequential
pivots are 10,10. Then

```math
2g(10)<g(1)+g(100),
```

so the local covariance lemma fails if logarithmic concavity is dropped.
This does not assert a counterexample to the optimized envelope theorem
for that g; it identifies the exact point at which this proof would fail.

## 6. Finite diagnostics

`computations/decisive_bridge_supersolution_audit_2026_09_06.py` generates
480 finite non-Gaussian examples in dimensions 2,3,4, with random child
channels, highly uneven atom weights, ordinary Hadamard and random
orthogonal mixing, nonlinear deterministic input dependencies, and
singular covariance cases. It computes actual conditional-variance
rewards and actual finite mutual informations, not Gaussian substitutions.

Seed: 202609067. Minimum reward slack: 0.002103268724010161.
Minimum information slack: -1.11e-15 (floating roundoff).
Minimum total slack: 0.01705508327669847. The full raw diagnostics are
in the same-prefix JSON output. These checks do not prove the theorem;
Sections 2--5 do.

## 7. Gaussian-boundary mechanism audit

I read the complete terminal-gap and unbounded-rigidity artifacts. The
essential chain is sound on reconstruction:

1. Compact Gaussian-kernel self-transport concavity is extended using
   locally uniformly bounded scaling measures. Gaussian tails permit
   locally uniform kernel convergence. The resulting potential has a
   quadratic upper bound and is integrable for every finite-second-moment
   source. The entropy comparison then identifies the actual optimizer.
2. Strict Gaussian-energy positivity gives strict source concavity of
   F_t. Combined with tensorization, equality at a dependent joint source
   is impossible: the product source is a maximum under fixed marginals,
   and strict concavity excludes a second maximum.
3. For an equal-marginal Bellman pair with mutual information I, the
   self-coupling entropy identity gives `F_t(pair)>=2F_t(nu)-I`.
   Tensorization after rotation then gives Bellman drift at least I/4.
   A zero-drift maximizing pair is therefore iid. Equality of output
   tensorization makes sum and difference independent, forcing a Gaussian
   source by the characteristic-function identity.
4. The logarithmic Gaussian cost bound controls rare high-energy source
   tails, giving weak continuity on bounded-second-moment sets despite
   possible moment escape. Compactness and zero-drift rigidity give a
   positive drift away from a small terminal-envelope gap.
5. The depth-r stopped-tree argument uses an exact nonnegative moment
   martingale. Its high-energy loss is bounded by
   `E_0 sup_(s>C) K(s)/s`, which vanishes as C grows. Its surviving-branch
   probability is at most `D_0/(r kappa)`. No uniform controlled CLT is
   inserted. The increasing lower iterates give a genuine fixed point H;
   their lower bound by the older E envelope is sufficient for stopping.
6. This establishes `B^r Phi -> H <- B^r G`. The finite-source T_g
   supersolution lies above each lower iterate; hence H<=T_g. Extending
   T_g to unbounded laws is unnecessary for a finite root and finite tree.

The compact primary input was subsequently checked directly: Feydy et al.,
https://proceedings.mlr.press/v89/feydy19a/feydy19a.pdf, Section 2.2,
Propositions 3--4, and the associated supplement Sections B.3--B.4,
https://proceedings.mlr.press/v89/feydy19a/feydy19a-supp.pdf. With their
epsilon=1 and C(x,y)=t|x-y|², their negentropy F_epsilon is our Phi_t,
not our positive F_t. On compact subsets this cost is Lipschitz and its
Gaussian kernel is positive universal, as required. Their representation
has the same auxiliary-measure energy and -1/2 normalization. Their strict
convexity proof uses strict Gaussian-kernel energy plus joint relative-
entropy convexity exactly as reconstructed. The unbounded extension is
still supplied by the repository argument, not by claiming the paper's
compact hypotheses cover unbounded supports. The separate rational root
certificate was not rerun here. No gap in the supersolution or boundary
mechanism was identified.
