# Conditional-variance kernels and arbitrary orthogonal gate factorization

Date: 2026-09-06. Structural consequences of the standalone reconstruction.
These are identities and inequalities for the source potential, not an
arbitrary-signing seed-transfer theorem.

For `t>0`, use

```math
g_t(v)=-tv(1-\rho)+\frac14\log(1-\rho^2),\quad
2tv=\frac\rho{1-\rho^2},\quad
T_t(\nu)=\sup_L\{\mathbb Eg_t(\operatorname{Var}(X\mid L))-I(X;L)\}.
```

Labels may be arbitrary standard-Borel random variables, and sources have
finite second moment. The finite-source definition is unchanged by requiring
finite labels. Write `K_t=-g_t`.

## 1. Exact variable-precision rate-distortion representation

The supporting lines of `g_t` are

```math
g_t(v)=\sup_{0<\lambda\le t}\{c_t(\lambda)-\lambda v\},\qquad
c_t(\lambda)=\frac14\log\frac{\lambda(2t-\lambda)}{t^2}.
\tag{1}
```

For `v>0` the optimizer is `lambda=t(1-rho)`; for `v=0` it is `t`.
Consequently

```math
-T_t(\nu)=\inf_{X\to(Y,\Lambda)}
\left\{I(X;Y,\Lambda)
+\mathbb E[\Lambda(X-Y)^2-c_t(\Lambda)]\right\},
\quad 0<\Lambda\le t.
\tag{2}
```

Proof in both directions: from a label `L`, choose its posterior mean as
`Y` and the optimizing value in (1) at its posterior variance as `Lambda`.
Replacing `L` by `(Y,Lambda)` decreases information, and the expected
distortion is exactly `E K_t(Var(X|L))`. Conversely, for any `(Y,Lambda)`
channel, conditional mean-square decomposition and (1) give

```math
\mathbb E[\Lambda(X-Y)^2-c_t(\Lambda)\mid Y,\Lambda]
\ge K_t(\operatorname{Var}(X\mid Y,\Lambda)).
```

Thus neither an equal-precision hypothesis nor a Gaussian posterior has
entered (2).

Eliminate the forward channel by its marginal reproduction distribution.
For `d(x;y,lambda)=lambda(x-y)^2-c_t(lambda)`, the entropy variational
identity gives

```math
T_t(\nu)=\sup_{Q\in\mathcal P(\mathbb R\times(0,t])}
\int\log\left[\int e^{c_t(\lambda)-\lambda(x-y)^2}
                         Q(dy,d\lambda)\right]\nu(dx).
\tag{3}
```

Indeed `I(X;Z)=inf_Q D(P_{XZ}||nu tensor Q)`, while for fixed `Q` the
pointwise Gibbs channel minimizes `D(P_{Z|x}||Q)+E[d(x;Z)|x]` at
`-log integral e^{-d(x;z)}Q(dz)`. Both operations are infima, so no minimax
exchange is required. Approximation handles infinite costs. The kernel is
positive and at most one, and the trivial choice `Y=0,Lambda=t` gives a
finite objective for every finite-second-moment source.

Formula (3) makes `T_t` source-convex directly. It also identifies the
difference from a common-temperature envelope: `Q` may mix different
precisions in one channel. This is a Gaussian-location/precision kernel
mixture, not a claim that the source or posterior is Gaussian.

For a source supported on an interval of diameter `D`, one may restrict
`Y` to its convex hull and `Lambda` to `[lambda_*,t]`, where `lambda_*`
is the optimizer of (1) at `v=D^2/4`. This follows from the first direction
of the proof of (2), since every posterior variance is at most `D^2/4`.

### A useful finite-source dual bound

Let `nu=sum_i p_i delta_{x_i}`. Any positive numbers `a_i` satisfy

```math
T_t(\nu)\le -\sum_i p_i\log a_i
+\log\sup_{y,\lambda}\sum_i p_i a_i
                    e^{c_t(\lambda)-\lambda(x_i-y)^2}.
\tag{4}
```

Proof: multiply the inner kernel in (3) by `a_i`, apply logarithmic Jensen
over `i`, and then bound its `Q` average by its pointwise supremum. The
parameter range can be restricted as above. Inequality (4) is an upper
certificate for the full posterior envelope, not a selected-channel lower
test. It may replace a two-dimensional posterior grid by a Gaussian-mixture
dual for a fixed source, but no stronger numerical bound is claimed here.

For the ternary symmetric source, choosing equal dual weights on its two
nonzero atoms reduces the inner supremum to

```math
e^{c_t(\lambda)-\lambda y^2}
\left[(1-p)a_0+p a_1e^{-\lambda/p}
          \cosh\left(\frac{2\lambda y}{\sqrt p}\right)\right].
```

Thus only two positive dual weights and a compact two-dimensional smooth
kernel maximum remain. Their certification is a separate computational
task; the existing exact root certificate remains the proved value.

## 2. Orthogonal mixing with heterogeneous parent marginals

Let `A=(A_1,...,A_d)` have finite second moment and finite total correlation

```math
\operatorname{TC}(A)=D(P_A\Vert\bigotimes_iP_{A_i}).
```

For every orthogonal `O`, put `U=OA`. Then

```math
\boxed{\ \sum_{j=1}^d T_t(P_{U_j})
\le\sum_{i=1}^d T_t(P_{A_i})+\operatorname{TC}(A).\ }
\tag{5}
```

No equal-parent-marginal assumption is present.

Choose independent child channels `M_j|U_j` conditional on `U`, and put
`L=(M_1,...,M_d)`. Label parent coordinate `A_i` by `(L,A_1,...,A_{i-1})`.
The information identity and product-channel inequality are

```math
\sum_i I(A_i;L,A_{<i})
=I(A;L)+\operatorname{TC}(A)
\le\sum_j I(U_j;M_j)+\operatorname{TC}(A).
\tag{6}
```

Condition on `L`. For its positive definite covariance `Sigma`, let
`delta_i` be the successive linear-regression innovation variances, namely
the squared diagonal entries of its Cholesky factor. They satisfy

```math
(\log\delta_i)_{i=1}^d\prec(\log\lambda_i(\Sigma))_{i=1}^d.
\tag{7}
```

An elementary proof uses exterior powers of the Cholesky factor `R`.
Any product of `k` diagonal entries of triangular `R` is an eigenvalue
of `wedge^k R`, hence its modulus is at most the product of the largest
`k` singular values of `R`. Square and choose the largest `k` pivots.
The total products agree by the determinant identity. These are exactly
the partial-sum inequalities in (7).

Conditional linear regression, Jensen, concavity of `g_t(exp(s))`, and
ordinary convexity of `g_t` now give

```math
\sum_i\mathbb E_{A_{<i}\mid L}g_t(\operatorname{Var}(A_i\mid A_{<i},L))
\ge\sum_i g_t(\delta_i)
\ge\sum_i g_t(\lambda_i(\Sigma))
\ge\sum_j g_t((O\Sigma O^T)_{jj}).
\tag{8}
```

For the final comparison, the diagonal vector is a doubly stochastic
average of the eigenvalue vector. Singular covariances follow by applying
the comparison to `Sigma+epsilon I` and sending `epsilon` to zero;
the limiting pivots are the linear-prediction residual variances.
Refining child labels from `M_j` to `L` can only increase the expected
`g_t` reward. Average (8), subtract (6), and optimize the child channels.
The resulting parent channels are admissible for the corresponding
`T_t(P_{A_i})`, proving (5).

## 3. Raw type Bellman supersolution for every orthogonal gate

The exact mixture-label inequality is

```math
dT_t(\nu)\ge\sum_{i=1}^d T_t(\nu_i)
-\sum_{i=1}^d D(\nu_i\Vert\nu),\qquad
\nu=\frac1d\sum_i\nu_i.
\tag{9}
```

To prove it, choose the mixture index uniformly and disclose it together
with an almost-optimal child label. Its information cost is exactly
`(1/d)sum_i D(nu_i||nu)` plus the averaged child information, while its
posterior rewards average without loss. This remains valid for unbounded
sources and is stronger than paying the entropy `log d` of the index.

For a symmetric source `nu`, let `pi` range over laws of a `d`-vector
whose average absolute marginal is the absolute law of `nu`. Let `nu_j'`
be the symmetrized laws of `(OA)_j`. Define

```math
(\mathcal B_O f)(\nu)=\sup_\pi\left\{
\frac1d\sum_jf(\nu_j')-\frac1dD(\pi\Vert\nu^{\otimes d})\right\}.
\tag{10}
```

Then

```math
\mathcal B_O T_t\le T_t
\quad\text{for every real orthogonal }O.
\tag{11}
```

Indeed simultaneous global reversal preserves each output absolute law
and decreases relative entropy, so assume that `pi` is globally symmetric.
Its individual signed marginals `nu_i` then average to `nu`; each signed
output is already symmetric. Use (5), the exact decomposition

```math
D(\pi\Vert\nu^{\otimes d})
=\operatorname{TC}(A)+\sum_iD(\nu_i\Vert\nu),
```

and (9). This proves (11) without averaging input permutations. Therefore
there is no hidden need for a transitive group of gate symmetries fixing
the individual child absolute laws.

For a flat orthogonal gate `O=H_d/sqrt(d)`, the block recursion
`diag(U_q^(1),...,U_q^(d))(O tensor I_q)g` again consists of normalized
Hadamards. A `d`-tuple table under a fresh signed permutation has logarithmic
probability `-s D(pi||nu^d)/d+O(log s)`. Thus (11) applies to exact
fixed-depth recursive ensembles built using any real Hadamard gate, not
only binary Walsh gates. For a general orthogonal gate without flat
entries, (11) remains analytic but does not produce a sign matrix.

## 4. The terminal drift and Gaussian rigidity also extend to Hadamard gates

Use the self-transport `F_t` and `Phi_t=-F_t/2` from the standalone proof.
The source-mixture upper bound can charge the actual information of the
mixture index, not its full entropy:

```math
F_t(\nu)\le\frac1d\sum_i F_t(\nu_i)
+\frac1d\sum_iD(\nu_i\Vert\nu),\qquad
\nu=\frac1d\sum_i\nu_i.
\tag{12}
```

Indeed conditionally self-couple given the mixture index `R`; the coupling
information is at most `I(X;R)+I(X;X'|R)`. This proves (12) directly.

For a vector self-coupling of `A` of finite relative entropy, the exact
multivariate information identity is

```math
I(A;A')-\sum_iI(A_i;A_i')
=\operatorname{TC}((A_i,A_i')_{i=1}^d)-2\operatorname{TC}(A)
\ge-\operatorname{TC}(A).
```

The last inequality is coordinatewise data processing of total correlation.
Adding costs and using (12) gives, for an admissible globally symmetric
raw type `pi`,

```math
F_t(\pi)\ge\sum_iF_t(\nu_i)-\operatorname{TC}(A)
\ge dF_t(\nu)-D(\pi\Vert\nu^{\otimes d}).
```

Orthogonal invariance and self-transport subadditivity over output
coordinates therefore give the quantitative policy inequality

```math
\Phi_t(\nu)-\left\{\frac1d\sum_j\Phi_t(\nu_j')
-\frac1dD(\pi\Vert\nu^{\otimes d})\right\}
\ge\frac1{2d}D(\pi\Vert\nu^{\otimes d}).
\tag{13}
```

For `d=2` this is the original quarter-information drift, now established
without swapping inputs or first forcing their individual marginals to
coincide. Global reversal alone was sufficient.

Suppose now `O=H_d/sqrt(d)` is a real Hadamard gate of order `d>1`.
Compactness of fixed-moment pair-policy sets, weak continuity of `Phi_t`,
and lower semicontinuity of relative entropy yield an attained maximizing
policy. Zero drift in (13) forces its inputs iid with law `nu`. Symmetry
of `nu` makes every output have the law of a normalized sum of `d` iid
copies. Equality then forces equality in transport subadditivity; strict
source concavity makes the outputs mutually independent.

Choose two distinct Hadamard rows. After multiplying each input by a fixed
sign, the first row is all positive and the second row is positive on
`d/2` coordinates and negative on `d/2`. Let `S,T` be the two independent
group sums. The corresponding outputs are `(S+T)/sqrt(d)` and
`(S-T)/sqrt(d)`, so the elementary independent sum/difference argument
shows that `S` is Gaussian. If `phi` is the real characteristic function
of the original symmetric source, this says

```math
\phi(z)^{d/2}=\exp(-dvz^2/4).
```

It is nonzero everywhere and equals one at zero. Continuity forces
`phi>0`, hence `phi(z)=exp(-vz^2/2)`. Conversely iid centered Gaussians
are fixed with zero cost. Thus the zero-drift set is again exactly the
centered Gaussians, including the zero point mass.

The direct-`T` stopped-tree proof consequently applies with `d` children
of equal branch weight. Their average second moment is still the parent
moment by orthogonality. The usual iid central-limit policy gives the
same Gaussian lower boundary. It follows that, for every fixed real
Hadamard gate of order greater than one,

```math
\lim_{r\to\infty}\mathcal B_O^r\Phi_t(\nu)\le T_t(\nu).
\tag{14}
```

Together with the `d`-tuple count following (11), arbitrary terminal
Hadamards, and the same fibre weave, the existing negative root
certificate gives the same strict all-order upper cap through any such
fixed gate. This is an architecture-independence statement for the proved
cap, not an improvement of its numerical constant.

## 5. What this removes, and what it does not

The potential no longer needs a privileged binary gate or equal parent
marginals. Its source dependence can be expressed through positive
Gaussian-mixture kernels and certified by (4). These are genuine
factorization and representation statements.

An arbitrary hollow signing divided by the square root of its order is
not generally orthogonal. Substituting it into (10) is not justified by
(11); the restricted-weave energy identity also needs orthogonal fibre
rows. Neither a small scalar quadratic cap nor a finite diagnostic implies
those hypotheses. A seed-transfer theorem still needs an additional
argument controlling this nonorthogonal defect or a different exact
construction that avoids it.
