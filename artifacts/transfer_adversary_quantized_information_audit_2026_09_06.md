# Independent audit: quantized Gaussian information and weighted resampling

Date: 2026-09-06, approximately 19:22 UTC. Status: **PASS**, with one
minor explicit-hypothesis clarification, for
`transfer_director_quantized_gaussian_information_2026_09_06.md`.
The pressure formula should state `beta>=0`, as is customary for the
inverse temperature and used in its displayed error bound. This was
reported to the owner. The finite information theorem permits either
sign of `s` under its stated norm restriction.

## 1. Finite quantizer operator and its exact matrix

Let one quantizer have `q` cells `C_a`, with standard Gaussian masses
`p_a>0`. Put `e_a=1_(C_a)/sqrt(p_a)` and
`v=(sqrt(p_1),...,sqrt(p_q))`. These cell indicators are orthonormal.
For a bivariate standard Gaussian pair of correlation `r`, define

```math
J_{ab}=\Pr\{G\in C_a,G'\in C_b\},\qquad
M_{ab}=\frac{J_{ab}}{\sqrt{p_ap_b}}-\sqrt{p_ap_b}.
```

The matrix `M` annihilates `v` and is positive semidefinite. On `v`'s
orthogonal complement, it represents the compression of `T_r` to the
centered quantizer subspace. The nonzero eigenvalues of
`T_sqrt(r)(P-Pi)T_sqrt(r)` equal those of this compression: apply the
standard equality of nonzero spectra of `UU*` and `U*U` to
`U=T_sqrt(r)(P-Pi)`. Therefore the source's `kappa` is exactly the
largest eigenvalue of `M`, with at most `q-1` nonzero eigenvalues.
It is not a Boolean response statistic of the ambient matrix.

This also verifies the required positivity, not merely an absolute
operator-norm bound. Constants are invariant and centered functions
stay centered. Thus each smoothed quantizer projection is bounded by
`Pi+kappa(I-Pi)` in positive-semidefinite operator order. This order
tensors as in the independently audited binary theorem. In fact
`kappa<=r<1` by the centered OU contraction, though the stated weaker
interval `[0,1]` is enough for the subset law.

## 2. Likelihood domain and random principal determinants

`|s|C<r` places every eigenvalue of `I+(s/r)B` in `(0,2)`, guaranteeing
both positive covariance and square-integrable Gaussian likelihood.
OU self-adjointness gives the density smoothing identity in the correct
direction. Conditioning that density on all quantizer cells produces
the actual output likelihood relative to its product marginals.

The hollow assumption matters: every Gaussian coordinate still has
variance one, so those actual marginals equal the standard Gaussian
cell laws. Tensor expansion of the bounding operators averages
orthogonal projections over independent selectors with probabilities
`kappa_i`. Integrating away unselected coordinates gives a principal
Gaussian marginal covariance, not a conditional Schur complement.
The exact squared likelihood norm is
`det(I-(s/r)^2 B_S^2)^(-1/2)`. Jensen under the output law gives the
KL/Renyi comparison with the additive `1+chi^2` in the correct place.

All of these claims hold for the nonflat symmetric hollow matrices in
the source. No eigenvector or exchangeability assumption is introduced.

## 3. Weighted fractional Holder step

Write `w_ij=B_ij^2` and `D_*=max_i sum_(j!=i) w_ij`. For `D_*>0`,
each positive edge has `theta_ij=w_ij/D_*` in `(0,1]`, and the sum
of incident exponents at every coordinate is at most one.

The sequential-integration proof of the displayed fractional Holder
inequality is valid. At a given coordinate, hold all other coordinates
fixed, apply Holder to the incident factors with exponents `1/theta_e`,
and insert a constant-one factor for any missing exponent mass.
Each incident function becomes its integral over that coordinate and
keeps the same outer power `theta_e`. The incidence bounds remain valid
at the next coordinate. After all coordinates are integrated, each
function has become its full product-measure expectation. Independence
of the selectors is essential here and is explicitly supplied.

The determinant logarithm is bounded by

```math
c\operatorname{tr}B_S^2
=2c\sum_{i<j}w_{ij}Z_iZ_j,
\qquad c=\frac{(s/r)^2}{2(1-(s/r)^2C^2)}.
```

Thus choosing `f_ij=exp(2cD_* Z_iZ_j)` and raising it to `theta_ij`
gives the correct original exponential. The factor two is necessary
and is present. Independence gives
`E f_ij=1+kappa_i kappa_j(exp(2cD_*)-1)`, proving the weighted log
bound. Finally `sum_(i<j)w_ij=tr B^2/2`, so the homogeneous-channel
coefficient `tr B^2/(2D_*)` is also correct. If `D_*=0`, symmetry and
hollowness force `B=0`; the separate zero-divergence case is valid.

## 4. Threshold dependence and pressure normalization

For binary thresholds, the one-dimensional channel operator has rank
one. The bivariate Gaussian density derivative yields exactly the
source's integral for `kappa`. Bounding its exponential by
`exp(z^2u)` and its denominator by `sqrt(1-r^2)` gives equation (5),
including the continuous interpretation at `z=0`.

For a flat signing, `D_*=(N-1)/N` and `tr B^2=N-1`. The weighted
entropy term is at most `N kappa^2(exp(2c)-1)/2`. This is a valid
slightly weaker form of the preceding weighted inequality.

After switching an actual maximizing spin and choosing its energy
sign, the threshold trial law has baseline energy `m^2 Q(A)`.
Parseval gives total centered threshold coefficient mass `v`, of which
the first Hermite contributes `a^2`. The remainder on an edge is
therefore at most `(v-a^2)s^2/N` in absolute value. For `beta>=0`,
multiplication by `beta/sqrt(N)` and summation over undirected edges
gives precisely

```math
\frac{\beta(v-a^2)s^2}{2}\frac{N-1}{\sqrt N}.
```

The linear contribution is `beta a^2s(N-1)/2`. Fixed marginals give
`H(X)=N h(delta)-D` exactly, and the chosen partition-function phase
is bounded by the sum of both phases. These observations establish
equation (6) with no missing trace term, phase factor, or factor two.
The physical parameters are fixed before the order limit.

The nonflat bounds depend explicitly on the threshold probability.
At fixed nonzero correlation, their displayed factor multiplying
`eta^2` diverges as the threshold becomes rare. They do not contradict
the earlier independent-pair obstruction to a coefficient uniform in
the rare probability. Neither theorem prohibits all threshold-dependent
nonflat information estimates.

## 5. Concrete guard against an unjustified cap/convergence implication

The information theorem and its deterministic pressure LOWER bound do
not, by themselves, imply any strict-subhalf cap upper bound. This
can be falsified on an actual eligible family, without replacing the
signing problem by an abstract scalar sequence.

Let `H=J_4-2I_4`, `W_k=H^(tensor k)`, `N=4^k`, and hollow it by
`A_k=W_k-(-1)^k I`. Then `W_k^2=N I`. The all-ones vector is a
Boolean eigenvector of eigenvalue `sqrt(N)`. A zero-sum Boolean
eigenvector in one tensor factor, tensored with all ones elsewhere,
has eigenvalue `-sqrt(N)`. Thus

```math
Q(A_k)=\frac{N\sqrt N+N}{2},\qquad
\|A_k\|_{op}=\sqrt N+1.
```

For `k>=5`, the normalized operator norm is at most `33/32`, so this
family satisfies even the narrower norm hypothesis of the seed's
factor-22 theorem. Both positive entropy/pressure credits apply, yet
its normalized cap is `1/2+1/(2sqrt N)` and tends to one half.

Accordingly, substituting an annealed partition upper budget as though
it held for every individual eligible signing would be invalid. An
actual upper certificate and its correctly scoped extraction are
essential. This example does not refute the legitimate extracted
strict-upper construction, and does not prove nonconvergence of the
original minima. It rules out the stronger, unjustified inference
that the new information inequality alone forces cap compression or
identifies a limiting optimal coefficient.

The theorem has no comparison of a given minimizing seed with a
different order. A positive original convergence mechanism remains a
separate obligation, as the source correctly states.

## Preservation

This audit used only read-only source inspection and symbolic
reconstruction. It generated no temporary outputs and changed no
source theorem. This audit file is the sole new file from the task.
