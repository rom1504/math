# Preserved live derivations and failures

2026-09-17. **Draft archive; claims are not promoted merely by appearing
here.** Canonical audited results are linked from the campaign record.

## A. Directional information lost by orbit statistics

For `mu_A=law((a_ij xi_i xi_j)_(i<j))`, coordinate multiplication by A
is an isometry for every coordinatewise absolute metric. Thus Cauchy
rate--distortion, ordinary Gaussian-channel information, and the fixed-law
Bernoulli/Gaussian widths cannot distinguish the signings. Gauge-equivariant
optimal coupling even identifies their common Bernoulli width with the
ordinary random-sign one-sided ground-state average (and with the absolute
average after including global coefficient reversal).

This defeats a proposed direct use of the Liu--Zadik rate integral to
recognize good seeds. It does not defeat its offset-sensitive theorem.
The required original direction is the fixed all-ones coefficient vector,
or an equivalent directional query. Replacing it by all rotationally or
coordinate-sign invariant data is an information loss, not a proof that
optimizers lack structure.

An attempted Cauchy area formula with edge-noise signs
`a_ij xi_i xi_j` would have endpoint Q(A), by switching invariance.
However the leave-one-coordinate-out proof uses independent coordinate
noise. Triangle constraints make these noise signs dependent. A vertex
rank-one Cauchy perturbation can even have a singular low-dimensional
support, revealing a planted finite input essentially exactly. No usable
Bayesian-risk comparison follows from the endpoint identity alone.

## B. Noise semigroup experiment

For the actual edge-cube cap function Q, let
`F_A(t)=E Q(A multiplied coordinatewise by independent signs of mean t)`.
Switching invariance forces Fourier support on Eulerian edge sets;
global reversal forces even edge cardinality. The first possible
nonconstant degree is four, with one common coefficient for 4-cycles.
This is an exact symmetry observation, not a sign or monotonicity theorem.

The director enumerated every switching class through order8, using an
integer Walsh transform and the unique Eulerian completion of a chord
character. All symmetry checks passed. The localization agent independently
enumerated full edge cubes through order7 and checked polynomial derivative
signs by rational root isolation. The data agree. Small complete-graph
minimizers have decreasing noise responses; this is only finite evidence.

The agent then found an exact counterexample on a sparse theta graph
with path lengths2,2,4. An actual minimizing signing has response
`13/2+t^4/2-t^6`, increasing near zero and decreasing near one.
Thus convexity, gauge symmetry, quadraticity and exact minimization alone
cannot prove fixed-branch monotonicity. Complete support would need to be
used. The optimized envelope `min_A F_A(t)` is monotone by the Markov
semigroup, but its cross-order convergence is not known: replacing M_n by
that envelope is not yet a reduction in difficulty.

## C. Orthogonal-dictionary attempt and its fluctuation failure

Suppose U is a Boolean Hadamard basis and select q of its rows as sign
columns of a bridge. A center from U has at most one nonzero field, of
size n. The bridge has operator norm sqrt(n), so a word at Hamming
distance rn from a center has bridge at most
`n+2sqrt(epsilon*r)n^(3/2)` when `q=epsilon n`.

Consequently a LOCAL energy barrier `deficit>=lambda*r*n^(3/2)` around
the complete near-extreme code would give incremental slope `1/lambda`.
This is an actual sign construction, but its cover/barrier hypothesis
has not been proved for any relevant minimizing sequence.

A tempting improvement averages all n dictionary columns and observes
that basis-center average absolute response is only one. This expectation
cannot be substituted for q sampled columns uniformly over exponentially
many words: leading fluctuations remain. Sampling without replacement
keeps the operator bound but does not give a uniform small restricted
norm on every linear-size sparse set. Partial-Hadamard RIP results from
the earlier archive do not supply that assertion at fixed linear density.
Do not silently replace a full dictionary by fractional column counts.

## D. Positive replacement now assigned for proof

Let `U=H_k tensor H_p`, n=kp, and use a sign column
`Z=h_L tensor g`, with L uniform among the k Hadamard rows and g a
uniform p-sign vector. It has EXACT mean0 and covarianceI, and Euclidean
subGaussian constant sqrt(k). For every center row of U,

```math
\mathbb E|f\cdot Z|=\mathbb E\left|\sum_{j=1}^p g_j\right|
\sim\sqrt{2/\pi}\sqrt{n/k}.
```

This reduces absolute response by a fixed factor without reducing its
variance. With q=epsilon*n independent columns, all n centers concentrate
simultaneously. A residual at Hamming radius r has proxy4krn per column.
The proposed shell bound is

```math
|x^TCy|/n^{3/2}
\le\sqrt{2/\pi}\,\epsilon/\sqrt{k}
 +\sqrt{8k\epsilon r[h(r)+\epsilon\log2]}+o_n(1).
```

If the FULL near-extreme code is covered by these centers with
`r(eta)h(r(eta))=o(eta)`, finite geometric shells appear to give extension
slope arbitrarily close to `sqrt(2/pi)/sqrt(k)`. The center code itself
has maximal factorization complexity n; the old low-factorization-center
test does not cover this situation. The Bernoulli researcher is proving
and falsifying this now. Neither the cover hypothesis nor selectable
optimizer applicability is claimed.

## E. Additional literature triage, not imported premises

The director inspected the introduction and precise normalization of
[Aronow--Lopatto2609.09103](https://arxiv.org/html/2609.09103), a recent
quantitative SK/Parisi paper. Its positive fractional moments, expected
ground state and upper-tail estimates concern random Gaussian disorder.
No transfer to minimum-over-signings or lower extreme disorder was
established. It is not a dependency of the current theorems.

Searches for a rigorous order-n-squared lower-tail law for the SK ground
state did not produce a theorem with verified applicable hypotheses.
Replica/numerical large-deviation reports were not used as premises.
