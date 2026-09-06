# Gaussian deformation stability and exclusion of finitely confined wells

Status: self-contained derivation, pending independent audit. This is a necessary
structural theorem for asymptotic minimizers, NOT a convergence theorem.

## Setup

Write `q_A(x)=sum_{i<j} A_ij x_i x_j`, `M(A)=max_x |q_A(x)|`,
`m_n=M_n/n^(3/2)`. Let `G` have independent standard Gaussian upper entries.
For a fixed signing define

```math
F_A(t)=n^{-3/2}\mathbb E\max_{s=\pm1,x\in\{\pm1\}^n}
s(q_A(x)+tq_G(x)).
```

## 1. Uniform Gaussian deformation inequality

For every signing A, every fixed `t>=0`,

```math
F_A(t)\ge\sqrt{1+t^2}(m_n-10n^{-1/6}).                 \tag{1}
```

Proof. Put `a=(1+t^2)^(-1/2)` and `sigma=sqrt(1-a^2)`.
Independently replace each A edge by a random sign with mean `a A_ij`.
Every realization is a signing, hence has norm at least M_n. Compare the
centered independent entries with `sigma G_ij` by Lindeberg replacement.
For

```math
\Phi(W)={1\over\beta n}\log\sum_{s,x}
\exp\left({\beta\over\sqrt n}s\sum_e W_e x_ix_j\right),
```

the third derivative in any edge is bounded by `8 beta^2/n^(5/2)`.
The centered sign has absolute third moment at most 8; the Gaussian third
moment is at most 2. Summing Taylor remainder bounds over fewer than n^2/2
edges gives error at most `(20/3) beta^2/sqrt(n)`. Softmax differs from the
normalized maximum by at most `2 log(2)/beta`. Choose `beta=n^(1/6)`.
Total error is less than `10 n^(-1/6)`. Divide the resulting inequality by a.
The assertion at t=0 follows directly. All estimates are uniform in A.

The inequality alone holds for arbitrary A, but is informative near minimizers:
if `M(A_n)/n^(3/2)-m_n -> 0`, the Gaussian gain must asymptotically be at least
`m_n(sqrt(1+t^2)-1)`.

## 2. A concrete forbidden landscape for near minimizers

Fix kappa>0. For each sign s choose centers `z_(s,j)`; let K_n be their total
number. Distances are modulo global spin reversal:
`d(x,z)=min(d_H(x,z),n-d_H(x,z))`.
Suppose every signed configuration satisfies

```math
s q_{A_n}(x)\le M(A_n)-\kappa\sqrt n\min_j d(x,z_{s,j}).       \tag{2}
```

Then an asymptotic minimizer sequence cannot have `log K_n=o(n)`.
More quantitatively, along a subsequence on which `m_n -> c>0` and
`log(K_n)/n -> eta`, one has, for every `0<t<=kappa/4`,

```math
c(\sqrt{1+t^2}-1)
\le3t\sqrt\eta+e\kappa\exp(-\kappa^2/(8t^2)).                \tag{3}
```

In particular, any t for which the left side exceeds the exponential term
gives the explicit positive entropy-rate lower bound

```math
\eta\ge\left[
{c(\sqrt{1+t^2}-1)-e\kappa e^{-\kappa^2/(8t^2)}\over3t}
\right]^2.                                                  \tag{4}
```

Proof. Assign x to a nearest center of the same sign. If its distance is r,
the Gaussian increment `s(q_G(x)-q_G(z))` has variance `4r(n-r)<=4nr`.
There are at most `2K_n binom(n,r)` such pairs. Therefore their expected
maximum is at most `sqrt(8nr log(2K_n binom(n,r)))`.
The expected maximum over center Gaussian values is at most
`sqrt(2 binom(n,2) log K_n)` (zero when K_n=1).

To interchange the expected maximum over r with the maximum of expectations,
each shell maximum is an n-Lipschitz Gaussian function because its individual
increment standard deviations are at most n. The usual Gaussian mgf bound
gives extra error at most `n sqrt(2 log(n+1))`. This concentration bound can
also be proved directly by Gaussian log-Sobolev; it is the sole standard
Gaussian concentration fact used here.

Writing rho=r/n, dividing by n^(3/2), and passing to the subsequence gives

```math
\limsup [F_{A_n}(t)-M(A_n)/n^{3/2}]
\le t\sqrt\eta+
\sup_{0\le\rho\le1/2}\{-\kappa\rho+t\sqrt{8\rho(H(\rho)+\eta)}\}
\le3t\sqrt\eta+
\sup_{0\le\rho\le1/2}\rho[-\kappa+\sqrt8t\sqrt{\log(e/\rho)}].
```

Here `H(rho)<=rho log(e/rho)` and `sqrt(8rho eta)<=2sqrt(eta)`.
The last supremum is at most `e kappa exp(-kappa^2/(8t^2))`:
positive values require `L=log(e/rho)>kappa^2/(8t^2)>=2`, and
`sqrt(8)t exp(1-L)sqrt(L)` decreases for L>=2. Compare with (1).
Since its exponential term is o(t^2), (3) rules out eta=0.

The positivity c>0 does not depend on recent reported bounds. Splitting
vertices into equal halves, choosing one side uniformly and optimizing the
other, and applying the elementary lower Khintchine bound gives
`M_n >= (1/4+o(1)) n^(3/2)`.

## Scope and next hard target

This uses global optimality through an exact signing perturbation, and controls
correlated Gaussian widths rather than raw layer cardinality. It does not say
that near-minimizers have exponentially many exact maximizers. Centers in (2)
must control the full deficit landscape, a substantially stronger demand.

Potential extension: replace (2) by multiscale coverings with sublinear deficit
moduli and obtain a quantitative lower bound on their Gaussian generic-chaining
functional. A convergence proof still needs a cross-order realization theorem;
this result supplies no such theorem and is not primary convergence progress.

## Initial independent proposals and archive collision

1. Weighted two-block completion using configuration-dependent slack: dropped
after finding extensive equivalent partition-function/entropy work in ledger
sections 10.8 onward.
2. Random-edge stochastic stability of exact minimizers: scalar entropy version
overlaps sparse-flip identities around ledger lines 1900 and stratified noise
section 10.30. The Gaussian correlated-width/confinement formulation above was
not located by targeted searches; novelty remains pending full archive audit.
3. Coding-covering concatenation: augmented cut-code equivalence is already
recorded and no new quantitative concatenation inequality was found.
