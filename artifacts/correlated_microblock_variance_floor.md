# Correlated microblock variance floor

## Status

This note audits whether negative dependence across macro edges can
improve the exact-sum traffic exponent in
`asymptotic_centered_width_recovery.md`.

The answer is sharply limited:

- changing the marginal law of one even-order block cannot improve the
  worst-state variance \(\kappa_s=s^3/(s+1)\);
- arbitrary correlations across macro edges can improve the uniform
  variance proxy by at most the factor
  \[
  \frac{\kappa_s}{s(s-1)}=\frac{s^2}{s^2-1};
  \]
- at the optimal nontrivial fibre size \(s=4\), the largest possible
  improvement is therefore only \(16/15\).

This does not exclude a useful *slack-dependent* anti-alignment.  It
does exclude a large gain from a generic low-rank orthogonal array or
from variance reduction alone.

## 1. Setup

For every macro edge \(e=ij\), let \(B_e\) be an \(s\times s\) sign
block with exact sum
\[
\mathbf1^\top B_e\mathbf1=a_es^{3/2}.
\tag{1}
\]
Put
\[
C_e=B_e-\frac{a_e}{\sqrt s}J_s.
\tag{2}
\]
Assume the randomized construction is entrywise unbiased:
\[
\mathbb E C_e=0.
\tag{3}
\]
This condition is needed for the coarse Hamiltonian to remain the
statewise mean.

Vectorize all \(C_e\), and let \(\Sigma\) be their joint covariance.
Because every realization satisfies \(\langle C_e,J_s\rangle=0\),
\(\Sigma\) is supported on
\[
\mathcal U=\bigoplus_e J_s^\perp.
\tag{4}
\]
Moreover, for every realization,
\[
\|C_e\|_F^2
=s^2+s-2s
=s(s-1).
\tag{5}
\]
Consequently
\[
\boxed{
\operatorname{tr}\Sigma
=\binom n2s(s-1).
}
\tag{6}
\]
Correlating the blocks can move this covariance, but cannot remove it.

For a microstate \(X=(x_1,\ldots,x_n)\), define the feature vector
\[
\psi_e(X)=\operatorname{vec}(x_ix_j^\top).
\tag{7}
\]
Then the residual variance is
\[
\operatorname{Var}_B R_B(X)
=\psi(X)^\top\Sigma\psi(X).
\tag{8}
\]

## 2. Global variance floor

Take the fibre spins \(x_i\) independently and uniformly from
\(\{\pm1\}^s\).  Distinct micro-edge characters are orthogonal, so
\[
\mathbb E_X\psi(X)\psi(X)^\top=I.
\tag{9}
\]
Combining (6), (8), and (9) gives
\[
\mathbb E_X\operatorname{Var}_B R_B(X)
=\operatorname{tr}\Sigma
=\binom n2s(s-1).
\tag{10}
\]
Hence
\[
\boxed{
\max_X\operatorname{Var}_B R_B(X)
\ge\binom n2s(s-1).
}
\tag{11}
\]

In particular, if a correlated ensemble is to have a uniform
subgaussian variance proxy \(\sigma_s^2\binom n2\), then necessarily
\[
\boxed{\sigma_s^2\ge s(s-1).}
\tag{12}
\]
The corresponding normalized endpoint rate coefficient can never
exceed
\[
\frac1{s(s-1)}.
\tag{13}
\]
For nontrivial exact square fibres \(s\ge4\), this is maximized at
\[
\frac1{4\cdot3}=\frac1{12}.
\tag{14}
\]

The independent exact-sum ensemble has worst-state proxy
\[
\kappa_s=\frac{s^3}{s+1}.
\tag{15}
\]
Therefore even an optimally correlated construction can improve that
proxy by at most
\[
\boxed{
\frac{\kappa_s}{s(s-1)}
=\frac{s^2}{s^2-1}.
}
\tag{16}
\]
At \(s=4\), this is \(16/15\).

## 3. Exact one-block minimax theorem

For even \(s\), the independent exact-sum marginal is already minimax
for one macro edge.

Let \(P_0=J_s/s\) and \(P_1=I-P_0\).  If \(x\) is uniform among
balanced sign vectors, then
\[
\mathbb E xx^\top=\frac{s}{s-1}P_1.
\tag{17}
\]
If \(x\) is uniform among the two constant sign vectors, then
\[
\mathbb E xx^\top=sP_0.
\tag{18}
\]

Choose a random pair \((x,y)\) by the following mixture:

- with probability \(1/(s+1)\), \(x\) is balanced and \(y\) is
  constant;
- with probability \(1/(s+1)\), \(x\) is constant and \(y\) is
  balanced;
- with probability \((s-1)/(s+1)\), both are independently balanced.

Writing \(w=\operatorname{vec}(xy^\top)\), equations (17)--(18) give
\[
\mathbb E ww^\top
=\frac{s^2}{s^2-1}
\left(
P_1\otimes P_0+
P_0\otimes P_1+
P_1\otimes P_1
\right)
=\frac{s^2}{s^2-1}P_{J^\perp}.
\tag{19}
\]

For any entrywise-unbiased exact-sum block ensemble, its covariance
\(\Sigma_e\) is supported on \(J^\perp\) and has trace \(s(s-1)\).
Thus (19) implies
\[
\max_{x,y}
\operatorname{Var}(x^\top B_ey)
\ge
\frac{s^2}{s^2-1}s(s-1)
=\frac{s^3}{s+1}
=\kappa_s.
\tag{20}
\]

The uniform exact-sum law has covariance
\[
\Sigma_e=\frac{s}{s+1}P_{J^\perp}.
\tag{21}
\]
Since
\[
\|P_{J^\perp}\operatorname{vec}(xy^\top)\|_2^2
=s^2(1-\mu_x^2\mu_y^2)\le s^2,
\tag{22}
\]
it attains the upper bound \(\kappa_s\), with equality for balanced
\(x,y\).  Hence
\[
\boxed{
\inf_{\text{unbiased exact-sum block laws}}
\max_{x,y}\operatorname{Var}(x^\top B_ey)
=\kappa_s
}
\tag{23}
\]
for even \(s\).

Thus no better traffic constant can come from changing the
single-block marginal.

## 4. What correlation could still do

The global lower bound (11) is slightly below the independent value,
so cross-edge negative dependence is not completely ruled out.
However, at \(s=4\) the entire numerical room is
\[
\frac{16}{15}-1=\frac1{15}.
\]
A low-rank orthogonal array does not evade (6): reducing covariance
rank while preserving its trace increases its nonzero eigenvalues.
It can help only if every relevant low-slack feature vector lies close
to the covariance kernel.

For a low-slack family \(\mathcal F\), define its empirical feature
Gram operator
\[
G_{\mathcal F}
=\frac1{|\mathcal F|}
\sum_{X\in\mathcal F}\psi(X)\psi(X)^\top.
\tag{24}
\]
Then every correlated ensemble obeys the exact identity
\[
\frac1{|\mathcal F|}
\sum_{X\in\mathcal F}
\operatorname{Var}_B R_B(X)
=\operatorname{tr}(\Sigma G_{\mathcal F}).
\tag{25}
\]
Therefore correlation can beat the independent traffic exponent only
through a genuine span defect of \(G_{\mathcal F}\), followed by a
covariance choice concentrated near that defect.  Large cardinality
alone does not imply the needed spectral lower bound.

The surviving target is consequently an entropy-versus-span
dichotomy:

1. either \(G_{\mathcal F}\) is quantitatively nondegenerate on
   \(\mathcal U\), forcing a variance floor through (25); or
2. its small-eigenvalue subspace yields a coherent edge-block
   direction that can be used in the global replacement identity to
   lower the centered width.

That is a structural problem about endpoint features, not a further
scalar variance optimization.
