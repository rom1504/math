# Independent audit of the exact-plateaued weave gap

Date: 2026-09-06. Full read and independent reconstruction of
`continued_director_weave_profile_entropy_2026_09_06.md`. The positive-gap
theorem passes in its declared exact row classes. It is not an
unrestricted Boolean cap bound.

## 1. Exact matrix and entropy normalization

For independent column-permuted Hadamard frames H_i and independent
symmetric off-diagonal signs S_ij, the weave

```math
K_{(i,a),(j,b)}=S_{ij}H_i(a,j)H_j(b,i)
```

is symmetric and K^2=m^2 I. In the square, summing over the inner
coordinate b first gives m delta_(i,k), and row orthogonality of H_i
finishes the calculation. Each block spin transform h_i has squared
norm m^2. Thus the director's defect identity

```math
2(m^3-\sigma x^TKx)
=\sum_{i,j}(h_i(j)-\sigma S_{ij}h_j(i))^2
```

has the correct factor two. Hollowing changes the quadratic value by
at most m^2. Fixed-point-free matching supports give Boolean eigenvectors
for both signs of the eigenvalue; those witnesses are not eliminated by
the theorem.

The row-entropy inequality on any undirected graph also checks. Each
edge appears in two row marginals, so the sum of row entropies is at
least twice the entropy of the full edge array. Consequently the number
of compatible edge arrays is at most the square root of the product
of the numbers of allowed row words. Conditioning on already exposed
diagonal or bad-edge symbols does not alter this fact.

## 2. Positive defect: complete counting reconstruction

An exact q-support plateaued row has nonzero magnitude sqrt(m/q),
support size t=qm, and magnitude-word count N=binomial(m,qm). Take
k=2 for q<1 and k=1 for q=1.

A failed off-diagonal equality contributes at least 2m/q to the
ordered defect sum: both orientations are included. Hence
|x^TKx|>=(1-eta)m^3 implies at most eta q m^2 failed unordered
edges, after choosing the appropriate sign sigma. Diagonal failures
may be ignored for this upper bound.

For a fixed set E of e failed edges, expose its two endpoint magnitude
symbols and all diagonal symbols. There are at most k^(2e+m) choices.
The remaining graph entropy bound gives at most product_i sqrt(N_i)
completions for each such choice. Independent row permutations are
uniform on their N_i magnitude words, so divide by product_i N_i.

At least (sum_i t_i-m-2e)/2 remaining edges have both endpoints nonzero.
Each independently constrains its random S_ij for a fixed sigma. Signs
attached to permuted spectral coordinates need not be independent:
conditioning on the entire permutations leaves the independent weave
edge signs available, which is all this step requires. The resulting
bound is

```math
k^m 2^{m/2}\prod_i(N_i2^{t_i})^{-1/2}(2k^2)^e.
```

Union over both sigma values, all row spins, and all E with
e<=eta q m^2. With binomial(m,2)=m^2/2+O(m), the normalized logarithm
of the latter edge-set cost is at most

```math
\frac12H_2(2\eta q)+\eta q\log_2(2k^2)+o(1),
```

provided 2 eta q<1/2. This reproduces exactly the director's sufficient
condition and exponentially small failure probability on the m^2 scale.

## 3. Primary theorem mapping and numerical margins

I fetched and read the definitions and Theorems 1(a),2 in
[Potapov, arXiv:2303.16547v3](https://arxiv.org/html/2303.16547v3#S5).
They use unnormalized Walsh coefficients, binary logarithms, and a fixed
plateau parameter s. Substituting m=2^d and q=2^-s gives c_1=11/32
and c_q=alpha/8+(H_2(q)+q)/4 for s>0, where
alpha=1+(3/8)log_2 6. The factor 1/8 comes from the fixed-shift central
Hamming ball in dimension d-2. The parity d+s even remains required.
This imports exact counting bounds, not an approximate-plateau stability
result or an independent reproof of the whole counting paper.

The resulting g_q values at q=1,1/2,1/4,1/8 are respectively

```text
0.1562500000, 0.1288298828, 0.0191494139, -0.0790290064.
```

For example eta=0.001 satisfies the sufficient positive-gap condition
for all three positive-margin classes: its entropy costs are respectively
0.0114070357, 0.0072038789, 0.0038520296. The q=1/8 class is not certified
by these imported counts. A union over the three allowed classes still
has vanishing failure probability when the appropriate parity is used.

## 4. Scope protection

The row classes are fixed independently of the random weave and are
EXACTLY plateaued. Neither a general sparse Walsh spectrum nor an
unspecified near-plateaued row is covered. Matching eigenvectors survive,
so the full hollowed weave still has normalized cap asymptotically 1/2.
An all-row stable profile count or an actual extraction mechanism remains
a separate obligation for the original convergence problem.
