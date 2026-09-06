# From typical Gaussian row profiles to a true annealed high-energy floor

Date: 2026-09-06. This argument combines the signed endpoint-count theorem
with a precise typical-profile hypothesis. The counting and energy bridge
below are proved here; the uniform retained-Hadamard profile theorem is a
separate dependency being established by the director. The conclusion is
an expectation lower bound, not a second-moment or typical-instance result.

## 1. The retained-Hadamard hypothesis needed

Let `k/m ->p in (0,1]`. In each fibre choose any real order-`m` Hadamard
`H_i` and any selector `T_i` of `k` rows. For a uniform row spin `x`, put

```math
v^{(i,x)}=H_i[T_i,:]^Tx/\sqrt k.
```

Its squared norm is exactly `m`. The hypothesis is that its absolute
empirical law tends in probability in `W_2` to the half-normal law `mu`,
uniformly over all bases and selectors.

For any fixed finite truncation/quantizer `Q:[0,infinity)->A subset[0,R]`
whose discontinuities have zero `mu` mass, this supplies good row-spin
sets of relative size `1-o(1)`, uniformly over fibres, such that their
quantized empirical frequencies tend uniformly to `nu=Q_*mu`, and

```math
\frac1m\sum_j(|v_j|-Q(|v_j|))^2
\le\mathbb E(|G|-Q(|G|))^2+o(1).
\tag{1}
```

These two conclusions, rather than any independence between row-spectrum
coordinates, are what the counting argument uses.

## 2. A finite signed Gaussian edge law

Fix `rho in (0,1)` and let `(G_1,G_2)` be centered unit Gaussians with
correlation `rho`. Define

```math
(A,B,S)=(Q(|G_1|),Q(|G_2|),\operatorname{sign}(G_1G_2)),
\qquad \theta=\mathcal L(A,B,S).
```

Its two endpoint marginals are `nu`, and `theta(a,b,s)=theta(b,a,s)`.
The Gaussian relative entropy and data processing give

```math
D(\theta\Vert\nu\otimes\nu\otimes\mathrm{fair})
\le-\frac12\log(1-\rho^2).
\tag{2}
```

To check the reference law and the sign factor, the same transformation
of two independent unit Gaussians yields `nu tensor nu tensor fair`.
Before quantization, the omitted common global sign is fair and independent
of the magnitudes and relative sign under both laws, so folding actually
preserves Gaussian relative entropy. Quantization can only decrease it.

Choose clipped increasingly fine quantizers for which
`e_Q^2=E(|G|-Q(|G|))^2 ->0`. Cauchy--Schwarz then gives

```math
\mathbb E_\theta[SAB]\longrightarrow\mathbb E G_1G_2=\rho.
\tag{3}
```

All uses of the endpoint-count theorem take `Q` and `rho` fixed before
letting `m` tend to infinity.

## 3. The actual randomization and its diagonal coordinates

Independently permute output columns of the fibre bases and choose each
off-diagonal `S_ij` independently fair. For any fixed tuple of row spins,
the signs of its original spectrum can be absorbed into `S_ij`; conditional
on every column permutation, the resulting effective signs on the
undirected edges remain independent and fair. Thus their absolute
quantized spectra fit the model of the signed counting theorem.

First require the spectrum coordinate selected for the diagonal position
of each fibre to have absolute value at most two. At least `3m/4` of its
coordinates do so, by the exact norm identity. This event has probability
at least `(3/4)^m`, which costs only `O(m)` in the logarithm. Conditional
on those chosen coordinates, all remaining row multisets are independently
uniform on their `m-1` off-diagonal positions. Removing one bounded
coordinate changes the frequency and squared-error conditions by `O(1/m)`.
The endpoint-count estimates are uniform over all these possible removed
coordinates, because their repair and entropy bounds depend only on the
uniform frequency discrepancy.

This diagonal conditioning is needed for a lower energy bound. The
within-fibre rank-one energy is `k sum_i S_ii v_i(i)^2`, not just the
physical diagonal of the signing. On the conditioned event its absolute
value is at most `4mk`; deleting the physical diagonal costs at most
another `mk` in the matrix quadratic form.

## 4. Quantization errors with both signs retained

View all off-diagonal endpoint magnitudes as a vector `a` on the `m(m-1)`
labeled endpoints, and their quantizations as `b`. Let `J` reverse the
two endpoints of every edge and let `D` apply its effective edge sign to
both endpoints. Then the signed oriented product sums are `a^TDJa` and
`b^TDJb`, and `DJ` is an isometry. Consequently

```math
|a^TDJa-b^TDJb|
\le\|a-b\|_2(\|a\|_2+\|b\|_2).
\tag{4}
```

The exact row energies give `||a||_2<=m`. Equation (1), summed over
the good rows, gives `||a-b||_2<=m(e_Q+o(1))`. Thus the normalized energy
error is `O(e_Q)+o(1)`, uniformly over every endpoint assignment and
every edge-sign pattern. In particular no false monotonicity assertion
for negative edge signs is used.

The signed counting theorem with `psi(a,b,s)=sab` and (2)--(3) now
shows: for every target `q<rho`, choose `Q` sufficiently accurate and
then let `m` grow. For every fixed tuple of good row spins,

```math
\liminf_{m\to\infty}\frac1{m^2}
\log\mathbb P\left\{\frac{x^TA_{\rm hollow}x}{m^2k}\ge q\right\}
\ge\frac14\log(1-\rho^2).
\tag{5}
```

Indeed the signed finite-color product average is within the chosen
margin of `rho`; (4) pays for the actual unbounded spectrum, and the
diagonal payments divided by `m^2k` vanish. The distinction between the
number of undirected edges `m(m-1)/2` and the oriented quadratic sum
only introduces the factor `(m-1)/m ->1`.

## 5. The true annealed spin-count floor

There are `product_i |G_i|=exp(p m^2 log2+o(m^2))` tuples of good row
spins. Sum (5) over them and then let `rho` decrease to the desired
threshold. For `0<=2c sqrt(p)<1`, if `N=mk` and

```math
Y_m(c)=\#\{x\in\{-1,1\}^N:x^TA_{\rm hollow}x/2\ge cN^{3/2}\},
```

then the typical-profile hypothesis implies

```math
\boxed{\quad
\liminf_{m\to\infty}\frac1{m^2}\log\mathbb E Y_m(c)
\ge p\log2+\frac14\log(1-4c^2p).
\quad}
\tag{6}
```

The expectation is over the actual independent output permutations and
fair off-diagonal weave signs. The bases and selectors may be arbitrary
deterministic choices; averaging over additional, even correlated, base
randomness preserves the same lower bound if the profile hypothesis is
uniform as stated.

At `c=sqrt(15)/8`, the right side is nonnegative for every `0<p<=1`, by
concavity in `p` and its zero values at `p=0,1`. It is strictly positive
for every fixed `p>0` when `c<sqrt(15)/8`. Thus, subject to the stated
profile theorem, the previous floor is not only a barrier to a Finner
upper bound: the true annealed number of high-energy spins itself has
that exponential lower floor. Improving the joint positive contraction
without changing this ensemble cannot make its first moment small.

No second moment is proved. An exponentially large expected spin count
does not by itself imply that a typical signing has a high-energy spin,
nor that every realization does. Nonuniform edge-sign laws, different
ensembles, or methods which exploit concentration beyond the first moment
are not excluded.
