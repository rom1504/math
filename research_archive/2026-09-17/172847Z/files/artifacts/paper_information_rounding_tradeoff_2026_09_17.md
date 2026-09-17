# Information, mean retuning, and dependence in sign recovery

2026-09-17. **PROVED below; independently reconstructed by the localization
researcher.** This is a scoped recovery theorem, not a proof of convergence
or an impossibility of globally changing a signing. No external novelty
claim is made. All logarithms are natural.

## 1. A uniform random-sign comparison

For an edge vector `w` on the complete graph of order `n`, write

```math
Q(w)=\max_{x\in\{\pm1\}^n}\left|\sum_{i<j}w_{ij}x_ix_j\right|,
\qquad \overline Q_n=\mathbb E Q(\varepsilon),
```

where the edges of epsilon are independent fair signs. If `R_r` has
independent sign edges of arbitrary means `r_e in [-1,1]`, then

```math
\boxed{\mathbb E Q(R_r)\ge \overline Q_n
 -\sqrt{2n\log2}\,\|r\|_2-8n^{4/3}.}                 \tag{1}
```

The numerical constant 8 is deliberately loose. No spectral hypothesis,
model of minimizers, Gaussian ground-state formula, or limit theorem is
needed.

**Proof.** Use the symmetric witness set
`{sigma*(x_i x_j): sigma=+-1, x_1=1}`, of cardinality at most `2^n`.
For `tau>0` let

```math
F_\tau(w)=\tau^{-1}\log\sum_{\sigma,x:x_1=1}
 \exp\{\tau\sigma\sum_e w_e x_i x_j\}.
```

Then `Q<=F_tau<=Q+n log2/tau`. Every coordinate of every witness is a
sign. Its Gibbs third central moment is `-2m(1-m^2)`, so each pure
third derivative of `F_tau` has absolute value at most `tau^2`.

Center each biased edge and replace it successively by an independent
Gaussian of the same variance `v_e=1-r_e^2`. Taylor expansion through
order two cancels the first two moments. Since the centered sign has
third absolute moment at most 2 and a variance-at-most-one Gaussian has
third absolute moment below 2, the replacement error is at most
`(2/3) tau^2 binomial(n,2)`. Thus

```math
\mathbb E Q(R_r)\ge
 \mathbb E Q(r+G_v)-n\log2/\tau
 -(2/3)\tau^2\binom n2.
```

The function Q is convex and even, and `G_v` is symmetric. Consequently
`E Q(r+G_v)>=E Q(G_v)`. Couple a standard edge Gaussian as
`G=G_v+(r_e g'_e)_e`, with independent `g'`. The triangle inequality and
the elementary finite Gaussian maximum bound give

```math
\mathbb E Q(G_v)\ge\mathbb E Q(G)
 -\sqrt{2n\log2\sum_e r_e^2}.
```

The same Lindeberg argument at zero bias compares `E Q(G)` with
`overline Q_n`, at an additional error no larger than
`n log2/tau+(2/3)tau^2 binomial(n,2)`. Take `tau=n^(-1/3)`.
The combined error is below `3n^(4/3)`, hence certainly below the
displayed conservative constant. All replacements concern independent
edge coordinates; no independence of the exponentially many witnesses
is assumed. This proves (1).

## 2. Information in a switching orbit of any signing

Fix **any** hollow full sign matrix A, possibly an exact minimizer.
Let `xi` be uniform on the vertex cube, and set
`X_ij=A_ij xi_i xi_j`. Then `Q(X)=Q(A)` pointwise. Let Y be any channel
whose input is X, and set

```math
b_e(Y)=\mathbb E[X_e\mid Y],\qquad I=I(X;Y).
```

One has the dimension-sensitive inequality

```math
\boxed{\mathbb E\sum_e b_e(Y)^2\le n I.}            \tag{2}
```

**Proof.** Conditional on Y, the lifted distribution nu of xi is
globally sign-symmetric, so every vertex marginal remains fair. Fix a
vertex i. Conditional entropy subadditivity implies

```math
\sum_{j\ne i}I_\nu(\xi_i;\xi_j)
\le n\log2-H(\nu)=D(\nu\|U_n).
```

For two fair signs with correlation c, their mutual information is
`[(1+c)log(1+c)+(1-c)log(1-c)]/2 >= c^2/2`.
Sum over i and divide the double count by two. This gives
`sum_(i<j)(E_nu xi_i xi_j)^2 <= n D(nu||U_n)` pointwise in Y.
Multiplication by `A_ij` does not change the square. Averaging gives
(2), since `I(xi;Y)=I(X;Y)` for a channel of X.

This proof is also recorded independently in
[the localization reconstruction](paper_localization_2026_09_17.md).
There is no claim that a channel breaks the global switching symmetry.

## 3. Three-resource theorem

Let Z be any full sign output conditional on Y; arbitrary dependence
between its edges is permitted. Define

```math
r_e(Y)=\mathbb E[Z_e\mid Y],\quad
R=\mathbb E\|r(Y)-b(Y)\|_2^2,
\quad D=\mathbb E D_{\rm KL}
 \left(P_{Z\mid Y}\middle\|\bigotimes_eP_{Z_e\mid Y}\right).
```

Then

```math
\boxed{\begin{aligned}
\mathbb E Q(Z)\ge\overline Q_n
&-\sqrt{2n\log2}\bigl(\sqrt{nI}+\sqrt R\bigr)\\
&-\sqrt{n(n-1)D}-8n^{4/3}.
\end{aligned}}                                             \tag{3}
```

**Product transport proof, including degeneracies.** For any law P on
`d` signs and the product P0 of its marginals, sequentially couple
`P(Z_i|Z_<i)` to the i-th marginal of P0 by optimal binary coupling.
The product-side conditional marginal is fixed, hence that side really
has product law. Chain rule and binary Pinsker, followed by
Cauchy--Schwarz, give

```math
\mathbb E d_H(Z,Z^0)
\le\sum_i\mathbb E\sqrt{D(P(Z_i\mid Z_{<i})\|P_i)/2}
\le\sqrt{dD(P\|P_0)/2}.
```

Coordinates with deterministic marginals cause zero error and can be
deleted; thus no positivity assumption is hidden. Changing one edge
changes Q by at most two. Apply the coupling conditionally on Y and
average to obtain a loss at most `sqrt(2dD)`.

Apply (1) to the conditional product with means r. Minkowski and (2)
give `sqrt(E||r||_2^2)<=sqrt(nI)+sqrt(R)`. Jensen and `2d=n(n-1)` now
prove (3).

In particular, if `I=o(n)`, `R=o(n^2)`, and `D=o(n)`, the output cannot
improve the random-sign expected cap by a fixed multiple of `n^(3/2)`.
If a target has normalized cap `c_n` below
`bar c_n=overline Q_n/n^(3/2)`, necessarily

```math
\sqrt{2\log2}\left(\sqrt{I/n}+\sqrt{R/n^2}\right)
 +\sqrt{D/n}\ge\bar c_n-c_n-O(n^{-1/6}).             \tag{4}
```

This remains true for channels of **actual exact minimizers**. The theorem
does not assume that they are flat, involutory, conference-like, or
random. It quantifies the information/retuning/dependence alternatives
instead of asserting that independent rounding is the only possible
recovery operation.

## 4. Scope, combination, and relation to the original problem

The proof combines a switching-specific entropy inequality, uniform
Lindeberg replacement with deterministic means retained, Gaussian
variance completion, and entropy transport. The elementary greedy bound
`overline Q_n >= (2/3)sqrt(2/pi)n^(3/2)-O(n^(4/3))` already makes the
tradeoff nonvacuous against cap `1/2` constructions; the theorem itself
retains the full random-sign benchmark, not merely that lower bound.

None of these steps alone supplies (3). Gaussian replacement alone does
not control posterior biases; information alone does not compare maxima;
product transport alone does not track mean retuning.

The theorem does **not** obstruct deterministic output of a selected
good seed: that operation can have `I=D=0` but `R=Theta(n^2)`. Nor does
it obstruct strongly dependent completion. Its implication for convergence
is a precise constraint on localization-plus-recovery proposals, not a
cross-order recurrence or an improved universal signing bound.

No external novelty claim is currently made: each imported mechanism is
classical, and the combined signing-specific statement is what has been
proved here. Further work must use, rather than merely rename, one of
the three resources to construct a good output.
