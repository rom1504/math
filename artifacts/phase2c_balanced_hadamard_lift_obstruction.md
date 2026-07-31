# Exact obstruction to a balanced-diagonal Hadamard lift family

Date: 2026-07-31. This is an agent-authored research report. It investigates
one explicit algebraically closed state; it does not modify user directives.

## 1. Candidate and why it was plausible

Let `A` be the saved order-14 symmetric conference signing with exact cap
`21`, so

```math
A^2=13I,\qquad {\operatorname{cap}(A)\over14^{3/2}}
=0.4008918628\ldots .                                  \tag{BH1}
```

For a balanced diagonal sign matrix `D` and the symmetric Sylvester matrix
`H_k`, define

```math
S_D(k)=A\mathbin\otimes H_k
       +D\mathbin\otimes(H_k-\operatorname{diag}H_k).  \tag{BH2}
```

Every off-diagonal entry of (BH2) is a sign and its diagonal is zero. The
state consists only of the fixed pair `(A,D)` and the Sylvester exponent, so
it has constant descriptive complexity and is closed under tensoring the
micro-coordinate by `H_4`.

This was a plausible repair of the earlier uniform Hadamard lift. Put
`Q=A+D` and `Delta_k=diag(H_k)`. Then

```math
S_D(k)=Q\mathbin\otimes H_k-D\mathbin\otimes\Delta_k. \tag{BH3}
```

For every Boolean spin `x`, the second quadratic form in (BH3) equals the
trace of its diagonal matrix, namely

```math
x^{\mathsf T}(D\mathbin\otimes\Delta_k)x
=\operatorname{tr}(D)\operatorname{tr}(\Delta_k)=0.   \tag{BH4}
```

Thus the balanced choice removes the `+n` shift that made the uniform
diagonal lift immediately unfavorable. Product spins retain the seed's
normalized Boolean energy, while the fixed-size operator norm of `Q` proves
that every member still has cap `O(N^(3/2))`.

## 2. Exact scalable obstruction

Take

```text
diag(D)=(1,-1,-1,-1,1,1,-1,1,1,1,-1,-1,1,-1).
```

The verifier records an explicit Boolean spin `z` of length 56 for which

```math
H_{S_D(4)}(z)=220.                                    \tag{BH5}
```

This is direct integer arithmetic, not a heuristic cap assertion.

Let `L=4^(r-1)`. The Sylvester identity is

```math
H_{4^r}=H_4\mathbin\otimes H_L.                       \tag{BH6}
```

The order-four matrix has a Boolean `+2` eigenvector
`v=(-1,-1,-1,1)`. Hence `y=v^(tensor(r-1))` is a Boolean `+sqrt(L)`
eigenvector of `H_L`, and

```math
y^{\mathsf T}H_Ly=L^{3/2}.                            \tag{BH7}
```

Use the product spin `z tensor y` in (BH3). Equations (BH4)--(BH7) give the
exact theorem

```math
\operatorname{cap}(S_D(4^r))\ge220L^{3/2}.            \tag{BH8}
```

Since its order is `N_r=56L`,

```math
\boxed{
{\operatorname{cap}(S_D(4^r))\over N_r^{3/2}}
\ge {220\over56^{3/2}}
=0.524977439470833\ldots >\frac12.}                  \tag{BH9}
```

The order-56 witness and materialized tensor checks through order 896 are in
`computations/results/balanced_hadamard_lift_obstruction.json`.

## 3. The obstruction is a linear landing gap

The canonical square-field Paley construction proves at every order

```math
M_N\le\left(\frac12+o(1)\right)N^{3/2}.               \tag{BH10}
```

Let `U_(N_r)=cap(S_D(4^r))`, `u_(N_r)=U_(N_r)^(2/3)`, and
`b_(N_r)=M_(N_r)^(2/3)`. Equations (BH9)--(BH10) imply

```math
u_{N_r}-b_{N_r}
\ge\left[
 \left({220\over56^{3/2}}\right)^{2/3}
 -2^{-2/3}-o(1)
\right]N_r
=(0.0208088299\ldots-o(1))N_r.                       \tag{BH11}
```

Thus this bounded-complexity algebraic family has a proved **linear** landing
gap in the project's `b` scale. This is a scalable obstruction, not an
isolated finite cap and not merely a failed certificate.

## 4. Order filling and research judgment

The native orders have ratio four. Taking a principal restriction of the
next native order can delete a fixed positive fraction of the vertices; the
ordinary operator estimate then loses a leading constant, so it is not a
summable order-filling operation. More importantly, (BH11) already rejects
the family on its native orders, before this order-filling defect matters.

The mechanism identifies the precise failure: balancing `D` protects all
separable product spins but does not control entangled Boolean spins across
macro and micro coordinates. Any new tensor/design family must therefore
prove an entangled-spin norm bound at the target scale. Closure of spectra,
conference identities, and preservation of the seed's product energies do
not suffice.

This inactivates the balanced-diagonal Sylvester lift as a near-optimal
structured family. It does not refute all association-scheme or code
families; it supplies a reusable falsification test for them: a single
finite entangled witness tensorizes to a leading obstruction whenever the
micro algebra has a Boolean extremal eigenvector.
