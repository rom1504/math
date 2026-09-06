# Wave 28 root route: random hashing around a planted low-row cut

## Status

The theorem below is proved from rectangular matrix Bernstein.  It shows that
every individual low-row cut can be embedded, without changing that cut, in a
balanced entropy-sized block coset whose **entire** support is low-row.  Thus
row regularization and exact local alignment are compatible even when the
diagonal is fixed by a planted witness.

This is only a local reduction.  It does not prove that every selector has a
low-row favorable completion, and it does not compress the cosets/cuts chosen
for different selectors into the number required by (10.808).  Hence it does
not prove (10.795) or convergence.

## 1. Fixed-diagonal random-partition theorem

Let `A` be any symmetric zero-diagonal sign matrix of order `n`, fix
`x in {+-1}^n`, and put `D=diag(x)`.  Hash the vertices independently
and uniformly into `k` labeled blocks, and let `P` be the resulting
`n x k` incidence matrix.  Define

```math
u=\log(8(n+k)),\qquad
L=\sqrt{(n-1)(1-1/k)},
```

```math
\nu=\max\left\{
(1-1/k)\lVert A\rVert_{\rm op}^2,\,
\frac{n(n-1)}k
\right\}.
\tag{R28.1}
```

If

```math
\frac nk\ge 12\log(16k),
\tag{R28.2}
```

then some realization is a nonempty balanced partition with
`n/(2k)<=|B_a|<=2n/k` for every block and

```math
\boxed{
\max_{z\in\{\pm1\}^k}\lVert ADPz\rVert_2^2
\le
2R_2(x)+8k\nu u+\frac{16}{9}kL^2u^2.
}
\tag{R28.3}
```

### Proof

Write `g(i)` for the random block, let
`\bar e=k^{-1}1_k`, and set

```math
v_i=x_iAe_i,\qquad
Z_i=v_i(e_{g(i)}-\bar e)^{\mathsf T}.
```

Then

```math
ADP=(Ax)\bar e^{\mathsf T}+\sum_iZ_i,
\qquad \mathbb E Z_i=0.
\tag{R28.4}
```

The summand bound and the two rectangular variances are exact:

```math
\begin{aligned}
\lVert Z_i\rVert_{\rm op}
&=\sqrt{(n-1)(1-1/k)}=L,\\
\sum_i\mathbb E Z_iZ_i^{\mathsf T}
&=(1-1/k)A^2,\\
\sum_i\mathbb E Z_i^{\mathsf T}Z_i
&=n(n-1)
\left(\frac1kI_k-\frac1{k^2}J_k\right).
\end{aligned}
\tag{R28.5}
```

Thus the rectangular matrix Bernstein inequality in
[Tropp](https://arxiv.org/abs/1004.4389) gives

```math
\Pr\left\{
\left\lVert\sum_iZ_i\right\rVert_{\rm op}
\ge \sqrt{2\nu u}+\frac23Lu
\right\}
\le(n+k)e^{-u}=\frac18.
\tag{R28.6}
```

Each block size is binomial with mean `n/k`.  Standard Chernoff bounds and
a union bound give

```math
\Pr\left\{\exists a:
|B_a|\notin[n/(2k),2n/k]\right\}
\le2k e^{-n/(12k)}\le\frac18.
\tag{R28.7}
```

Hence the two events intersect.  On their intersection,
`\|(Ax)\bar e^T\|=sqrt(R_2(x)/k)`, so

```math
\begin{aligned}
\max_z\lVert ADPz\rVert_2^2
&\le k\lVert ADP\rVert_{\rm op}^2\\
&\le2R_2(x)+2k
\left(\sqrt{2\nu u}+\frac23Lu\right)^2\\
&\le2R_2(x)+8k\nu u+\frac{16}{9}kL^2u^2.
\end{aligned}
```

This proves (R28.3).

## 2. Exact-minimizer scale

Fix `0<c<1/4` and take

```math
k=\left\lfloor
\frac{n^{3/4-c}}{C_0\log n}
\right\rfloor.
\tag{R28.8}
```

For an exact minimizer,
`\|A\|_{\rm op}^2<=2q_n=O(n^(3/2))`.  Also

```math
\frac{n^2}{k}=O(n^{5/4+c}\log n)=o(n^{3/2}).
```

Therefore (R28.3) becomes

```math
\boxed{
\max_zR_2(DPz)
\le2R_2(x)+O(n^{9/4-c}).
}
\tag{R28.9}
```

The partition has maximum block size
`O(n^{1/4+c}\log n)=o(sqrt n)`, and the coset has at most
`2^k=exp(O(n^{3/4-c}))` oriented words.

## 3. Consequence and remaining gap

If a selector `S` has a favorable full cut
`d=(sigma,x)` satisfying `R_2(x)=O(n^{9/4-c})`, choose
`D=diag(x)`.  The all-ones block word reproduces `x`, so the partition
from (R28.3) obeys

```math
Q(P^{\mathsf T}DH_SDP)
\ge |x^{\mathsf T}H_Sx|
\ge X_d(S).
\tag{R28.10}
```

In particular, a low-row completion of an exact oriented child ground gives
perfect local alignment, because (10.825) has
`X_d(S)>Y_A(S)`.  The whole surrounding coset is row-good at the target
scale by (R28.9).

This proves a clean equivalence at the local level: once the low-row
favorable witness exists, forcing it into a row-regular block architecture
costs only the already allowed entropy and row exponent.  It does **not**
produce that witness uniformly for every selector.  Even if it did, assigning
one planted coset per selector could still use `exp(Theta(n))` different
cosets.  A shared-completion collision theorem, a high-excess selector count,
or a distribution giving the hit bound (10.838) is still necessary.

## 4. Falsification criterion for this reduction

The theorem itself is unconditional.  Its use as a route to (10.808) would be
blocked along a target sequence if either:

1. some high-excess selectors have no favorable completion with
   `R_2=O(n^{9/4-c})`; or
2. such completions exist but every choice map needs
   `exp(omega(n^{3/4-c}))` planted cosets (equivalently no one coset has
   the coverage required by (10.795)).

The current ledger has neither an asymptotic example nor a theorem deciding
these alternatives.
