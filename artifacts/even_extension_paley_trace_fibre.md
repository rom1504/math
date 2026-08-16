# A trace-fibre Boolean state in even-degree Paley matrices

## Status

**Verified selected-family theorem.**  This gives a scalable arithmetic state
inside even-extension Paley conference matrices.  It is not a lower bound for
arbitrary signings and therefore is not a nonconvergence mechanism.
The construction uses the classical
[Paley conference matrix](https://doi.org/10.1002/sapm1933121311) and the
standard quadratic Gauss-sum diagonalization (a modern statement of the
extension formula is [Zheng, Theorem 1.1](https://doi.org/10.1016/j.jnt.2017.04.005)).

## Theorem

Let `p` be an odd prime, let `q=p^(2r)`, and let `B_q` be the Paley core on
`F_q`,

```math
(B_q)_{a,b}=\chi_q(a-b),\qquad (B_q)_{a,a}=0,
```

where `chi_q` is the quadratic character.  Choose nonzero `t in F_q` and a
Boolean function `phi:F_p->{+-1}` satisfying

```math
\sum_{u\in F_p}\phi(u)=1.
```

Then the Boolean vector

```math
x_a=\phi\bigl({\rm Tr}_{F_q/F_p}(ta)\bigr)
```

satisfies

```math
\boxed{
|x^TB_qx|=\left(1-{1\over p^2}\right)q^{3/2}.}       \tag{1}
```

For the bordered symmetric Paley conference matrix

```math
C_q=\begin{pmatrix}0&\mathbf1^T\\ \mathbf1&B_q\end{pmatrix}
```

of order `q+1`, the one-copy quadratic cap obeys

```math
\boxed{
Q(C_q)\ge {1\over2}\left(1-{1\over p^2}\right)q^{3/2}
             +{q\over p}.}                         \tag{2}
```

In particular, for `p=5`,

```math
Q(C_{5^{2r}})\ge {12\over25}q^{3/2}+{q\over5}.
```

## Proof

Additive characters diagonalize `B_q`.  For a nontrivial character indexed
by `s in F_q^*`, its eigenvalue is a quadratic Gauss sum of magnitude
`sqrt(q)`, with sign depending on `chi_q(s)`.  Because the extension degree is
even,

```math
\chi_q(k)=\chi_p(k)^{2r}=1
\qquad(k\in F_p^*).
```

The Fourier transform of `x` is supported on the `p` additive characters
indexed by `kt`, `k in F_p`.  All `p-1` nonconstant characters therefore lie
in one eigenspace of `B_q`, with one common eigenvalue of absolute value
`sqrt(q)`.  The constant character has eigenvalue zero because
`B_q\mathbf1=0`.

The vector has squared norm `q` and mean `1/p`, so its constant projection has
squared norm `q/p^2`.  Its nonconstant projection consequently has squared
norm `q(1-1/p^2)`.  The spectral decomposition gives (1).

Also

```math
\mathbf1^Tx={q\over p}.
```

Choose the border spin to align the sign of the border contribution with
`x^TB_qx`.  The matrix quadratic value is then the right side of (1) plus
`2q/p`.  Since the project convention is
`Q(C)=max_z|z^TCz|/2`, equation (2) follows.

## Boundary

For odd extension degree, `chi_q` restricts to `chi_p` rather than the
constant character on `F_p^*`; the nonzero trace-line Fourier modes split
between the two Paley eigenspaces, so this particular identity disappears.
That parity effect concerns the selected Paley family only.  A genuine
nonconvergence proof would additionally need an all-signings lower theorem on
one order sequence and a strictly smaller construction on another.
