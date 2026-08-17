# Independent audit: automatic tensor-prefix phase law

Audit target: `automatic_tensor_prefix_phase.md`.

Verdict: **REPAIR (minor statement repairs; the mathematical argument passes).**

The theorem is correct for a nontrivial order `h>1`, with tensor coordinates
put in the displayed lexicographic order.  No trace-zero hypothesis is needed.
The two repairs are:

1. state `h>1` explicitly (the formal order-one Hadamard satisfies the written
   hypotheses, but `[1,h]`, `log_h`, and the mantissa parametrization then
   degenerate); and
2. make the uniform assertion precise by introducing

   ```math
   G_r(t)=
   {Q(A_{\lfloor th^r\rfloor})\over\lfloor th^r\rfloor^{3/2}}.
   ```

   What is proved is `sup_(1<=t<=h)|G_r(t)-L(t)| -> 0`.  This follows from
   uniform convergence of `F_r` and the uniform floor estimate, but is only
   described informally as "uniformly after the displayed rescaling."

It would also improve the statement to define a base-`h` rational as
`p/h^k` with integer `p` and `h^k<=p<=h^(k+1)`, and to declare the recursive
Kronecker order `H_(r+1)=H tensor H_r`.  These are clarifications, not extra
mathematical assumptions.

## 1. Coherent prefix identity

Take `H_0=(1)` and order the tensor coordinates lexicographically so that

```math
H_(r+1)=H\mathbin\otimes H_r.
```

The first `h^r` coordinates correspond to the first outer coordinate.  Since
`H_(1,1)=1`, their principal block is exactly `H_r`.  Hence the finite powers
are compatible and define the claimed infinite symmetric sign matrix.

For `t=p/h^k`, `r>=k`, and `R_r=h^r`, one has exactly (with no floor error)

```math
\lfloor tR_r\rfloor=p h^(r-k).
```

Associativity in the same coordinate order gives

```math
H_(r+1)=H_(k+1)\mathbin\otimes H_(r-k).
```

Its leading `p h^(r-k)` principal block is therefore
`B_(p,k) tensor H_(r-k)`.  The endpoints `p=h^k` and
`p=h^(k+1)` are both admissible; the latter merely gives the whole outer
matrix.

## 2. Fixed-template convergence and normalization

Write `s=r-k`, `B=B_(p,k)`, and `N=p h^s=tR_r`.  The regular-Hadamard theorem
applied to the fixed `p by p` symmetric matrix `B` gives convergence of

```math
{1\over2N^(3/2)}\max_x |x^T(B\otimes H_s)x|.
```

The diagonal of `B tensor H_s` consists of signs.  If `E` is that diagonal,
then for every Boolean `x`, `x^TEx=tr(E)`, and in any event

```math
{1\over2N^(3/2)}|x^TEx|<= {1\over2\sqrt N}.
```

Thus hollowing changes the normalized absolute response by at most
`1/(2sqrt(N))`; no trace condition is required.  Finally

```math
F_r(t)=t^(3/2){Q(A_N)\over N^(3/2)},
```

so the fixed-template limit proves convergence of `F_r(t)` on the stated
dense set.  The positive Boolean eigenvalue is exactly what the cited
amplification theorem needs.  It also implies that `sqrt(h)` is integral
(each coordinate of `Hu` is an integer), but that arithmetic consequence is
not otherwise used.

## 3. Principal deletion and the operator estimates

For the block decomposition at orders `n<=m`, choose a Boolean maximizer `x`
for `A_n` and a sign `sigma` such that

```math
sigma (x^TA_nx)/2=Q(A_n).
```

If the missing Boolean coordinates `y` are independent and unbiased, then
the cross term and the hollow internal quadratic have mean zero.  Therefore

```math
E_y[ sigma (x,y)^TA_m(x,y)/2 ]=Q(A_n),
```

and one extension realizes at least that signed value.  This proves
`Q(A_n)<=Q(A_m)` for every hollow symmetric matrix, as used in the draft.

The rectangular cross block is a coordinate compression of `H_(r+1)`, so

```math
||C||<=||H_(r+1)||=sqrt(hR_r).
```

If `D_0` is the corresponding unhollowed principal block, then
`||D_0||<=sqrt(hR_r)`.  Hollowing subtracts a diagonal sign matrix of norm
one, hence `||D||<=sqrt(hR_r)+1`.  Consequently, for Boolean blocks,

```math
|x^TCy|<=sqrt(n)sqrt(d)sqrt(hR_r),

|y^TDy|/2<=d(sqrt(hR_r)+1)/2,
```

which is exactly (ATP.8).  Pointwise comparison of the old and new energies,
together with principal deletion, gives (ATP.9) with all factors of two
correct.

## 4. Asymptotic equicontinuity and uniform convergence

The floor difference satisfies

```math
d/R_r<=s-t+1/R_r,
```

while `n/R_r<=h`.  Thus (ATP.9) supplies, uniformly in `r`, an estimate of
the form

```math
0<=F_r(s)-F_r(t)
<=C_h(sqrt(s-t)+s-t)+eta_r,
\qquad eta_r->0.
```

The spectral estimate also gives a uniform bound, for example

```math
F_r(t)<= {1\over2}(n/R_r)(sqrt(h)+R_r^(-1/2))<=O_h(1).
```

Given `epsilon`, choose a finite base-`h` rational net fine enough for the
common modulus and then choose the level so that the finitely many net
values are Cauchy and `eta_r` is small.  Comparing arbitrary points to
neighboring net points proves that `(F_r)` is uniformly Cauchy.  Its limit
is continuous.  Each `F_r` is nondecreasing by principal deletion, so the
limit is nondecreasing as well.  The fact that the `F_r` themselves are
step functions causes no problem because the jump error is included in
`eta_r`.

For the normalized functions, put

```math
a_r(t)={\lfloor tR_r\rfloor\over R_r}.
```

Uniformly on `[1,h]`, `|a_r(t)-t|<=R_r^(-1)` and both quantities are bounded
away from zero.  Hence

```math
G_r(t)={F_r(t)\over a_r(t)^(3/2)}
       \longrightarrow {F(t)\over t^(3/2)}=L(t)
```

uniformly.  For an arbitrary integer `n`, taking
`r=floor(log_h n)` and `t_n=n/h^r` gives
`floor(t_nh^r)=n` exactly, so the all-order mantissa conclusion follows.

## 5. Endpoint value and absence of a trace hypothesis

At `N=h^r`, the unhollowed prefix is `H_r`, with
`||H_r||=sqrt(N)`, and `u^(tensor r)` is a Boolean eigenvector of eigenvalue
`sqrt(N)`.  If `E_r=diag(H_r)`, then `||E_r||=1` and

```math
A_N=H_r-E_r.
```

The spectral upper bound and the Boolean eigenvector lower bound yield

```math
{1\over2}-{1\over2sqrt(N)}
<= {Q(A_N)\over N^(3/2)}
<= {1\over2}+{1\over2sqrt(N)}.
```

This proves the left endpoint `L(1)=1/2` without assuming `tr(H)=0`.
At `t=h`, the same expression is the next geometric order
`hR_r=h^(r+1)`, so division by `h^(3/2)` in the definition of `L(h)` gives
`L(h)=1/2` as well.

Finally, uniform all-order approximation shows that the normalized prefix
sequence converges when `L` is constant.  Conversely, if that sequence
converges, each fixed phase `floor(th^r)` has the same limit, so `L` must be
constant.  The stated if-and-only-if is therefore valid.

## Audit conclusion

No hidden trace-zero, conference-order, or exact-hollowness assumption was
found.  The proof uses only:

- a nontrivial symmetric Hadamard generator;
- its positive regular Boolean eigenvector;
- top-left entry one and the fixed lexicographic Kronecker order; and
- standard operator compression and diagonal-deletion bounds.

After adding `h>1` and making the uniform normalized convergence explicit,
the draft is ready for canonicalization.

## Final recheck after repairs

Verdict: **PASS.**

The revised draft now states `h>1`, fixes the lexicographic recursion
`H_(r+1)=H tensor H_r`, specifies the admissible integer range for `p`, and
defines the normalized functions `G_r`.  Its final division uses exactly

```math
(\lfloor th^r\rfloor/h^r)^(3/2) -> t^(3/2)
```

uniformly on `[1,h]`, so the asserted uniform convergence of `G_r` and the
arbitrary-order mantissa conclusion both follow.  No remaining mathematical
or normalization issue was found.
