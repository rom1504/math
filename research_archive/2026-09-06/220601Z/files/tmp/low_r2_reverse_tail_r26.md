# Wave 26 memo: centered sign mass does not align the adaptive envelope

## Status

The algebraic reduction and the transitive-subspace sign lemma below are
**proved**.  The exact `A_9` obstruction is checked by
`tmp/reverse_sign_r26.py`.  The low-row approximate-envelope criterion is a
restatement of the required coverage and remains **open**; no codebook
construction was proved.

Put

```math
Y_A(S)=Q(A[S])-p^{3/2}q_n,
\qquad
X_d(S)=c_A(S,d)-p_2\langle A,d\rangle.
```

Then the retained-deficit definition (10.792) simplifies pointwise to

```math
\boxed{\widehat\ell(S,d)=Y_A(S)-X_d(S).}
\tag{R26.1}
```

For each fixed cut `d`, `X_d` is a degree-two polynomial on the Johnson
slice and

```math
\mathbb E_{U_m}X_d=0.
\tag{R26.2}
```

Thus the problem is not a scalar reverse tail for one centered polynomial;
it is alignment of that polynomial with the adaptive envelope `Y_A`.

## A polynomial sign-mass lemma

Let `Omega` be a finite transitive space with uniform measure and let `V` be
an invariant `D`-dimensional real function space.  Its reproducing kernel has
constant diagonal `D`, so every `f` in `V` satisfies

```math
\lVert f\rVert_\infty\le\sqrt D\,\lVert f\rVert_2.
\tag{R26.3}
```

If `f` is nonzero and centered, write `M=||f||_infty`.  Since
`f^2<=M|f|` and `E|f|=2E f_+`,

```math
\Pr\{f>0\}
\ge{\mathbb E f_+\over M}
\ge{\lVert f\rVert_2^2\over2M^2}
\ge{1\over2D},
```

and the same argument applied to `-f` gives

```math
\boxed{
\Pr\{f>0\},\Pr\{f<0\}\ge{1\over2D}.}
\tag{R26.4}
```

On the `m`-slice, degree-at-most-two functions lie in a transitive space of
dimension at most

```math
D\le1+n+\binom n2=O(n^2).
```

Consequently every nonzero `X_d` has both signs on a polynomial fraction of
the slice.  This is a genuine one-sided reverse-sign theorem, requiring no
lattice estimate.

It is far below (10.795): that event asks for

```math
X_d(S)\ge Y_A(S)-O(n^{3/2-c}),
```

and (R26.4) gives no correlation between `X_d` and the unrelated nonlinear
envelope `Y_A`.

## Exact finite wall

For the exact order-nine minimizer `A_9`, there are cuts with the minimum
row square `R_2=16` having the following selector statistics:

| `m` | slice size | `# {X_d>0}` | `# {hat ell<=0}` | `max X_d` | range of `Y_A` |
|---:|---:|---:|---:|---:|:---|
| 8 | 9 | 5 | 0 | `16/9` | constant `3.886740446...` |
| 7 | 36 | 22 | 0 | `16/3` | `[1.537547398...,5.537547398...]` |

Thus even majority positive centered payoff, at the smallest available row
cost, need not produce a single zero-residual effective-loss selector.  The
obstruction at `m=7` is genuinely one of alignment: `max X_d` exceeds
`min Y_A`, but not on the same selectors.

## Exact envelope reformulation

A sufficient codebook theorem is immediate but currently unproved.  If a
family `C` of full cuts satisfies

```math
|C|\le\exp\{O(n^{3/4-c})\},
\qquad
\max_{d\in C}R_2(d)=O(n^{9/4-c}),
```

and, for every selector in the active slice,

```math
\boxed{
\max_{d\in C}X_d(S)
\ge Y_A(S)-O(n^{3/2-c}),}
\tag{R26.5}
```

then the coverage sets of the codewords cover the whole slice.  Pigeonhole
gives one `d` with coverage at least `1/|C|`, exactly (10.795).

Equation (R26.5) isolates a low-row approximate-envelope problem.  It is not
a result: no subexponential codebook with this property is known.  A slice
mgf controls deviations of each already chosen centered polynomial in the
opposite direction and does not construct or compress the adaptive envelope.

## Conclusion

Degree-two transitivity supplies polynomial mass on both sides of the mean,
but retained-deficit coverage is an alignment problem rather than a sign
problem.  Further reverse-hypercontractive work is useful only if it couples
`X_d` to `Y_A`, or if it constructs the low-row envelope (R26.5); a scalar
reverse tail for `X_d` alone should not remain an active route.
