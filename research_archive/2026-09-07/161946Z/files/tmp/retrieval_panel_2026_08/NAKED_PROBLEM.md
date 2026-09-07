# Naked problem packet

This packet intentionally omits the project archive and prior route history.

For a symmetric zero-diagonal matrix `A=(a_ij)` with `a_ij in {+1,-1}` for
`i != j`, define

```math
H_A(x)=\sum_{i<j}a_{ij}x_ix_j,
\qquad
Q(A)=\max_{x\in\{\pm1\}^n}|H_A(x)|,
\qquad
M_n=\min_AQ(A).
```

Determine whether

```math
{M_n\over n^{3/2}}
```

converges as `n` tends to infinity.  Convergence to any constant is success.
A proof of nonconvergence must establish strict `liminf < limsup`, not merely
an unusual finite sequence or a counterexample to one proof method.

The rigorous asymptotic frontier supplied to the panel is

```math
0.336493364431\ldots
\le \liminf_{n\to\infty}{M_n\over n^{3/2}}
\le \limsup_{n\to\infty}{M_n\over n^{3/2}}
\le {1\over2}.
```

## Exact code formulation

Let `E=binom(n,2)` and

```math
\mathcal C_n^+
=\{(\sigma x_ix_j)_{i<j}:\sigma,x_i\in\{\pm1\}\}.
```

This is the antipodal augmentation of the cut/coboundary code of `K_n`.
Under the sign-to-bit map,

```math
Q(a)=E-2d(a,\mathcal C_n^+),
\qquad
M_n=E-2\rho(\mathcal C_n^+).
```

The size of `C_n^+` is `2^n` for the relevant orders.  A theorem about
one-sided maximum cut, frustration, packing, or minimum distance is not
automatically a theorem about this absolute maximum or the largest coset
distance; check the map exactly.

## Finite anchors

The recorded exact values for orders `3` through `14` are

```math
(M_3,\ldots,M_{14})=(3,4,4,5,9,10,12,13,17,18,20,21).
```

Use these only for falsification or normalization checks.

## Conference-scale examples

A symmetric conference matrix `C` of order `n` is hollow, has off-diagonal
entries in `{+1,-1}`, and satisfies

```math
C^2=(n-1)I.
```

Its spectral bound gives

```math
Q(C)\le {n\over2}\sqrt{n-1},
```

so conference and dense restriction constructions naturally live at
normalized scale `1/2`.  They are examples, not known universal minimizers.

## Required output from a domain researcher

Use the separate report schema.  In particular, freeze no more than three
mathematically distinct architectures and give an exact implication through
one fully quantified missing lemma.  Prove why the lemma contains strictly
less information than full minimization or the full coset-weight histogram,
and give a decisive falsifier.  Do not ask for the archive until the proposal
is frozen.
