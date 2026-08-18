# Adversarial audit: actual-child inverse-mixture barrier

**Object audited:**
[`drafts/actual_child_inverse_mixture_barrier.md`](../drafts/actual_child_inverse_mixture_barrier.md)

**Verdict:** **PASS.**  The directions of Jensen, KL/Renyi monotonicity, and
the Hellinger-to-total-variation bound are all correct.  IM.12 is only a
lower bound on the cost of the displayed Jensen domination certificate; the
draft does not promote it to a distance theorem.  IM.18--IM.20 are a genuine
distance theorem in their stated large-`beta`, `lambda=1` regime.

## 1. Reversed channel and component complexity

For one bridge bit,

```math
(1+\tanh(t)Qb)^{-\lambda}
=(\cosh t)^\lambda e^{-\lambda tQb}.
```

Averaging over the fair bit gives
`c_lambda(t)=(cosh t)^lambda cosh(lambda t)`.  Dividing by this constant
therefore produces exactly the reversed product channel IM.5, of mean
`-tanh(lambda t)Q`.  Its row Renyi-two divergence is

```math
n\log(1+\tanh^2(\lambda t)),
```

so IM.7 has the correct normalization and fixed-parameter scale.

Convexity of `z->z^(-lambda)` gives

```math
(E k_Q)^{-\lambda}\le E k_Q^{-\lambda},
```

which is the direction in IM.8.  After normalization it yields
`q<=e^J r`, hence `D(q||r)<=J`.  The rejection sampler in IM.11 has
acceptance exactly `e^{-J}` for that displayed envelope.  This proves only
`D_infinity(q||r)<=J`, not equality or optimality, exactly as the draft
warns.

## 2. IM.12 Jensen lower and upper bounds

Central symmetry decomposes the child latent law into pairs `{Q,-Q}`.  For
each pair,

```math
p_Q^{pair}(B)={\cosh(tS_Q(B))\over(\cosh t)^L},
```

and under the fair bridge `S_Q` is a sum of `L` fair signs.  Therefore every
pair has the same inverse moment

```math
Z_pair=(\cosh t)^{\lambda L}E\cosh(tS_L)^{-\lambda}.
```

Applying convexity once more to the mixture of pair densities gives
`Z_actual<=Z_pair` (the potentially confusing direction).  Thus

```math
J\ge L\log\cosh(\lambda t)
      -\log E\cosh(tS_L)^{-\lambda}
 \ge L\log\cosh(\lambda t),
```

because the expectation is at most one.  Conversely `E_U p=1` and Jensen
give `Z_actual=E p^{-lambda}>=1`, hence

```math
J\le L\{\lambda\log\cosh t+\log\cosh(\lambda t)\}.
```

This verifies both sides of IM.12 and the coefficients in IM.13.  Neither
inequality lower-bounds `D(q||r)`.

## 3. Rank-one support and IM.18

At finite temperature every child-spin state has positive weight, and the
words `Q=tau xy^T` range over all rank-one sign matrices.  The simultaneous
sign redundancy leaves exactly `2^(m+n-1)` distinct words.

For fixed `Q`, `S_Q(B)` is centered `L`-subgaussian under a fair bridge.
The exponential maximum bound therefore gives

```math
E_U\max_QS_Q\le\sqrt{2L\log K}.
```

Since

```math
\log p(B)\le t\max_QS_Q(B)-L\log\cosh t,
```

negating and averaging proves IM.18 with the stated direction.  The bound
may be negative below its threshold; the draft does not claim otherwise.

## 4. KL, Hellinger, and TV directions in IM.19--IM.20

At `lambda=1`, `tilde k_Q=k_{-Q}`.  Central symmetry therefore gives
`r_1=p` exactly.  With `Z=E_U p^{-1}` and `q=p^{-1}/Z`, the Hellinger
affinity is

```math
A(q,p)=E_U\sqrt{qp}=Z^{-1/2}.
```

Jensen gives

```math
\log Z\ge E_U[-\log p]=D(U||p).
```

Renyi monotonicity gives

```math
D(q||p)\ge D_{1/2}(q||p)=-2\log A=\log Z,
```

which verifies the first inequality of IM.20.  Finally

```math
1-TV(q,p)=E_U\min(q,p)\le E_U\sqrt{qp}=A,
```

so

```math
TV(q,p)\ge1-e^{-D(U||p)/2}.
```

All KL and TV directions are therefore correct.

Substituting `L=theta(1-theta)N^2+o(N^2)` and
`log K=(N-1)log 2` gives

```math
\gamma={\beta^2\over2}\theta(1-\theta)
       -\beta\sqrt{2\theta(1-\theta)\log2},
```

positive precisely for
`beta>sqrt(8log2/[theta(1-theta)])`.  In that regime the KL is extensive
and the TV tends to one exponentially.

## 5. Scope

IM.2 separates the actual negative escort only from the **canonical**
reversed child-spin proposal.  It does not lower-bound the reverse
projection onto all row products or rule out a different tight latent
mixture.  The Jensen gap IM.1 remains extensive even in hypothetical cases
where proposal and target coincide, so it cannot be used as a proxy for
distance.  These limitations are all stated correctly in the draft.
