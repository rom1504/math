# The spectral-roof assumption gap in the proposed near-minimizer selector lemma

Date: 2026-08-17.

Status: proof draft for director review.  This checks whether Candidate P1
in `nearmin_response_sufficiency_targets.md` is genuinely weaker than the
original objective.  The answer is negative for exact-sign cores: its
vanishing-defect form already forces the conjectural constant `1/2`.

## Proposition SR.1

Let `A` and `B` be hollow exact signings of order `n`.  Suppose

```math
d_\square(A,B)=\max_x|H_A(x)-H_B(x)|\le h.                   \tag{SR.1}
```

Let `r>=||B||_(2->2)`, let `W` be a Boolean port frame, and let `tau` be an
odd Boolean selector.  In the notation of Candidate P1, suppose for at least
one endpoint `epsilon` that

```math
(a^epsilon)^T(G-R)a^epsilon\le\beta,                         \tag{SR.2}
```

where

```math
G={Z^TZ\over n},
\qquad R={Z^TBZ\over rn},
```

and `Za^epsilon=x_epsilon` is the Boolean selector spin.  Then

```math
\boxed{
Q(A)\ge {1-\beta\over2}n\sqrt{n-1}-h.}                       \tag{SR.3}
```

Consequently, if every exact minimizer admits such an exact-sign core with

```math
\beta=o(1),
\qquad h=o(n^{3/2}),                                         \tag{SR.4}
```

then

```math
{M_n\over n^{3/2}}\longrightarrow {1\over2}.                \tag{SR.5}
```

Thus the exact-sign version of the proposed implication

```text
near-minimality -> vanishing spectral-roof selector defect
```

is not currently a strict reduction of convergence.  It proves the stronger
conjectural value.

### Proof

Because `x_epsilon=Za^epsilon` is Boolean,

```math
(a^epsilon)^TGa^epsilon
={||x_epsilon||_2^2\over n}=1.                                \tag{SR.6}
```

Equation (SR.2) therefore gives

```math
{x_epsilon^TBx_epsilon\over rn}
=(a^epsilon)^TRa^epsilon\ge1-\beta.                           \tag{SR.7}
```

Hence

```math
Q(B)\ge {1\over2}x_epsilon^TBx_epsilon
      \ge {1-\beta\over2}rn.                                 \tag{SR.8}
```

The exact-sign Frobenius identity is

```math
||B||_F^2=n(n-1).
```

Therefore `r>=||B||op>=||B||F/sqrt(n)=sqrt(n-1)`.  Finally
`|Q(A)-Q(B)|<=d_square(A,B)`, which proves (SR.3).

Under (SR.4), (SR.3) gives

```math
\liminf_n{M_n\over n^{3/2}}\ge {1\over2}.
```

The known conference/randomized construction bound gives the matching
`limsup<=1/2`, proving (SR.5). `square`

## Boundary of the obstruction

The proposition does not rule out three materially different repairs.

1. A **cap-relative selector defect** may compare selector energy with
   `Q(B)`, not the spectral roof `rn/2`.  It avoids forcing `1/2`, but the
   tensor/composition law of Theorem 21.54 no longer follows from positivity
   of `G-R`.
2. A weighted core can evade the Frobenius lower bound on `r`.  To be useful,
   however, its all-Boolean response distance from the exact signing must be
   `o(n^(3/2))`, and exact-sign all-order realization becomes a new missing
   lemma.
3. The semantic contracting-fibre candidate P2 does not compare energy with
   a spectral roof and is untouched by SR.1.

The director conclusion is therefore to demote exact-sign P1 from the top
near-minimality lemmas.  A repaired P1 must first prove a cap-relative
composition theorem; otherwise `near-minimality -> P1` is at least as
ambitious as proving the conjectural constant, not merely convergence.
