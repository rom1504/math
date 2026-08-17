# Independent audit: the selector spectral-roof assumption gap

Date: 2026-08-17.

Status: independent algebraic and repository audit.  No canonical file is
modified.

## Verdict

**REPAIR.**  Proposition SR.1 and its exact-sign consequence are correct.
For the exact-sign Candidate P1 currently stated in
`nearmin_response_sufficiency_targets.md`, there is no normalization loophole:
a vanishing joint selector defect produces a Boolean vector at the positive
spectral roof, and exact off-diagonal sign mass puts that roof at least at
`(1/2+o(1))n^(3/2)`.  Hence a uniform theorem

```text
near-minimality => exact-sign P1 with h=o(n^(3/2)), beta=o(1)
```

would prove `M_n/n^(3/2) -> 1/2`, not merely convergence to an unknown
constant.  In this proof-obligation sense P1 is not a strict reduction of the
original objective.  This does not prove a logical equivalence, or prove that
P1 is false; it proves that establishing it settles the stronger conjectural
value.

Two repairs are needed around that correct core conclusion.

1. The parenthetical claim in Candidate P1 that "a rational weighted core
   would give the same implication" is safe for the one-step response
   decoder, but not for the SR.1 convergence consequence unless an additional
   Frobenius-mass hypothesis is imposed.  Weighted or otherwise nonexact cores
   are a genuine escape from the spectral-roof no-go.
2. A cap-relative scalar selector defect has a valid one-step response bound,
   but no scalar cap-relative tensor/composition theorem was found in the
   repository.  The available generic composition theorem propagates the
   full upper response roof under a closed bi-affine feature algebra; it does
   not propagate one cap-relative selector number.

The recommended promoted pair is therefore:

1. **P2, edit-thick contracting fibre with sparse recurrent cohomology.**
2. **P3, edit-thick switched finite-type response replacement.**

Exact-sign spectral-roof P1 should be demoted from the top two.  A
cap-relative P1 should remain a research prerequisite only after a genuinely
compressed composition law is proved.

## 1. Normalization audit

Let `B` be symmetric and hollow, with `b_ij in {+-1}` for `i ne j`.  With the
repository convention

```math
H_B(x)=\sum_{i<j}b_{ij}x_ix_j={1\over2}x^TBx,
\qquad Q(B)=\max_x|H_B(x)|,                         \tag{A.1}
```

there is no missing factor of two in SR.1.

Let `Z` contain the separately labelled product columns in the nonzero
normalized Walsh support of the Boolean selector `tau`, and put

```math
a_S^\epsilon=\widehat\tau(S)\prod_{i\in S}\epsilon_i.
```

Fourier inversion is a rowwise identity, even if two product columns happen
to coincide:

```math
x_\epsilon:=Za^\epsilon,
\qquad
x_\epsilon(j)=
\tau(\epsilon_1w_1(j),\ldots,\epsilon_pw_p(j))\in\{+-1\}.
                                                               \tag{A.2}
```

Consequently

```math
(a^\epsilon)^TGa^\epsilon
={\|Za^\epsilon\|_2^2\over n}=1.                  \tag{A.3}
```

This identity uses Fourier inversion plus the Boolean range of `tau`; it does
not follow merely from Parseval.  Parseval separately gives
`||a^epsilon||_2=1`.

For `r>=||B||_(2->2)>0`,

```math
D=G-R={1\over n}Z^T(I-B/r)Z\succeq0,              \tag{A.4}
```

because `B<=rI`.  In fact the precise selector identity is

```math
(a^\epsilon)^TD a^\epsilon
=1-{x_\epsilon^TBx_\epsilon\over rn}.             \tag{A.5}
```

Thus `(a^epsilon)^TD a^epsilon<=beta` gives

```math
x_\epsilon^TBx_\epsilon\ge(1-\beta)rn.            \tag{A.6}
```

There is no sign error: this is the positive spectral channel.  If the defect
is exactly zero, (A.4) gives more than Rayleigh equality.  Namely,

```math
0=(x_\epsilon)^T(I-B/r)x_\epsilon
```

with `I-B/r` positive semidefinite, so

```math
Bx_\epsilon=rx_\epsilon.                           \tag{A.7}
```

Since `x_epsilon` is nonzero and `r>=||B||op`, exact zero defect also forces
`r=||B||op` and makes `x_epsilon` a Boolean positive-top eigenvector.  For a
small nonzero defect, (A.5) gives a near-top Boolean Rayleigh vector.  It does
not, without a spectral gap, imply norm closeness to a particular top
eigenspace; SR.1 neither states nor needs that stronger conclusion.

Candidate P1 defines `Delta_tau` as a maximum over endpoints.  Therefore
`Delta_tau<=beta` supplies (A.6) for every endpoint, while SR.1 only needs one.
The declared family must be nonempty; the nonvacuity paragraph in the
candidate supplies this condition.

## 2. Exact-sign Frobenius lower bound

For a hollow exact signing,

```math
\|B\|_F^2=\sum_{i\ne j}b_{ij}^2=n(n-1).            \tag{A.8}
```

If `s_1,...,s_n` are its singular values, then

```math
n(n-1)=\sum_js_j^2\le n\|B\|_{op}^2,
```

and hence

```math
r\ge\|B\|_{op}\ge\sqrt{n-1}.                     \tag{A.9}
```

Symmetry is not needed for this singular-value inequality, although it is
needed for the quadratic interpretation.  It is immaterial whether the
largest absolute eigenvalue of `B` is initially positive or negative: the
small selector defect itself supplies the positive near-`r` Rayleigh value.

Combining (A.1), (A.6), and (A.9) gives

```math
Q(B)\ge {1\over2}|x_\epsilon^TBx_\epsilon|
      \ge {1-\beta\over2}n\sqrt{n-1},              \tag{A.10}
```

whenever `beta<=1`.  For `beta>1`, the displayed SR.1 lower bound is
nonpositive and follows trivially; the asymptotic application has
`beta=o(1)`.

## 3. Transfer through `d_square`

For two matrices on the same labelled vertex set,

```math
d_\square(A,B)=\max_x|H_A(x)-H_B(x)|.
```

For every `x`, reverse triangle inequality gives

```math
\big||H_A(x)|-|H_B(x)|\big|
\le |H_A(x)-H_B(x)|.
```

Taking maxima in both directions yields the exact Lipschitz estimate

```math
|Q(A)-Q(B)|\le d_\square(A,B).                     \tag{A.11}
```

Thus `d_square(A,B)<=h` and (A.10) prove

```math
Q(A)\ge {1-\beta\over2}n\sqrt{n-1}-h,              \tag{A.12}
```

exactly as claimed.  This step would need a separate compiler theorem for a
smaller-order core or a different spin space; the current P1 uses a same-order
core, so no such issue arises.

## 4. Consequence for Candidate P1 and the order of limits

Suppose first that for every exact minimizer `A_n`, one has exact-sign cores
with `beta_n=o(1)` and `h_n=o(n^(3/2))`.  Since `Q(A_n)=M_n`, (A.12) gives

```math
\liminf_n {M_n\over n^{3/2}}\ge {1\over2}.         \tag{A.13}
```

The repository's rigorous frontier records

```math
\limsup_n {M_n\over n^{3/2}}\le {1\over2},         \tag{A.14}
```

from conference matrices and dense principal restrictions of orders
`n+o(n)`.  Principal restriction is legitimate because `Q` is hereditary:
fix spins on the retained vertices, extend the omitted spins randomly, and
average.  Equations (A.13)--(A.14) prove convergence to `1/2`.

The reference in SR.1 to a "conference/randomized construction bound" should
preferably be replaced by this precise repository-supported formulation.
The elementary random-sign union bound recorded elsewhere in the repository
only gives the larger constant `sqrt(log 2)`; it is not by itself (A.14).

Candidate (2.13) has an outer `epsilon` as well.  Its consequence is still
the same.  For any fixed `epsilon>0`, every exact minimizer lies in
`N_n(epsilon)`.  Choose any frame in the required nonempty declared family and
write `beta<=mathfrak C_tau`.  If the remainder in (2.13) is uniform over the
near-minimizer class, (A.12) gives

```math
\liminf_n {M_n\over n^{3/2}}
\ge {1\over2}-f(\epsilon).                         \tag{A.15}
```

Here `sqrt(1-1/n)` supplies only `o(1)`, and
`beta/2+h/n^(3/2)<=mathfrak C_tau+h/n^(3/2)`.
Sending `epsilon` to zero proves (A.13).  No diagonal choice of `epsilon_n`
is needed.  The final candidate theorem should explicitly state that its
`o_n(1)` is uniform over `A in N_n(epsilon)` and over the presentations it
asserts; otherwise the class-wide implication is formally ambiguous.

Accordingly, the exact-sign implication is demonstrably as hard in
consequence as the conjectural `1/2` lower bound.  Claims that P1 retains
fewer response data or that `beta=0` does not determine `Q(B)` do not alter
this conclusion: a predicate may forget most of a landscape and still be an
extremely strong universal property of minimizers.

## 5. Weighted and nonexact cores: the genuine loophole

For an arbitrary real weighted `B`, the same selector calculation remains
valid, but (A.8) disappears.  The general singular-value bound only gives

```math
r\ge {\|B\|_F\over\sqrt n},
```

and hence

```math
Q(A)\ge {1-\beta\over2}\sqrt n\,\|B\|_F-h.         \tag{A.16}
```

Thus the `1/2` conclusion survives for a weighted core only if, for example,

```math
\|B\|_F=(1-o(1))\sqrt{n(n-1)}.                    \tag{A.17}
```

Rationality alone says nothing about (A.17).

Nor does the allowed response distance immediately restore it.  For
`C=A-B`, Walsh orthogonality gives

```math
\mathbb E_x H_C(x)^2=\sum_{i<j}c_{ij}^2={1\over2}\|C\|_F^2,
```

so `d_square(A,B)<=h` implies only `||A-B||F<=sqrt(2)h`.  At the permitted
scale `h=o(n^(3/2))`, this is far too weak to imply
`||A-B||F=o(n)` or (A.17).  A weighted core may therefore evade the
Frobenius roof.  Whether near-minimizers actually admit useful low-mass
weighted cores with `o(n^(3/2))` all-Boolean response error is a new theorem,
not a consequence already present in P1.

The same warning applies to sparse, rectangular, smaller-order, or
compiler-mediated cores.  They can escape SR.1 only by changing a stated
hypothesis, after which same-order exact-sign realization and response
transfer must be proved separately.

## 6. Cap-relative replacement and repository composition laws

A one-step cap-relative repair is elementary.  For the selector from (A.2),
define, after orienting the positive quadratic channel,

```math
\delta_{cap}(B,W,\epsilon)
=1-{H_B(x_\epsilon)\over Q(B)}.                    \tag{A.18}
```

Since every Boolean quadratic value is at most `Q(B)`, the proof of Theorem
21.54 with the Boolean cap in place of the spectral roof gives

```math
0\le Q(B)+m\|W\epsilon\|_1-
\mathcal B_\epsilon(B;W,m)
\le Q(B)\delta_{cap}(B,W,\epsilon).                \tag{A.19}
```

This avoids every Frobenius lower bound because its baseline is the unknown
actual cap itself.

It does not inherit Theorem 21.54's tensor law.  Spectral normalization gives
the positive semidefinite matrix

```math
Z^T(I-B/r)Z/n,
```

and positivity is exactly what makes the Schur-product subtraction in the
tensor proof harmless.  The formally analogous cap-relative matrix would be

```math
D_{cap}=G-{Z^TBZ\over2Q(B)}
={1\over n}Z^T\left(I-{nB\over2Q(B)}\right)Z.       \tag{A.20}
```

Boolean cap control holds only on Boolean vectors and does not make (A.20)
positive semidefinite on the whole product-column span.  Moreover Boolean
caps do not have the simple multiplicative tensor behavior needed to
normalize a scalar recursion.  The repository's explicit "tensor
submultiplicativity" counterexample in
`artifacts/quadratic_signing_limit_research_log.md`, Section 3.8, is direct
evidence of this obstruction (with that file's full-quadratic convention).

There is one existing but materially larger fallback.  Theorem 4.1 and
equations (4.7)--(4.9) of `drafts/sufficiency_axioms_report.md` give exact
bi-affine composition and additive error propagation for the **full upper
response roofs**, provided the parent feature algebra closes in the tensor
span.  This can carry an actual-cap-relative response because it carries all
declared response directions.  It is not a composition theorem for the
single scalar (A.18), and it supplies no sublinear state bound by itself; with
a full-spin feature it can simply recode the landscape.

Therefore the repository currently contains:

- a correct one-step cap-relative estimate, (A.19);
- a generic full-roof bi-affine composition theorem; but
- no compressed cap-relative selector composition theorem replacing
  Theorem 21.54's spectral-defect subadditivity.

That missing theorem must precede promotion of a repaired P1.

## 7. Recommended revisions and ranking

The two source drafts should be interpreted or revised as follows.

1. Keep SR.1, but cite the all-order upper bound specifically as conference
   matrices plus dense principal restrictions, and state uniformity of any
   class-wide `o_n(1)` remainder.
2. Qualify the rational-weighted-core parenthesis: it gives the same response
   estimate, but gives the same `1/2` consequence only under an asymptotic
   Frobenius-mass condition such as (A.17).
3. Replace the original ranking `P1 > P2 > P3` by `P2 > P3`, with
   cap-relative P1 unranked pending a compressed composition theorem.

P2 ranks first because its near-minimizer hypothesis does not compare a
Boolean state with a universal spectral roof and it can retain low-rank
persistent coherence rather than forcing it to vanish.  P3 ranks second
because its exact-sign response transfer and bi-affine composition are
already closed, although the near-minimality-to-finite-type theorem and
cross-order compatibility remain substantial.  These are rankings of proof
targets, not claims that either implication is presently established.
