# Independent audit: stopped asymmetric powered-aggregation obstruction

Date: 2026-09-06. Audited source:
`transfer_fresh_subhalf_asymmetric_power_obstruction_2026_09_06.md`,
especially its strengthened Section 7. Verdict: PASS, within its stated
asymmetric pointwise scope. No convergence/nonconvergence implication
for the original minimum values has been obtained.

The audit reconstructed the sampling import and stopped argument
independently; it did not use an earlier agent verdict as a premise.

## 1. Primary sampling theorem and the actual Boolean norm

[Rudelson--Vershynin, arXiv:math/0503442v3, Theorem 1.5, p. 4](https://arxiv.org/pdf/math/0503442)
applies to a common Bernoulli selector for rows and columns. Its
expected cut-norm bound has terms `p^2||A-diag A||_C`,
`p||diag A||_C`, and `p^(3/2)(||A||_Col+||A^T||_Col)`, up to an
absolute multiplicative constant. Here the column norm is the sum of
Euclidean column lengths. This primary statement has the required
uniformity and is not a theorem about independent row/column selectors.

The remaining norm conversion is elementary. For a symmetric hollow
matrix and `beta(A)=max_(x,y in {−1,1}^n)|x^T A y|`,

```math
Q(A)\le\beta(A)/2\le2\|A\|_C,
\qquad \|A\|_C\le\beta(A)\le4Q(A).
```

The first cut comparison expands signs into their positive and negative
rectangles. The converse allows indicator vectors inside the two cubes.
For the final inequality use `u=(x+y)/2`, `v=(x-y)/2` and
`x^T A y=u^T A u-v^T A v`. Each quadratic form has absolute value at
most `2Q(A)` by multilinearity on the cube. Hollow diagonals are needed
for that last extension. For a signing the column norm is exactly
`n sqrt(n-1)`. Thus Bernoulli sampling gives the asserted bound in Q,
with changed absolute constants and no spectral-norm assumption.

For exact-size sampling take `p=2h/n`, with `4<=h<=n/2`.
Conditional on `E={|R|>=h}`, a uniformly chosen `h`-subset of `R` is
uniform among all `h`-subsets of `[n]`. This follows by permutation
invariance, including after conditioning on E. Principal-cap monotonicity
then gives

```math
\Pr(E)\,E_{|T|=h}Q(A_T)\le E Q(A_R).
```

Even Chebyshev suffices for the probability factor: the binomial mean
is `2h`, variance at most `2h`, and `Pr(|R|<h)<=2/h<=1/2`.
Consequently the exact-size lemma is uniform in all such n,h.
It is essential to keep this probability factor; no equality-size
conditioning of an expectation bound is being silently invoked.

Principal monotonicity itself follows by fixing an optimizing spin on
the smaller block and independently averaging every omitted spin:
all omitted and crossing edges have mean zero. Thus this step is in
the actual hollow Boolean cap, not a regularized norm.

## 2. Fixed rational retention and the liminf starting sequence

Write `q=1-s=u/v` in lowest terms. The integer
`r=ceil(log(1/rho)/(-log q))+1` is FIXED once s is fixed, and

```math
\rho q^2\le q^r<\rho.
```

Deleting at most `v^r-1` vertices from an exact minimizing liminf
sequence to obtain order divisible by `v^r` is legitimate bounded
deletion. Principal monotonicity supplies the upper normalized limit
ell; `Q(A_N)>=M_N` and the definition of liminf supply the matching
lower limit. The retained starting parents need not remain exact
minimizers, and the theorem does not claim that they do.

All r nested q-restriction orders are integral and tend to infinity.
Their terminal marginal is exactly uniform in the initial vertex set.
The already banked uniform small-fixed-retention theorem with starting
cap bound `1/2` and target `1/2<2/pi` makes its terminal cap exceed
`1/2` with probability tending to one. No assertion about a typical
selector at every retention is required.

For EVERY possible selected signing at these finitely many orders,
the normalized cap is at least `ell-o_j(1)`. Hence
`u_i=C_i^(2/3)>=lambda-o_j(1)`, where `lambda=ell^(2/3)`.
This is a lower bound uniform over the finite process states, not a
claim that selected matrices minimize their orders.

## 3. Stopping, discounting, and lower expectation

Let tau be first crossing of `b>ell`, or r. Direct cancellation gives

```math
q^tau(u_tau-lambda)-(u_0-lambda)
=sum_(i<r) q^i 1_(tau>i) [q u_(i+1)+s lambda-u_i].
```

The exact `s=1-q` term is required. On a hit, the terminal term is
at least `q^r(b^(2/3)-lambda)`; on no hit it is at least `-o_j(1)`.
The starting difference tends to zero. Therefore the liminf expected
left side is at least the fixed positive constant

```math
G=(rho/4)(b^(2/3)-lambda)
```

for q>=1/2. This lower estimate does not assume monotonicity of
normalized caps, and it does not discard a potentially negative no-hit
term without the liminf lower bound.

## 4. Complement repair and bad-event control

At each pre-stop state the parent has normalized cap below b. Its
uniform complement has size `h=s n_i`. The exact-size lemma and the
all-order upper bound on M_h give, uniformly in the state,

```math
E[R_i | history] <= C_R s^(3/2)+o_j(1),
R_i=[Q(old complement)+M_h]/n_i^(3/2).
```

C_R can be independent of sufficiently small s: the extra parent term
is `O(s^2)`, hence `O(s^(3/2))`. The order threshold may depend on
the already FIXED s. Choosing an exact minimizing replacement D_h is
permitted at every complement order, with its labels transported to
that complement. All bridge entries and the large principal block are
unchanged.

The triangle inequality gives `c(P)<=C_i+R_i`. Concavity of the power
function, with its derivative bounded on the eventual universal
positive lower-cap interval, then yields

```math
J_i:=q u_(i+1)+s lambda-u_i
<=F_i+L_0 R_i+s o_j(1),
```

where F_i is the repaired powered defect divided by n_i. The use of
`M_h^(2/3)/h>=lambda-o_j(1)` has the CORRECT lower-bound direction.
No upper closeness of M_h to ell is required.

On good repairs `R_i<=a_1`, the repaired parent is below `b+a_1<1/2`.
The large child is below `b q^(-3/2)<1/2` by principal monotonicity;
the exact small child is below the fixed all-order upper constant
`U+o_j(1)<1/2`. These give a common positive margin for fixed s.

On bad repairs the particularly important exact inequality is

```math
q u_(i+1)=Q(B)^(2/3)/n_i<=Q(A_i)^(2/3)/n_i=u_i.
```

Therefore `J_i<=s lambda` regardless of the repaired cap. Markov bounds
the bad probability by `C_R s^(3/2)/a_1`, after absorbing fixed-s order
errors. Its possible positive contribution is only `O(s^(5/2))`.
There is no uncontrolled high-cap tail, and negative contributions
only improve this upper estimate.

## 5. Extraction and the nonexchange of limits

Let K_j be the largest positive repaired defect F_i over all good
repairs at all possible pre-stop states of the finite process.
This deterministic finite maximum can be used in every conditional
expectation. Since `sum_(i<r)q^i<=1/s`, the conditional estimates imply

```math
G <= liminf_j K_j/s + L_0 C_R sqrt(s)
                        +(C_R lambda/a_1)s^(3/2).
```

The liminf direction is correct: a pointwise upper bound for the
stopped expectation, combined with its liminf lower bound, gives the
displayed lower bound on liminf K_j. First choose fixed rational s
so the two errors are at most G/2; only then take the order limit.
This gives `liminf K_j>=Gs/2`, hence an actual good repair with
`F_i>=Gs/4` at every sufficiently large starting-sequence index.
Its order is at least `q^r N_j` and therefore tends to infinity.

This use of a finite maximum is an existence extraction, not a proposed
efficient sampling algorithm. Growing r with the starting order would
require a new proof and is not used.

## 6. Exact boundary of the theorem

The small child D is an EXACT minimizer. The large child B and repaired
parent P have uniformly strict-subhalf actual caps, but are NOT shown
to be exact or asymptotic minimizers of their own orders. In particular,

```math
Q(B)^(2/3)+M_h^(2/3)-Q(P)^(2/3) >= kappa_s n
```

does not remain a proved lower bound when Q(B) is replaced by M_m.
The decrease could consume the entire defect. The assertion that
Q(P)>=M_n does not repair this missing comparison. Nor does the
construction treat comparable splits: it requires sufficiently small
fixed s.

Thus the strengthened result genuinely defeats the indicated
pointwise powered-aggregation principle even with an exactly minimizing
small block. It does NOT defeat a recurrence solely for actual minimum
values, and does not settle the original convergence question. No
further selector or numerical extension was pursued in this audit.
