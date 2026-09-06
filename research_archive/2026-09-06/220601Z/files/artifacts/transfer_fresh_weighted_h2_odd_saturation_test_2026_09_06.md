# Odd-Walsh exact saturation test for the weighted H2 seed

Date: 2026-09-06. Bounded follow-up to the stabilized-norm literature audit.
This file proves a finite Boolean feasibility equivalence and identifies
small-dimension obstructions. It proves neither `R=T` nor a uniform strict
`R<T` gap, and does not provide an original full-sign seed below one half.

## 1. Scope and frozen test

Use the precise Walsh-only stabilization

```math
R(B)=\frac12\sup_{k\ge0}
 \frac{\beta(W_{4^k}\otimes B)}{(4^k)^{3/2}},\qquad
B=\begin{pmatrix}-1&2\\2&4\end{pmatrix}.
```

Independent row and column sign/permutation equivalence identifies this
with the prescribed regular order-four Hadamard powers. The already-audited
majorant is `T(B)=5/sqrt(2)`; see
`transfer_adversary_weighted_h2_reciprocal_2026_09_06.md` for its exact
primal and cut-covariance dual certificates.

The concrete test frozen on 2026-09-06 at approximately 20:29 UTC was:

> For some integer r>=1 and N=2^(2r+1), do four Boolean arrays f,h,g1,g2
> satisfy the following two integer linear equations?

```math
 W_Nf=2^r(-g_1+2g_2),\qquad
 W_Nh=2^{r-1}(g_1+2g_2).                         (L)
```

This is a sufficient finite certificate for `R(B)=T(B)`, despite the odd
outer dimension not itself belonging to the prescribed catalyst family.
It is not asserted necessary for asymptotic equality.

## 2. Exact equivalence and transfer to allowed even dimensions

Write `U=W_N/sqrt(N)` and use normalized counting inner products. The
known exact defect identity is

```math
 D=\mathbb E\left(Uf-\frac{-g_1+2g_2}{\sqrt2}\right)^2
 +\mathbb E\left(2Uh-\frac{g_1+2g_2}{\sqrt2}\right)^2,
\qquad
 \frac{\beta(W_N\otimes B)}{2N^{3/2}}
 =\frac5{\sqrt2}-\frac1{2\sqrt2}\min D.             (1)
```

Since `sqrt(N)=2^r sqrt(2)`, the condition `D=0` is exactly (L), with no
rounding or omitted scaling. Parseval then forces

```math
 \langle f,h\rangle=\langle g_1,g_2\rangle=\frac34,
\qquad d_H(f,h)=d_H(g_1,g_2)=N/8.                  (2)
```

For example, the first squared norm in (L) gives
`1=(5-4<g1,g2>)/2`, and the pointwise transformed product is `3/4`.

Suppose (L) holds at one fixed odd dimension d. For any odd dimension e,
the product of optimizing bilinear sign vectors gives

```math
 \frac{\beta(W_{2^{d+e}}\otimes B)}{2(2^{d+e})^{3/2}}
 \ge \frac5{\sqrt2}
       \frac{\beta(W_{2^e})}{(2^e)^{3/2}}.
```

The scalar ratio tends to one along odd e. Indeed Schmidt proves that
there exist Boolean z with `||W_{2^e}z||_infinity/sqrt(2^e)<=mu_e`,
where `mu_e->1`; Parseval yields `E|Uz|>=1/mu_e`, while the spectral
upper bound is one. Since d+e is even, these are allowed outer orders.
The already-proved `R<=T` therefore gives `R(B)=T(B)`.

This tensor argument uses the primary scalar theorem only in its actual
scope; it does not replace a prescribed Walsh matrix by an arbitrary
orthogonal optimizer. Primary source: Kai-Uwe Schmidt,
[Asymptotically optimal Boolean functions, Theorem 1 and Section 2](https://math.uni-paderborn.de/fileadmin-eim/mathematik/AG-Diskrete_Mathematik/Publications-schmidt/pw.pdf).

## 3. Small odd dimensions are rigorously excluded

### Elementary obstructions at N=8 and N=32

At N=8 the second equation requires odd Walsh coefficients, whereas every
Walsh coefficient of an even-length Boolean array is even.

At N=32 the first spectrum consists of `+/-4,+/-12`, so the number of
negative entries of f is even. The second consists of `+/-2,+/-6`, so
that number for h is odd: for dimension at least two, all Walsh
coefficients have the same residue modulo four, equal to minus twice
the Boolean weight. But (2) requires Hamming distance four, which makes
the two weights have the same parity. Contradiction.

### N=128 is excluded by an existing exact covering-radius result

Schmidt's introduction records `mu_d=sqrt(2)` for odd d=1,3,5,7, where
`mu_d=min_h ||U h||_infinity`. Equation (L), however, requires

```math
 |Uh|\in\left\{\frac1{2\sqrt2},\frac3{2\sqrt2}\right\},
```

whose maximum is strictly less than `sqrt(2)`. This excludes d=7 as well
as the two smaller dimensions. The d=7 exclusion imports the primary
paper's cited exact covering-radius result; it is not claimed to follow
from the elementary parity argument.

More quantitatively, for every Boolean quadruple at any of these odd
dimensions, a coordinate attaining `|Uh|>=sqrt(2)` contributes at least
`4(1/(2sqrt(2)))^2/N=1/(2N)` to D. Thus (1) gives the finite gap

```math
 \frac{\beta(W_N\otimes B)}{2N^{3/2}}
 \le\frac5{\sqrt2}-\frac1{4\sqrt2 N}
 \quad (N=8,32,128).                              (3)
```

This bound vanishes with N and is not an asymptotic separation.

### The next dimension already demands a difficult Boolean spectrum

At d=9, N=512, equation (L) requires
`W_N h in {+/-8,+/-24}`. Such h has nonlinearity exactly
`(512-24)/2=244` (the larger level occurs, by Parseval). The primary
construction of Kavut and Yucel gives nonlinearity 242:
[9-variable Boolean Functions with Nonlinearity 242 in the Generalized Rotation Class](https://arxiv.org/abs/0808.0684).
No claim that 242 remains the current record is needed here. This
comparison shows that the frozen exact test would require a substantial
Boolean-spectrum construction plus the compatible partner f, not merely
an easily available bent function. No N=512 infeasibility is asserted.

## 4. Reproducible integer checker and actual run statuses

The script `computations/transfer_fresh_weighted_h2_odd_exact.py` encodes
(L) with OR-Tools CP-SAT. It adds the implied two Hamming constraints (2).
An agreeing input coordinate exists; translation and common reversal
therefore permit the symmetry choice `f(0)=h(0)=1`. If a feasible answer
is returned, every Walsh equation and correlation is independently
recomputed using integers, and the complete witness is printed.

Recorded diagnostic calls:

```text
.venv/bin/python computations/transfer_fresh_weighted_h2_odd_exact.py --r 1 --seconds 10 --workers 1
INFEASIBLE in presolve; 0 branches; 0 conflicts; 0.000299273 seconds.

.venv/bin/python computations/transfer_fresh_weighted_h2_odd_exact.py --r 2 --seconds 30 --workers 2
UNKNOWN; 454337 branches; 25735 conflicts; 30.030 seconds.

.venv/bin/python computations/transfer_fresh_weighted_h2_odd_exact.py --r 3 --seconds 180 --workers 2
Manually interrupted after the d=7 primary obstruction was noticed.
UNKNOWN; 19976 branches; 2764 conflicts; 134.093 seconds.
```

Neither UNKNOWN result is evidence for nonexistence. The exact exclusions
above do not depend on either computation. No N=512 run was launched.

One bounded direct lower-witness check reused the existing script:

```text
.venv/bin/python computations/transfer_adversary_weighted_walsh_search_2026_09_06.py --order 128 --batch 128 --cycles 100 --seed 260907
```

Its best integer bilinear score was 9592, normalized value
`9592/(2*128^(3/2))=3.3118009009479312`, input correlation `0.640625`,
and reciprocal-product L1 defect `0.4755859375`. The script's final
integer recomputation passed. This is weaker than the known stabilized
lower bound 7/2 and supplies no evidence for an upper bound; no new
extremal witness is claimed. The seed and command reproduce the arrays.

## 5. What remains open after this test

At every finite allowed even dimension the target coordinates in (1) are
irrational, whereas normalized Walsh coefficients are rational, so exact
attainment is impossible there. This also does not yield a uniform gap.
The asymptotic alternative remains exactly whether the infimum of D over
all Boolean quadruples and all allowed even Walsh dimensions is zero,
or is bounded away from zero.

An exact odd-dimensional solution of (L) would settle equality for this
one weighted seed. Nonexistence of (L) at every finite odd dimension
would still not settle its asymptotic equality. Even a strict gap for
this weighted seed would not by itself produce an original full-sign
seed with normalized stabilized cap below one half.

## 6. Exact inverse-spectrum diagnostic for all nonzero correlations

The separate script
`computations/transfer_fresh_perfect_dyadic_inverse_classify.py` checks
all Boolean f up to common reversal at orders 4,8,16. For each positive
integer c<=N divisible by four, it uses the unique candidate spectrum
`Wg=c/(Wf)`. A zero coefficient of Wf excludes a nonzero perfect pair;
otherwise divisibility and the inverse Walsh computation are integer
exact. The test `W(Wg) in {+/-N}^N` is necessary and sufficient.

The restriction `4|c` is exhaustive here: all Walsh coefficients of a
Boolean array of even order are even, so their constant nonzero product
is divisible by four. Negative c follows by reversing g.

The reproducible command

```text
.venv/bin/python computations/transfer_fresh_perfect_dyadic_inverse_classify.py
```

returned these exact counts (f is quotiented only by common reversal):

| N | f examined | f with no zero Walsh coefficient | Perfect pairs c<N | Perfect pairs c=N |
| --- | ---: | ---: | ---: | ---: |
| 4 | 8 | 4 | 0 | 4 |
| 8 | 128 | 64 | 0 | 0 |
| 16 | 32768 | 16832 | 0 | 448 |

The c=N cases necessarily have f=g and f bent. Thus no intermediate
nonzero perfect correlation exists at these orders. This is a finite
classification, not an all-order rigidity theorem.

One elementary all-order necessary congruence is slightly stronger than
`4|c`: for N divisible by four, a nonzero perfect pair cannot have
`v_2(c)=3`. Indeed `4|c` makes `d_H(f,g)=(N-c)/2` even, so the Boolean
weights of f and g have the same parity. If both are odd, every Walsh
coefficient is 2 modulo four, so `v_2(c)=2`; if both are even, all Walsh
coefficients are divisible by four, so `v_2(c)>=4`. This recovers the
N=32, c=24 exclusion without using the particular two-level palette.
It does not exclude c=3N/4 at arbitrary growing odd dimensions.

## 7. Exact compression to a full-sign seed, and its normalization

The weighted test also has a precise full-sign realization, namely

```math
 S=\begin{pmatrix}-1&1&1\\1&1&1\\1&1&1\end{pmatrix}.
```

For every real outer matrix H,
`beta(H tensor S)=beta(H tensor B)`: two clone coordinates have identical
rows and columns, so in bilinear maximization their maximizing signs can
be taken equal on each outer fibre. Grouping their multiplicities gives
exactly `B=diag(1,2) [[-1,1],[1,1]] diag(1,2)`.
In particular `R(S)=R(B)`.

The majorants agree too, with an explicit certificate. Let

```math
 L=\begin{pmatrix}1&0\\0&1\\0&1\end{pmatrix},\quad
 D_0=L^TL=\operatorname{diag}(1,2),\quad
 B_0=\begin{pmatrix}-1&1\\1&1\end{pmatrix},\quad
 P_S=\sqrt2 LL^T,\quad
 C_S=L\begin{pmatrix}1&3/4\\3/4&1\end{pmatrix}L^T.
```

Here `S=L B0 L^T`, and `sqrt(2) I>=+/-B0` proves `P_S>=+/-S`.
Its maximal Boolean quadratic value is `5sqrt(2)`. The matrix C_S is a
cut covariance: take the last two signs equal and their correlation with
the first equal to 3/4. Direct multiplication gives
`S C_S S=P_S C_S P_S`; hence the cut-covariance nuclear dual yields
`2T(S)=tr(P_S C_S)=5sqrt(2)`. Thus `T(S)=T(B)` exactly.

This does not put the candidate near the original one-half threshold:
`R(S)>=q(S)=7/2`, so its normalized stabilized cap is at least
`7/(6sqrt(3))>1/2`. Accordingly a strict R/T gap for this seed would
disprove the universal equality conjecture on full sign matrices, but
would still not provide the full-sign below-half seed requested by the
positive transfer route. S is full, not hollow.

## 8. Primary endpoint-factorization check

Theorems 1.3 and 6.1 of Hinrichs,
[Hilbert space factorization and Fourier type of operators](https://www.impan.pl/shop/en/publication/transaction/download/product/89984),
were rechecked directly in this pass. The nonuniform comparison uses
L2-based amplifications, varying operators `T_n:ell_1^n->ell_infinity^n`,
and an outer matrix at the prescribed same gradation n. It neither
proves equality nor provides a counterexample for the present endpoint
`ell_infinity`-to-`ell_1` Walsh supremum with one seed fixed while outer
dimension tends to infinity. No new norm terminology or unproved endpoint
identification is used in this artifact.
