# Independent audit: equal-square amplification and local response variances

Date: 2026-09-05. This is a conditional structural lemma plus a replay
audit of the exact finite certificates supplied by the algebra and
variational tracks. It concerns general bounded-operator-norm signings,
not near-minimizers.

## 1. Exact amplification lemma

Let `C,D` be symmetric full sign matrices of the same fixed order `n`,
with `C^2=D^2`. Suppose a rational symmetric matrix `P` satisfies
`P >= C`, `P >= -C` and

```
t = (1/2) max_(x in {+/-1}^n) x^T P x < beta(D)/2 = b.
```

Here `beta(D)=max_(f,g signs) f^T D g`; the full quadratic convention
is `q_full(K)=(1/2)max_x |x^T K x|`, including the diagonal.
Put `H4=J4-2I4`, `H_s=H4^(tensor a)`, `s=4^a`, and `N=ns`.
For all `a>=1` the two full parents `H_s tensor C` and `H_s tensor D`
have identical squares `s I_s tensor C^2`, and

```
q_full(H_s tensor C) <= s^(3/2) t,
q_full(H_s tensor D) >= s^(3/2) b.
```

For the upper bound, diagonalize the symmetric orthogonal matrix
`H_s/sqrt(s)`. The hypotheses on `P` imply
`sqrt(s) I_s tensor P >= +/- (H_s tensor C)`. Maximizing its quadratic
form block by block gives the displayed upper bound.

For the lower bound choose Boolean `f,g` with `f^T D g=beta(D)`.
The exact same-spin lift `z=(f,f,g,g)` satisfies

```
z^T (H4 tensor D) z = 8 f^T D g,
q_witness(H4 tensor D) = 4 beta(D).
```

Regularity `H_(s/4) 1=sqrt(s/4) 1` then propagates this witness and
gives the lower bound. In particular, using `q_full(D)=beta(D)/2`
without separately checking it is unnecessary and can be false.

Delete the actual diagonal of each full parent, obtaining hollow sign
matrices `A_N^C,A_N^D`. Since every deleted diagonal entry has modulus
one, each full/hollow quadratic optimum differs by at most `N/2`.
Consequently

```
q(A_N^D)/N^(3/2) - q(A_N^C)/N^(3/2)
  >= (b-t)/n^(3/2) - 1/sqrt(N).
```

Thus a strictly positive asymptotic cap gap is proved along common
orders, without asserting convergence of either individual sequence.

## 2. Their normalized squares become indistinguishable

Use full normalization `B_N^C=(H_s tensor C)/sqrt(N)` and the analogous
`B_N^D`. Their squares agree exactly and their operator norms are the
same fixed number `L0=||C||op/sqrt(n)`.

The actual hollow normalization is
`Bhat_N^C=A_N^C/sqrt(N-1)`. Directly,

```
||Bhat_N^C-B_N^C||op
 <= 1/sqrt(N-1) + L0*(sqrt(N/(N-1))-1) = O(N^(-1/2)).
```

The same estimate holds for `D`. Therefore the hollow normalized
operator norms stay bounded and

```
||(Bhat_N^C)^2-(Bhat_N^D)^2||op = O(N^(-1/2)).
```

A continuous rule recovering the quadratic optimum from the square
operator alone would contradict the positive gap in Section 1. This
does not contradict a rule specifically restricted to near-minimizers:
the finite examples here have normalized caps far above that regime.
It also does not give an upper-preserving all-order enlargement rule.

## 3. Exact equality of every local Schur-channel variance

The current finite examples have the form

```
C = [[P0, J_(k,l)], [J_(l,k), E]],
D = [[-P0,J_(k,l)], [J_(l,k), E]],     P0 1=0,
```

where `P0,E` are symmetric sign matrices (in the current examples
`E=-J2`). Hence `C^2=D^2=:Q0`, and the cross block of `Q0` is
`J_(k,l) E`, whose rows are identical. For every entrywise integer power
`r>=1`, the cross block of `R0=Q0^(circ r)` likewise has identical rows,
so it is annihilated on the left by `P0`.

Write `C=K+Ptilde`, `D=K-Ptilde`, where `Ptilde` is supported on the
upper-left block and `K` has zero upper-left block. Then

```
C R0 C - D R0 D = 2(Ptilde R0 K + K R0 Ptilde).
```

The lower-right diagonal blocks of both summands are zero. Their
upper-left diagonal blocks contain `P0` times a matrix with identical
rows, or its transpose, and vanish because `P0 1=0`. Thus

```
diag[C (C^2)^(circ r) C] = diag[D (D^2)^(circ r) D]   (r>=1).
```

This equality persists under scalar normalization and under any common
symmetric Hadamard amplification. Indeed its normalized square is
`I_s tensor (C^2/n)`, and its transported Schur power is
`I_s tensor [ (C/sqrt(n)) (C^2/n)^(circ r) (C/sqrt(n)) ]`.
In particular, all odd Hermite response variances agree root by root,
and so do all convex mixtures used in joint variance normalization.

After actual diagonal removal the equality holds up to a uniform
`O_r(N^(-1/2))` error. To check this without an entrywise-to-operator
shortcut, all full and hollow squares are correlation matrices. For
two correlation matrices `Q,Q'`, the Schur telescoping identity and
operator contractivity of Schur multiplication by a correlation matrix
give

```
||Q^(circ r)-(Q')^(circ r)||op <= r ||Q-Q'||op.
```

Combine this with Section 2 and the fixed operator cap in
`B Q^(circ r) B`. Taking a diagonal entry cannot increase operator norm.

Scope: equality is asserted for the diagonal/local variances, not for
the entire transported matrices `B Q^(circ r) B`. Their cross blocks
can distinguish the examples. Thus the result identifies a limitation
of square data even supplemented by the present local variances, not
a limitation of every possible unmarked response or joint-root state.

More explicitly, under standard Gaussian input, distinct normalized
Hermite degrees have exactly zero covariance, since
`E[h_r(G_j)h_t(G_k)]=1{r=t} Q_jk^r`. Thus, for any fixed finite list
of odd degrees, the entire covariance matrix of the unmarked vector
at a single root is diagonal with the matching entries `v_(r,i)`.
Any two finite polynomial odd responses `h=sum a_r h_r` and
`h'=sum b_r h_r` have matching local covariance
`sum a_r b_r v_(r,i)` as well. This proves joint local covariance
equality, not merely equality for one selected mixture.

For the actual hollow amplified sequences, the bounded-op finite-chaos
theorems already proved in the campaign turn these matching local
covariances into matching limiting one-root joint Gaussian laws, also
jointly with each fixed old marked-tree family. The Gaussian statement
is asymptotic; orthogonality of distinct degrees by itself would not
establish finite-order independence. No equality of covariances at two
different roots is asserted for the transported unmarked fields.

## 4. Finite replay status

The order-fourteen candidate is in
`computations/fresh_limit_cosquare_scalable_certificate.py`. Its initial
assertion equating the high seed's quadratic cap with half its bilinear
cap was caught independently during this audit; the exact `H4` lift
above is the correct replacement. The intended certified numbers are
`beta(C)=72`, `beta(D)=80`, and `t=39.9659395<40`.
The patched replay was independently rerun successfully. It enumerates
all `16384` seed spin inputs, checks both rational positive-definiteness
claims by exact Fraction LDL elimination, and verifies the explicit
56-coordinate lift. Its exact output is

```
t=79931879/2000000,
b-t=68121/2000000,
q_full(D)=34, beta(D)=80, q_witness(H4 tensor D)=320.
```

Consequently the asymptotic normalized cap gap is at least
`68121/(2000000*14^(3/2))>0`. The proof uses no optimizer output as a
certificate and no bilinear-to-same-spin identification beyond the
explicit checked lift.

## 5. Stronger order-twelve certificate

The separate, newer order-twelve pair in
`computations/fresh_cosquare12_strict_certificate.py` is not the earlier
order-twelve diagnostic pair in `fresh_limit_cosquare_certificate.py`.
The newer replay was independently rerun successfully during this
audit. Exact integer enumeration and exact Fraction LDL verify

```
q_full(C)=26, beta(C)=52,
q_full(D)=30, beta(D)=60,
T(C)<=59563/2000=29.7815,
30-T(C)>=437/2000>0.
```

The explicit high-seed quadratic witness has signed energy `-60`, so
this pair does not need the `H4` bilinear lift: regular all-ones outer
replication already supplies the lower value `30` at every outer order.
Sections 1--3 therefore apply, with strict normalized gap at least
`437/(2000*12^(3/2))`, equal full normalized squares, asymptotically equal
hollow normalized squares, and matching local joint odd-channel
covariances/variances. The parent artifact
`fresh_cosquare12_scalable_gap_2026_09_05.md` was read in full; its
normalizations, diagonal deletion bound, and near-minimizer scope pass.

Both sequences lie outside the minimizing regime: even regular
replication of the low seed gives normalized cap at least
`26/12^(3/2)>1/2`. This is therefore a rigorous general square-state
counterexample, not a counterexample to original convergence or to
near-minimizer square-state recovery.
