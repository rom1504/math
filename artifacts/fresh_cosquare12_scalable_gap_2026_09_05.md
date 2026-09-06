# A scalable Boolean-cap gap with identical squared operators

Date: 2026-09-05. Status: exact certificate, independently replayed by
the algebra agent. This is a structural counterexample, not an improved
upper construction for the original minimum and not a nonconvergence
proof.

## 1. The exact full symmetric signing pair

The executable certificate
`computations/fresh_cosquare12_strict_certificate.py` contains two full
symmetric order-12 signings. The first is

`C=[[P,J_(10,2)],[J_(2,10),-J_2]]`,

where `P` is the displayed symmetric balanced signing (`P1=0`). The
second, `C'`, reverses the entire upper-left 10 by 10 block and leaves
all other entries unchanged. Direct block multiplication gives

`C^2=(C')^2`.

All 4096 sign vectors are checked by exact integer arithmetic. For
`Q_full(C)=max_x |x^T Cx|/2` and
`beta(C)=max_(x,y) |x^T Cy|`, the exact values are

`Q_full(C)=26`, `beta(C)=52`,

`Q_full(C')=30`, `beta(C')=60`.

An explicit quadratic witness for `C'` is

`x=(-1,1,-1,-1,1,-1,-1,-1,1,1,1,1)`,

with `x^T C' x=-60`. In this pair the larger value is genuinely a
same-spin quadratic witness; it is not a bilinear value relabeled as a
quadratic one.

## 2. An exact rational upper certificate for C

The same script gives a symmetric rational matrix `D`, with denominator
1000, and proves

`D-C>0`, `D+C>0`

by exact `Fraction` LDL elimination, with all pivots strictly positive.
It also enumerates the entire Boolean cube to verify

`max_x x^T D x/2 = 59563/2000 = 29.7815 <30`.                 (1)

Thus the absolute-PSD majorant norm satisfies

`T(C)<=59563/2000`,

where `T(C)=(1/2)inf_(D>=C,D>=-C) max_x x^T D x`. No numerical
optimizer is used in the replay. A diagnostic SDP suggested the matrix;
the certificate is the fixed rational data and exact checks alone.

## 3. Hadamard amplification preserves a strict gap

Let `H_4=J_4-2I_4`, and for `s=4^k` let `H_s=H_4^(tensor k)`.
Then `H_s` is symmetric, `H_s^2=s I`, and
`H_s 1=sqrt(s) 1`. Define the full symmetric signings

`A_s=H_s tensor C`, `A'_s=H_s tensor C'`, of order `N=12s`.

Their squares are identical exactly:

`A_s^2=(A'_s)^2=s I_s tensor C^2`.                            (2)

The explicit Boolean vector `1_s tensor x` gives

`Q_full(A'_s)>=30 s^(3/2)`.                                  (3)

For the other signing, put `U_s=H_s/sqrt(s)`. Since its eigenvalues are
`+-1` and `D>=+-C`,

`I_s tensor D >= +- U_s tensor C`.

One can check positivity by decomposing into the two eigenspaces of
`U_s`: the two resulting blocks are `D-C` and `D+C`. For an arbitrary
Boolean vector with seed blocks `z_a`, this gives

`|z^T(H_s tensor C)z|/2`
` <= sqrt(s) sum_(a=1)^s z_a^T D z_a/2`
` <= (59563/2000) s^(3/2)`.                                  (4)

Equations (3)--(4) yield a strict order-`N^(3/2)` gap:

`Q_full(A'_s)-Q_full(A_s) >= (437/2000) s^(3/2)`.

In particular the regularized norms themselves satisfy

`R(C')>=30>59563/2000>=R(C)`.

The normalized full operators `B_s=A_s/sqrt(N)` and
`B'_s=A'_s/sqrt(N)` have identical squares

`B_s^2=(B'_s)^2=(I_s tensor C^2)/12`,                         (5)

and uniformly bounded operator norms. Nevertheless their normalized
Boolean caps are separated by at least

`437/(2000*12^(3/2)) >0`.                                   (6)

Thus the squared operator does not determine the regularized Boolean
cap, even within full symmetric flat sign matrices and even with a fixed
operator-norm bound.

## 4. Hollow matrices retain the asymptotic gap

Delete the diagonals, setting

`Abar_s=A_s-diag(A_s)`, `Abar'_s=A'_s-diag(A'_s)`.

These are legitimate symmetric hollow signings of the original problem.
A diagonal changes every quadratic energy by a constant of absolute
value at most `N/2`, so (3)--(4) imply

`limsup Q(Abar_s)/N^(3/2) <= (59563/2000)/12^(3/2)`,

`liminf Q(Abar'_s)/N^(3/2) >= 30/12^(3/2)`.

The strictly positive separation (6) therefore persists. Their squared
operators are no longer exactly equal at finite order, but they become
equal in the strong operator-norm sense:

`|| [Abar_s/sqrt(N-1)]^2-[Abar'_s/sqrt(N-1)]^2 ||op ->0`.     (7)

Indeed (2) cancels the full squares. Each diagonal has operator norm
one, its square is the identity, and
`||A_s||op=||A'_s||op=sqrt(s)||C||op`. Expanding both hollow squares
bounds their difference before normalization by
`4 sqrt(s)||C||op`, which is `O(sqrt(N))`; division by `N-1` proves (7).

This explicitly refutes a proposed continuous determination of the
normalized Boolean cap by the squared operator on bounded-operator
hollow signings. It does not refute a richer state that retains signed
operator action or higher noncommutative data.

## 5. Scope and reproduction

Both normalized caps in this construction are above the original
all-order upper constant `1/2`; in particular
`30/12^(3/2)=5/(4sqrt(3))`. For the lower-cap sequence, regular replication
of its seed witness also gives the lower bound `26/12^(3/2)>1/2`.
Neither sequence is claimed to minimize
the original problem. Consequently this example does not show
`liminf c_n<limsup c_n`, nor does it improve the original upper bound.

The obstruction is specific: bounded operator norm and complete squared-
operator information are insufficient to recover the Boolean cap, even
after allowing Hadamard regularization. The new Gaussian lower-bound
identities remain valid on both sequences; they never asserted such
completeness.

Reproduce the exact finite certificate with

`.venv/bin/python -B computations/fresh_cosquare12_strict_certificate.py`.

The algebra agent independently read the entire script and reran all
4096 quadratic/bilinear/certificate checks and both rational LDL proofs.
The director and literature agent separately repeated that exact replay and
checked the tensor and diagonal-removal arguments. The accompanying
`fresh_cosquare_asymptotic_and_variance_audit_2026_09_05.md` proves that
all fixed odd-channel one-root Gaussian covariance data also agree, without
claiming equality of cross-root transported operators or full action limits.
