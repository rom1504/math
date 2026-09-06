# Transitive majorant floor and a scalable weighted Cayley obstruction

Date: 2026-09-06. These two elementary statements distinguish the
majorant realization gap from the actual sign-entry constraint. They
do not prove an amplification identity or original convergence.

## 1. Majorant floor for every transitive real symmetric seed

For a real symmetric k-by-k matrix B, define

\[
T(B)=\inf_{P-B\succeq0,\ P+B\succeq0}
\frac12\max_{x\in\{\pm1\}^k}x^TPx.
\tag{1}
\]

The constraints imply `P>=0`. Suppose a group of permutation matrices
acts transitively on the coordinates and preserves B by conjugation.
Then

\[
\boxed{\qquad T(B)\ge\frac{\sqrt k}{2}\|B\|_F.\qquad}
\tag{2}
\]

No sign-entry or hollow hypothesis is needed. Average a feasible P over
the transitive group. Feasibility is preserved, its Boolean maximum
cannot increase by convexity, and its diagonal becomes a constant p.
If p=0 then P=0 and B=0; otherwise the matrix `P/p` is a Gaussian
covariance with unit diagonal. Taking the signs of such a Gaussian
vector and using the elementary arcsine covariance identity gives

\[
\begin{aligned}
\max_x x^TPx
&\ge kp+\frac{2}{\pi p}\sum_{i\ne j}P_{ij}^2\\
&=kp\left[1+\frac2\pi(r-1)\right]
\ge kp\sqrt r=\sqrt k\|P\|_F,
\end{aligned}
\tag{3}
\]

where `r=||P||_F^2/(kp^2)>=1`. The first inequality follows since
`u arcsin(u)>=u^2`; the last follows from `2/pi>1/2` and
`(r+1)/2>=sqrt(r)`. The arcsine identity itself is the two-dimensional
rotational-symmetry formula for the signs of jointly normal variables;
the same rounding identity is used in the campaign's elementary
Grothendieck proof.

Finally,

\[
\|P\|_F^2-\|B\|_F^2
=\operatorname{Tr}[(P-B)(P+B)]\ge0,
\tag{4}
\]

since both factors are positive semidefinite. Equations (3)--(4) prove
(2). In particular, every full sign transitive seed has
`T(B)>=k^(3/2)/2`.

The equality case is also rigid: for B nonzero, equality in (2) holds
if and only if `B^2=(||B||_F^2/k)I`. Indeed, an optimum exists by
bounding `Tr P` on a minimizing level set and compactness. The strict
inequality in (3) unless `r=1` forces an averaged optimal P to be `pI`.
Then `|lambda(B)|<=p`, while equality of Frobenius norms forces every
eigenvalue to have magnitude p. Conversely this flat spectrum admits
the majorant `pI` and attains (2).

## 2. Weighted Cayley symmetry and Frobenius mass alone cannot reach one-half

Let A_9 be the Paley core on the additive group `F_3^2`. An explicit
description not requiring a field multiplication convention is as
follows. Give its four projective frequency lines eigenvalues
`3,3,-3,-3`, and give the zero frequency eigenvalue zero. Fourier
inversion gives a real even kernel with value zero at zero and values
`+1,-1` everywhere else. Exactly,

    A_9^2=9I-J,    ||A_9||_F^2=72,    Q(A_9)=12.

The last assertion is independently replayed by all 256 projective
Boolean spins in
`computations/resumed_convergence_weighted_paley9_barrier_verify_2026_09_06.py`.

For any `N=3^r`, set

\[
C=A_9\otimes I_N,
\qquad q=9N.
\tag{5}
\]

This is a hollow real symmetric Cayley matrix on `F_3^(r+2)`.
Its kernel is the A_9 kernel on the first factor times the delta
function at zero on the second. It is a block sum of N copies of A_9,
so exactly

\[
\|C\|_F^2=72N=8q,
\qquad Q(C)=12N=\frac43q,
\qquad
\boxed{\quad
\frac{Q(C)}{\sqrt q\|C\|_F}=\frac{\sqrt2}{3}
<\frac12.\quad}
\tag{6}
\]

This is a scalable obstruction to proving a one-half Frobenius lower
bound from hollow weighted Cayley structure alone, even with all
finite Boolean quotient tests. Rescaling C to the Frobenius mass of
a signing preserves this strict ratio, but creates large nonzero
entries and many zeros. It is NOT an actual hollow sign matrix of the
original problem. The real-space sign/flat-entry condition would have
to reenter any proof of one-half for all Cayley signings.

In particular, this example does NOT falsify a statement that also
requires bounded entries and near-full Frobenius mass:
`||C||_max<=1` together with `||C||_F^2=(1-o(1))q^2`.
The unscaled example has Frobenius mass only `8q`; the rescaling to
mass `q(q-1)` violates the entry bound. The obstruction is to the
unrestricted homogeneous weighted statement, not to that dense bounded
version, which would retain essentially the sign-entry information.

At the same time, (2) guarantees the one-half floor for the majorant of
every member of this transitive family. Already at order nine,

    Q(A_9)=12 < 9sqrt(2) <= T(A_9).

Thus the failure occurs at actual Boolean-cap realization of the
majorant, not at its universal lower mass. Nothing here proves that an
available Hadamard amplification realizes T or composes its modules.

## 3. What the weighted extension DOES transfer to nearby actual signings

The characteristic-three incidence proof remains valid for any hollow
real Cayley C after replacing its signing Frobenius identity by the
actual `||C||_F^2`. With the exact finite constant `c_(3,q)` in that
theorem, it gives

    Q(C)>=c_(3,q) q ||C||_F/sqrt(q-1).

If B is an arbitrary hollow symmetric real matrix, its translation
average P(B) is Cayley, `Q(P(B))<=Q(B)`, and P is an orthogonal
projection in Frobenius inner product. For any characteristic-three
Cayley signing A, Cauchy--Schwarz therefore gives the valid lower test

\[
\frac{Q(B)}{q^{3/2}}
\ge c_{3,q}\frac{|\langle A,B\rangle_F|}{q(q-1)}.
\tag{7}
\]

In particular, d coefficient flips from A give the factor
`c_(3,q)|1-4d/[q(q-1)]|`; `o(q^2)` flips preserve the asymptotic
constant `2080/(9sqrt(269441))`. The full independently audited proof
is in
`resumed_bound_audit_characteristic_three_projection_stability_2026_09_06.md`.
Unlike the exact Paley radial ground-state law, (7) is a projection
inequality and does not assert a covariance formula for A's ground
states.
