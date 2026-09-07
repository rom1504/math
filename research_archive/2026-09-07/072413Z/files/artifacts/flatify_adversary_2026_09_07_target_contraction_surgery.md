# Actual sign recovery by contracting the target before rounding

Date: 2026-09-07. Status: **proved; independently reconstructed by root
and the constructive agent**. This supersedes the leverage-mask/Schur-bias
implementation for nuclear-budget perturbations. It needs no Grothendieck
inequality, no incoherence, and no a priori cap bound on the input signing.
It is a same-order operation, not a completed cross-order recurrence.

## 1. General theorem

Let A be any hollow full signing of order n. Let Delta be ANY real
symmetric matrix and put L=||Delta||_*=Tr|Delta|. There exists an actual
hollow full signing A' such that, with a=(n+2)log2,

    Q(A')<=Q(A-offdiag(Delta))+2sqrt(n L a)+4a/3.       (1)

Moreover A' can be chosen to change at most 2nL edges of A simultaneously
with (1). The edit bound is harmless if it exceeds the total edge count.
For L=0 use A'=A directly. In particular L=o(n) gives o(n^(3/2)) cap
error. No rank or operator-norm bound on Delta is required.

More generally, the SAME statement and bounds hold whenever supplied
nonnegative numbers ell_i satisfy |Delta_ij|<=sqrt(ell_i ell_j) for i!=j,
with L=sum_i ell_i. The nuclear choice ell_i=|Delta|_ii is one option;
the entrywise choice ell_i=max_(j!=k)|Delta_jk| is another. This unifies
the entrywise-small and nuclear-budget operations without paying bias.
The general envelope budget always obeys
||offdiag(Delta)||_F<=L, so it does not evade the flatification scale gap.

### Coefficient feasibility by diagonal target contraction

For the nuclear version set ell_i=|Delta|_ii. The spectral absolute value satisfies

    |Delta_ij|<=sqrt(ell_i ell_j).

Define s_i=(1+ell_i)^(-1/2), S=diag(s_i), and the hollow target
F=A-offdiag(Delta). The mean coefficient matrix B=SFS obeys

    |B_ij|<=s_i s_j(1+|Delta_ij|)<=1,

because

    (1+ell_i)(1+ell_j)>=(1+sqrt(ell_i ell_j))^2.

Thus B can be rounded independently, entry by entry, to actual signs
with exactly that mean. Crucially it is the TARGET F that was contracted.
Separate affinity on the cube proves Q(B)<=Q(F), with no bias payment.

### Variance and full-cube error

Let A'_ij be independent signs of means B_ij for i<j, and mirror them.
Their centered coefficients are bounded by 2. For each spin vector the
variance sum is the same V=sum_(i<j)(1-B_ij^2). Direct expansion yields

    1-B_ij^2
      <=(1-s_i^2s_j^2)+2s_i^2s_j^2 A_ij Delta_ij
      <=ell_i+ell_j+2|Delta_ij|.

Since sum_i ell_i=L and (sum_i sqrt(ell_i))^2<=nL,

    V<= (n-1)L + 2sum_(i<j)|Delta_ij|
      <=2(n-1)L<=2nL.                                (2)

Bernstein with threshold 2sqrt(nLa)+4a/3 bounds the two tails for each
spin by 2exp(-a). The union bound over 2^n spins has failure probability
at most 1/2. Some actual signing therefore has Q(A'-B) below this
threshold; triangle inequality and Q(B)<=Q(F) prove (1).

### Actual number of changed edges

The probability of changing edge ij is (1-A_ij B_ij)/2. Also

    1-s_i s_j <= (ell_i+ell_j)/2,

by applying 1-(1+x)^(-1/2)<=x/2 separately to the two factors.
Consequently

    E[number changed]
      <=(1/2)sum_(i<j)[1-s_i s_j+|Delta_ij|]
      <=(n-1)L/2<=nL/2.

Markov gives probability at most 1/4 of exceeding 2nL edits. Intersecting
with the rounding event has positive probability, proving the joint claim.

## 2. Finite-rank and nuclear-budget stationarity, with improved error

Take Delta=theta sqrt(n) G and tau=||G||_*. Equation (1) gives normalized
error

    e_n=2sqrt(theta tau a)n^(-3/4)+4a/[3n^(3/2)].       (3)

For fixed theta,tau this is O(n^(-1/4)), improving the earlier masked
construction's O(n^(-1/6)). For fixed theta, tau=o(sqrt(n)) still suffices
for e_n=o(1). In particular G=UTU^T with ||T||op<=1 has tau<=r and the
bound is uniform over every rank-r feature space without incoherence.

If Q(A)<=M_n+delta n^(3/2), actual sign optimality implies

    Q(A-theta sqrt(n) offdiag(G))
      >=Q(A)-(delta+e_n)n^(3/2).                      (4)

For a fixed rank-r U, minimizing over the operator-norm unit ball of T
and applying the previously audited minimax argument gives a law on
(sigma,x) with p=U^Tx/sqrt(n), h=H_A/n^(3/2), q=Q(A)/n^(3/2), such that

    E[q-sigma h]+(theta/2)||E[sigma p p^T]||_*
      <=delta+e_n+theta r/(2n).                       (5)

Alternatively use the FULL symmetric nuclear-norm unit ball of G. Its
duality gives ONE full-space law satisfying

    E[q-sigma h]+(theta/2)||E[sigma xx^T/n]||op
      <=delta+e_n+theta/(2n),                        (6)

with tau=1 in (3). This is operator norm, not nuclear norm, and does not
force the two polarity masses or the energy midpoint to balance.

## 3. Weighted extension

Let C_ij=a_ij A_ij have nonzero amplitudes in [b,B0], 0<b<=B0. Require
Delta_ij=0 on every forbidden zero edge. Put ell_i=|Delta|_ii and
s_i=(1+ell_i/b)^(-1/2). Then the contracted target

    B_ij=s_i s_j(C_ij-Delta_ij)

lies in [-a_ij,a_ij], and Q(B)<=Q(C-offdiag(Delta)). Round independently
to the SAME permitted amplitudes +/-a_ij. Centered coefficients are
bounded by 2B0, and their variance sum is at most

    (n-1)L(B0^2/b+B0).                               (7)

This follows from the identical expansion, bounding the first term by
(B0^2/b)(ell_i+ell_j) and the cross term by 2B0|Delta_ij|.
Bernstein therefore gives cap error O_(b,B0)(n sqrt(L)+n). Thus the
weighted-profile stationarity extension also holds without masking or
Schur bias, under its explicit zero-pattern restriction.

The weighted proof also works with any supplied off-diagonal magnitude
envelope ell, not only the spectral choice. The constructive agent
independently reconstructed this weighted extension, including the exact
(n-1) factor in (7).

## 4. What this does not solve

The smallness condition is nuclear perturbation budget L=o(n), not
merely a small number of named modes or a formal variance interpolation.
For a direct comparison of a uniform full signing with the row-regular
weighted two-child target, every fixed full-sign comparison has a
Frobenius discrepancy of order n, so its nuclear budget cannot be o(n).
The theorem by itself therefore does not provide favorable flatification.

It is instead a rigorously implementable same-order sign operation, and
an optimizer-specific stationarity theorem derived from actual sign
optimality. Any further composition application must exhibit the needed
small nuclear correction or another separately paid mechanism.
