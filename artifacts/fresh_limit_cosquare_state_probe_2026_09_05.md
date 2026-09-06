# Cosquare sign seeds: a finite falsifier and an exact amplification target

Date: 2026-09-05. This probes whether squared-operator data could supply
the missing recovery theorem for original convergence. Sections 1--3
record the first finite candidate. The stronger construction in
Section 4 now proves an asymptotic cap separation.

## 1. Exact finite counterexample

The fixed matrices `C,Cprime` in
`computations/fresh_limit_cosquare_certificate.py` are symmetric full
sign matrices of order twelve. In block notation,

`C=[[P,J_(10,2)],[J_(2,10),D]]`,

`Cprime=[[-P,J_(10,2)],[J_(2,10),D]]`,

where `P1=0` and `D=[[-1,-1],[-1,1]]`. Therefore their squares agree
exactly: the diagonal blocks are unchanged, and `PJ=0` makes the
off-diagonal blocks unchanged. Exhausting all 4096 Boolean inputs gives

`Q_full(C)=27`, `beta(C)=54`,

`Q_full(Cprime)=29`, `beta(Cprime)=58`.

The convention is `Q_full(A)=max_x |x^T A x|/2`, including the diagonal.
Thus neither the full quadratic optimum nor the bilinear optimum is a
function of the square of a symmetric full sign matrix.

Earlier exact enumeration through orders five and six found no such
bilinear difference: 16,384 and 1,048,576 globally signed-normalized
matrices yielded 8,224 and 617,424 square classes, respectively. Those
failed smaller searches are recorded separately; they are not used in
the proof. The order-twelve witness was discovered at trial two of the
fixed-seed regular-block search `(k,m,seed)=(10,2,20260905)`.

## 2. Exact tensor-stable upper certificate

Let `T(C)=min_(R>=C,R>=-C) max_x x^T R x/2`. The fixed rational matrix
`R=P_SCALED/1000000` in the replay script satisfies both positive
definiteness assertions by exact rational LDL elimination. Exhausting
the Boolean cube also proves

`T(C) <= 60806691/2000000 = 30.4033455`.

No numerical optimizer is needed to check these claims. Numerical SDP
output was used only to suggest the rational certificate.

For every symmetric Hadamard outer `H_s`, diagonalizing `H_s/sqrt(s)`
shows

`sqrt(s) I_s tensor R >= +/- (H_s tensor C)`.

Hence the proved upper bound is tensor-stable:

`Q_full(H_s tensor C) <= s^(3/2) * 30.4033455`.

## 3. The remaining exact witness target

An explicit Boolean vector of length 192 with

`Q_witness(H4 tensor H4 tensor Cprime)>=1946`

would suffice for a genuine amplified separation, since

`1946/64 = 30.40625 > 30.4033455`.

Here `H4=J4-2I4`. The difference of the two displayed rational values is
strict and exact. A numerical search below that target would not be an
upper bound and would not establish or refute the separation.

For clarity, the tempting order-48 target `Q>=244` is not the practical
target: a diagnostic SDP for `Cprime` has value about `30.48518`, below
`244/8=30.5`. Only the upper certificate for `C` has been made rigorous
in the current replay script.

If the order-192 witness is obtained, regular outer replication of it
gives a positive normalized cap gap at arbitrarily large common orders.
The two full parents have exactly the same squares. Removing their
diagonals changes their normalized quadratic caps by `O(N^(-1/2))`;
their normalized operator norms remain bounded, so the difference of
their normalized squared operators also tends to zero in operator norm.
That would disprove a general continuous square-state closure theorem.
The resulting caps would be well above the minimizer regime, however,
so it would not itself rule out a closure theorem restricted to
near-minimizers, nor settle convergence.

The stated order-192 witness was not banked. A stronger seed below
bypasses that witness target entirely.

## 4. Successful stronger seed: a scalable exact falsifier

Scanning the same zero-row-sum block construction further found a
different order-twelve pair at trial four. Its exact matrix and rational
certificate are saved in
`computations/fresh_cosquare12_strict_certificate.py`, written by the
variational agent and independently read and replayed by the algebra
agent. Every assertion below uses integer cube enumeration and exact
Fraction LDL, not an optimization status.

Call this new pair `D,Dprime` to distinguish it from Sections 1--3.
Then

`D^2=Dprime^2`, `Q_full(D)=26`, `Q_full(Dprime)=30`,

`beta(D)=52`, `beta(Dprime)=60`,

`T(D)<=59563/2000=29.7815<30`.

Here the high seed's displayed quadratic witness has signed energy
`-30`, so even the one-step bilinear lift is unnecessary. If `H_s` is
any regular symmetric Hadamard outer, with `H_s 1=sqrt(s) 1`, tensoring
that witness with the all-one vector proves

`Q_full(H_s tensor Dprime)>=30 s^(3/2)`.

The rational majorant proves simultaneously

`Q_full(H_s tensor D)<= (59563/2000) s^(3/2)`.

The full parents have exactly equal squares `s I_s tensor D^2`.
At common order `N=12s`, their normalized cap gap is at least

`delta=(437/2000)/12^(3/2)>0`.

Removing the two parent diagonals costs at most `N/2` in either cap.
Their normalized operator norms are bounded uniformly in `s`; hence
the normalized squared operators of the hollow parents differ by
`O(N^(-1/2))` in operator norm. The normalized cap gap stays at least
`delta-o(1)`. Thus even bounded-operator symmetric hollow signings can
have asymptotically equal squared operators but separated Boolean caps.

There is a further equality relevant to the current Gaussian lower
modules. For every fixed positive odd integer `r`, the profiles

`diag(B (B^2)^(circ r) B)`

are equal for the two full parents, and differ uniformly by `o(1)`
after diagonal deletion. In the finite block construction, every
sign-sensitive diagonal term contains the zero product `P1`; the
remaining diagonal terms involve `P` twice. The literature agent is
writing the independent detailed asymptotic audit of this observation.

The separating upper/lower thresholds are approximately `.716` and `.722`,
not claimed limits of the two caps. The low seed's replicated witness
also gives a lower bound `26/12^(3/2)>.5`. Both families are above the
known minimizer upper bound `.5`. This is therefore a
falsifier of a **general** square-state completeness claim, not of a
claim restricted to near-minimizers. It supplies no nonconvergent
sequence of minimizer values and no upper-preserving all-order
enlargement operation.

The completed independent audit is
`fresh_cosquare_asymptotic_and_variance_audit_2026_09_05.md`; the exact
successful order-12 construction is presented separately in
`fresh_cosquare12_scalable_gap_2026_09_05.md`.
