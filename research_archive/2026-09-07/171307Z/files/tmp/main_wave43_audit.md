# Main-agent Wave 43 derivations

## A. Incidence-support extraction (independent derivation)

Let `U` be the full uniform selector slice and let `P(S,D)` be any joint
law supported on `S in I_D`, where `D` is a full parent projective cut.  Put

```
Htot = D(P_S || U) + I(S;D),   Rbar = E R_2(D).
```

The KL chain rule gives

```
Htot = D(P_{S,D} || U x P_D)
     = E_D D(P_{S|D} || U).
```

For each output `d`, support on `I_d` and the KL chain rule relative to
`U(.|I_d)` give

```
D(P_{S|d} || U)
 = D(P_{S|d} || U(.|I_d)) + log(1/U(I_d))
 >= log(1/U(I_d)).
```

Thus `E_D log(1/U(I_D)) <= Htot`.  Applying the elementary averaging
principle to the two nonnegative normalized costs shows that some output in
the support has

```
log(1/U(I_d)) <= 2 Htot,   R_2(d) <= 2 Rbar,
```

with the zero-denominator cases handled separately.  Hence if
`Htot=O(n^(3/4-c))` and `Rbar=O(n^(9/4-c))`, that one cut has project row
and incidence mass `exp{-O(n^(3/4-c))}`.  In the high-ratio window,
(10.1045) says all of `I_d` is favorable, so (10.795) holds with `t=0`.

Conversely, a single row-good `d` with `U(I_d)=u` yields such a kernel by
taking `P_S=U(.|I_d)` and `D=d`, with `Htot=log(1/u)`.  Therefore the first
two clauses of (10.1089) are already equivalent, at the target scale, to a
single row-good complement incidence column.  Its conflict clause and the
colored-codegree implementation are logically unnecessary for the bare-tail
recurrence as currently stated.

Audit points: `D` must really be the full parent cut used in `I_d`; if it is
only a local label, the argument does not apply.  The definitions in
(10.1087)--(10.1090) do use a full projective word.  Anchoring costs only a
constant density and is already included in `D(P_S||U)`.

## B. Aggregate strict-pressure extraction

For uniform selectors `S` and uniform full spin centers `z`, put

```
d(S,z)=Q(A[S])-|z_S^T A[S] z_S| >= 0,
q(z)=E_S exp(-theta d(S,z)),
Z=E_z q(z).
```

Tilt the joint uniform law by `exp(-theta d)/Z`.  Its relative entropy is

```
D(Ptilt || U_S x U_z) = -theta E_tilt d - log Z <= -log Z.
```

Data processing gives the same upper bound for the `z` marginal.  Applying
the verified row entropy transport (10.870) with `K=-log Z` gives

```
E_tilt R_2(z)
 <= n(n-1)+O(n^(7/4)sqrt(K)+n^(3/2)K).
```

Thus `K=O(L0)`, `L0=n^(3/4-c)`, implies project-scale average row for
`0<c<1/4`.  Markov gives a row-good set carrying at least half of tilted
mass, and the identity `Ptilt_z(z)=U_z(z)q(z)/Z` extracts a row-good center
with `q(z)>=Z/2`.  If `theta H=bL0` and `Z>=exp(-aL0)` with `a<b`, the bounded
soft-to-hard argument yields

```
U_m{d(S,z)<=H} >= exp(-(a+o(1))L0).
```

This is a valid scalar sufficient lemma, but it is stronger than the
one-exceptional-center target because averaging over all centers can lose
order `n` entropy.  Exact minimality has not been shown to establish the
scalar premise.

## C. Integrated harmonic context-cost bound

For the vertex coordinates, let

```
L_s(D)=sum_i k_{i,D_-i}(s)
```

be the local conditional KL load from (10.1081).  On the endpoint selector
posterior define the symmetric context cost

```
C_i(e)=-(1/2) log(r_i(d^+) r_i(d^-)) >= 0,
Ctot(D)=sum_i C_i(D_-i).
```

The reveal/affinity/cross-level factorization (10.1144) decomposes `C_i`.
The exclusion-ratio identity gives

```
|chi_i(e)| = |log r_i(d^+) - log r_i(d^-)| <= 2 C_i(e).
```

For a binary exponential tilt, its time-`s` KL obeys

```
0 <= k_i(s) <= s |chi_i| <= 2s C_i.
```

Consequently `0<=L_s<=2s Ctot` pointwise, but this raw-cost majorant has a
fatal baseline discussed below.  From the exact covariance collapse
(10.1081), Cauchy--Schwarz gives the sharper intrinsic bound

```
A := [-int_0^1 Cov_mu_s(g,L_s) ds]_+
 <= sqrt(H) J_L,
H = int_0^1 s Var_mu_s(g) ds = Ent_nu(f),
J_L^2 = int_0^1 Var_mu_s(L_s)/s ds.
```

Writing

```
X(D)=sum_i |chi_i(D_-i)|,
```

one has `L_s<=s X` and hence the baseline-free majorant

```
J_L^2 <= int_0^1 s E_mu_s X(D)^2 ds =: J_chi^2.
```

Thus `J_L^2=O(a_n)` (or the stronger `J_chi^2=O(a_n)`), together with
endpoint vertex/orientation cost
`O(a_n)` and the existing restoring input, closes the same quadratic
bootstrap for `H=O(a_n)`.  Unlike a pointwise score bound, this weights large
reveal costs by the actual interpolation context mass; the A9 bad edge is
therefore exponentially suppressed.  This is a sufficient target, not a
proof of its project-scale estimate.

The tempting raw reveal-cost target

```
int_0^1 s E_mu_s Ctot(D)^2 ds = O(a_n)
```

is falsified even by the constant likelihood `f=1`: then
`r_i=(n-m)/n`, so every `chi_i`, `k_i`, and `L_s` vanishes, but
`C_i=-log((n-m)/n)=Theta(1)` and `Ctot=Theta(n)`.  Therefore the three-factor
cost can only be used after removing this common omission baseline, or via
its actual endpoint variation `|chi_i|`; its uncentered sum is not a viable
project target.

Audit points: the orientation coordinate is not covered by omitted-vertex
probabilities and remains separate.  At finite temperature every completion
weight is positive, so every `r_i` is positive and `C_i` is finite.  The
intrinsic load variance avoids both an avoidable sum of coordinatewise L2
norms and the uncentered omission baseline.
