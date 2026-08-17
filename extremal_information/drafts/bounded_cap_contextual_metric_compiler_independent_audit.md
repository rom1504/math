# Independent audit: bounded-cap contextual metric compiler

**Verdict:** PASS WITH ONE MINOR REPAIR.  The construction, all normalization
constants, absolute-cap estimates, projective factors, edge ownership, and
linear-order scale conversion are correct.  The only formal defect is that
the proof of BCX.1 samples `floor(exp(gamma n))` words while claiming at
least `exp(gamma n)` words.  Use `ceil(exp(gamma n))`, or sample at exponent
`2gamma` and state the lower exponent `gamma`.

I ran
`experiments/verify_bounded_cap_contextual_metric_compiler.py`.  It passes
all exact Walsh, parent-cap, channel, metric, and ownership regressions at
orders `n=4,16` (seven metric pairs).

## 1. Regular Walsh conventions

Let `W` be the raw Walsh matrix of order `n=2^(2m)` and
`b(u,v)=(-1)^(u dot v)`.  The elementary Walsh transform gives
`Wb=q b`, where `q=sqrt n`.  Therefore

```math
mathcal H=D_bWD_b
```

is symmetric, satisfies `mathcal H^2=nI`, and obeys
`mathcal H 1=q1`.  Diagonal conjugation preserves trace, while

```math
tr W=(tr [[1,1],[1,-1]])^(2m)=0.
```

Hence hollowing introduces no Boolean calibration:

```math
H_A(x)=x^Tmathcal Hx/2.
```

The spectral bound gives `Q(A_s)<=qn/2`, and `x=s` switches to the positive
eigenvector `1` and attains equality.  BCX.1--BCX.3 are correct.

## 2. Hanson--Wright code

For uniform Boolean `Z`, the quadratic form has mean
`E Z^Tmathcal HZ=tr mathcal H=0`.  The relevant norms are exactly

```math
||mathcal H||_F=n,
\qquad ||mathcal H||_(2->2)=q.
```

At threshold `qn/4`, the two Hanson--Wright parameters are

```math
{q^2n^2\over16n^2}={n\over16},
\qquad {qn\over4q}={n\over4}.
```

Thus the failure probability is `2exp(-c_1 n)`.  Pair products of two
independent uniform switches are uniform, so a union bound over
`O(exp(2gamma n))` pairs succeeds for sufficiently small `gamma`.  The
Rayleigh condition excludes equal and antipodal words because both give
absolute Rayleigh value `qn`.

**Required repair:** `M=floor(exp(gamma n))` is strictly less than the lower
bound claimed in BCX.4.  Replacing `floor` by `ceil` changes the union bound
only by a constant factor and fixes the theorem.  No asymptotic constant or
later argument changes.

## 3. Exact parent and channel reduction

The parent is genuinely hollow and complete:

* `A_s` supplies every old--old sign;
* `B_t=t1_q^T` supplies every old--new sign;
* `J_q-I_q` supplies every new--new sign.

The ownership is modular: these blocks depend respectively only on the
child, only on the query, and on public data.  There is no joint `(s,t)`
edge.

After omitting the clique, for fixed old spin the absolute value is convex
in `1 dot y`, so its maximum occurs at `+-q`.  Switching variables by `s`
and using the evenness of the child gives exactly the two channels

```math
R_\pm(w)=\max_u\{\pm u^Tmathcal Hu/2+q w^Tu\}.
```

Thus BCX.14 does account for both the endpoint sign and the outer absolute
value; no channel is missing.

## 4. Complete-square constants

For the positive channel,

```math
K_-=2qI-mathcal H,
\qquad K_-^{-1}=(2qI+mathcal H)/(3q^2),
```

because `mathcal H^2=q^2I`.  Completing the square on the sphere containing
the Boolean cube yields

```math
R_+(w)
\le qn+{1\over2}(qw)^TK_-^{-1}(qw)
=qn\left(1+{2+rho(w)\over6}\right).
```

The negative channel has `2-rho` in the numerator.  If
`|rho|<=1/4`, both are at most

```math
qn(1+(2+1/4)/6)=11qn/8.
```

The factor `11/8` is therefore exact for this relaxation.

For the positive clique,

```math
H_C(y)=((1\cdot y)^2-q)/2.
```

Since `q` is even, its most negative value is `-q/2`, while its positive
endpoint is `E_q=q(q-1)/2`; hence `Q(C)=E_q` (including `q=2`).  Adding it
costs at most `E_q` in absolute cap.

On the diagonal, `w=1`, the child contributes `qn/2` and the rank-one cross
block contributes at most `qn`; `x=s,y=1` attains both simultaneously and
also attains `E_q` in the clique.  This proves the exact diagonal value
`3qn/2+E_q`.  Off diagonal, the preceding upper bound proves
`11qn/8+E_q`.  The gap is exactly `qn/8=n^(3/2)/8`.

## 5. Projective metric and original distance

At query `s`, the response difference `F_s-F_t` is at least `qn/8`; at
query `t` it is at most `-qn/8`.  Its oscillation is therefore at least
`qn/4`, and the factor `1/2` in projective distance gives

```math
d_C(s,t)>=qn/8.
```

Changing the child changes only its old block.  Absolute maxima are
nonexpansive under a pointwise perturbation, so

```math
|F_s(r)-F_t(r)|<=Q(A_s-A_t)=d_0(s,t).
```

All response-coordinate differences lie in `[-d_0,d_0]`; half their
oscillation is at most `d_0`.  Finally
`d_0<=Q(A_s)+Q(A_t)=qn`, so the absolute lower gap implies
`d_C>=d_0/8`.  This verifies both sides of BCX.23 and all projective factors.

The direct lower bound is also correct.  At `x=s`,

```math
H_(A_s)(s)-H_(A_t)(s)
={qn-w^Tmathcal Hw\over2}
={1-rho(w)\over2}qn
\ge3qn/8.
```

## 6. Cap and order scaling

All parents satisfy

```math
Q(P_(s,t))<=3n^{3/2}/2+O(n).
```

Their order is `N=n+sqrt n`, so `N/n->1`; the cap remains `O(N^(3/2))`,
the fixed `n^(3/2)/8` response gap remains a fixed `N^(3/2)` gap, and
`log|mathcal S|=Omega(n)=Omega(N)`.  The construction therefore has genuine
linear response rate at the total parent scale.

## 7. Scope

This is a positive theorem at the contextual **metric** level.  It does not
recover the individual negative-clone coordinate
`Q(A_s-A_t)`, nor does it produce a pointwise overlay compiler.  It applies
to a selected linear-rate subcode of one cap-`1/2` regular-Hadamard switching
orbit, not to the alternating-form Gram family, lower-cap near-minimizers, or
arbitrary bounded-cap signings.  Off-diagonal optimizers are deliberately
unconstrained and may depend jointly on child and query even though edge
ownership does not.

These limitations are accurately stated.  After the one-line cardinality
repair in BCX.1, the theorem is rigorous and is a substantive positive
bounded-cap exact-sign/disjoint compiler result.
