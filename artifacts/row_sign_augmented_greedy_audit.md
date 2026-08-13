# Row-sign recoupling by augmented one-flip ascent

## Status

This note proves a polynomial-time, samplewise same-spin certificate for the
concrete asymmetric law

~~~math
X\sim\operatorname{Unif}\{\pm1\}^n,
\qquad Y=\operatorname{sign}(AX),
~~~

and audits it only after fixing the theorem.  The certificate keeps the
anchored cross field from the exact recoupling theorem and applies a
deterministic one-flip ascent to the resulting weighted shore.

The finite outcome is substantially stronger than the spectral surrogate.
On the stored exact cases through order 14 the greedy certificate often
equals or closely approaches the exact augmented-shore certificate.  Its
sampled normalized defect on Paley conference matrices is `0.00128` at order
30 and falls below `0.000010` at orders 90 and 98.  Fixed held-out random
signings have zero defect in every sampled realization.  This is a
reproducible computational scaling law that isolates a precise open uniform
lemma, not an asymptotic proof.  In contrast, the weighted projector and
uncorrected `sign(h)` witnesses have positive normalized defects and miss the
required coefficient on the tested conference orders; they are not retained
as leading routes.

All energies below use doubled normalization

~~~math
Q(A)=\max_{z\in\{\pm1\}^n}|z^{\mathsf T}Az|.
~~~

At a zero row field the audit uses the switching-equivariant convention
`Y_i=X_i`.  The response identity below is independent of this tie rule.

## 1. The response is exact and signing-independent

Put

~~~math
V_A(X)=X^{\mathsf T}AY=\sum_i|(AX)_i|.
~~~

For every fixed row, `(AX)_i` is a sum of `n-1` independent Rademacher
variables.  Therefore every signing satisfies the exact identity

~~~math
\boxed{
\mathbb E_XV_A(X)
=n\,\mathbb E|S_{n-1}|
=n(n-1)
 {\binom{n-2}{\lfloor(n-2)/2\rfloor}\over2^{n-2}},}       \tag{1}
~~~

where `S_m` is a length-`m` Rademacher sum.  In particular,

~~~math
\mathbb E_XV_A(X)
=\left(\sqrt{2/\pi}+o(1)\right)n^{3/2}.                 \tag{2}
~~~

Thus improving the current doubled constant
`c_*=0.672986728863...` reduces to retaining all but

~~~math
\sqrt{2/\pi}-c_*=0.1248978\ldots                         \tag{3}
~~~

of normalized response during same-spin recoupling.

## 2. Exact anchored-shore reduction and clipped defect

For one realization, let

~~~math
I=\{i:X_i=Y_i\},\qquad J=\{i:X_i=-Y_i\},
\qquad p=X_I,\quad q=X_J,
~~~

and

~~~math
P=p^{\mathsf T}A[I]p,\qquad R=q^{\mathsf T}A[J]q.
~~~

Then `V_A(X)=P-R`.  If `PR>=0`, ordinary recoupling is lossless.  If
`PR<0`, positivity of `V_A(X)` forces `P>0>R`.  Define

~~~math
h_J=A[J,I]p,
\qquad
E_J(p)=
\begin{pmatrix}A[J]&h_J\\h_J^{\mathsf T}&0\end{pmatrix},              \tag{4}
~~~

and define `E_I(q)` symmetrically.  The exact anchored-shore theorem gives

~~~math
Q(A)\ge
\max\left\{
 |P|+C_{\operatorname{sgn}P}(E_J(p)),
 |R|+C_{\operatorname{sgn}R}(E_I(q))
\right\}.                                                   \tag{5}
~~~

Here `C_sigma(E)=max_w sigma w^T E w`.  Each branch optimizes only a free
shore plus one collapsed anchor coordinate, but exact evaluation still costs
roughly `2^(n/2)` on balanced shores.

For any computable lower witnesses `L_J<=C_sgn(P)(E_J)` and
`L_I<=C_sgn(R)(E_I)`, put

~~~math
K(X)=\max\{|P|+L_J,|R|+L_I\},
\qquad
\Delta(X)=[V_A(X)-K(X)]_+.                              \tag{6}
~~~

The positive part is essential: a recoupled spin can exceed the asymmetric
response.  Samplewise and after averaging,

~~~math
\boxed{Q(A)\ge V_A(X)-\Delta(X),\qquad
Q(A)\ge\mathbb EV_A-\mathbb E\Delta.}                  \tag{7}
~~~

Every “certificate” in the tables is `E V - E Delta`, never an unclipped
average of `K`.

## 3. Deterministic augmented coordinate ascent

Consider either weighted augmented shore `E` in (4), let `sigma` be the
anchor-energy sign, and write its last column as `(h,0)`.  Define

~~~math
w^{(0)}=(\operatorname{sign}(\sigma h),1),              \tag{8}
~~~

At a zero entry of `h`, use the corresponding original free-shore spin as
the sign in (8).  This deterministic tie rule makes the initialization
covariant under switching and global negation.  Starting at `w^(0)`,
repeatedly flip the coordinate with largest strictly positive gain in
`Phi(w)=sigma w^T E w`.  Gain ties are resolved by the smallest coordinate
index.  The exact flip gain is

~~~math
\boxed{
\Phi(w^{\oplus u})-\Phi(w)
=-4\sigma w_u(Ew)_u.}                                  \tag{9}
~~~

Stop when every gain is nonpositive and return

~~~math
G_\sigma(E)=\max\{0,\Phi(w^{\rm term})\}.               \tag{10}
~~~

This is a rigorous lower bound on `C_sigma(E)`, independent of whether the
local optimum is global.  Consequently (6)--(7) with

~~~math
L_J=G_{\operatorname{sgn}P}(E_J(p)),\qquad
L_I=G_{\operatorname{sgn}R}(E_I(q))                    \tag{11}
~~~

is a polynomial-time same-spin certificate.

The separately reported one-shot witness is related but not identical to the
initial aligned value.  With `r=sign(h_J)` (using the same free-spin tie at
zero), the two orientations of the
collapsed anchor give the exact elementary certificate

~~~math
Q(A)\ge
|P+r^{\mathsf T}A[J]r|+2\|h_J\|_1.                    \tag{11a}
~~~

The aligned initialization instead uses `r=sign(sigma h_J)` and lower-bounds
the sign-specific augmented cap before ascent.  The checker evaluates these
as distinct columns.  Whenever the exact augmented cap is enumerated, the
checker asserts samplewise that it dominates the initial-aligned,
weighted-projector, and greedy witnesses.  The absolute one-shot certificate
need not be dominated by that sign-aligned branch; where `Q(A)` is known, it
is instead checked directly against `Q(A)` samplewise.

### Complexity

Let

~~~math
L(E)=\sum_{u<v}|e_{uv}|.
~~~

The objective lies in `[-2L(E),2L(E)]`.  Because `E` is integral, every
strict gain in (9) is at least four, so ascent takes at most `L(E)` flips.
For a collapsed shore with `i+j=n`,

~~~math
L(E_J)
\le {j\choose2}+\sum_{v\in J}|(h_J)_v|
\le {j\choose2}+ij=O(n^2).                              \tag{12}
~~~

Maintaining the field vector `Ew` makes one flip cost `O(n)`.  Initialization
costs `O(n^2)`, so each anchored run costs `O(n^3)` in the worst case and
uses `O(n^2)` input/state space.  The implementation uses exactly this
incremental update; there is no hidden Boolean enumeration in the greedy
column.

## 4. Exact terminal condition and the missing uniform lemma

At termination, (9) says

~~~math
\sigma w_u(Ew)_u\ge0\quad\text{for every }u.
~~~

Therefore its energy has the exact local-field form

~~~math
\boxed{
G_\sigma(E)=\sum_u|(Ew)_u|.}                            \tag{13}
~~~

For the branch anchored on `I`, writing `w=(r,t)`, this is

~~~math
T_J(r,t)=
\sum_{v\in J}|(A[J]r+t h_J)_v|+|h_J^{\mathsf T}r|.     \tag{14}
~~~

The terminal spin also obeys the original-block stability conditions

~~~math
\sigma r_v(A[J]r+t h_J)_v\ge0\ (v\in J),
\qquad
\sigma t h_J^{\mathsf T}r\ge0.                         \tag{15}
~~~

In the hard branch `P>0>R`, the greedy defect is exactly bounded by

~~~math
\boxed{
\Delta_{\rm greedy}(X)
\le\min\left\{
 [|R|-T_J(r,t)]_+,
 [|P|-T_I(s,u)]_+
\right\}.}                                             \tag{16}
~~~

Thus the concrete sufficient terminal condition for an `O(n)` samplewise
defect is

~~~math
T_J(r,t)\ge |R|-O(n)
\quad\text{or}\quad
T_I(s,u)\ge |P|-O(n).                                  \tag{17}
~~~

The weaker condition actually needed for a constant improvement is that the
expectation of the minimum in (16) have coefficient strictly below (3).
The strongest clean next lemma, suggested by the data, is

> **Row-sign local-repair lemma.** Uniformly over project-scale signings,
> deterministic ascent (8)--(10) on the two anchored shores satisfies
> `E_X Delta_greedy(A,X)=o(n^(3/2))`.

This is genuinely simpler than full parent minimization: for each `X` it runs
two deterministic polynomial local searches and checks only the terminal
fields (14)--(15).

### Conference proof attempt

Switching by `X` makes `p=q=1` on their shores.  If `C^2=(n-1)I`, put
`D=C[J]`, `h=C[J,I]1`, and `ell=D1+h`.  Row-sign selection gives
`ell_v<0` on `J`, while the conference block equations give

~~~math
\|h\|_2^2=(n-1)|I|-\|C[I]1\|_2^2.                     \tag{18}
~~~

These identities do not yet prove (17).  They control the total squared
cross field but not its signs relative to the terminal vector.  The trivial
integrality/parity floor for the terminal local-field sum is only `O(n)` in
the worst case, whereas `|R|` is on the `n^(3/2)` scale.  Monotonicity from the initialization
also does not suffice: the uncorrected aligned `sign(h)` defects approach a
positive normalized constant in the finite data.  The numerical gain comes
from about `O(sqrt(n))` correcting flips, each harvesting correlated local
field, and a proof must use the row-selected sign pattern—not conference
singular values alone—to lower-bound the accumulated gains or the terminal
sum (14).  No conference proof is claimed here.

## 5. Exact finite audit

The checker exhausts every projective `X` for all stored minimizer orbit
representatives through order 8, the exact order-10 matrix, and witnesses of
orders 11--14.  It also exhausts the conference cases through order 14.
For every such case it independently recomputes `Q`, the exact shore caps,
the exact augmented caps, and all polynomial witnesses.  Agreement-shore
spectra are cached by subset; augmented data are cached by agreement subset
and the restricted anchor spin.

Every number below is a resulting clipped certificate divided by `n^(3/2)`.
`old-exact`, `old-proj`, and `nuclear` discard the cross field; `aug-exact`
keeps it but enumerates the free shore; `aug-proj` is the weighted spectral
surrogate; and `greedy` is (8)--(11).

| case | response | old-exact | old-proj | nuclear | aug-exact | aug-proj | greedy | one-shot `sign(h)` |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| min `n=3` | 0.577350 | 0.577350 | 0.577350 | 0.577350 | 0.577350 | 0.577350 | 0.577350 | 0.577350 |
| min `n=4` | 0.750000 | 0.750000 | 0.696213 | 0.687500 | 0.750000 | 0.701717 | 0.750000 | 0.625000 |
| min `n=5` | 0.670820 | 0.670820 | 0.574604 | 0.559017 | 0.670820 | 0.616651 | 0.670820 | 0.670820 |
| min `n=6` | 0.765466 | 0.595362 | 0.534025 | 0.518097 | 0.680414 | 0.579939 | 0.680414 | 0.510310 |
| min `n=7`, orbit 0 | 0.708683 | 0.678311 | 0.582934 | 0.568663 | 0.708683 | 0.621146 | 0.708683 | 0.698559 |
| min `n=7`, orbit 1 | 0.708683 | 0.708683 | 0.619281 | 0.604544 | 0.708683 | 0.667921 | 0.708683 | 0.708683 |
| min `n=7`, orbit 2 | 0.708683 | 0.688435 | 0.590023 | 0.575716 | 0.708683 | 0.629191 | 0.708683 | 0.708683 |
| min `n=8`, orbit 0 | 0.773398 | 0.745777 | 0.570467 | 0.552517 | 0.759587 | 0.603063 | 0.759587 | 0.563476 |
| min `n=8`, orbit 1 | 0.773398 | 0.767874 | 0.618430 | 0.599630 | 0.767874 | 0.648548 | 0.767874 | 0.707107 |
| exact `n=10` | 0.778217 | 0.736712 | 0.601906 | 0.587489 | 0.751782 | 0.644728 | 0.751782 | 0.692984 |
| witness `n=11` | 0.742001 | 0.725405 | 0.571990 | 0.556086 | 0.740395 | 0.612344 | 0.731293 | 0.707738 |
| witness `n=12` | 0.781453 | 0.760779 | 0.575273 | 0.558196 | 0.774405 | 0.608704 | 0.770388 | 0.660795 |
| witness `n=13` | 0.750795 | 0.744837 | 0.553590 | 0.537261 | 0.750795 | 0.589095 | 0.747045 | 0.689042 |
| witness `n=14` | 0.783775 | 0.760024 | 0.554909 | 0.539973 | 0.771900 | 0.588967 | 0.771667 | 0.676897 |

The greedy and exact augmented columns coincide on all minimizer-orbit cases
through order 10.  They separate on the later witnesses, proving that the
greedy algorithm is not merely hiding exact shore optimization.  The greedy
certificate is above `c_*` in every listed case from order 6 onward; orders 3
and 5 are finite small-order exceptions.

## 6. Conference scaling and falsifiers

Order 18 uses an exact exhaustive recomputation of `Q=66`, the analytic
exact response (1), and `32768` sampled `X` values for shore defects.  Orders
26 and above use reproducible Monte Carlo for defects; their spectral and
greedy values are rigorous samplewise witnesses, but their averages are
statistical estimates.  Intermediate Paley matrices are generated exactly
and checked against `C^2=(n-1)I`.

| `n` | defect mode | response | aug-exact defect | aug-proj defect | greedy defect | greedy certificate | one-shot defect | one-shot certificate |
|---:|---|---:|---:|---:|---:|---:|---:|---:|
| 6 | exhaustive | 0.765466 | 0.085052 | 0.185527 | 0.085052 | 0.680414 | 0.255155 | 0.510310 |
| 10 | exhaustive | 0.778217 | 0 | 0.206590 | 0.002965 | 0.775252 | 0.377991 | 0.400226 |
| 14 | exhaustive | 0.783775 | 0.011875 | 0.194808 | 0.012015 | 0.771760 | 0.106878 | 0.676897 |
| 18 | sample `32768` | 0.786885 | 0.000170 | 0.197557 | 0.000380 | 0.786505 | 0.135632 | 0.651253 |
| 26 | sample `16384` | 0.790251 | -- | 0.206169 | 0.005894 | 0.784357 | 0.118843 | 0.671409 |
| 30 | sample `4096` | 0.791264 | -- | 0.213070 | 0.001278 | 0.789987 | 0.119835 | 0.671429 |
| 38 | sample `4096` | 0.792653 | -- | 0.221547 | 0.001240 | 0.791413 | 0.121274 | 0.671379 |
| 42 | sample `4096` | 0.793150 | -- | 0.223552 | 0.000423 | 0.792726 | 0.122910 | 0.670240 |
| 54 | sample `4096` | 0.794199 | -- | 0.230997 | 0.000124 | 0.794075 | 0.127286 | 0.666913 |
| 62 | sample `4096` | 0.794674 | -- | 0.235809 | 0.000133 | 0.794541 | 0.130804 | 0.663870 |
| 74 | sample `4096` | 0.795194 | -- | 0.240600 | 0.000034 | 0.795160 | 0.134049 | 0.661144 |
| 90 | sample `4096` | 0.795671 | -- | 0.245592 | 0.000009 | 0.795662 | 0.137743 | 0.657928 |
| 98 | sample `2048` | 0.795852 | -- | 0.246547 | 0 | 0.795852 | 0.137277 | 0.658574 |

The final JSON records standard errors, positive-defect frequencies,
quantiles, maxima, cache sizes, and flip counts.  The weighted projector is
far below `c_*`.  The one-shot witness is already below `c_*` at order 26
and is more than seven estimated standard errors below it at order 98.  Those
two implementations are therefore stopped as leading mechanisms.  The
greedy repair is the surviving lead.

For the conference scaling rows, the exact Monte Carlo details for the
greedy defect are:

| `n` | seed | samples | normalized defect | standard error | positive fraction | largest raw defect | mean cached flips |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 18 | 26081318 | 32768 | 0.00037963 | 0.00003776 | 0.00348 | 14 | 2.312 |
| 26 | 26081326 | 16384 | 0.00589383 | 0.00015235 | 0.11987 | 24 | 2.771 |
| 30 | 26081330 | 4096 | 0.00127778 | 0.00012293 | 0.03882 | 20 | 3.454 |
| 38 | 26081338 | 4096 | 0.00124026 | 0.00010133 | 0.05029 | 18 | 4.278 |
| 42 | 26081342 | 4096 | 0.00042336 | 0.00005220 | 0.02295 | 20 | 4.842 |
| 54 | 26081354 | 4096 | 0.00012428 | 0.00002530 | 0.01001 | 18 | 6.125 |
| 62 | 26081362 | 4096 | 0.00013303 | 0.00002472 | 0.00952 | 18 | 7.130 |
| 74 | 26081374 | 4096 | 0.00003375 | 0.00001273 | 0.00244 | 20 | 8.575 |
| 90 | 26081390 | 4096 | 0.00000915 | 0.00000657 | 0.00098 | 22 | 10.540 |
| 98 | 26081398 | 2048 | 0 | 0 | 0 | 0 | 11.420 |

The exact augmented-state audit also exposes its exponential cost.  The
exact order-10 minimizer used `799` cached anchor states and `50912` Boolean
candidates across them.  The exact order-14 conference used `11622` cached
anchor states and `4172896` Boolean candidates.  The sampled order-18 run
already used `52534` anchor states and `143544352` exact Boolean candidates.
The greedy calculation uses the same `(P,h,E)` inputs but no candidate
enumeration.  The largest observed cached flip count was 26 for the order-98
conference case and 32 across all audited order-98 cases.

Two small post-selected random stress matrices were also audited exactly:
the worst seeds in the declared sweeps `9300:9308` at order 12 and
`9400:9412` at order 14.  Their normalized greedy defects were respectively
about `0.000963` and `0.000107`.  Fixed, non-post-selected random signings at
orders 30, 62, and 98 are included as held-out Monte Carlo cases in the
result file.  Their matrix seeds are respectively `271830`, `271862`, and
`271898`; their `X`-sample seeds are those values plus one million.  Each had
zero greedy defect in every sampled realization (`4096`, `4096`, and `2048`
samples respectively).

## Reproduction

~~~text
.venv/bin/python computations/audit_row_sign_recoupling_law.py \
  --output computations/results/row_sign_recoupling_law_audit.json
~~~

The result file records source matrices and hashes, exact versus sampled
classification, seeds and sample counts, exact rational responses, exact cap
checks where feasible, every clipped defect/certificate, state-cache sizes,
greedy flip counts, and a canonical payload hash.  Numerical eigensolver
outputs and Monte Carlo averages are diagnostics; the response identity,
recoupling theorem, and greedy certificate are proved independently of them.

The final canonical payload SHA-256 is
`591acb823f980cd88ad1689081b1032e9071a02f757b8f23d0f1ad0a7f570cb7`.
