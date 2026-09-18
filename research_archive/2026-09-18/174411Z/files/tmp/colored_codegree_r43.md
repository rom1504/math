# Wave 43: information--incidence collapse and the redundancy of conflict

## Verdict

The weighted colored-codegree package (10.1141) is sufficient but
unnecessarily strong.  In fact, the first two clauses of (10.1089)—total
selector/certificate information and parent row—already extract one full
parent cut with the required favorable selector mass.  Anchored conflict,
plurality decoding, Gram transport, shared-priority consensus, and the extra
shared-latent code length are not needed for this implication.

The underlying posterior-support extraction is not new: it is precisely the
general theorem (10.862)--(10.863), specialized to the complement incidence
relation.  What was missed in the later strategic interpretation is that
complement incidence satisfies that old theorem's hypothesis *and* already
implies favorability by (10.1045).  Consequently the complement-kernel
implementation is equivalent up to fixed extraction constants to the existing
bare row-good incidence-column target.  This is a correction to the strategic
interpretation of (10.1089), not a convergence proof: constructing such a
low-information, low-row law/column remains open.

The additional exact scalarization below makes the same collapse explicit for
every linear information--row objective.  Random kernels can convexify the
two-cost Pareto frontier, but cannot improve a scalarization.

All logarithms below are natural, consistently with (10.1140).

## 1. Exact conditional-information identity

Let `U` be the reference law on the relevant selector slice (the full
`m`-slice, or first the anchored slice).  Let `D` be the finite set of full
oriented projective parent cuts.  For each cut define its actual complement
incidence column

```math
\mathcal I_d
=\{S:\langle B_S,d\rangle\ge q_n\},
\qquad
\alpha_d=U(\mathcal I_d),
\qquad
h_d=-\log\alpha_d.
\tag{R43.1}
```

Take any joint law `P(S,D)` supported on actual incidences:

```math
P\{S\in\mathcal I_D\}=1.
\tag{R43.2}
```

The chain rule gives exactly

```math
\boxed{
\mathsf H_{\rm tot}
:=D(P_S\Vert U)+I(S;D)
=\mathbb E_DD(P_{S\mid D}\Vert U).}
\tag{R43.3}
```

For any probability `mu` supported on a set `I`,

```math
D(\mu\Vert U)
=D(\mu\Vert U(\cdot\mid I))+\log\frac1{U(I)}
\ge\log\frac1{U(I)}.
```

Applying this conditionally in (R43.3) yields the exact support bound

```math
\boxed{
\mathbb E_D h_D\le\mathsf H_{\rm tot}.}
\tag{R43.4}
```

This is the posterior-support argument already recorded as (10.862).  Its
previously overlooked application here uses neither conflict nor any geometric
property of the labels.

## 2. Simultaneous row extraction

Put

```math
\overline R=\mathbb E R_2(D).
```

When both denominators are positive, (R43.4) gives

```math
\mathbb E\left[
\frac{h_D}{\mathsf H_{\rm tot}}
+\frac{R_2(D)}{\overline R}
\right]\le2.
```

Therefore some full oriented parent cut `d` satisfies simultaneously

```math
\boxed{
h_d\le2\mathsf H_{\rm tot},
\qquad
R_2(d)\le2\overline R,
\qquad
U(\mathcal I_d)\ge e^{-2\mathsf H_{\rm tot}}.}
\tag{R43.5}
```

The zero-denominator cases are easier: zero total information forces the
conditional selector law to equal `U`, while zero mean row forces row zero
almost surely.  Thus (R43.5) remains valid with the evident convention.

In (10.1089), `D` is indeed a full projective parent word `z^D`, with its
response orientation included.  Its row square is the parent quantity
`R_2(z^D)` and is independent of response orientation.  Hence there is no
type mismatch in (R43.5).

If `P_S` is uniform on an anchored family of relative density `beta` inside
the anchored slice, then relative to the full slice

```math
D(P_S\Vert U_m)=\log\beta^{-1}+\log(n/m).
```

The additional term is `O(1)` in the fixed high-ratio window.  Alternatively
one may apply (R43.3) first with the anchored reference law and then multiply
the resulting incidence mass by `m/n`.  Either way anchoring changes only a
fixed factor and does not obstruct extraction.

## 3. Complement incidence already gives favorability

Write the selected cut as `d=(sigma,x)`, with parent energy

```math
E_d=\langle A,d\rangle=q_n-\Delta_d.
```

For every `S in I_d`, the exact complement identity is

```math
\langle B_S,d\rangle=2c_S(d)-E_d\ge q_n,
```

so

```math
c_S(d)\ge q_n-\frac{\Delta_d}{2}.
\tag{R43.6}
```

Principal monotonicity `Q(A[S])<=q_n` and `p_2>=1/2` then give exactly the
already verified estimate (10.1045):

```math
\boxed{
h_d^{\rm fav}(S)
\le Q(A[S])-q_n+(1/2-p_2)\Delta_d-B_{n,m}
\le-B_{n,m}\le0.}
\tag{R43.7}
```

Thus every selector counted by `U(I_d)` is favorable for this **same** full
parent cut.  At project scales,

```math
\mathsf H_{\rm tot}=O(n^{3/4-c}),
\qquad
\overline R=O(n^{9/4-c}),
```

(R43.5)--(R43.7) supply the bare arbitrary-cut tail (10.795), with only fixed
constant changes.  The established inverse-tail recurrence and landing
argument then prove convergence.  No plurality center or agreement estimate
is used.

## 4. Exact scalarization of the old collapse

Specializing (10.862)--(10.863) to actual complement columns also gives the
following exact scalarization, for every `lambda>=0`:

```math
\boxed{
\begin{aligned}
&\inf_{P:\,S\in\mathcal I_D\ {m a.s.}}
\left\{
D(P_S\Vert U)+I(S;D)
+\lambda\mathbb E R_2(D)
\right\}\\
&\hspace{35mm}=
\min_{d:\alpha_d>0}
\left\{-\log U(\mathcal I_d)+\lambda R_2(d)\right\}.
\end{aligned}}
\tag{R43.8}
```

Proof of the lower bound: use (R43.3) and the support inequality conditionally,
then average

```math
\mathsf H_{\rm tot}+\lambda\mathbb ER_2(D)
\ge\mathbb E[h_D+\lambda R_2(D)]
\ge\min_d[h_d+\lambda R_2(d)].
```

Proof of equality: choose one minimizing `d`, take `D=d` deterministically,
and take `S` uniformly from `I_d`.  Then `I(S;D)=0`,
`D(P_S||U)=h_d`, and the row is `R_2(d)`.

Equivalently, the achievable `(total information, mean row)` region has lower
boundary equal to the convex hull of the column points

```math
(h_d,R_2(d)).
```

Indeed, for fixed `P_D`, replacing every `P_(S|d)` by
`U(.|I_d)` makes (R43.4) an equality without changing row.  Finite minimax
therefore reduces the honest two-price problem to

```math
\boxed{
\begin{aligned}
\Theta_{H,R}
&=\min_{\rho\in\Delta(\mathcal D)}
\max\left\{
\frac{\mathbb E_\rho h_D}{H_*},
\frac{\mathbb E_\rho R_2(D)}{R_*}
\right\}\\
&=\max_{0\le t\le1}\min_d
\left\{
t\frac{h_d}{H_*}+(1-t)\frac{R_2(d)}{R_*}
\right\}.
\end{aligned}}
\tag{R43.9}
```

Any feasible convex mixture yields a single column within factor two by
(R43.5).  Conversely a single column supplies the constant-certificate
kernel exactly.  This is the correct reduced separation dual; the conflict
price in (10.1137) is unnecessary for reaching the bare tail.

## 5. Consequences for weighted races and fractional covers

For weighted MinHash, the marginal incidence channel has

```math
\pi_S(d)=\frac{w_d\mathbf1_{d\in\mathcal I_S}}{Z_S}.
```

If its marginal information and row already meet project budgets, apply
(R43.5) directly to this channel.  There is no reason to retain the shared
clock `J`, pay `I(S;J|D)`, prove a colored codegree estimate, or fix one race
realization.  The marginal channel's independent conflict is irrelevant.

Likewise, a full-slice row-good fractional cover of total weight `W` gives
the incidence channel in (10.1090), with `I(S;D)<=log W`; (R43.5) extracts a
row-good column of mass at least `W^{-O(1)}`.  If every cover state is already
row-good, the simpler averaging identity

```math
1\le\sum_dw_dU(\mathcal I_d)
```

directly gives a column of mass at least `1/W`.

The converse direction must be scoped carefully: one high-degree column is
enough for the bare tail and gives a kernel on its conditioned selector
family, but need not cover the entire slice.  Therefore a full-slice
fractional cover remains a stronger construction, whereas the affordable
family formulation of (10.1089) is exactly equivalent, up to constants, to
one row-good column.

## 6. Random-map entropy converse

The full-code-length issue from Wave 42 has a separate exact formulation,
although it is no longer needed for the convergence reduction.  Let `J` be
independent shared randomness and `D=f_J(S)` an incidence-preserving map.
For

```math
\deg(d)=|\{S:d\in\mathcal I_S\}|,
```

conditioning on `(D=d,J=j)` leaves at most `deg(d)` possible selectors.
Therefore

```math
\boxed{
\mathbb E_J I(S;D\mid J)
\ge H(S)-\mathbb E\log\deg(D).}
\tag{R43.10}
```

This holds for randomized maps, not only deterministic extractions.  Combined
with

```math
\mathbb E_J I(S;D\mid J)=I(S;D)+I(S;J\mid D),
```

it quantifies the shared-latent price exactly.

There is also a row-truncated list version.  If
`epsilon=Pr{R_2(D)>R_0}` and

```math
\deg_{R_0}^{\max}
=\max_{d:R_2(d)\le R_0}\deg(d),
```

then

```math
\boxed{
\mathbb E_J I(S;D\mid J)
\ge(1-\epsilon)
\log\frac{|\mathcal G|}{\deg_{R_0}^{\max}}
-h_2(\epsilon).}
\tag{R43.11}
```

This follows by revealing the high-row indicator and applying the list-size
bound on each branch.  It recovers the row-sensitive incidence-degree barrier
behind (10.1047), in nats.

## 7. Exact finite priority audit: useful but now secondary

Exhaustion on actual incidence fibers gives a sharp warning about the shared
priority ansatz itself.

For `A5,m=4`, among all `10^4` incidence maps:

```text
57 have zero conflict;
27 zero-conflict maps are induced by a common total priority;
the best zero-conflict map has row 16 and entropy log 2.
```

One such map uses the same anchored word `+++--` throughout and the two
response orientations on two selectors each.  Every selected complement
energy is exactly `q_5=8`.  Hence the finite point `(row,information,conflict)
=(16,log 2,0)` preserves actual incidence exactly and is priority-realizable.

For the exact conference minimizer `A6,m=5`, among all `12^5=248832`
incidence maps:

```text
22 have zero conflict;
all have row 30 and minimum entropy log 5;
none is induced by any common total priority;
every common-priority map has conflict at least 2/5,
and exactly 40 maps attain 2/5.
```

The displayed zero-conflict selection uses five distinct states, all with
complement energy exactly `q_6=10` and row exactly `30=n(n-1)`.  Its
precedence graph already has the two-cycle `2<4<2`: the first selector must
prefer state `2` to state `4`, while the second must prefer `4` to `2`.
Exhaustion verifies that every zero-conflict incidence map has some precedence
cycle.

Any realization of arbitrary shared continuous state priorities induces a
total order.  Thus the `2/5` lower bound applies not just to exponential
MinHash rates but to **every distribution over common priorities**.  The full
global-map LP can attain zero, so there is no universal conflict-preserving
rounding from the honest incidence-map simplex to shared priorities.

This is a finite exact-minimizer obstruction to the race implementation, not
a scalable asymptotic obstruction and not a falsifier for the bare incidence
column.  In light of (R43.3)--(R43.9), that distinction is now strategic:
priority consensus is unnecessary.

## 8. What remains open

The surviving target is simply to prove, uniformly in the active high-ratio
window, that some exact-minimizer parent cut satisfies

```math
\boxed{
R_2(d)=O(n^{9/4-c}),
\qquad
-\log U_m(\mathcal I_d)=O(n^{3/4-c}).}
\tag{R43.12}
```

By (R43.7), this is the bare favorable-tail target already known to imply
convergence.  Equivalently one may construct an incidence-supported law with
the corresponding total-information and mean-row bounds, but (R43.8) shows
that this is a proof device rather than a genuinely weaker structured target.

No result in this memo proves (R43.12), and no scalable actual-minimizer family
falsifying it is known.  The correction removes an artificial conflict
obstruction; it does not solve the remaining row-versus-incidence-congestion
problem.

## Verification

`tmp/colored_codegree_r43_check.py` independently verifies:

- the chain-rule and shared-latent information identities on a nontrivial
  randomized map;
- conditional-KL support extraction and scalar column minima on
  `A5,A6,A8,A9`;
- the oriented complement-incidence inequality (10.1044) on every audited
  fiber, with `p_2>=1/2` and principal cap at most `q_n`;
- the complete `A5` and `A6` global-map and priority-rationalizability counts;
- the exact `A6` priority lower bound `2/5` and its explicit precedence
  two-cycle.

It ends with

```text
PASS colored_codegree_r43_check
```
