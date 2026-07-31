# Walsh-shell structured landing: definition and quantifier audit

Status: a strongest noncircular family suggested by the second-order shell
theorem is defined below.  Its internal composition implication is rigorous
and its landing statement is exact.  No verified theorem places an arbitrary
exact minimizer, or an asymptotically equivalent signing, in the family.

## 1. Fixed data, independent of the unknown optimum

Fix once and for all:

1. a small leaf order `n_0`;
2. a power `delta>0`;
3. for every balanced pair `m=floor(N/2)`, `n=ceil(N/2)`, an explicit finite
   bridge catalogue `W_(m,n)`.

The first catalogue element is the `m` by `n` northwest rectangle of the
Sylvester matrix of order

```math
q=2^{\lceil\log_2\max(m,n)\rceil}.
```

One may add a fixed finite list of transposes, deterministic row/column
permutations, or other explicitly generated Hadamard rectangles.  The list
is fixed before seeing a signing and does not contain a bridge obtained by
minimizing the parent cap.

For child signings `A,D` and a catalogue bridge `W`, put

```math
T_\delta(A,D)
=\left(\operatorname{cap}(A)^{2/3}
       +\operatorname{cap}(D)^{2/3}\right)^{3/2}
 +(m+n)^{3/2-\delta},                                \tag{WS1}
```

and

```math
g_\delta(A,D,W)
=\operatorname{cap}(A)+\operatorname{cap}(D)
 +L_W-T_\delta(A,D),                                \tag{WS2}
```

where `L_W=max_(u,v)|u^T Wv|`.  Define the child and bridge deficit
histograms `N_eta,N_zeta` exactly as in the second-order amalgamation theorem.
Call the triple `(A,D,W)` **shell-admissible** if `g_delta<=0`, or if

```math
\sum_{e+z<g_\delta}N_\eta(e)N_\zeta(z)<2^{m+n-2}.   \tag{WS3}
```

This is one exact shell-rate inequality.  It can equivalently be replaced by
a finite certificate consisting of high-moment bounds through order
`kappa(m+n)` that implies (WS3), but fixed-order moments do not suffice.
Neither formulation mentions `M_N`.

## 2. The recursively certified family

Let `R_N(delta)` be defined recursively.

- For `N<=n_0`, use a fixed finite list of leaf signings.
- For `N>n_0`, split `N=m+n` as above.  Choose
  `A in R_m(delta)`, `D in R_n(delta)`, and `W in W_(m,n)` for which the
  triple is shell-admissible.  The histogram switching theorem guarantees row
  and column switches `p,q` with amalgamation gain at least `g_delta` when
  `g_delta>0`; when `g_delta<=0`, the independent bound already gives (WS5).
  Put the corresponding block signing in `R_N(delta)`.  A lexicographically
  first successful switching makes this a completely definite construction.

To keep the family nonempty even if no recursive certificate exists at an
order, define

```math
\mathcal F_N^{\rm WS}
=\mathcal R_N(\delta)\cup\{P_N\},                   \tag{WS4}
```

where `P_N` is the canonical square-field Paley restriction already defined
in `constructive_family_phase2_report.md`.  Every ingredient of (WS4) is
specified without the unknown optimum.  Exact membership may be expensive,
because shell enumeration is a counting problem, but it is not logically
circular.

The family is deliberately stronger than a fixed-moment family.  Banica's
fixed-moment Gaussian glow theorem for a Hadamard bridge controls only the
central scale.  Condition (WS3) asks for the extreme glow and child shells at
exponential accuracy, and includes exactly that missing information.

## 3. Verified implication inside the family

For every shell-admissible node, the second-order theorem and (WS1)--(WS3)
give

```math
\operatorname{cap}(S)
\le T_\delta(A,D).                                  \tag{WS5}
```

Since every signing has cap `Theta(N^(3/2))` on the relevant low-cap scale,
the mean-value theorem converts (WS5) to

```math
\operatorname{cap}(S)^{2/3}
\le \operatorname{cap}(A)^{2/3}
   +\operatorname{cap}(D)^{2/3}
   +O(N^{1-\delta}).                                \tag{WS6}
```

The error in (WS6) is geometrically summable.  This is the rigorous reason
for using the recurrence-calibrated target (WS1).  Requiring a raw residual
`o(N^(3/2))` at every recursive node would incorrectly force an impossible
zero limiting constant; a valid shell target must retain the leading amount

```math
\left(a^{2/3}+d^{2/3}\right)^{3/2}-a-d
```

allowed by `2/3`-power composition.

## 4. Exact landing statement

Define

```math
u_N=\min_{S\in\mathcal F_N^{\rm WS}}
       \operatorname{cap}(S)^{2/3},
\qquad
b_N=M_N^{2/3}.                                      \tag{WS7}
```

Because `F_N^WS` is a subclass of all signings,

```math
0\le u_N-b_N.                                       \tag{WS8}
```

The required landing theorem is exactly

```math
\boxed{0\le u_N-b_N=o(N).}                         \tag{WS9}
```

On the known `Theta(N^(3/2))` cap scale, (WS9) is equivalent to

```math
\min_{S\in\mathcal F_N^{\rm WS}}\operatorname{cap}(S)
=M_N+o(N^{3/2}).                                    \tag{WS10}
```

Equations (WS6) and (WS9) would respectively supply the composition and
landing halves of the steering package.  Only the first is presently
available, and only conditionally at nodes satisfying (WS3).

## 5. Audit from arbitrary exact minimizers

Let `A_N` be an arbitrary exact minimizer.  The following verified operations
do **not** place a member within `o(N^(3/2))` cap loss in (WS4).

### 5.1 Switching, permutation, and global negation

These operations preserve the cap and the complete energy histogram.  They
cannot change a generic balanced cross block into the switching orbit of a
fixed Walsh rectangle.  The number of bridges in that orbit is only
exponential in `N`, while the number of rectangular sign matrices is
exponential in `N^2`; no rigidity theorem for minimizers closes this gap.

### 5.2 Principal restriction

A restriction of `A_N` has cap at most `M_N`, but this is an order-`m`
statement with the wrong normalization.  There is no verified balanced
partition for which the two restriction caps satisfy

```math
\left(\operatorname{cap}(A_N[U])^{2/3}
     +\operatorname{cap}(A_N[V])^{2/3}\right)^{3/2}
\le M_N+o(N^{3/2}).                                 \tag{WS11}
```

Nor is there a verified inverse completion that restores order `N` with a
Walsh bridge and preserves cap to lower order.  This is the cross-order
restriction obligation, not a consequence of monotonicity.

### 5.3 Universal tail bounds

Exact optimality gives the endpoint `|H_(A_N)(x)|<=M_N`; it gives no bound on
how many spins lie near that endpoint.  Bonami and general Hanson--Wright
bounds have only `exp(-Theta(sqrt(N)))` project-scale strength, while (WS3)
requires rate exceeding the full `Theta(N)` spin entropy.  Therefore no
verified implication puts the child shells of an exact minimizer in the
shell-admissible class.

### 5.4 Existing conference and completion maps

Conference completion, shallow principal restriction, and the explicit
Paley fallback produce cap constant `1/2+o(1)`.  Their landing statement is
equivalent to the unproved assertion
`M_N=(1/2+o(1))N^(3/2)`.  They do not transform an arbitrary minimizer and do
not prove (WS9).

Consequently there is currently **no verified implication** from an
arbitrary exact minimizer to a member of `F_N^WS` with lower-order cap loss.

## 6. Circularity and quantifier audit

The valid quantifier order is

```math
\text{define }\mathcal F_N^{\rm WS}
\quad\longrightarrow\quad
\text{prove }\exists S_N\in\mathcal F_N^{\rm WS}
\text{ with }
\operatorname{cap}(S_N)\le M_N+o(N^{3/2}).          \tag{WS12}
```

The following apparent repairs reverse or hide these quantifiers.

1. Defining the family as constructions from "near-optimal children" uses
   `M_m,M_n` in membership and assumes landing at smaller orders.
2. Centering (WS3) at `b_N` or choosing its threshold from `M_N` puts the
   unknown target directly into the definition.
3. Choosing `W` by minimizing the full parent cap is exactly full bridge
   optimization, not a fixed structured catalogue.
4. Taking one canonical representative of every exact energy histogram gives
   tautological landing, because the endpoint of the histogram is already
   the cap.  It supplies neither a smaller optimization problem nor a
   composition rule.
5. Saying that exact minimizers have thin extreme shells because their cap is
   small confuses an endpoint bound with a shell-count bound.  The latter is
   the missing theorem.

The shell histogram is genuinely smaller than the full response alignment
for composition, but it is not automatically a simpler state for landing:
its support endpoint encodes the objective exactly, and its extreme masses
are not controlled by exact optimality.

## 7. Precise falsifiable missing theorem

A noncircular theorem strong enough to establish landing is the following.

> **Minimizer-to-Walsh-shell transfer.** There exist fixed `delta>0`, a
> finite explicit bridge catalogue as in Section 1, and `epsilon_N->0` such
> that every signing `A` with
> `cap(A)<=N^(3/2)` admits a recursively balanced Walsh-shell reconstruction
> `T_N(A) in F_N^WS` satisfying
>
> ```math
> \operatorname{cap}(T_N(A))
> \le\operatorname{cap}(A)+\epsilon_NN^{3/2}.        \tag{WS13}
> ```
>
> At every internal node the reconstruction must exhibit either the exact
> convolution certificate (WS3) at the recurrence target, or a stated
> high-moment certificate through order `Theta(|Q|)` that implies it.

The hypothesis uses only the cap of the supplied signing and a fixed explicit
upper scale; it does not refer to `M_N`.  Applying (WS13) to an exact minimizer
proves (WS9).

This theorem is falsifiable in two independent ways:

- find a scalable low-cap family for which every bridge in the fixed
  catalogue fails (WS3) by an exponential shell-rate margin at some balanced
  node; or
- prove that every recursively certified output has cap at least
  `M_N+cN^(3/2)` along an infinite sequence.

No current result proves either outcome.  The first missing implication is
already visible at the root: exact optimality supplies neither the balanced
restriction budget (WS11) nor the exponential shell inequality (WS3).
Without one of those genuinely new ingredients, (WS13) is a named landing
hypothesis rather than progress toward (WS9).
