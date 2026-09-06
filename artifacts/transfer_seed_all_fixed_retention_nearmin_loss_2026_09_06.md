# Every fixed random retention loses on some liminf-realizing near-minimizers

Date: 2026-09-06. Seed-transfer track. This upgrades the small-retention
random-selector obstruction to EVERY fixed retention `q<1`, with a
necessary quantifier: the offending parents are actual asymptotically
near-minimizing signings, not asserted to be exact minimizers.

## Theorem

Write `ell=liminf_n M_n/n^(3/2)`. The verified strict upper gives
`ell<1/2`. For EVERY fixed `q` in `(0,1)`, there exist positive constants
`epsilon_q,eta_q` and actual hollow signings `B_j` of orders `d_j->infinity`
such that

```math
\frac{Q(B_j)}{d_j^{3/2}}\longrightarrow\ell,
\qquad
\frac{Q(B_j)-M_{d_j}}{d_j^{3/2}}\longrightarrow0,           (1)
```

but for a uniform principal subset `T_j` of exact size `floor(q d_j)`,

```math
\Pr\left\{
\frac{Q((B_j)_{T_j})}{|T_j|^{3/2}}
\ge\frac{Q(B_j)}{d_j^{3/2}}+\epsilon_q\right\}
\ge\eta_q                                                   (2)
```

for every sufficiently large `j`.

The parents `B_j` can be chosen as principal descendants, at ONE fixed
depth of nested `q`-sampling, of exact minimizers along a liminf order
sequence. Their own orders necessarily also realize the original liminf.

To make the probability constant quantitative, take `rho>0` from the
uniform small-fixed-retention theorem for parent cap bound `C0=1/2`
and target child coefficient `1/2`. Choose any integer `r>=1` with
`q^r<rho`. The proof gives `eta_q=1/(2r)`. It does not give an explicit
positive value of `epsilon_q`; this gap is obtained by compactness.

Thus a uniform random `q`-restriction cannot preserve cap with vanishing
normalized loss on ALL liminf-realizing near-minimizer sequences, even
when `q` is arbitrarily close to one. This does NOT exclude a special
well-chosen seed sequence with a successful one-step restriction rule.

## 1. The finite sampling process and its terminal obstruction

Choose exact minimizers `A_j` at orders `D_j` with
`Q(A_j)/D_j^(3/2)->ell`. Discard finitely many terms so their normalized
caps are at most `1/2`. Put

```math
d_{j,0}=D_j,\qquad d_{j,s+1}=\lfloor q d_{j,s}\rfloor,
\quad 0\le s<r.
```

Starting from the entire vertex set, choose each `T_{j,s+1}` uniformly
inside `T_{j,s}` at size `d_{j,s+1}`. Every marginal `T_{j,s}` is exactly
uniform among all `d_{j,s}` subsets of `[D_j]`: this follows by counting
the equal number of parent subsets containing each possible child, or
by permutation symmetry. In particular the terminal selector is not
a biased or adaptive selector law.

Define the random normalized caps

```math
C_{j,s}=Q((A_j)_{T_{j,s}})/d_{j,s}^{3/2}.
```

For fixed `r`, every size tends to infinity and
`d_{j,s}/D_j->q^s`. The director theorem
`transfer_director_exponential_selector_cost_2026_09_06.md`
therefore gives

```math
\Pr\{C_{j,r}\le1/2\}\longrightarrow0,                   (3)
```

since `d_{j,r}<=rho D_j` for all large `j`.

There are also deterministic lower and upper bounds

```math
C_{j,s}\ge M_{d_{j,s}}/d_{j,s}^{3/2}\ge\ell-o(1),
\qquad
C_{j,s}\le\frac{Q(A_j)}{D_j^{3/2}}
                  (D_j/d_{j,s})^{3/2}.                   (4)
```

The second inequality is exact principal monotonicity. Since `s<=r`
is fixed, all coordinates are bounded by a fixed constant independent
of `j`. This supplies genuine compactness; no cap convergence at all
orders is assumed.

## 2. A first-exit compactness argument gives an actual near-minimal parent

Pass to a subsequence on which the joint law of
`(C_{j,0},...,C_{j,r})` converges weakly to `(C_0,...,C_r)`.
By (3)--(4), almost surely

```math
C_0=\ell,\quad C_s\ge\ell\ (0\le s\le r),\quad C_r\ge1/2>\ell.
```

The disjoint events

```math
E_s=\{C_0=\cdots=C_{s-1}=\ell,\ C_s>\ell\},
\qquad1\le s\le r,
```

partition the probability space up to a null set. Hence some fixed
`s` satisfies `Pr(E_s)>=1/r`. There is an `epsilon>0` such that

```math
\Pr\{E_s,\ C_s>\ell+4\epsilon\}>\frac{3}{4r}.            (5)
```

This follows by increasing the events as `epsilon` decreases to zero;
the limiting event is `E_s`.

For every fixed `delta>0`, the open event
`{C_(s-1)<ell+delta, C_s>ell+4 epsilon}` contains the event in (5).
Portmanteau therefore implies its corresponding prelimit probability
has liminf greater than `3/(4r)`. A diagonal choice gives deterministic
`delta_j->0`, with `delta_j<epsilon` eventually, for which

```math
\Pr\{C_{j,s-1}<\ell+\delta_j,
                 C_{j,s}>\ell+4\epsilon\}>\frac{2}{3r}   (6)
```

for all sufficiently large `j`. The diagonal choice concerns only this
fixed finite process; no iteration depth grows with order.

Condition now on the entire selected parent at step `s-1`. Its next
selector is uniform at exact retention `floor(q d_(j,s-1))`. Averaging
the conditional failure probabilities in (6) shows that at least one
deterministic parent realization `B_j` has

```math
Q(B_j)/d_{j,s-1}^{3/2}<\ell+\delta_j
```

and conditional child probability at least `1/(2r)` of exceeding
`ell+4 epsilon`. Otherwise the left side of (6) would be at most
`1/(2r)`. Select one such parent at every order.

The lower bound in (4) proves its normalized cap tends to `ell`.
Moreover

```math
0\le\frac{Q(B_j)-M_{d_{j,s-1}}}{d_{j,s-1}^{3/2}}
\le\delta_j+\ell-M_{d_{j,s-1}}/d_{j,s-1}^{3/2}\longrightarrow0.
```

On the child event its normalized cap exceeds its parent's by at
least `3 epsilon` eventually. Taking `epsilon_q=epsilon` and
`eta_q=1/(2r)` proves (1)--(2), with room in both constants.

## 3. A finite quantitative version at a fixed near-minimality tolerance

The compactness extraction above obtains a vanishing parent tolerance
but leaves the positive loss size ineffective. There is also a useful
fully finite first-exit statement.

For any deterministic starting parent of normalized cap `b`, run `r`
nested fixed-size restrictions. If the terminal normalized cap is
greater than `b+r epsilon` with probability at least `1-zeta`, then
some step has an actual deterministic parent whose normalized cap is
at most `b+(s-1)epsilon` and whose next uniform restriction increases
that cap by more than `epsilon` with probability at least
`(1-zeta)/r`.

Indeed use the first violation of the thresholds `b+s epsilon`.
These first-exit events are disjoint, their union contains the terminal
event, and there are `r` possibilities. Condition on the parent at the
step with largest first-exit probability. Its conditional child event
is a uniform restriction and its increase is more than `epsilon`.
The parent selected in this finite argument is only within a FIXED
near-minimality tolerance; Section 2 is needed for the stronger (1).

## Exact scope of the conclusion

The universal quantifier over near-minimizer sequences matters.
An exact minimizer restricted once need not be an exact minimizer, so
this argument does not show that every `q` has an offending sequence
of EXACT optimizing parents. It does produce actual signings with
vanishing normalized optimality gap, on liminf-realizing orders.

Nor does failure for some near-minimizer sequences rule out the
sufficient convergence target, which only needs one well-chosen seed
sequence and a suitable realization mechanism. The theorem rules out
uniform random-restriction transfer as a black-box principle for all
near-minimizers. Exceptional or specially chosen selectors remain open.

The proof uses exact sampling composition, elementary compactness and
conditioning, and the already established uniform small-fixed-retention
theorem. There is no numerical step requiring a new finite checker.
