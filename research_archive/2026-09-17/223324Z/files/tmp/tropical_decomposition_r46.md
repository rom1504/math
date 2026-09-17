# Wave 46C: matched tropical normalization and tight principal decomposition

## Status

- **Verified:** the actual matched selector/completion endpoint has tropical
  score

  ```math
  E_d-q_*-a_d,
  \qquad
  a_d=\min_{|S|=m} k_S(d[S]),
  ```

  where `k_S(y)` is the maximum external completion energy.  The proposed
  expression retaining `q_S` inside the maximum is not the normalization of
  the actual endpoint: the selector factor `q_S-q_*` cancels exactly against
  component normalization.
- **Verified:** its interpolation envelope is affine exactly when a parent
  ground is a best external lift of a ground child on a selector attaining
  `q_*`.  External optimality is automatic once the parent and child are both
  ground states.
- **Exact finite evidence:** the existential tight-decomposition property
  holds at every selector size for every switching-normalized exact minimizer
  through order eight (all `7473` representatives), and for the stored
  order-nine minimizer.  The stronger claim that every maximal child ground
  lifts is false already at small order and on `A_8,m=4,5`.
- **Open minimizer conjecture:** no deduction of the existential property from
  exact signing minimality was found.  Section 5 gives the exact block-game
  inequality a proof must establish and the exact falsifier.
- **Verified exclusion link:** a crossing score is exactly a log ratio of
  posterior selector-exclusion probabilities.  This bounds the incidence at
  the charged (low-base) endpoints, but it does not bound reverse fan-in at a
  high-base active state and therefore does not yet control the overlap `D`.

The checker and output are
`tmp/tropical_decomposition_r46_check.py` and
`tmp/tropical_decomposition_r46.out`.

## 1. Oriented completion notation

Let an oriented cut be `d=(sigma,x)` and put

```math
E_d=\sigma x^{\mathsf T}Ax,
\qquad
c_S(d[S])=\sigma x_S^{\mathsf T}A[S]x_S.
```

For an oriented child state `y` define its optimal external energy

```math
k_S(y)=\max_{d:d[S]=y}\{E_d-c_S(y)\}.
\tag{R46C.1}
```

Thus the actual outside partition in the matched construction is

```math
K_{\beta,S}(y)
=\sum_{d:d[S]=y}
  \exp\{\beta(E_d-c_S(y))\},
\qquad
\beta^{-1}\log K_{\beta,S}(y)\longrightarrow k_S(y).
\tag{R46C.2}
```

Write

```math
q_n=Q(A),\qquad q_S=Q(A[S]),\qquad
q_*=\max_{|S|=m}q_S,
\qquad
a_d=\min_{|S|=m}k_S(d[S]).
\tag{R46C.3}
```

All energies here use the full symmetric-matrix convention, so flipping one
internal edge changes a child energy by `4`, not `2`.

## 2. Exact cancellation and the corrected endpoint exponent

Retain the canonically normalized construction (10.938):

```math
b_S(d)=e^{-F_{\beta,S}(d[S])},\qquad
F_{\beta,S}(y)=\log\frac{K_{\beta,S}(y)}{2^{n-m}},
```

```math
z_S=\mathbb E_{\nu_\beta}b_S
=\frac{2^{n-m}Z_S(\beta)}{Z_A(\beta)},
\qquad
q(S)=\frac{z_S}{\sum_Tz_T},
\qquad h_S=\frac{b_S}{z_S}.
```

The endpoint likelihood therefore satisfies the finite-temperature identity

```math
\boxed{
f_\beta(d)=\sum_Sq(S)h_S(d)
=\frac{\sum_S b_S(d)}{\sum_Tz_T}.}
\tag{R46C.4}
```

In particular, the factors `z_S` cancel before any asymptotics are taken.
Now

```math
\beta^{-1}\log z_S\longrightarrow q_S-q_n,
\qquad
\beta^{-1}\log\sum_Tz_T\longrightarrow q_*-q_n,
\qquad
\beta^{-1}\log b_S(d)\longrightarrow-k_S(d[S]).
```

Consequently

```math
\boxed{
\beta^{-1}\log f_\beta(d)
\longrightarrow q_n-q_*-a_d.}
\tag{R46C.5}
```

Since `nu_beta(d)` has normalized exponent `E_d-q_n`, the endpoint law has
normalized exponent

```math
\boxed{E_d-a_d-q_*.}
\tag{R46C.6}
```

Thus `E_d+max_S{q_S-k_S(d[S])}-q_*` double-counts the selector exponent in
this construction.  It would correspond to weighting the unnormalized
components by the low-temperature selector weight without dividing each
component by `z_S`.  The checker verifies (R46C.4) directly at finite
temperature with maximum error `3.56e-15`.

For `mu_s proportional nu_beta f_beta^s`, all state-independent terms drop
and the tropical interpolation lines are exactly

```math
L_d(s)=E_d-sa_d,
\qquad
\Lambda(s)=\max_d L_d(s).
\tag{R46C.7}
```

## 3. Endpoint maximum and tight-decomposition equivalence

For every `d,S`, the current completion is admissible in (R46C.1), so

```math
E_d-k_S(d[S])\le c_S(d[S])\le q_S\le q_*.
\tag{R46C.8}
```

It follows that `E_d-a_d<=q_*`.  Conversely, take a selector with
`q_S=q_*`, a positive oriented child ground `y`, and an external optimizer
`d` in (R46C.1).  Then

```math
E_d-k_S(y)=c_S(y)=q_*,
```

and hence

```math
\boxed{\max_d(E_d-a_d)=q_*.}
\tag{R46C.9}
```

Equality for a particular `d` forces equality throughout (R46C.8).  Thus
`d` is endpoint-active if and only if there is an `S` such that

```math
q_S=q_*,\qquad c_S(d[S])=q_*,\qquad
E_d=c_S(d[S])+k_S(d[S]).
\tag{R46C.10}
```

The last equality says precisely that `d` is an optimal external completion.
It follows that a common base/endpoint maximizer exists if and only if some
`d,S` satisfy

```math
\boxed{
E_d=q_n,\quad q_S=c_S(d[S])=q_*,\quad
q_n=q_*+k_S(d[S]).}
\tag{R46C.11}
```

The external-optimality clause in (R46C.11) is actually automatic from the
first two clauses: every completion has total energy at most `q_n`, while
the displayed parent ground supplies external energy `q_n-q_*`.

For a finite collection of lines, a common endpoint maximizer is equivalent
to affine interpolation:

```math
\boxed{
\Lambda(s)=(1-s)q_n+s q_*\quad(0\le s\le1)
\quad\Longleftrightarrow\quad\text{(R46C.11).}}
\tag{R46C.12}
```

Indeed every line is below the chord joining the two endpoint maxima.  A
common maximizer attains the chord.  Conversely, equality at one interior
`s` forces any attaining line to attain both endpoint maxima separately.

## 4. Parent-ground face versus maximal-child cylinders

Let

```math
\mathcal G_n(A)=\{d:E_d=q_n\},
```

and let the maximal-child cylinder union be

```math
\mathcal C_m^*(A)=
\bigcup_{\substack{|S|=m\\q_S=q_*}}
\{d:c_S(d[S])=q_*\}.
\tag{R46C.13}
```

Then the exact affine/tight property is simply

```math
\boxed{\mathcal G_n(A)\cap\mathcal C_m^*(A)\ne\varnothing.}
\tag{R46C.14}
```

This formulation contains no temperature, posterior, or harmonic notation.
It is the appropriate **tight principal decomposition conjecture**:

> For every exact order-`n` minimizer and every `1<=m<n`, (R46C.14) holds.

Its exact falsifier is an exact minimizer with

```math
\boxed{
\delta_m(A):=q_*-\max_{\substack{d\in\mathcal G_n(A)\\
|S|=m,\ q_S=q_*}}c_S(d[S])>0.}
\tag{R46C.15}
```

Such a matrix has no common base/endpoint state, so its tropical envelope is
strictly below the endpoint chord at every interior time.  The conjecture is
existential.  The stronger assertion that every maximal child ground has a
parent-ground lift is already false.

## 5. What exact signing minimality gives, and the missing separation

Fix one selector `S` and freeze all edges outside its internal block.  The
external profile `k_S(y)` is then fixed.  Replacing `A[S]` by any complete
order-`m` signing `B` gives the exact norm

```math
\boxed{
Q(A^{S\leftarrow B})
=\max_y\{k_S(y)+c_B(y)\}.}
\tag{R46C.16}
```

Therefore global exact minimality implies the finite block-game inequalities

```math
\max_y\{k_S(y)+c_B(y)\}\ge q_n
\qquad\text{for every signing }B.
\tag{R46C.17}
```

This is the cleanest exact route from minimality, but it does not by itself
prove (R46C.14).  Let

```math
s_S(y)=q_n-[k_S(y)+c_{A[S]}(y)]\ge0.
```

If `B` is obtained by flipping an internal edge set `F`, then

```math
c_B(y)=c_{A[S]}(y)-4\sum_{ij\in F}a_{ij}y_{ij}.
```

To contradict minimality one must find one `F` and one `epsilon>0` satisfying

```math
\boxed{
4\sum_{ij\in F}a_{ij}y_{ij}
\ge \epsilon-s_S(y)
\qquad\text{for every oriented child state }y.}
\tag{R46C.18}
```

Empty intersection in (R46C.14) says only that the zero-slack states of the
external-profile game are not maximal child grounds.  It does not supply the
simultaneous signed inequalities (R46C.18); active restrictions can point in
incompatible edge directions, and nonactive restrictions can rise under the
same flip.  No valid separation argument deriving (R46C.18) from empty
intersection was found.  This is why the finite evidence below is recorded
as evidence for a new minimizer conjecture, not as a proof from minimality.

## 6. Exact finite exhaustion

Switching invariance lets us set every edge from vertex zero to `+1`, giving
one labeled representative per switching orbit.  Exhaustion produced:

| `n` | `q_n` | exact switching representatives | failures of (R46C.14), all `m` |
|---:|---:|---:|---:|
| 2 | 2 | 1 | 0 |
| 3 | 6 | 2 | 0 |
| 4 | 8 | 6 | 0 |
| 5 | 8 | 12 | 0 |
| 6 | 10 | 12 | 0 |
| 7 | 18 | 3240 | 0 |
| 8 | 20 | 4200 | 0 |

Thus all `7473` exact representatives through order eight pass at every
selector size.  The stored `A_9` also passes for `m=1,...,8`.

This does not support the stronger all-pairs property.  For example, on the
stored `A_8,m=4`, there are ten maximal child-ground pairs but only eight
have a parent-ground lift; at `m=5`, the counts are `120` and `104`.
Affineness needs only one common state, and all eight parent grounds of this
`A_8` are common at both sizes.

## 7. Exact posterior-exclusion link and its limitation

For a vertex coordinate `i`, let `x,y=x^i` be its two parent states and put

```math
B_i(x_{-i})=\mathbb E_S[\mathbf1_{\{i\notin S\}}b_S(x)],
\qquad
r_i(d)=\frac{B_i(d_{-i})}{U(d)}
=\Pr_{\pi_d}\{i\notin S\}.
\tag{R46C.19}
```

The numerator is identical at the two endpoints.  Since `f` is proportional
to `U`, this gives the finite-temperature identity

```math
\boxed{
\chi_i(x,y):=\log\frac{f(y)}{f(x)}
=\log\frac{U(y)}{U(x)}
=\log\frac{r_i(x)}{r_i(y)},
\qquad
\sum_i r_i(d)=n-m.}
\tag{R46C.20}
```

Define the best omitted-selector completion level

```math
a_i^0(e)=\min_{S:i\notin S}k_S(e[S]),
\qquad
\delta_i(d)=a_i^0(d_{-i})-a_d\ge0.
```

Then

```math
\boxed{
-\beta^{-1}\log r_i(d)\longrightarrow\delta_i(d),
\qquad
a_x-a_y=\delta_i(y)-\delta_i(x).}
\tag{R46C.21}
```

In particular, `delta_i(d)>0` means every completion-minimizing selector at
`d` contains `i`.

This locates every positive leading crossing edge:

- if both endpoints are active and their base energies differ by `h>0`, the
  lower-base endpoint `l` has `delta_i(l)>=h`;
- if the active endpoint has a neighbor whose base energy is higher by
  `h>0` but which is endpoint-inactive (the reversal edge of (10.1180)), the
  active endpoint has `delta_i>h`.

Thus each crossing has a distinguished charged endpoint with
`r_i<=exp{-beta h+o(beta)}`.  At any one state, all charged coordinate labels
belong to every leading selector, so their raw number is at most `m`.

This does **not** give a uniform weighted-overlap constant `D`.  An
active--active edge is charged at its lower-base endpoint, while a high-base
active state may be the uncharged endpoint of many such edges.  The tropical
star obstruction is exactly this reverse-fan-in geometry: fixed-size
selector exclusion controls the leaves but not the high center, where the
loads add before squaring.  A project theorem still needs a genuinely new
bound on weighted high-side fan-in (or a replacement for nodewise squaring),
not merely the identity `sum_i r_i=n-m`.

## 8. Reproduction

Run

```bash
.venv/bin/python tmp/tropical_decomposition_r46_check.py
```

The script verifies the exact cancellation numerically at finite temperature,
checks the tropical endpoint maximum and common-line equivalence directly,
exhausts every switching-normalized exact minimizer through order eight at
every selector size, audits the stored `A_8,A_9`, and records explicit
failures of the stronger all-child-ground lift claim.
