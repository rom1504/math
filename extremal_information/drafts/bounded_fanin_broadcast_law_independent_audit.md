# Independent audit of `bounded_fanin_broadcast_law.md`

**Verdict: PASS, with scope qualifications for the finite verifier and the
public description complexity.**

The incidence theorem, neighbour sharpness construction, support discrepancy
constant, quadratic packing construction, and Gram-coordinate double count
all survive independent reconstruction.  In particular, the packing uses a
fixed (although exponentially large) continuation language; it does not use a
future that observes the hidden child.  The draft also correctly keeps the
unrestricted packing separate from simultaneous spectral flatness.

No change to a displayed theorem or constant is required before
canonicalization.  The finite script is a wind tunnel for some local
ingredients, not a computational verification of the two probabilistic
existence arguments in BF.4; this limitation should remain explicit wherever
the script is cited.

## 1. BF.1: contextual cancellation and incidence constants

For a fixed public continuation `C_theta`, uniform-norm nonexpansiveness gives

```math
|R_\theta(z)-R_\theta(z+e_i)|
\le \|H_z-H_{z+e_i}\|_\infty
\le \sum_{e:i\in I_e}|c_e(z)-c_e(z+e_i)|.
```

Taking the supremum over `theta`, summing over `i`, and using
`|c_e(z)-c_e(z+e_i)|<=omega_e` gives exactly

```math
\sum_i d(z,z+e_i)\le\sum_e|I_e|\omega_e.
```

There is no extra factor from the absolute-value response.  If
`|c_e|<=B`, then `omega_e<=2B`, so

```math
h\epsilon S\le2BtE,
\qquad
h\le {2BtE\over\epsilon S}.
```

For `E<=N(N-1)/2` and `S=N^(3/2)`, this becomes

```math
h\le {Bt(N-1)\over\epsilon\sqrt N}
   < {Bt\over\epsilon}\sqrt N,
```

as stated.  The coordinate-dependent version follows by replacing
`h epsilon` with `sum_i epsilon_i`.

The cancellation requires only that the **language** of continuations be
fixed before the realized state is supplied.  A language containing all
`-H_T` is legal: for a particular pair one may choose the already declared
coordinate `T=z` when taking the supremum.  This is not a hidden-dependent
future.  A mechanism that generated a new continuation after observing `z`
would lie outside BF.1, as the draft says.

## 2. BF.2: neighbour separation and its information content

If bit `i` belongs to group `G_j`, its flip changes every coefficient in
`E_j` from `sigma_j` to `-sigma_j` and changes no other coefficient.  Hence

```math
\|H_z-H_{z+e_i}\|_\infty
=2\max_x\left|\sum_{\{u,v\}\in E_j}x_ux_v\right|
=2d;
```

the all-one spin attains the last maximum.  In context `-H_z`, one response
is zero and the other is `2d`, proving the reverse inequality.  Thus every
neighbour distance is exactly `2d`, and

```math
\sum_i d(z,z+e_i)=2gtd
=\sum_e|I_e|\omega_e.
```

The information warning is essential and correct.  Neighbour-separated
hypercube vertices need not be pairwise separated: a response cell of
diameter below the neighbour gap is only forced to be an independent set of
the cube, whose maximum size is `2^(h-1)`.  Therefore two cells (one bit),
not `2^h` cells, are the universal consequence.  BF.2 itself has only the
`g` parity bits as semantic state when `t>1`.

## 3. BF.3: support discrepancy and the Khintchine factor

Let a vertex bipartition retain `m_cross>=m/2` supported edges, and let
`d_u` be the supported cross degree of a left vertex.  Since every nonzero
coefficient has magnitude two, sharp real Khintchine at `p=1` gives

```math
\mathbb E_y\left|\sum_vD_{uv}y_v\right|
\ge {1\over\sqrt2}(4d_u)^{1/2}
=\sqrt{2d_u}.
```

After choosing the left signs rowwise,

```math
\max_{x,y}|x^TD_{LR}y|
\ge\sqrt2\sum_u\sqrt{d_u}
\ge{\sqrt2\over\sqrt N}\sum_ud_u
\ge {m\over\sqrt{2N}}.
```

Here `d_u<=N` is sufficient for the middle inequality.  Negating all left
spins reverses the cross term and preserves the two within-part terms, so
one of the two complete quadratic values has absolute value at least the
cross value.  Thus BF.13 has the stated constant; there is neither a missing
factor of two nor a bilinear-to-quadratic polarization loss.

## 4. BF.4: construction and all parameter regimes

### Spin library

For two independent uniform spins, Hamming distance is `Bin(N,1/2)`, and

```math
\Pr\{|d_H-N/2|>N/4\}\le2e^{-N/8}.
```

There are fewer than `q^2/2` pairs with `q=2^t`.  Thus a pairwise balanced
library exists uniformly for every `t<=cN` once, for example,
`2c log 2<1/8` (with harmless room for the leading factor two).  This also
covers `t=1`; no large-alphabet assumption is used here.

### Balanced cell partition

If `v=u_au_b` has Hamming weight `s in [N/4,3N/4]`, then the number of edges
crossing its cut is

```math
s(N-s)\ge3N^2/16,
```

which is a fraction at least `3/8` of `binom(N,2)`.  In a uniformly random
balanced partition, each size-`d` cell is marginally a hypergeometric
sample.  Its mean crossing count is at least `3d/8`; the multiplicative
lower-tail estimate at threshold `d/4` is at most `e^(-d/48)`.  A union
bound over the `g q^2` cell/pair choices is valid without independence
between cells.  Since `d=Theta(N^(3/2))`, it succeeds throughout the more
restrictive spin-library regime `t<=cN`.

### Fan-in and neighbour response

An edge in cell `j` reads only the `t` bits of symbol `a_j`; the shared
lookup table does not increase fan-in.  If one bit changes, the two symbols
are distinct, at least `d/4` cell coefficients flip, and at the old switching
spin each nonzero coefficient difference contributes with the same sign.
Consequently the context `-H_z` gives response at least
`2(d/4)=d/2=Theta(N^(3/2))`.  The public language contains this context for
every possible `z` before any one child is queried.

### Outer code and pairwise packing

For `q>=4`, greedy deletion of radius-`g/4` balls, together with

```math
|B(g,g/4)|\le2^gq^(g/4),
```

gives size at least `q^(g/4)` already in the worst case `q=4`.  For `q=2`,
the binary Gilbert bound at relative distance `1/4` has rate
`1-H_2(1/4)>0`.  Hence one absolute positive rate works for every
`q=2^t`.

Two codewords differ in at least `g/4` disjoint cells and therefore flip at
least `gd/16=Theta(N^2)` distinct edges.  BF.3 gives
`Theta(N^(3/2))` response in the already-declared context `-H_z`, while the
code has

```math
\log_2|C|=Omega(g\log_2q)=Omega(t\sqrt N).
```

All floors are harmless for sufficiently large `N`, including the endpoint
`t=1` and the upper regime `t=floor(cN)`.  The same absolute theorem constant
can be chosen as the minimum of the library, code-rate, neighbour-gap, and
pair-gap constants.

## 5. Gram-basis double count

Write an alternating form in any linear basis `B_1,...,B_h`.  The coordinate
vector `w_e` of edge evaluation satisfies

```math
(w_e)_i=B_i(p_u,p_v).
```

Every basis vector `B_i` is a nonzero alternating form.  The sampler's
minimum-distance property therefore gives

```math
#\{e:(w_e)_i=1\}\ge E/4
```

for every `i`.  Double counting incidences yields

```math
\sum_e|supp(w_e)|
=\sum_i#\{e:(w_e)_i=1\}
\ge hE/4.
```

Thus average edge fan-in is at least `h/4` and worst-case fan-in is at least
`h/4`; the trivial upper bound is `h`.  With `k=64r^2` and
`h=r(r-1)/2`, this is `Theta(k)` in every basis.  Likewise every coordinate
influences at least `E/4=Theta(k^2)` coefficients.  There is no basis-choice
loophole and no duplicated counting of ordered edges: `E=binom(k,2)` and all
counts are over unordered pairs.

## 6. Scope and flatness audit

BF.4 proves sharpness for unrestricted exact hollow signings, but its
children are not asserted to have cap `O(N^(3/2))`; indeed the draft
explicitly notes that they need not.  The Gram family imports the separate
simultaneous-flatness theorem and escapes BF.1's local scale by linear
fan-in and quadratic coordinate influence.  These two roles are not
conflated.

Two additional resource qualifications are worth keeping visible:

1. BF.4 uses a public spin library of size `2^t` and an exponentially large
   fixed continuation language.  It proves a sharp **fan-in/incidence** law,
   not a low-description or small-query-algebra theorem.
2. A context `-H_w` is a same-support additive overlay.  The child and the
   context are separately signings, while their sum has coefficients in
   `{0,+-2}`.  The result is not closure of exact signings under an appended
   disjoint composition.

These are limitations of interpretation, not defects in BF.1--BF.4.

## 7. Verifier audit

Running

```text
source .venv/bin/activate
python extremal_information/experiments/verify_bounded_fanin_broadcast_law.py
```

completed successfully.  It exhaustively checks BF.2 on its small parity
instance, tests BF.3 on random supports through order nine, and checks the
coherent switching witness on a small local cell.  Its Gram panel verifies
the incidence double count for the displayed deterministic sample lists.

It does **not** construct the asymptotic BF.4 spin library, balanced edge
partition, or outer code.  In the local switching check the acceptance
threshold is two flipped edges for a cell of size twelve, rather than the
theorem's `d/4=3`; that check validates coherent signs, not the asymptotic
cell-balance constant.  The Gram panel uses lists of length `4r^2`, not the
probabilistic theorem's `64r^2`, and checks basis-column degrees rather than
the minimum distance over all nonzero alternating forms.  None of these
limitations affects the analytic proofs, but the script should not be
described as a finite construction or verification of the probabilistic
existence assertions.

