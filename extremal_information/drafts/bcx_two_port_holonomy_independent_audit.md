# Independent audit: BCX two-port orientation holonomy

**Verdict:** PASS, with one scope clarification.  The one-port orientation
bound, exact block-energy normalization, spectral cap estimates, gauge
statement, and narrow-future ceiling are all correct.  No mathematical
repair is required.  The words `serial` and `reusable quotient` should be
read only for common quadratic continuations that flatten onto one auxiliary
vertex set whose **total** width satisfies the stated bound; they do not prove
closure under adaptive re-encoding or indefinitely many stages.  I ran
`experiments/verify_bcx_two_port_holonomy.py`; all exact and spectral checks
pass.

## 1. One-port blindness

After removing the auxiliary clique, write

```math
F_\sigma^0(t)=\max_{x,y}|\sigma H_A(x)+L_t(x,y)|.
```

The change `y -> -y` reverses `L_t`, and multiplication of the expression
inside absolute value by `-1` gives

```math
F_-^0(t)=F_+^0(t)
```

for every query, with no optimizer assumption.  Adding `H_C` perturbs each
cap by at most `Q(C)`, so

```math
|F_+(t)-F_-(t)|<=2Q(C)=q(q-1)<n.
```

A coordinatewise bound by `2Q(C)` gives the same sup bound; its oscillation
is at most `4Q(C)`, so half-oscillation is at most `2Q(C)`.  Both normalized
one-port distances therefore vanish at the `n^(3/2)` scale.

## 2. Exact parent energies and the absolute cap

Trace zero is essential and is used correctly.  Since
`A=mathcal H-diag(mathcal H)` and `tr mathcal H=0`, Boolean energies satisfy

```math
H_A(x)=x^Tmathcal Hx/2,
\qquad H_(-A)(x)=-x^Tmathcal Hx/2.
```

Thus the complete exact-sign parents in BH.9 are precisely one half of the
quadratic forms of

```math
M_+=T_+\otimes mathcal H,
\qquad M_-=T_-\otimes mathcal H.
```

The actual parent diagonals are zero; inserting the diagonals of
`+-mathcal H` in this calculation adds zero total Boolean calibration because
each block has trace zero.  Every genuine off-diagonal coefficient is an
exact sign.

The block spectra are

```math
||T_+||=2,
\qquad T_-^2=2I,
\qquad ||T_-||=sqrt2.
```

Since `||mathcal H||=q`, the two full matrix norms are `2q` and
`sqrt2 q`.  For a Boolean vector of length `2n`,

```math
|z^TMz|/2<=||M||\,||z||_2^2/2=n||M||.
```

This is an upper bound on the **absolute** Boolean cap, not merely its
one-sided maximum.  It yields exactly

```math
Q(P_+)<=2qn,
\qquad Q(P_-)<=sqrt2 qn.
```

At `(1,1)`, regularity gives two internal contributions `qn/2` and one
cross contribution `qn`, attaining `2qn`.  Therefore the first cap is exact
and the fixed gap `(2-sqrt2)qn` is rigorous without knowing the second cap.
Dividing by `(2n)^(3/2)` gives the stated total-order constant.

## 3. Gauge covariance and logical congruence claim

For a fixed switch `s`, conjugating both shores by `D_s` sends

```math
(A_s,+-A_s,mathcal H_s)
\longmapsto(A,+-A,mathcal H).
```

Hence the block calculation is identical.  The one-port blindness proof is
matrix-independent and also applies to each orientation pair
`A_s,-A_s`.  The continuation exposing a given pair is fixed before choosing
between its two orientations; no coefficient depends jointly on that hidden
orientation.

The conclusion is correctly scoped.  It disproves congruence of the narrow
one-port BCX response metric under enlargement to general dense block
composition.  It does not claim that `-A_s` belongs to the selected BCX
switching code, nor that every pair of distinct switching states is close.

## 4. Narrow-future ceiling

For an arbitrary old--new block `B`, deleting the auxiliary internal energy
again makes the two orientations exactly equivalent by the global change
`y -> -y`.  Reintroducing a hollow exact signing `C` changes each cap by at
most `Q(C)`, so

```math
|R_+-R_-|<=2Q(C).
```

Every energy of an `m`-vertex exact signing is a sum of
`binom(m,2)` unit terms, hence

```math
2Q(C)<=2{m\choose2}=m(m-1).
```

This proves the claimed ceiling independently of the size, signs, or norm
of `B`.  If `m=o(n^(3/4))`, its square is `o(n^(3/2))`.  A bounded number of
BCX-width ports has total width `O(sqrt n)` and remains well inside this
regime, even when all interactions among appended spins are absorbed into
`C`.

This also fixes the precise scope of the word `serial`.  A fixed sequence of
ordinary exact-sign graph attachments, with all variables retained until the
final maximization, flattens to BH.20: collect every newly appended spin in
`y`, put every original--new edge in `B`, and put every edge among new spins
(including cross-stage edges) in `C`.  The proposition therefore covers such
serial attachments whenever their total appended width is `o(n^(3/4))`.  It
does **not** by itself cover a procedure that, between stages, applies an
orientation-dependent map, discards and re-encodes a maximized value, or
changes the allowed future language.  Nor does it give a finite quotient
reusable for unbounded depth: after enough stages the total width can leave
the regime of BH.21.

For canonical use, the safest replacement for the concluding transition
sentence is:

> Orientation is invisible at the `n^(3/2)` scale to every common quadratic
> continuation that can be flattened into `o(n^(3/4))` appended vertices.
> This includes any fixed number of retained-variable BCX-width stages, with
> total appended width `O(sqrt n)`.

The proposition assumes the auxiliary internal coefficients form an exact
signing for the final edge-count bound; the orientation identity itself is
valid for arbitrary `B`.  This matches the statement.

## 5. Evidence and scope

The verifier confirms:

* every one-port query at `n=4` and the full explicit MM query source at
  `n=16`;
* exact two-shore caps at `n=4`;
* both block spectral norms at `n=4,16`;
* 24 independent narrow-future instances.

The theorem is a scalable exact falsifier, not a finite numerical inference.
It isolates one relative orientation bit as composition-created information.
It does not establish incompressibility for arbitrary bridges or say
anything about near-minimizers.  Subject to these explicit limitations, the
draft is rigorous and ready for canonical use.
