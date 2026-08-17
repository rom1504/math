# Independent audit of the planted-ground cap falsifier

Date: 2026-08-17.

**Verdict: PASS after two wording repairs.**  The planting identity, exact cap
equality, uniform hitting argument, constant `32` in the projective-radius
bound, and subexponential support conclusion are correct.  The two repairs
are asymptotic bookkeeping rather than mathematical changes:

1. Theorem PG.1 should either be stated for sequences `(a_n,r_n)`, or its
   finite statement should explicitly assume `1<=r<=N_D`; expressions such
   as `n=o(r)` do not apply to one fixed signing.
2. In Corollary PG.2, “for every `n`” should read “for all sufficiently large
   `n`,” because `r_n<=N_D` is obtained asymptotically.

No canonical file was changed in this audit.

## 1. Exact planting calculation

Let `z_0` be an oriented ground word for `a`, let

```math
D=\{e:a_e(z_0)_e=-1\},\qquad |D|=N_D={E-Q(a)\over2},
```

and obtain `b` by reversing `a_e` on an `r`-set `F subseteq D`.  For every
oriented augmented cut `z`, put

```math
T_z=\{e:z_e\ne(z_0)_e\},\qquad
d_a(z)=Q(a)-\langle a,z\rangle .
```

Since `a_e=-(z_0)_e` on `F`,

```math
\langle b,z\rangle
=\langle a,z\rangle+2\sum_{e\in F}(z_0)_ez_e.
```

The sum is `r-2|F cap T_z|`.  Therefore

```math
Q(a)+2r-\langle b,z\rangle
=d_a(z)+4|F\cap T_z|.                                  \tag{A.1}
```

Both terms on the right are nonnegative.  At `z=z_0` they vanish, so no
upper-bound relaxation is involved:

```math
Q(b)=Q(a)+2r.                                           \tag{A.2}
```

Writing `m_z=|D cap T_z|` and `t_z=|T_z|`, the original deficit is

```math
d_a(z)=2\sum_{e\in T_z}a_e(z_0)_e
=2(t_z-2m_z),
```

and hence

```math
t_z=2m_z+{d_a(z)\over2}.                               \tag{A.3}
```

Thus PG.9, PG.10, and the exact cap assertion PG.3 all pass.

## 2. Sampling domain and simultaneous hitting

The sample space is precisely the uniform law on the `r`-subsets of `D`.
It is available whenever `1<=r<=N_D`.  Under the bounded-cap hypothesis,

```math
N_D={E-Q(a)\over2}=\Theta(n^2),
```

so `r=o(n^2)` implies `r<=N_D` for all sufficiently large `n`.

For fixed `G subseteq D`, the variable `X=|F cap G|` is hypergeometric with
mean

```math
\mu={r|G|\over N_D}.
```

The standard sampling-without-replacement Chernoff bound gives

```math
\Pr\{X<\mu/2\}\le \exp(-\mu/8).                       \tag{A.4}
```

If `|G|>=16nN_D/r`, then `mu>=16n`, so the failure probability is at most
`e^(-2n)`.  The oriented augmented-cut family has at most

```math
2^{n-1}\cdot2=2^n
```

distinct words: vertex assignments are quotiented by global sign and there
are two energy orientations.  Consequently

```math
2^n e^{-2n}=e^{-(2-\log2)n}<1,
```

and a single `F` hits every large set `D cap T_z` as asserted.  Duplicate
sets only decrease the union-bound cost.

## 3. Constants in the global radius

For `z` in the positive oriented `Delta`-shell of `b`, (A.1) gives

```math
d_a(z)\le\Delta,\qquad |F\cap T_z|\le\Delta/4.         \tag{A.5}
```

The hitting dichotomy says either

```math
m_z<{16nN_D\over r},
```

or

```math
|F\cap T_z|\ge {rm_z\over2N_D},
```

in which case (A.5) gives `m_z<=Delta N_D/(2r)`.  A bound valid in either
case is therefore

```math
m_z\le {16nN_D\over r}+{\Delta N_D\over2r}.            \tag{A.6}
```

Combining (A.3), (A.5), (A.6), and `N_D<=E` yields exactly

```math
{t_z\over E}
\le {32n\over r}+{\Delta\over r}+{\Delta\over2E}.     \tag{A.7}
```

Thus neither the factor `32` nor either `Delta` term is missing a factor of
two.

## 4. Orientation, projective geometry, and entropy

Write `z/z_0=tau c(y)`, where `tau in {+-1}`.  If `tau=+1` and `y` differs
from the all-one vertex word in `k` coordinates, then

```math
d_E(z,z_0)=k(n-k).                                     \tag{A.8}
```

If `tau=-1`, then

```math
d_E(z,z_0)=E-k(n-k)
\ge E-\lfloor n^2/4\rfloor
=(1/2-o(1))E.                                         \tag{A.9}
```

Hence an `o(E)` oriented cap cannot mix the two orientations.  After
projectivizing the vertex word, take `k<=n/2`.  Equation (A.8) and
`d_E(z,z_0)<=eta_n E`, with `eta_n=o(1)`, imply `k=o(n)`.  The number of
possible projective words is consequently at most

```math
\sum_{k\le o(n)}{n\choose k}=\exp(o(n)).               \tag{A.10}
```

Moreover, if two shell words are each within `eta_n E` of `z_0`, their
mutual Hamming distance is at most `2eta_nE`, and therefore

```math
{\langle z,z'\rangle\over E}\ge1-4\eta_n.
```

Taking an absolute value only weakens this inequality.  PG.8a and the
failure of fixed-scale projective packing are correct.

## 5. Quantifier audit

For any prescribed `Delta_n=o(n^(3/2))`, put

```math
h_n=\max\{n,\Delta_n\},\qquad
r_n=\left\lceil\sqrt{n^{3/2}h_n}\right\rceil.
```

Because `h_n=o(n^(3/2))`, one has

```math
r_n=o(n^{3/2}).                                        \tag{A.11}
```

If `h_n=n`, then `r_n` is of order `n^(5/4)`, so both `n/r_n` and
`Delta_n/r_n` vanish.  If `h_n=Delta_n>n`, then

```math
{\Delta_n\over r_n}
=\sqrt{{\Delta_n\over n^{3/2}}}\longrightarrow0,
```

and `n/r_n` vanishes as well.  Thus

```math
\max\{n,\Delta_n\}=o(r_n),\qquad
r_n=o(n^{3/2}),                                        \tag{A.12}
```

which is the desired choice
`max(n,Delta_n) << r_n << n^(3/2)`.  Starting from an exact minimizer gives

```math
Q(b_n)=M_n+2r_n=M_n+o(n^{3/2}),
```

while (A.7) tends to zero uniformly over the prescribed shell.

The more parameterized statement is also correct: if `epsilon_n->0` and
`epsilon_n sqrt(n)->infinity`, choosing
`r_n` proportional to `epsilon_n n^(3/2)` gives `n=o(r_n)` and vanishing
normalized excess, and confines every shell with
`Delta_n=o(epsilon_n n^(3/2))`.  For fixed `epsilon>0`, the choice
`r=floor(epsilon n^(3/2)/4)` puts the perturbed signing inside the declared
`epsilon`-halo and confines every `o(n^(3/2))` shell.  Floors, ceilings, and
the condition `r>=1` matter only at finitely many orders.

## 6. Exact falsification scope

The theorem rigorously falsifies:

- a statement that **every** vanishing near-minimizer has a fixed-scale
  projective packing in its exact active shell;
- for each preassigned `Delta_n=o(n^(3/2))`, a statement that every
  vanishing near-minimizer has such a packing in its positive
  `Delta_n`-shell;
- the analogous universal fixed-`epsilon` halo statement when the queried
  shell width is `o(n^(3/2))`.

The quantifiers matter.  The construction does **not** falsify:

- a theorem restricted to exact minimizers;
- a halo theorem whose allowed shell width is comparable to or larger than
  the parent's planted excess `2r`;
- one assertion that there is a shell width chosen adaptively after seeing
  the parent, unless that width is required to be `o(Q(b)-M_n)`;
- the existence of one universal perturbed sequence defeating every
  possible subleading shell schedule simultaneously.

In particular, the result says that exact optimality may be structurally
discontinuous at the `o(n^(3/2))` coefficient-edit scale.  It does not say
that exact minimizers themselves have concentrated shells.

## 7. Archive collision and genuine increment

The identity `Q(b)=M_n+2r` and its frozen-coordinate consequence already
appear as the geodesic planted-face mechanism in
`nearmin_deterministic_inequalities.md`, equations (3.3)--(3.5).  Those
results show that thin-shell barycentric balance can be destroyed by
planting negative ground edges.

The genuine increment in PG.1 is narrower but real: the hypergeometric
hitting set is chosen simultaneously for the complete augmented-cut family.
Together with the exact deficit identity (A.3), this upgrades coordinate
freezing to a uniform geometric theorem for the whole shell, and then to
subexponential projective support.  The planting mechanism itself should
not be counted as new.

## 8. Finite identity checks

No dedicated verifier accompanied the draft at audit time.  As an
independent arithmetic check, I enumerated all oriented augmented cuts for
random signings at every order `3<=n<=7`, selected random nonempty subsets
`F subseteq D`, and checked:

```math
Q(b)=Q(a)+2|F|,
```

```math
Q(a)+2|F|-\langle b,z\rangle
=d_a(z)+4|F\cap T_z|,
```

and

```math
|T_z|=2|D\cap T_z|+d_a(z)/2
```

for every oriented augmented cut in each trial.  All checks passed.  The
finite test validates the signs and factors in the exact identities; the
simultaneous asymptotic conclusion is supplied by the proof above, not by
the experiment.

## 9. Final judgment

After the two wording repairs, the draft is suitable for canonicalization
as a scalable negative result.  Its correct headline is:

> Arbitrarily small normalized near-minimality slack can support a planted
> shell whose entire preassigned sub-slack-width positive face lies in one
> vanishing projective cap.

The exact-minimizer projective-shell question remains open and is precisely
the case not addressed by this construction.
