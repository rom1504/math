# Independent audit of the balanced affine atlas

Date: 2026-08-17.

Audited files:

- `nearmin_balanced_affine_atlas.md`;
- `experiments/verify_nearmin_balanced_affine_atlas.py`;
- the cited FB.1--FB.2, MP.1--MP.3, SC.1, and MC.4 inputs.

## Verdict

**REPAIR (scope and documentation only).**  The mathematical statements
AA.1--AA.4, including all displayed constants and the substitution
`kappa=n^(-1/6)`, `K=ceil(n^(1/3))`, `q=ceil(n^(1/6))`, are correct.  I
found no hidden reconstruction of the full shell in the *stored atlas* and
no circular use of a target response value.  Four claims should nevertheless
be repaired before canonicalization:

1. invoke explicitly the known lower bound `M_n>=c n^(3/2)` when using the
   absolute-bias clause FB.6 in AA.4;
2. call the result an **existential, signing-dependent designed-interface
   certificate**, not an efficient encoder or a response quotient over one
   common exogenous query class;
3. describe the Walsh law as *optimally/asymptotically balanced*, not
   literally zero-bias or “perfectly balanced”; the three bad Walsh ports do
   not by themselves form the balanced atlas--they sit inside the full
   balanced Walsh family;
4. narrow the verifier description: it checks AA.1--AA.2 and the
   selector-energy ingredient of AA.13, but does not instantiate a shell law
   and hence does not computationally check AA.11 or maximize the response in
   AA.13.

These repairs do not alter a theorem, constant, or exponent.  The source
draft should retain its final warning: this is not a mixed-chart theorem,
fixed-ratio recurrence, or reusable cross-order congruence.

## 1. Constant audit for AA.1

Set `P=Q-d` and let `J={i:ell_i<0}`.  If every vertex of `J` is selected
independently with probability `p`, the exact flip identity gives

```math
E[H_D(1^R)-P]
=2pL_-+4p^2\sum_{\{i,j\}\subseteq J}d_{ij}.
```

Randomly completing the spin outside `J` shows

```math
\sum_{\{i,j\}\subseteq J}d_{ij}\ge-Q_-(D[J])\ge-Q(D),
```

and switching/orienting preserves `Q(D)=Q(A)=Q`.  Since every realization
has energy at most `Q`,

```math
d=Q-P\ge2pL_--4p^2Q.                              \tag{A.1}
```

At `p=1`, (A.1) gives

```math
L_-\le(d+4Q)/2\le3Q,
```

because `0<=d<=2Q`.  Thus `p=L_-/(4Q)<=3/4` is admissible, and substituting
it in (A.1) gives exactly

```math
d\ge L_-^2/(4Q),
\qquad L_-\le2\sqrt{Qd}.
```

There is no lost factor of two.  For a nonempty signing `Q>0`, so the
division is harmless.

## 2. Constant audit for AA.2

Since `sum_i ell_i=2P=L_+-L_-`, AA.1 gives

```math
L_+=2P+L_-
\le2Q+2\sqrt{Qd}.                                  \tag{A.2}
```

For a balanced `q`-cell partition and

```math
C_b=2\sum_{i\in J_b}(ell_i)_+ +4Q_-(D[J_b]),
```

MP.1 and (A.2) give

```math
\sum_bC_b
\le2L_++4Q
\le8Q+4\sqrt{Qd}.                                  \tag{A.3}
```

Every cell has at least `k=floor(n/q)` vertices.  Passing from the cheapest
cell to a `k`-subset decreases the positive-field sum, while random
completion proves monotonicity of the one-sided principal cap
`Q_-(D[I])<=Q_-(D[J_b])`.  Finally,

```math
P-H_D(1^S)
=2\sum_{i\in S}ell_i-4\sum_{\{i,j\}\subseteq S}d_{ij}
\le2\sum_{i\in I}(ell_i)_++4Q_-(D[I]).
```

Combining this with (A.3) proves

```math
sigma H_A(x^S)
\ge Q-d-(8Q+4\sqrt{Qd})/q.
```

Thus AA.3's constant is correct.  For a general atom this is a statement
with one fixed **signed orientation**; it only guarantees that the actual
energy has that sign when the displayed total deficit is less than `Q`.
That nontriviality does hold in AA.4 for all sufficiently large `n`.

Odd products of either the spins or their fixed-`sigma` signed-cut words
XOR the masks, so the closure assertion is exact.

## 3. AA.3 and the balance calculation

Uniformly averaging the masks on a `k`-set kills exactly the edges with at
least one endpoint in that set.  Their number is

```math
{n\choose2}-{n-k\choose2}
=kn-{k(k+1)\over2}.
```

Consequently

```math
theta_(n,k)
={kn-k(k+1)/2\over {n\choose2}}
\le {2k\over n-1},
```

so AA.7 is exact.  For `K` independent shell samples, coordinatewise
Jensen gives

```math
E\left|K^{-1}\sum_{r=1}^K Z_(r,e)-E Z_e\right|
\le K^{-1/2}.
```

Only the edge-average is required; no union bound over edges is hidden.
The triangle inequality between the original mean, empirical centre mean,
and chart mean proves AA.11 with precisely

```math
delta+K^{-1/2}+2\lfloor n/q\rfloor/(n-1).
```

Duplicate sampled centres merely reduce the number of distinct charts, so
“at most `K`” is correct.

The proof is existential.  It invokes the minimax measure on the full shell
and principal one-sided cap optimizations to establish that a short atlas
exists.  The final object stores only its selected centres and supports; it
does **not** store or reconstruct the full shell.  On the other hand, the
proof supplies no polynomial-time procedure for finding the measure,
centres, or cheapest support.  Information compression and computational
compression must therefore not be conflated.

## 4. One-sided response claim

Let `I'` be even, `W=(x,(x^{\{i\}})_(i in I'))`, and let `p=|I'|+1`.
Because `p` is odd, every coordinate of `W epsilon` is nonzero.  Relative
to the sign of the common outside-row sum, its sign selector differs from
`x` only inside `I'`; hence it is projectively in the chart.  It realizes
the full field norm.  Therefore

```math
Q-D_(q,d_0)+m\|W epsilon\|_1
\le R_(A,sigma)(mW epsilon)
\le Q+m\|W epsilon\|_1,
```

for every real `m>=0` (and hence for every integer `m` if fields are
restricted to integers).  AA.13 follows with no missing factor of `m` or
two.

The columns of `W` are linearly independent when `I'` is a proper subset
of `[n]`: the outside row forces the total coefficient sum to vanish and
the exceptional rows then force each singleton coefficient to vanish.
Thus the `2^p` sign vectors yield `2^p` distinct labelled fields.  The
claim of `exp(Theta(n/q))` declared endpoints per chart is genuine.

This is nevertheless an `A`-dependent query bank: the centre `x` and
support `I'` are chosen from `A`.  It is a valid compressed answer table for
that designed bank, but not by itself a common contextual pseudometric
quotient in which all landscapes are tested against the same predeclared
external contexts.

## 5. AA.4 substitution and state accounting

For an exact minimizer, put `epsilon=0` and `kappa=n^(-1/6)` in FB.1.  Then

```math
eta_n(kappa)
=C(n^(-1/3)+n^(-1/2))=O(n^(-1/3)).
```

The known rigorous positive lower bound `M_n>=c n^(3/2)` makes
`Q>2kappa n^(3/2)` for large `n`, so FB.6--not merely FB.5--gives

```math
d_0=2kappa n^(3/2)=2n^(4/3),
\qquad delta=O(eta_n(kappa)/kappa)=O(n^(-1/6)).
```

This lower bound needs to be stated in the AA.4 derivation.  The random-sign
union bound supplies the separate upper bound `Q=M_n=O(n^(3/2))`.

With `K=ceil(n^(1/3))` and `q=ceil(n^(1/6))`,

```math
k=floor(n/q)=Theta(n^(5/6)),
K^(-1/2)=O(n^(-1/6)),
2k/(n-1)=O(n^(-1/6)),
```

and

```math
D_(q,d_0)
=2n^(4/3)+O(n^(3/2-1/6))+O(n^(17/12-1/6))
=O(n^(4/3)).
```

All three lines of AA.15 follow.  In particular `D/Q=o(1)`, using both
the upper and positive lower bounds on `M_n`, so the common orientations
are genuinely one-sided rather than merely formal.

A direct labelled description of one chart costs

```math
n+log_2 {n\choose k}+O(log n)
```

bits for `x`, `I`, `sigma`, and scalar parameters.  For `K` charts this is
at most `O(Kn)` bits, and in AA.4 it is `O(n^(4/3))=o(n^2)`.  This is less
than an arbitrary edge signing or complete response table.  It is not a
sublinear number of bits: “sublinear atlas” can only mean a sublinear
number of charts, while its total presentation is superlinear in `n`.

## 6. Stress tests and verifier scope

- **Local ascent:** PASS.  AA.1 controls total negative local-field mass
  without moving the centre, so its signed-edge barycentre is not silently
  changed.
- **Singleton exact-active shell:** PASS.  The input is the FB.1 thick
  shell.  No exact-active multiplicity is claimed.
- **Sparse-flip/planted face:** PASS as a compatibility statement.  Mask
  averaging pays exactly its incident-edge fraction, while inherited
  centre balance is paid separately.  This is not a theorem that a thinner
  planted-face shell becomes unfrozen.
- **Walsh:** PASS after wording repair.  SC.1 supplies a full optimally
  balanced exact-shell law containing three centres whose mixed majority
  has zero energy.  The three centres alone are not asserted or proved to
  have the balance of the full law.  Hence SC.1 kills cross-chart closure,
  not any within-chart statement.
- **Physical composition:** PASS as a negative scope conclusion.  AA.4
  supplies no cancellation theorem for the `Theta(n^2/q)` aligned port
  scale and no reusable cross-order transition.

The executable verifier returned

```text
{'status': 'PASS', 'signings': 1146,
 'oriented_atoms': 82064, 'star_frames': 104448}
```

It exhausts all signings through `n=5` and samples fixed-seed signings for
`n=6,7,8`.  It checks AA.1, the actual minimum partition-cell construction
in AA.2, every mask in the resulting chart, and the selector-energy
ingredient for the largest even star frame.  It does **not** build an
arbitrary measure `mu`, sample its empirical atlas, verify AA.11, or compute
the maximum in AA.13.  Those statements are proved analytically above, but
the script's docstring should not describe itself as a complete verifier of
AA.3.

## 7. Canonical status

After the four wording/scope repairs in the verdict, AA.1--AA.4 may be
recorded as theorem-level results.  Their exact implication is

```text
FB fractional shell balance
 + direct in-place affine thickening
    => a short existential atlas of large coherent local query banks
    != a common all-context quotient
    != mixed-chart coherence
    != a fixed-ratio or cross-order recurrence.
```

The atlas is a real strict one-block certificate.  Its missing information
is cross-chart orientation/transition coherence, not further enumeration
of the atoms already generated inside a chart.
