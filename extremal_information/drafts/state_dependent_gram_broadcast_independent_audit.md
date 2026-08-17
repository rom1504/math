# Independent audit of `state_dependent_gram_broadcast.md`

**Verdict: REPAIR, with the core theorem and all displayed constants passing.**

The finite packing theorem at order `k=64r^2` is rigorous as stated.  The
sampler, discrepancy, simultaneous spectral-flatness, contextual packing,
and information-scaling arguments all survive independent reconstruction.
The construction also genuinely lies outside Theorem 21.23: the hidden form
changes a dense set of cross coefficients, rather than entering through
bounded onsite children.  It does not use a hidden-state-dependent future.

Before canonicalization, however, the draft should narrow three claims about
the block lift and the continuation language.  These are scope repairs, not
repairs to the finite theorem GB.4.

## 1. Probability and constant audit

### Lemma GB.1

For a nonzero alternating form `B` of rank `rho>=2`, uniform independent
`p,q in F_2^r` give

```math
Pr(B(p,q)=1)
=Pr(p\notin\operatorname{rad}B)/2
=(1-2^{-\rho})/2\ge3/8.
```

Changing one sampled label changes at most `k-1` unordered-pair indicators.
For

```math
t={1\over8}{k\choose2}={k(k-1)\over16},
```

McDiarmid gives

```math
\exp\left(-{2t^2\over k(k-1)^2}\right)=\exp(-k/128).
```

With `h=r(r-1)/2` and `k=64r^2`, the union exponent is

```math
h\log2-k/128
={\log2\over2}r(r-1)-{r^2\over2}<0.
```

Thus the distance-`1/4` evaluation code exists.  Repeated labels are allowed
and cause no issue, since the proved positive minimum support itself implies
injectivity.

**Minor notation repair.**  Equation GB.9 uses `c_B(i,j)` for both orders and
for the diagonal, whereas GB.2 only defines it for `i<j`.  Define the natural
symmetric extension `c_B(i,j)=B(p_i,p_j)` and `c_B(i,i)=0` before GB.9.

### Lemma GB.2

If `m>=k^2/16` unordered coefficients are nonzero, some vertex bipartition
has at least `e>=m/2>=k^2/32` nonzero cross entries.  For a row with `d_i`
coefficients of magnitude two, sharp real Khintchine at `p=1` yields

```math
E_y\left|\sum_jD_{ij}y_j\right|
\ge {1\over\sqrt2}(4d_i)^{1/2}=\sqrt{2d_i}.
```

Since `d_i<=k`,

```math
\sum_i\sqrt{d_i}\ge {e\over\sqrt k}.
```

After choosing the left signs rowwise, this gives the displayed
`(sqrt(2)/32)k^(3/2)` cross value.  Negating all left signs preserves both
within-part terms and reverses the cross term, so one of the two complete
spin vectors has full absolute quadratic value at least the cross value.
There is no missing factor of two.

### Lemma GB.3

For a unit vector `z`, the Rademacher coefficients of
`z^T A_B z` are `2z_iz_j`, and

```math
\sum_{i<j}(2z_iz_j)^2
=2\left(1-\sum_i z_i^4\right)\le2.
```

Hoeffding therefore gives `2 exp(-t^2/4)`.  A `1/4`-net has at most `9^k`
points, and the symmetric quadratic-form net inequality costs exactly a
factor two.  Taking `t=4sqrt(k)` gives operator norm at most `8sqrt(k)`.
The total failure probability is

```math
2\exp(k\log9+h\log2-4k)<1,
```

because `h/k<1/128` (indeed the exponent is strongly negative already at
`r=2`).  Multiplying the random base signs by a fixed character preserves
entrywise independence for each fixed `B`; independence between different
`B` is neither asserted nor needed for the union bound.

### Theorem GB.4

If `B!=T`, then `B+T` is nonzero, so at least
`k(k-1)/8>=k^2/16` unordered coefficients of `A_B-A_T` have magnitude two.
Lemma GB.2 gives the lower bound.  Also

```math
Q(A_B-A_T)
\le {k\over2}\|A_B-A_T\|_{op}
\le {k\over2}(16\sqrt k)=8k^{3/2}.
```

For distinct `B,C`, the predeclared coordinate `T=B` has
`R_B(B)=0` and `R_B(C)>=delta k^(3/2)`, where
`delta=sqrt(2)/32`.  Hence the response vectors are pairwise separated.
A deterministic summary with uniform error strictly below `delta k^(3/2)/2`
cannot identify two states, so its range has at least `2^h` elements.  Finally

```math
h={r(r-1)\over2}\ge {r^2\over4}={k\over256}
```

for `r>=2`.  The information claim is therefore `h` bits, not merely `h`
states.

## 2. Query legality and information content

The context family is fixed before the hidden state is selected:

```math
\mathcal T=\{-H_T:T\in Alt(V)\}.
```

For any one query `T`, its coefficients depend on `T`, the public labels,
and the public base signing, but not on the queried child state `B`.  Choosing
`T=B` only after selecting a pair to distinguish is the ordinary operation
of taking a supremum over a fixed contextual language.  It is **not** a
state-dependent context.

Nor does a query store an arbitrary `Theta(k^2)` coefficient table: after
the shared nonuniform public data `(A,P)` are fixed, both a child and a query
are generated from `h=Theta(k)` bits.  The parity identities of the
alternating-form evaluation code are real constraints.  The result is thus a
conditional response-packing theorem for one fixed public base family.

There are two important limitations which should be stated alongside that
conclusion.

1. The public base signing itself has `Theta(k^2)` nonuniform description and
   is obtained only by probabilistic existence.  The theorem compresses the
   **hidden variation conditional on this base**; it is not an explicit
   uniform encoder for arbitrary dense signings.
2. The separator context is a negative clone on the same variable and edge
   set.  Although the child `A_B` and context `-A_T` are each exact hollow
   sign matrices, their sum has coefficients in `{0,+-2}`.  Thus GB.4 is
   valid for the declared additive/overlay response language, but it does
   not by itself give a parent which remains in the exact-sign class after
   composition, nor an appended disjoint-future theorem.  The draft should
   say this explicitly rather than letting “genuine dense block bridges” be
   read as closure of exact signings under the query operation.

The term “symplectic Gram form” is also slightly too strong for a possibly
degenerate `B`; “alternating (presymplectic) Gram form” is accurate.

## 3. Block-lift audit

Let `D=A_B-A_T` and let `W_n` be symmetric Hadamard with a Boolean regular
vector `u`, `W_nu=sqrt(n)u`.  On block spins `x_i=s_i u`,

```math
\sum_{i<j}D_{ij}x_i^TW_nx_j
=n^{3/2}\sum_{i<j}D_{ij}s_is_j.
```

Since `\|D\otimes W_n\|=\|D\|\sqrt n`, the lower and upper response bounds
indeed multiply by `n^(3/2)` and become the displayed constants times
`N^(3/2)` for `N=kn`.  Symmetry, off-diagonal signs, and hollowness are exact;
common hollow sign diagonal blocks can be added and canceled by the negative
context.

Two qualifications are required.

1. If the draft claims that the **individual lifted children** retain an
   `O(N^(3/2))` Boolean cap uniformly when both `k` and `n` vary, the common
   diagonal blocks must themselves be chosen spectrally flat
   (`O(sqrt(n))` operator norm), or that individual-child claim must be
   dropped.  Arbitrary common sign diagonal blocks cancel in the contrast
   but need not be spectrally flat.
2. The hidden information is `h=Theta(k)=Theta(N/n)` bits.  It is a positive
   rate per total Boolean variable only when `n` stays bounded (taking
   `n=1`, or any fixed regular-Hadamard order, suffices).  If `n` grows, the
   rate is `Theta(1/n)` and vanishes.  The response gap remains a fixed
   multiple of `N^(3/2)` either way.

These do not affect unlifted Theorem GB.4, where `n=1` and every child is an
exact hollow signing with cap at most `4k^(3/2)`.

## 4. Collision with existing no-go theorems

Theorem 21.23 assumes that hidden state enters through bounded state-local
children while every cross interaction is state-independent.  Here every
nonzero change of the hidden alternating form changes at least
`k^2/16` cross-edge coefficients.  At total order `k`, this exceeds the
`Omega(k^(3/2))` atom-influence threshold of (21.118) by
`Omega(sqrt(k))`.  Therefore GB.4 realizes the explicitly advertised escape
case and does not contradict or evade a hypothesis covertly.

It also does not collide with the connected Walsh-flux packing: that theorem
keeps hidden bits in bounded-size local gadgets and consequently loses
`h^(-3/2)` after total normalization.  GB.4 instead broadcasts every
nonzero hidden difference over a constant density of the common edge set.
The price is a tailored exponentially large cancellation language and a
state-dependent dense coefficient family.

The construction uses neither scalar finite-fibre decomposition nor
independently paid channels.  The lower bound is a direct Boolean
discrepancy estimate on the whole coefficient contrast.  On the other hand,
it should not be advertised as a low-rank bridge: the random base dressing
generically destroys low rank, and no rank bound is proved or used.

## 5. Verifier audit

Running

```text
source .venv/bin/activate
python extremal_information/experiments/verify_state_dependent_gram_broadcast.py
```

completed successfully.  It checks:

- the three probability inequalities for `2<=r<=49`;
- exhaustive bilinearity/alternation panels;
- deterministic samplers of size `4r^2` for `2<=r<=5`;
- complete absolute response tables for the full label sets at `r=2,3`;
- the order-four scaling arithmetic for the smallest lift.

The output has minimum normalized off-diagonal responses `0.75` at `r=2`
and `1/sqrt(2)` at `r=3`.

The script intentionally does **not** construct the simultaneous random base
`A`, check all its spectral norms, construct a regular Hadamard lift, or test
the common diagonal blocks.  Those steps rest on the analytic proofs above.
The current module docstring already says that the probabilistic existence
claims are not computationally constructed, so this is acceptable; the
script should not be cited as empirical verification of Lemma GB.3.

## 6. Required repairs before canonicalization

1. Define the symmetric extension of `c_B` used in GB.9 and replace
   “symplectic” by “alternating/presymplectic” where degeneracy is allowed.
2. State that the contextual result uses same-support additive overlay; the
   child and future are signings separately, but their combined queried
   landscape is not an exact signing.
3. In the block lift, distinguish the fixed normalized response gap from the
   information rate `Theta(1/n)` per total variable, and either require flat
   common diagonal blocks or refrain from claiming uniform individual-child
   cap control for varying `n`.
4. State that the `h`-bit compression/packing is conditional on a shared
   nonuniform `Theta(k^2)`-bit public base signing.

With those repairs, the core result is safe to canonicalize as a rigorous
structured counterexample to any total-scale ceiling based only on hidden
coordinate count or individual spectral flatness.

## Repair resolution

Before canonicalization, the draft defined the symmetric evaluation word,
used presymplectic language, made same-support additive overlay explicit,
charged the quadratic shared public base, and separated fixed-inner-order
positive rate from the `Theta(1/n)` lifted rate.  It also conditions any
varying-`n` individual cap claim on spectrally flat common diagonal blocks.
These changes implement all four required repairs.
