# Wave 31 root route: partial mesoscopic collisions

## Status

Verified finite-set-system theorem.  It strictly weakens the Wave 30
common-coset sufficient event from all `r` sampled selectors to any fixed
positive fraction of them.  It does not prove that the resulting event has
the required probability for signing minimizers.

## 1. Exact partial-collision sandwich

Let `Omega` be finite, let `H_a subseteq Omega` (`a in mathcal A`) be an
arbitrary family of `M` hit sets, and put

```math
delta_* = inf_{w in Delta(Omega)} max_a w(H_a).
```

For integers `1<=s<=r`, let `J_(r,s)` be the event that some candidate hit
set contains at least `s` entries of an iid `r`-sample, counted with
multiplicity.  Write `h_w=max_a w(H_a)`.  Then

```math
h_w^s <= Pr_(w^r)(J_(r,s))
       <= M binom(r,s) h_w^s.                         (P1)
```

For the lower bound, retain a candidate attaining `h_w` and require the first
`s` samples to fall in it; the other samples are unrestricted.  For the
upper bound, if a candidate hits at least `s` positions, some `s`-subset of
positions is entirely in that hit set.  Union-bound over the candidate and
the `binom(r,s)` position subsets.

Taking `s`-th roots and infima gives the exact analogue of (10.890):

```math
[M binom(r,s)]^(-1/s)
 inf_w Pr_(w^r)(J_(r,s))^(1/s)
 <= delta_*
 <= inf_w Pr_(w^r)(J_(r,s))^(1/s).                    (P2)
```

No exchange of an infimum with a nonlinear operation is hidden here:
`x -> x^(1/s)` is increasing, so `inf_w h_w^s=(inf_w h_w)^s`.

## 2. Target exponent

Use the Wave 30 candidate bound `M<=(2k)^n`, put `L=n^(3/4-c)`, and take

```math
r=ceil(n log(2k)/L).
```

Fix `alpha in (0,1]` independent of `n` and set `s=ceil(alpha r)`.  Then

```math
log M/s <= L/alpha,
log binom(r,s)/s <= h(alpha)/alpha+o(1)=O_alpha(1).
```

Consequently the following **partial mesoscopic collision lemma** is already
sufficient for the fractional-cover target:

> Uniformly for every selector law `w`, an iid `r`-batch has probability at
> least `exp(-K rL)` that one eligible row-good coset hits at least `alpha r`
> of its selectors.

Indeed (P2) gives

```math
delta_* >= exp[-(K/alpha+1/alpha+o(1))L],
tau_*   <= exp[(K/alpha+1/alpha+o(1))L].                (P3)
```

Conversely, `delta_*>=exp(-K L)` implies
`Pr(J_(r,s))>=delta_*^s>=exp(-K sL)`, so for fixed `alpha` the partial-event
criterion is exponent-equivalent to the fractional target, just as the
all-`r` event was.

## 3. Bounded batch covers are enough

Let `C_T` be the event that all `r` sampled selectors can be covered by the
union of at most `T` candidate hit sets.  On `C_T`, one of those candidates
hits at least `ceil(r/T)` sample positions.  Therefore

```math
Pr(J_(r,ceil(r/T))) >= Pr(C_T).                          (P4)
```

For every fixed `T`, a bound `Pr(C_T)>=exp(-K rL)` proves
`tau_*=exp(O_(K,T)(L))`.  In particular, a deterministic theorem partitioning
every mesoscopic batch into a bounded number of common-coset subbatches is
sufficient.  One globally common coset is not required.

There is also a power-saving tradeoff when `T` grows.  Start with a structural
parameter `c_0>0`, so `L_0=n^(3/4-c_0)`, and suppose
`T<=n^eta` for some fixed `eta<c_0`.  Taking `s=ceil(r/T)` in (P2), the
candidate and event-probability losses are both `O(TL_0)`, while the binomial
loss is only `O(log T)` after the `s`-th root.  Hence

```math
tau_* <= exp[O(TL_0)]
       <= exp[O(n^(3/4-(c_0-eta)))].                     (P5)
```

The original row and tolerance bounds with exponent `c_0` are stronger than
those required with `c'=c_0-eta>0`.  Thus a decomposition into `n^eta`
common-coset groups still proves convergence whenever `eta<c_0`; it merely
spends part of the power saving.  The same conclusion holds for a partial hit
fraction `alpha_n>=n^(-eta)`.

This does not make pointwise planting sufficient.  Planting covers a batch
by `T=r` candidates.  Then `s=1`, and (P2) loses the full candidate entropy
`log M=Theta(n log n)`, recovering no improvement over self-incidence.  To
retain any power saving, `T` must grow strictly slower than `n^c` (or a
stronger event-probability estimate must compensate for its growth).

## 4. Completion/signature consequence

Combine the event with (10.894)--(10.896).  It now suffices that, with
probability `exp(-O(rL))` against every selector law, **some subbatch of at
least `alpha r` selectors** admits favorable completions with

```math
J=O(L/log n),\qquad C_sig=O(n^(9/4-c)).                   (P6)
```

The constrained-hashing theorem builds one eligible coset hitting that
subbatch, and (P3) completes the fractional exponent.  More generally, it is
enough to split the full batch into `O(1)` completion families satisfying
(P6) separately.  More generally, `n^eta` such families are enough under the
preceding `eta<c_0` tradeoff.

Both structural statistics are monotone under discarding witnesses.  The
coordinate-signature partition of a subfamily is a coarsening of the full
partition, so its class count cannot increase.  Its coarsest signature coset
is a subset of the old coarsest coset (after regauging by any retained base
word), hence its `C_sig` cannot increase either.  Positive-fraction pruning
can therefore remove exceptional completions without creating a new row-cap
cost.

This change matters for agreement methods: a theorem that decodes one global
object explaining a constant fraction of local views has the right *form* of
conclusion.  The hard premises remain—the present project has not proved an
overlap-test success probability, adversarial-law uniformity, signature
entropy, or `C_sig` control.

## 5. A low-VC signature successor

The positive-fraction formulation gives a concrete higher-order alternative
to spanning-tree Hamming variation.  For a chosen subbatch of witness words,
view each coordinate signature as a subset of the subbatch indices.  If this
family of coordinate subsets has VC dimension at most `d`, Sauer--Shelah gives

```math
J <= sum_(j=0)^d binom(s,j)=O_d(r^d).                    (P7)
```

At the project scales, this meets `J=O(L/log n)` whenever

```math
d(1/4+c)<3/4-c.                                         (P8)
```

In particular, VC dimension at most two suffices for every fixed `c<1/12`,
and VC dimension at most one suffices for every fixed `c<1/4`.  Thus a
possible successor to pairwise agreement is: with the required probability,
find a positive-fraction favorable subbatch whose coordinate signatures do
not shatter any triple (for `c<1/12`) and whose `C_sig` has the target bound.
This is only a reformulation of the signature-count input, not a proof that
principal minimizers have low VC dimension.  Generic binary word families
need not satisfy it.

## 6. Falsification criterion

For fixed `alpha>0`, a sequence of selector laws with

```math
Pr_(w^r)(J_(r,ceil(alpha r)))=exp(-omega(rL))             (P9)
```

falsifies the fractional coset implementation by the upper half of (P2).
The Wave 30 abstract independent-label construction has this stronger
property for every fixed `alpha`, but remains outside the exact-minimizer
class.  A minimizer-specific instance of (P9) is still unknown.
