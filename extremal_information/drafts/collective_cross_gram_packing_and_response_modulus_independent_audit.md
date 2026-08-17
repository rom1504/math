# Independent audit: collective cross-Gram packing and response modulus

**Audited files.**
`collective_cross_gram_packing_and_response_modulus.md` and
`verify_collective_cross_gram_packing.py`.

**Verdict.**  CP.1--CP.18 and the labelled-table separation have the stated
constants.  The square-root continuity proof, including replacement of the
trust parameter `t`, is valid.  The hard pair is genuinely realizable by one
involution with fixed one-port data.  Four scope qualifications should be
made explicit:

1. `d_q` is a pseudometric on raw matrix pairs and a metric only after its
   null directions are quotiented (or on a fixed-diagonal slice).
2. CP.11 uses SA.19 in the regime `m/r=kappa/p`, which violates SA.3's stated
   auxiliary assumption `2m>r` when `p` grows.  The trust-region derivation
   itself does not use that assumption, so the formula remains correct, but
   this extension must be stated rather than silently invoked.
3. The `Theta(kappa sqrt(c))` wording for CP.19 needs a uniform margin from
   `kappa^2(1-c)=8`, or simply fixed `kappa<sqrt(8)` as `c` tends to zero.
4. The packing is separated in the **labelled spherical response table**.
   It is not separated after the endpoint label is optimized away, nor is a
   Boolean-cap separation proved.

An independent verifier supplements the canonical checks with exact
projective-ball counts, a pointwise SA.19/CP.11 identity check, and an
explicit common-involution realization.

## 1. CP.1 and its exact factors

Let `Delta K^+=(Delta G+Delta R)/2` and
`Delta K^-=(Delta G-Delta R)/2`.  The two values of the outer sign give

```math
\Delta G+\Delta R=2\Delta K^+,
\qquad
\Delta G-\Delta R=2\Delta K^-.                       \tag{ACP.1}
```

Consequently

```math
d_q={2\over p^2}\max_{\epsilon}max_{\tau\in\{+,-\}}
 |\epsilon^T\Delta K^\tau\epsilon|,                  \tag{ACP.2}
```

so CP.2 has the correct factor two.

For the rank-one states, the negative sectors vanish and
`Delta K^+=ss^T-tt^T`.  After gauging by `s`, write `A` and `D` for the
agreement and disagreement coordinates.  Then

```math
(\epsilon^Ts)^2-(\epsilon^Tt)^2
=4\left(\sum_Au_i\right)\left(\sum_Du_i\right).       \tag{ACP.3}
```

The signs on `A` and `D` are independently selectable, so the maximum is
exactly `4h(p-h)`.  Multiplication by the factor `2/p^2` in (ACP.2) yields

```math
d_q={8h(p-h)\over p^2}
=2\left(1-(s^Tt/p)^2\right).                         \tag{ACP.4}
```

The canonical verifier checks this on its finite code; the independent
verifier checks every projective pair through `p=6`.

Strictly speaking, `d_q` is not a metric on all raw pairs.  Boolean
quadratic probes recover every off-diagonal entry of a symmetric matrix but
only the sum of its diagonal entries.  For example,
`diag(1,-1)` is invisible to every `epsilon^T(.)epsilon`.  Thus `d_q` is a
response pseudometric.  It becomes a genuine metric on its response quotient
and, in particular, on a slice where every sector diagonal is fixed.  This
does not affect CP.2, because the trust response uses exactly the quadratic
queries retained by the pseudometric.

## 2. Projective Hamming packing

There are `2^(p-1)` antipodal classes.  For `r<=p/2`, a projective ball of
radius `r-1` has **exactly**

```math
\sum_{j=0}^{r-1}\binom pj                              \tag{ACP.5}
```

classes: every class at distance below `p/2` has a unique representative
near the center.  Greedy deletion of these balls gives CP.9.  With
`r=ceil(eta p)` and fixed `eta<1/2`, the usual binomial-volume estimate is

```math
\sum_{j<r}\binom pj=2^{(H_2(\eta)+o(1))p},             \tag{ACP.6}
```

which proves rate `1-H_2(eta)-o(1)`.  Since
`h(p-h)` is increasing on `[0,p/2]`, minimum projective distance
`eta p` gives the separation `8eta(1-eta)`.  All constants in
CP.6--CP.9 are correct.  The independent verifier also enumerates the exact
ball formula through `p=6`.

## 3. CP.11 from SA.19

For a fixed outer sign, write `k_d` for the dangerous sector and `k_s` for
the safe sector.  Since

```math
g=k_d+k_s,
\qquad
\sigma h=k_d-k_s,                                    \tag{ACP.7}
```

the numerator in SA.19 is

```math
2\alpha g+\sigma h
=(2\alpha+1)k_d+(2\alpha-1)k_s.                      \tag{ACP.8}
```

Put `t=2alpha-1`; then `2alpha+1=t+2` and
`4alpha^2-1=t(t+2)`.  Therefore SA.19 divided by `rn` is

```math
{1+t\over2}+{(m/r)^2\over2}
 \left({k_d\over t}+{k_s\over t+2}\right).           \tag{ACP.9}
```

Substitution of `k_d=p^2a`, `k_s=p^2b`, and
`kappa=pm/r` gives CP.11 exactly.  The independent verifier checks this
pointwise for 600 choices of `p`, outer sign, sectors, `kappa`, and `t`.

There is a formal hypothesis mismatch worth recording.  SA.3 was stated
after “retain `2m>r`,” inherited from the anti-pin compiler.  Bounded
`kappa=pm/r` gives `m/r=kappa/p`, hence `2m>r` fails for large `p`.
Inspection of SA.21 and the trust-region dual shows that the formula itself
is valid for every `m>=0`; the inequality `2m>r` is never used in that
derivation.  CP should either restate this unrestricted trust-region lemma
or explicitly explain the extension.  It cannot formally cite SA.3 as
currently stated without doing so.

At this budget the number of repeated auxiliary variables is
`pm=kappa r=Theta(sqrt(n))`; the total order is `n+O(sqrt(n))`, and `rn` is
the correct `Theta(n^(3/2))` normalization.  This is a different regime from
the original anti-pin choice `m=r` when `p` grows.

## 4. Square-root continuity and the `t` replacement

Equation ACP.2 gives

```math
|a-a'|,|b-b'|\le\delta/2.                            \tag{ACP.10}
```

Every objective inside the infimum is coordinatewise increasing, hence so
is `Psi_kappa`.  Suppose `a` is raised by `e` and choose an arbitrarily
accurate old minimizer `t>0`.

- If `t>=kappa sqrt(e)`, retaining `t` costs at most
  `kappa^2e/(2t)<=kappa sqrt(e)/2`.
- If `t<kappa sqrt(e)`, put `t'=kappa sqrt(e)`.  The linear term rises by at
  most `kappa sqrt(e)/2`.  Since `t'>t`, both the old dangerous term
  `a/t` and safe term `b/(t+2)` decrease.  The new increment contributes
  exactly `kappa^2e/(2t')=kappa sqrt(e)/2`.

Taking the approximation error to zero proves CP.15.  The cases `e=0` or
`kappa=0` are immediate and may be separated before defining `t'`.  Raising
`b` by `e` costs at most `kappa^2e/4`, because `t+2>=2`.

For unordered pairs, compare each one to their coordinatewise maximum.
Both original values lie below that maximum response and within the same
one-sided bound, so their mutual gap is bounded by that bound, not twice the
bound.  Setting `e=delta/2` gives exactly

```math
\kappa\sqrt{\delta/2}+\kappa^2\delta/8.              \tag{ACP.11}
```

Thus the monotonicity proof and its constants are sound.

## 5. Hard case and realizability

At `a=0`, differentiation of the remaining one-variable objective shows
that its infimum is at the boundary `t downarrow0` when `kappa^2b<4`, giving
CP.17.  The elementary inequality

```math
{1\over t+2}\ge {1\over2}-{t\over4}                  \tag{ACP.12}
```

turns the difference from CP.17 into

```math
{1\over2}\left(1-{\kappa^2b\over4}\right)t
 +{\kappa^2a\over2t}.                                \tag{ACP.13}
```

Minimizing this expression proves the lower half of CP.18.  Testing
`t=kappa sqrt(a)` makes the dangerous-plus-linear cost `kappa sqrt(a)`
and only decreases the safe term, proving the upper half.  Both constants
are correct.

The matrices in CP.19 are genuinely realizable.  For example, take a
four-dimensional involution with one positive unit vector `u` and two
negative orthonormal vectors `v_1,v_2`, and put

```math
w_i^{\mathcal A}=\sqrt{nc}\,s_i u+\sqrt{n(1-c)}\,v_i,
\qquad
w_i^{\mathcal B}=\sqrt{nc}\,t_i u+\sqrt{n(1-c)}\,v_i. \tag{ACP.14}
```

Their squared norms are `n`; their two sector Gram matrices are precisely
CP.19.  They have common diagonal data, and direct evaluation gives
`d_q=2c` and CP.21.  The independent verifier constructs these vectors and
checks all identities numerically.

The stated `Theta` conclusion needs a margin: the lower constant is

```math
\sqrt{1-\kappa^2(1-c)/8}.                             \tag{ACP.15}
```

Strict positivity alone is not a uniform `Theta` constant if the left side
approaches zero.  For the sharp-modulus limit `c downarrow0`, it suffices to
fix any `kappa<sqrt(8)`, as the canonical verifier does with `kappa=1`.
Alternatively, replacing the common safe sector `(1-c)I` by
`(1-c)ss^T` makes `b=0` at query `epsilon=t`; then the gap is exactly
`kappa sqrt(c)` for every `kappa`.  The independent verifier checks this
stronger admissible variant too.

This is Euclidean Gram--Rayleigh realizability, not Boolean-port
realizability.  The draft says “Euclidean-realizable,” so its formal claim is
honest; square-root sharpness has not been proved on the Boolean slice.

## 6. Labelled response-table packing

With `b=0`, direct minimization gives

```math
\Psi_\kappa(a,0)={1\over2}+\kappa\sqrt a.             \tag{ACP.16}
```

For state `s` and query label `epsilon`, this is CP.22.  Given states
`s,t` at projective Hamming distance `h<=p/2`, the common labelled query
`epsilon=s` produces values differing by

```math
\kappa\left(1-{|p-2h|\over p}\right)
={2\kappa h\over p}.                                  \tag{ACP.17}
```

Thus the projective code really is an `exp(Omega(p))` packing in sup norm
on the labelled spherical response tables, at fixed separation when
`kappa` is bounded below.

This claim has two deliberate limits.  First, maximizing over `epsilon`
gives `1/2+kappa` for every rank-one state, so the unlabelled scalar response
collapses.  Second, `Psi` is the spherical certificate, not the exact
Boolean response.  To call the table packing contextual for a concrete
exact-sign model, the declared future-context family must expose or pin the
label `epsilon`; SA.17 itself internally maximizes that endpoint sign.  The
draft acknowledges both collapses, but “contextually visible” should be
read conditional on this labelled query interface rather than as an
exact-Boolean continuation theorem.

## 7. Verifier assessment

The canonical verifier correctly checks:

- finite rank-one Boolean realizations and CP.5;
- random instances of CP.13;
- CP.18 at `kappa=1`;
- one labelled table gap.

It does not check the SA.19 normalization or explicitly build CP.19.  The
independent verifier adds 1,373 exact metric/projective-ball checks, 600
pointwise normalization identities, four explicit common-involution hard
realizations, and nine checks of the unrestricted-`kappa` sharp variant.

Run both with

```bash
./.venv/bin/python \
  extremal_information/experiments/verify_collective_cross_gram_packing.py

./.venv/bin/python \
  extremal_information/experiments/verify_collective_cross_gram_packing_independent_audit.py
```

Both pass.
