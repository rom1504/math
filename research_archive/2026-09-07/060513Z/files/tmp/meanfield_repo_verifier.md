# Independent repository audit: Theorem 16.13

**Scope.** I audited `Theorem 16.13`, Example 29, Section 17 of
`axioms.md`, `drafts/heterogeneous_meanfield_response.md`, and
`experiments/verify_meanfield_response_state.py`. I rederived the exact
contextual quotient, metric identity, histogram count/error, bilinear roof
congruence, strict quotient, and sharp curvature threshold. I also ran the
experiment and performed a separate exhaustive small-instance test of roof
congruence under appended blocks.

## Verdict

**The mathematical theorem is sound for the explicitly declared context
family and maximization convention.** I found no circular use of the desired
state, no false sharpness assertion, and no counterexample to the composition
law. The integration has incorporated the substantive qualifications from
the earlier theorem audit: full chemical-potential range, anchored baseline,
fixed pair normalization, fixed-mass versus size-uniform thresholds, and the
distinction between exact-grid and macroscopic metric entropy.

Three precision/verification issues remain. The first is worth fixing before
calling the statement fully theorem-grade; the other two are documentation
and test-coverage improvements rather than mathematical errors.

## Concrete remaining issues

### 1. Equation (16.98) has no quantified parameter range

The phrase “whenever the displayed scale is nontrivial” is not a precise
hypothesis for the claimed `Omega` bound. The construction works, but the
statement should say, for example:

> There are universal constants `c,C>0` such that, for `B>0`,
> `epsilon>0`, and
> `L=min{n,sqrt(B/epsilon)}>=C`, one has
> `log_2 K_(epsilon n)>=cL`.

Indeed choose

```math
q=floor(L/16),\qquad s=floor(n/(4q)),\qquad d=B/(2q).
```

For `L` larger than an absolute constant, `q<=n/16`,
`q^2<=B/(256 epsilon)`, and `s>=n/(8q)`. Put `q` disjoint response tents of
half-width `d` inside the field range and allocate `2s` sites to each.
Every tent has height

```math
sd>=Bn/(16q^2)>2 epsilon n,
```

so the `2^q` bit choices are pairwise separated at decoding radius
`epsilon n`. This supplies the omitted quantifiers without changing the
claimed order.

### 2. The executable does not test the central roof congruence

`check_quadratic_roofs` checks the two-site strict example, the endpoint
threshold, and **raw** associative quadratic merges. It does not compute
least concave majorants or test (16.100), despite its name. Thus the script is
not presently executable verification of the most novel assertion.

I independently exhaustively grouped all field multisets of masses at most
three over `{-1,0,1}` by their upper concave hull, for `J=-2,...,3`, appended
every block in the same finite family, and compared terminal responses on a
half-integer field grid. All 27,300 same-roof contextual comparisons passed.
It would be useful to incorporate a rational upper-hull version of this
check into the committed script.

### 3. Formula (16.99) should declare domains

For complete formal clarity, state that if the child masses are `n_A,n_C`,
then `u in [0,n_A]`, `v in [0,n_C]`, `t in [0,n_A+n_C]`, with the roofs
piecewise-linearly interpolated on their intervals. The draft makes this
interpretation recoverable, but the theorem formula itself leaves the
optimization range implicit.

There is also one mildly ambiguous summary sentence in `README.md`: the
histogram has simultaneously sublinear bits and sublinear error only after a
choice such as `eta_N=B/sqrt(N)` (which the proof draft explicitly gives).
Adding that choice to the summary would prevent “fixed eta” from being read
as giving sublinear absolute error.

## Mathematical checks

### Exact local-field quotient

The declared future family is not a hidden lookup table. It consists of one
uniform scalar field plus anonymous appending. Exchange gives
`p_A(k)=sum_(j<=k)a_j`, and every `k` is supported by a scalar field in
`[-B,B]`. Hence biconjugacy recovers the whole concave profile. For equal
masses, max-plus conjugacy and biconjugacy are inverse nonexpansive maps in
the sup norm, proving the exact metric identity (16.94). Empty future queries
give necessity; sorted multiset union gives sufficiency and composition.

### Histogram rate

The endpoint grid has actual spacing
`Delta=2B/(M-1)<=eta`, so nearest rounding costs at most `eta/2` per occupied
site. Every weak `M`-bin composition of mass `n` is attainable, giving the
exact binomial count. Distinct histograms have sorted lists whose first
differing entry changes a prefix sum by at least `Delta`, so strict error
below `Delta/2` indeed needs all grid states. No site is requantized after a
merge, which is exactly why error depends on total mass and not depth.

### Bilinear roof congruence

The key identity can also be seen directly:

```math
L_(A sqcup C)(lambda)
=max_l {q_C(l)+lambda l+L_A(lambda+Jl)}.
```

Thus equality of the linear response of `A` is preserved after appending
any `C`; symmetry gives a two-sided congruence. Equivalently, independent
occupancy mixtures prove the lower comparison in the concave-envelope
formula, while pure occupancy pairs plus outer concavification prove the
reverse comparison. The cocycle

```math
Juv+J(u+v)z=Jvz+Ju(v+z)
```

proves associativity. This does not reconstruct the raw hidden fibres.

### Strict quotient and curvature threshold

For `0<a<min(B,J/2)`, `(0,0,J)` and `(0,a,J)` have the same endpoint chord
and terminal response `max{0,J+2lambda}`. The fixed-mass bound

```math
q_A(k)-k q_A(n)/n
<= k(n-k)(2B/n-J/2)
```

proves sufficiency of `J>=4B/n`. Taking `k` fields at `B` and the rest at
`-B` attains the heterogeneity term and proves necessity uniformly over the
class. Maximizing over `n>=2` gives the stated size-uniform sufficient
threshold `J>=2B`; equality may tie, as correctly noted.

## Executable result

The committed script completed successfully:

```text
linear-field biconjugacy coordinates: 18900
contextual metric identities: 3150
max-plus/sorted-union identities: 3150
quantized response bounds: 87000
histogram merge identities: 3480
strict quadratic roof collapses: 81
sharp-threshold chord inequalities: 13000
below-threshold chord obstructions: 8
quadratic associative merges: 2400
```

`python3 -m py_compile`, control-character scanning, and `git diff --check`
also passed.

## Bottom line

I recommend accepting Theorem 16.13 after quantifying (16.98). Adding the
roof-congruence finite test and explicit domains in (16.99) would make the
verification trail match the strength of the theorem, but neither changes
the mathematical conclusion.
