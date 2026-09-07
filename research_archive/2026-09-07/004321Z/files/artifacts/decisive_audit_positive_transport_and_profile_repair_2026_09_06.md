# Positive transport: independent audit and a uniform profile-repair lemma

Date: 2026-09-06. Independent audit PASS for
`decisive_director_positive_transport_limit_2026_09_06.md` and
`decisive_bridge_positive_homogeneous_pressure_2026_09_06.md`.
The new finite continuity lemma below concerns the actual nonnegative
row-permutation partition function, not a replacement Hilbert norm.

## 1. Audit of the exact finite-profile limit

The normalization is the product of `n` row multinomial counts, whose
logarithm divided by `n^2` tends to `sum_a pi_a H(nu_a)`. The numerator
has one variable per UNDIRECTED edge, so its entropy/energy contribution
has factor `1/2`. Averaging within an ordered pair of vertex classes
gives `gamma_ba=gamma_ab^T`; symmetry of the physical kernel makes this
orientation change harmless. Only the NEIGHBOR-AVERAGED row marginal
is constrained to `nu_a`. No separate fixed marginal per class pair is
licensed by the prescribed row multisets.

For the lower bound, independent edge-pair draws from a fixed feasible
family have every row count within `o(n)+O(sqrt(n log n))` of its target,
uniformly over vertices. Changing the color at ONE endpoint changes
only that endpoint's row count, so all rows can be repaired independently.
The total number of changed endpoint incidences is `o(n^2)`, and the
preimage count is `exp(o(n^2))`. Fixed strict positivity of `K` bounds
the log-weight loss by this number times `osc(log K)`. Information-density
and energy concentration for independent edges give the correct number
and weight of typical preimages. Zero entries of the optimal coupling
cause no problem because they are never sampled; the nonzero entries
are fixed before taking the order limit.

This proves the director's formula at all orders admitting the indicated
integer row counts and class proportions. No transitivity of vertex
classes or homogeneous-profile optimality is assumed. The homogeneous
formula has the equivalent `-D(gamma||nu tensor nu)` normalization.

The physical quantization estimate in the companion note is also correct:
on the radius-`B` ball each endpoint gradient of
`log K_R(x,y)=-t(||x||^2+||y||^2)+log cosh(2t x^T R y)` has norm at most
`4tB`. Two endpoint displacements at most `eta` therefore change the
log kernel by at most `8tB eta`, uniformly in orthogonal `R`.

## 2. Exact coupling of two prescribed profile systems

Fix a positive finite kernel `K` and let `L=osc(log K)`. At vertex `i`
let `k_i,l_i` be two color-count vectors, both summing to `n-1`.
Let `Z(k),Z(l)` denote their independent uniform-row-permutation
partition functions with the same kernel. Then

```math
 |\log Z(k)-\log Z(l)|
 \le \frac L2\sum_i\|k_i-l_i\|_1.                       (1)
```

To prove it, form at each row a deterministic list of `n-1` pairs of
colors. Include `min(k_ic,l_ic)` copies of `(c,c)` for each color, and
match the remaining surplus first colors to deficit second colors
arbitrarily. Exactly `||k_i-l_i||_1/2` pairs have unequal colors.
Apply a uniform random permutation to this PAIR list and use its first
and second coordinates as the coupled row arrangements. Each marginal
is exactly the appropriate uniform multiset arrangement. Take these
couplings independently over rows.

The resulting two whole arrays differ at exactly
`s=sum_i ||k_i-l_i||_1/2` directed incidences. At most `s` undirected
kernel factors change, and each log factor changes by at most `L`.
Their products compare pointwise by `e^(-Ls)` and `e^(Ls)`.
Expectation proves (1). No entropy or preimage bound is needed here.

Writing the profiles as `p_i=k_i/(n-1)` and `q_i=l_i/(n-1)`, (1) gives

```math
 n^{-2}|\log Z(p)-\log Z(q)|
 \le L\frac{n-1}{n}\frac1n\sum_i\|p_i-q_i\|_{\rm TV}.    (2)
```

Thus arbitrarily varying row-profile families can first be rounded to a
FIXED finite net of the color simplex, at uniformly vanishing pressure
cost. Integer rounding adds only `O_D(n)` changed incidences. The limiting
frequencies of the finitely many rounded profiles may be extracted on a
subsequence, after which the director's theorem applies. This removes
the finite-number-of-profiles restriction as a TECHNICAL compactness
issue. It does not supply the exponential number of actual spin tuples
having a given joint profile.

Equivalently the pressure admits a unique continuous extension from
finite profile distributions to arbitrary probability distributions on
the finite-color simplex. Its modulus is at most `L` times Wasserstein-1
distance with profile total variation as the ground metric. To see this,
match finitely many profile classes with a coupling, realize rational
class proportions, and apply (2) after relabeling vertices. Density of
finite-support distributions and completeness give the extension.

## 3. What is and is not retained about the seed

In the grouped actual weave the reflection `R_A` still appears inside
every positive edge kernel. Unlike a local Hilbert-norm estimate, the
transport theorem does not erase it algebraically. However the theorem
is conditional on the full outgoing JOINT physical profile. The exact
spin-summed interface remains the source enumerator `N_H(p)` in the
companion note. Neither the entropy objective nor (1) estimates that
enumerator from the original Boolean cap `Q(A)`.

In particular, unrestricted coordinate-sign gauge invariance cannot be
used to assert that the ACTUAL joint Hadamard source is seed blind:
the required entrywise gauge need not preserve a Boolean spectrum.
The separate gauge note gives an explicit four-dimensional failure.
Conversely, merely retaining `R_A` symbolically does not prove a
cap-sensitive upper bound. A useful seed-transfer statement must still
bound this full positive transport-plus-source-enumerator expression by
an expression involving the actual seed cap, with controlled all-order
loss.
