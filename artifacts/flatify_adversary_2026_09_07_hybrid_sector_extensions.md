# Reserved-row removal and all-order hybrid exact-support sectors

Scope update: the fixed-center energy-window application below is now
dominated by `flatify_adversary_2026_09_07_dephased_hadamard_energy_window_audit.md`.
These narrower extensions remain valid but are no longer needed for that
application.

2026-09-07. This extends
`flatify_director_hybrid_exact_support_bridge_2026_09_07.md`, whose exact
formulas were independently read and reconstructed. The formulas for
a1,a2 and D1 agree with the independent pointwise GP/PP derivation in
`flatify_adversary_2026_09_07_mixed_pinned_sector_audit.md`. In particular
a2=(1/2)log(4/3), and the type-weighted norm error is bounded directly by
||Lambda||op*||u-g||*(||u||+||g||), not by assuming unsigned cancellation
of total squared-norm errors.

These are sector extensions only. No child edges or child energies are
modified in their statements or proofs.

## 1. Removing the four reserved-row profile conditions

Keep the reserved pin destinations and the construction exactly as in the
source. Before drawing bridge permutations and signs, choose in each
reserved fibre a fixed physical target word which is G-typical at noise
level 1/10 about its prescribed center. Such words exist for all large k:
independent biased signs in a normalized Hadamard basis have the required
two-dimensional pair marginals and four-dimensional two-pair convergence,
while their pinned mean and partner errors have expected squared size O(1).
The argument is uniform over normalized Hadamard bases and centers.

For ANY candidate pair, replace its spins on the two reserved fibres per
side by these fixed prototype words. Leave every other fibre unchanged.
The modified pair meets the reserved-G conditions; the average error on
the other fibres is unchanged up to the harmless normalization, and the
prototype error tends to zero. Apply the source theorem with a vanishing
tolerance enlarged to include that prototype error.

Each side changes at most 2k spins, so the norm of its change is at most
2sqrt(2k). Since ||C||op=sqrt(n), the uniform bridge comparison costs at
most 4n sqrt(2k)=O(n^(5/4)). Therefore one actual bridge satisfies

    |x^T C y| <= (.72+o(1))*n^(3/2)

for all original pairs whose NONRESERVED fibres admit the G/arbitrary-P
declarations, with NO profile or support requirement on the reserved
fibres themselves. The exact center bridge remains zero. The overwrite
is solely an auxiliary comparison of bridge bilinear forms. Use the
ORIGINAL values H_A(x),H_D(y) in the parent identity; no bound on the
effect of overwriting those spins inside a child is invoked or needed.
The fixed strict saving at c>=.47 absorbs the O(n^(5/4)) bridge loss.

## 2. Hadamard rather than Walsh pairs suffice

For an arbitrary normalized Hadamard H_h of order h, use
H_k=H2 tensor H_h, k=2h. Physical coordinates have a first H2 bit and an
h-index. The fixed P constraint is equality of the relative physical
word at the two values of that first bit, for each h-index. Pair the two
frequency columns having the same H_h column and different H2 bits,
listing the constant H2 bit first.

P's inactive frequency is exactly zero. Its active coefficients are
arbitrary real numbers, which the hybrid proof treats pointwise. The
center coefficient is sqrt(k) in the first column and zero elsewhere,
since H_h is normalized. G-typical prototypes and profile approximation
follow from the same bounded-coefficient fixed-dimensional replacement
argument for ANY Hadamard, not from Walsh algebra. Thus every step of the
hybrid theorem, including signs, pair permutations, pins, exact center
zero, and entropy k/2 for P words, remains valid for this H2 factor.
Other Walsh translations are not asserted at these general orders.

## 3. Dense orders and arbitrary original sizes

Use Hadamard orders k=2^a12^b with a>=1. Their multiplicative gaps tend
to one by the elementary irrational-rotation argument, and each has the
required H2 factor. For every original size n choose k>=sqrt(2n) with
k/sqrt(2n)->1, put m=k/2 and n_+=mk. Then r=n_+-n=o(n).

Place the original n coordinates first in the full fibre ordering. Call
all incomplete or deleted fibres, plus the four reserved fibres across
the two sides, exceptional. Per side their full coordinate count L is
at most r+3k, hence o(n). On every fully retained nonexceptional fibre
define G and P using the actual center and H2 tensor H_h as above. Average
their profile errors with full-size denominator n_+ (equivalently any
asymptotically equal retained-size denominator). Require that average to
vanish. Impose NO support, mean, or profile condition on original spins
in the exceptional fibres.

Pad the centers arbitrarily to n_+. Build the full hybrid bridge using
fixed G prototypes on every exceptional fibre when comparing target
vectors, and declare all these fibres G. Its pin reservoir is among the
fixed reserved fibres. The source theorem already unions over every
allowed declaration on the other fibres. The prototypes need not match
the original candidates on the exceptional coordinates.

Let C be the original n-by-n restriction of the resulting full bridge
C_+. Its entries are signs. Embed the original target x by zeros on
deleted coordinates, calling it x^0; let x~ be its full modified sign
word. Then ||x^0-x~||<=2sqrt(L), and the same holds on the other side.
The exact operator bound gives, uniformly in all original targets,

    |x^T C y-x~^T C_+ y~| <= 4 n_+ sqrt(L)=o(n^(3/2)).

For the center, comparison only deletes r coordinates and costs at most
2 n_+ sqrt(r), so its restricted normalized bridge value tends to zero.
Consequently, at EVERY sufficiently large original order n there exists
an actual full sign bridge satisfying

    x0^T C y0=o(n^(3/2)),
    |x^T C y| <= (.72+o(1))*n^(3/2)

simultaneously throughout the inherited G/arbitrary-P sector on all fully
retained nonexceptional fibres. There is no condition on the O(r+k)
exceptional physical coordinates. Their effect is paid uniformly by the
bridge norm, rather than by an unjustified entropy deletion.

## 4. Original-value statement and limitations

Both actual optimal children remain at the ORIGINAL order n with all
their edges unchanged. For original candidates satisfying the usual
energy-window condition

    |H_A(x)|,|H_D(y)| <= (.64*c+o(1))*n^(3/2),

the parent expression is at most (.72+1.28*c+o(1))*n^(3/2), strictly below
2sqrt(2)c*n^(3/2) for c>=.47. This conclusion uses no child-energy estimate
for the overwritten or padded comparison words.

The sector is defined explicitly on the fully retained Hadamard fibres;
it is not falsely called an unpadded Walsh profile. Typical iid G noise
and translation-paired noise on those fibres give examples in it, while
arbitrary P spectra are allowed outright. The theorem still does not
cover arbitrary non-paired exceptional rows, arbitrary energy windows,
different noise levels, or a free choice of pairing after the target is
seen. It supplies no whole-parent recurrence and no convergence proof.
