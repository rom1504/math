# Director audit: actual augmentation and forced soft collective excursions

2026-09-07. **PASS.** Independently reconstructed from the complete
`principle_invent_2026_09_07_balanced_valley_augmentation.md` and
`principle_construct_2026_09_07_critical_barrier_augmentation.md`.

The primary Gram--Schmidt Walk Theorem1.4 has all-directions MGF parameter
sqrt(40), not variance40 per coordinate only. The supplied factorization
G=L V, column norms at most1 and center row norms at mostgamma, is inserted
as actual vectors (aV_i,b e_i). These have norm at most1 when a^2+b^2=1.
Testing with (L_g/a,(x-g)/b) gives the exact old-spin linear form v dot x.
There is no missing correlation between feature and identity components:
they belong to orthogonal summands of the test space.

With gamma^2/n=g_n->0, the choice a^2=sqrt(g_n), b^2=1-sqrt(g_n) gives
variance proxy40n sqrt(g_n)+160d/(1-sqrt(g_n)). If a denominator vanishes
at an early finite order, omit that order; asymptotic hypotheses ensure
the displayed positive denominators eventually. A nonempty sign center
family has gamma>=1, so a is not zero.

VC dimension at mostgamma^2 follows by shattering and a random-sign sum
of the unit feature columns. Sauer's elementary recursion then gives
log|G|=o(n). Union over every distance shell and all2^q new spin words
therefore yields the uniform bridge envelope

    sqrt(320 epsilon r[h(r)+epsilon log2])+o_n(1).

Columns are actual full signs from independent GS walks. All old edges
are retained. The bridge is paid by its absolute value, so the proof
respects the child-reversal identity exactly. A random full-sign new child
has cap at mostq^(3/2); this supplies a real parent at every n+floor(epsilon n).

For a barrier f(r), the actual increment is at mostepsilon^(3/2) plus the
supremum of envelope minusf. If f(r)>=kappa' r^2 log(e/r) near zero,
write w=r sqrt(log(e/r)). The first envelope component is at most
sqrt(320epsilon)w, whose excess overkappa'w^2 is at most80epsilon/kappa'.
The other component contributes at mostO(epsilon sqrt(delta)) on r<=delta.
A positive barrier away from zero pays every r>=delta. The limit order
is n first, thenepsilon, thendelta andkappa'; it gives the claimed slope80/kappa.

If seeds realize the original liminf c_inf, a slope strictly below1.5c_inf
would lower the normalized cap at a fixed slightly larger order and violate
the definition of liminf. Hence kappa<=160/(3c_inf) under the stated
low-complexity-center and macroscopic-isolation hypotheses. This is not an
assumption that the whole sequence M_n/n^(3/2) converges.

The double-limit stiffness variant has the same correct quantifiers:
choosing a coefficient strictly below the alleged liminf stiffness supplies
an eventual uniform lower bound for each fixed small delta. There is no
illegitimate uniformity across growing n-dependent delta.

If G contains ALL exact absolute grounds, every noncenter excitation has
energy deficit at least2, since full-sign quadratic values have a common
parity. A soft excitation with bounded deficit/[r^2 log(e/r)] therefore
has d^2 log(en/d)>=c sqrt(n), hence d>=c'n^(1/4)/sqrt(log(en)). This lower
size is conditional on that branch and on including all exact grounds.

The hypotheses are nonempty for actual signs. In the affine-plane example
on F_q^2, choosing(q+3)/2 directions positive gives eigenvalues2q-2 on
the constant vector andq-2 or-q-2 elsewhere. Forq>=9 the two absolute
grounds are exactly+-1, and their deficit is bounded below by a constant
times normalized Hamming distance. This is a bounded-cap but SUBOPTIMAL
family, not evidence that actual minima have the same geometry.

The important gain over a conditional reformulation is the explicit
old-edge-preserving sign extension and its quantitative slope. The missing
convergence obligation has NOT disappeared: near-minimizers are forced to
avoid the easy rigid/low-complexity case, and the construction does not yet
control their remaining soft, high-complexity landscape.
