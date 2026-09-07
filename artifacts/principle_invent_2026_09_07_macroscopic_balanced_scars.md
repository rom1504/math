# Macroscopic balanced scars: a strict-face certificate is not a variance bound

2026-09-07. **Actual balanced-block construction, with a fixed positive
minority density.** Unlike a microscopic scar, this falsifier survives
an additive o(N^(3/2)) error in a proposed variance-linear response.
It does not refute the existence of a favorable unscarred ensemble.

## 1. Hypotheses and exact switch operation

Partition N=mq physical vertices into m fibres of even size q. Let B be
symmetric, zero within each fibre, and FULL signs between distinct
fibres. Suppose every cross-fibre block is a rank-one sign matrix with
exactly balanced rows and columns. Thus B1=0. Assume its balanced-face cap

    Q_bal(B)=max_(x Boolean, sum_(a in fibre i) x_a=0 all i)
                       |x^T B x|/2
                <=(b_0+o(1)) N^(3/2),             (1)

where b_0<1/2. Such blocks are exactly the output of the balanced-column
repair/compiler: signs are outer products of balanced transform columns.

Fix an integer 1<=ell<=q/4. In each fibre choose ell vertices and label
them 1,...,ell; denote their set by S_i. In each cross-fibre block force
every equal-label entry to be positive, using the following exact
balanced switch whenever that entry is initially negative.

If its row and column are a,c, choose a companion row u outside S_i and
a companion column v outside S_j. In the original rank-one block choose
u with row sign opposite a and v with column sign opposite c. Also choose
both companions unused by previous switches in this block. There are at
least q/2-ell-(ell-1)>0 choices on each side. Selected target rows and
columns are distinct and never serve as companions; new companions have
not occurred in any previous switch. Thus the four entries have not been
altered and are exactly

    [B_ac B_av; B_uc B_uv]=[-1 +1; +1 -1].

Flip all four signs. Every affected row and column retains its sum, the
target becomes positive, and no other entry with BOTH endpoints selected
is changed. The transpose block is changed correspondingly, preserving
the symmetry of the whole physical matrix. At most ell switches occur
per cross-fibre block. The resulting B' still has FULL signs between
fibres and exact row/column balance, although its blocks need no longer
be rank one.

## 2. The strict balanced-face cap survives small fixed density

For a single switch and Boolean endpoint words, its half-Hamiltonian
increment is

    2(x_a-x_u)(y_c-y_v),

whose absolute value is at most8. Therefore the safe deterministic bound
is

    |Q_bal(B')-Q_bal(B)|<=8 ell binom(m,2).         (2)

There is no assumption that only half of the targets were negative.
If q/m->p>0 and ell/q->r in (0,1/4), (2) gives

    limsup Q_bal(B')/N^(3/2)<=b_0+4r/sqrt(p).      (3)

In particular choose any fixed

    0<r<min {1/4, (1/2-b_0) sqrt(p)/4}.

Then the balanced-face cap of B' remains strictly below1/2.

## 3. A macroscopic biased word with coefficient above one half

Choose the ordered selected sets independently and uniformly in the
different fibres. Let z be the indicator of their union, and put
x=1-2z, so every fibre has the SAME mean a=1-2ell/q.
Write H_B(z)=z^T Bz/2. Exact block balance gives

    E H_B(z)=0.

For each equal-label target, its original expected sign is also zero:
the two selected locations are independent uniform physical coordinates,
and the corresponding block has total sum zero. Since the only changed
selected-selected entry is the target, the switch procedure gives exactly

    H_(B')(z)=H_B(z)+sum_(equal-label targets)(1-B_target).

Consequently

    E H_(B')(z)=ell binom(m,2).

Some ordered selection therefore yields H_(B')(z) at least that value.
Because B'1=0, the corresponding ACTUAL Boolean word obeys

    H_(B')(x)=4 H_(B')(z)>=4 ell binom(m,2).        (4)

Its total centered variance is exactly

    sum_i ||P_i x_i||^2=4m ell(1-ell/q).          (5)

Combining (4)--(5), for q/m->p and ell/q->r,

    liminf H_(B')(x)/(sqrt(N) sum_i||P_i x_i||^2)
                 >=1/[2 sqrt(p)(1-r)].           (6)

For p<1 this is strictly above1/2. The variance density in (5) tends to
4r(1-r)>0, so no additive o(N^(3/2)) term can restore a strict-subhalf
variance-linear bound for this sequence.

## 4. Scope and consequence

Equations (3) and (6) coexist in the same actual balanced signing: strict
subhalf control of EVERY fibre-balanced Boolean word and zero energy on
the constant-fibre mode do not imply a strict-subhalf response linear in
the centered variance. The example has only O(r N^(3/2)) sign edits and
preserves every cross-block row and column sum exactly.

Within-fibre zeros describe the BULK operator, as in the balanced compiler.
They can be filled by J_q-I_q to give an actual full hollow signing; this
adds only -N/2 to balanced words. The response assertion (6) remains a
statement about its explicitly separated balanced bulk, not a false claim
that the entire filled signing has small unrestricted cap.

This does not show that the original randomized repaired weave necessarily
contains such scars, nor that a favorable choice of all its signs cannot
avoid them. It proves that the two endpoint certificates and exact balance
alone are insufficient; an additional quantitative mixed-profile property
of the chosen construction really is needed.
