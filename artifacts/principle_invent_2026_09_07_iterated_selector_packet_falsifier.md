# Product-packet obstruction to naive iterated selector compilation

2026-09-07. **Derived finite block-certificate falsifiers; Sections1--5
independently reconstructed PASS by the construction agent, including
both distinct normalizations. Not actual signing lower bounds.** This
tests the natural tensor/Kronecker realization of
several one-hole selector stages. It does not rule out a differently
specified marked-spine compiler with intermediate randomization.

## 1. The concrete product model

Take J dephased Hadamards H_(L_j), all L_j>=4, and the total block

    H=H_(L_1) tensor ... tensor H_(L_J),
    L=product_j L_j,  k=product_j(L_j-1),  p=k/L.

Selectors may be chosen by independently retaining L_j-1 coordinates
in each factor and taking their Cartesian product. More generally the
argument applies to any fixed finite selector distribution which gives
positive probability to the product selectors used below. The full
block-word pressure counts each rare physical packet once and applies
the scalar E_t separately to the L transformed output channels, exactly
as in the stratified-selector artifact.

For a permitted nonempty packet S of size s, write C_j=(H1_S)_j.
Put a=1-2r, mix a rare occurrence of this packet into the all-plus baseline
selector law, and choose its mass theta=k r/s. Its probability cost and
selector likelihood cost outside the rare Bernoulli entropy are O(r),
because the finite selector law is fixed before r tends to zero. The
rare-information lemma therefore gives, at t=tau log(1/r), the packet's
exact leading row contribution

    [theta log(1/r)/L]
          [1-sum_j min(4tau C_j^2/k,1)].           (1)

The coefficient is a LOWER witness for the supremum defining that block
pressure. There is no assumption that the rare packet is typically drawn.

## 2. Two explicit packet types

First take a singleton. Then s=1 and all L values C_j^2 equal1. In the
unsaturated range tau<k/4, the defect improvement represented by (1),
normalized by t*4r, is

    delta_single(tau)=1-p/(4tau).                 (2)

For a strict-subhalf variance certificate one needs delta>1-sqrt(p),
so necessarily tau>sqrt(p)/4. This necessary lower bound on tau remains
valid also when tau>=k/4, trivially because k>=sqrt(p).

For each factor take the positive support of a non-DC Hadamard row. It
has L_j/2 entries and can be contained in an (L_j-1)-selector. Let S be
their Cartesian product. Then

    d=2^J,  s=L/d,

and H1_S has exactly d nonzero coordinates, each of absolute value s.
Once tau>=p d^2/(4L), all these coordinates are saturated in (1), and
the normalized improvement is

    delta_product(tau)=p d(d-1)/(4tau L).         (3)

Hence a strict-subhalf certificate also requires

    tau<p d(d-1)/[4L(1-sqrt(p))].                (4)

When the singleton lower threshold already lies in the product packet's
saturated range, (2)--(4) can be compatible only if

    L(1-sqrt(p))<sqrt(p) d(d-1).                  (5)

Thus increasing the tensor depth or multiplying a very large new group
order need not increase the available common-temperature bandwidth.
The new sparse-frequency product packet creates an opposing constraint.

## 3. Exact two-stage example

Take J=2 and L_1=L_2=32. Then

    L=1024, k=961, p=(31/32)^2, d=4,
    sqrt(p)=31/32.

The singleton requires tau>31/128. Its threshold exceeds the saturation
threshold p*16/(4*1024)=p/256 of the product packet. But the two sides
of the necessary compatibility inequality (5) are exactly

    L(1-sqrt(p))=32,
    sqrt(p) d(d-1)=93/8=11.625.

They have the WRONG strict order. Therefore NO common logarithmic
temperature tau makes this natural two-stage block pressure certify a
strict-subhalf variance coefficient on every sufficiently rare biased
slice. This remains true although each separate one-hole stage has a
proved strict-subhalf rare-band certificate.

The mechanism is explicit: a singleton needs enough temperature to pay
one rare event; a product of half-column supports has only four spectral
atoms despite occupying one quarter of the full product group, and at
that temperature its entropy gain per unit variance is too small.

## 4. Scope

This invalidates a multiplicative-bandwidth argument for the specified
tensor block functional. It does not prove that the actual product
signing has a large biased cap, nor that every hierarchical selector
distribution admits these packets. In particular, a marked-spine-only
construction with fresh intermediate ordinary bases needs its own exact
selector/packet law. No impossibility is asserted for that unspecified
variant. The positive one-stage polynomial-band theorem remains intact.

## 5. A separate test of the natural marked-spine-only iteration

There is also an explicit test for a different, less destructive model.
At level1 apply H_(L_1) independently inside each smallest group. Send
its DC coordinates into a second L_2-way Hadamard stage, while its other
coordinates enter independent ordinary deep children. Continue only the
DC branch through J stages. The final pre-child transform F is an
orthogonal multiscale transform, NOT the full Kronecker Hadamard above.
Choose nested one-hole selectors by retaining L_j-1 subgroups at level j;
within retained subgroups apply the earlier-level choices. Product
selectors occur with positive probability in this finite model.

The same rare-packet calculation applies, with clipping sum

    W_lambda(S)=sum_b min(lambda |(F1_S)_b|^2,1),
    lambda=4tau/p.

The normalized improvement represented by a pure packet is

    delta_S(lambda)=[W_lambda(S)-1]/(lambda |S|).  (6)

For the product of positive half-column supports used above,

    s=product_j L_j/2^J,
    d_star=sum_(ell=1)^(J-1)
                  [product_(h>ell) L_h/2^(J-ell)]+2.       (7)

Here at stage ell there is one nonzero ordinary character, repeated on
the remaining half-support indices. Its coefficient magnitude is
sqrt(product_(h<=ell)L_h)/2^ell. The final DC coefficient has that same
last-stage magnitude. This proves (7) and shows that all nonzero squared
magnitudes are at least L_1/4 when every L_j>=4. Therefore for lambda>=1
all of them are clipped and W_lambda=d_star.

For a singleton packet, the squared F coefficients sum to1 and are at
most1/L_1, so W_lambda=lambda for lambda<=L_1. In the high-retention
regime 1/sqrt(p)<=L_1, strict-subhalf control requires

    lambda>1/sqrt(p),
    lambda<(d_star-1)/[s(1-sqrt(p))].             (8)

Consequently the common-temperature LOG-DENSITY bandwidth of this
particular iterated block certificate is bounded above by

    sqrt(p)(d_star-1)/[s(1-sqrt(p))].             (9)

For two stages of order32, one has d_star=18, s=256, and sqrt(p)=31/32.
The bandwidth upper bound is exactly527/256, approximately2.05859375.
Thus this natural iteration does NOT multiply the one-stage bandwidth;
in this example it makes its upper limit smaller.

If later group orders grow rapidly while the first one is large and
fixed, d_star/s is dominated by2/L_1, and p is dominated by1-1/L_1.
The upper bound (9) remains near4, not a product of factors near4.
This is a precise test of the defined nested/DC-spine selector law.
It does not cover arbitrary feedback-dependent selector laws or a
different distribution of fresh intermediate transformations.
