# A common-temperature mixed-profile consequence of the balanced certificate

2026-09-07. **Positive theorem; director independent audit PASS after the
entropy inequality correction recorded below.** Uses the actual stratified marked ensemble
and its proved one-row/Finner compiler. This is not endpoint interpolation.

Fix an actual ensemble of the canonical stratified construction at p=k/L,
and a fixed temperature t>0 for which the BALANCED one-row exponent is
at most A. Assume that for a chosen b<1/2,

    d=1-2b sqrt(p),       sigma=-(A+t d)>0.              (1)

The certificate A must already include the exact selector entropy fee
and valid same-law inheritance. The safe Sylvester construction is one
such setting when its numerical certificate is supplied. Fix any
0<xi<sigma, and choose its finite child depths so the balanced exponent
error is less than xi/2.

For a physical word, let B be its exactly balanced fibres. All other
fibres may have arbitrary means, with minority fractions r_i in[0,1/2].
Define the explicit nonnegative entropy cost

    R(r)=p h(r)+4t d r(1-r).                            (2)

Then the SAME actual construction, with asymptotically high probability,
satisfies uniformly over every word obeying

    sum_(i not in B) R(r_i) <= (sigma-xi)|B|,            (3)

the variance response

    |H_bulk(x)| <= b sqrt(N) sum_i ||P_i x_i||^2
                                                +o(N^(3/2)). (4)

Every fibre may be constant or have either majority polarity. In
particular there is no lower bound on positive r_i and no restriction
on the number of distinct density scales in (3).

## Proof and quantifiers

For ANY bias r, the finite-depth row contribution has the exact crude
upper exponent p h(r), since every Gaussian orbital norm is at most1.
Equivalently the selector-averaged count of q-coordinate sign words of
that bias is binom(q,qr). In the finite group formula one can use

    H(T,x)-log binom(L,k) <= k h(r),

by conditioning on T and subadditivity/concavity of binary entropy.
There is no approximate-envelope error on these nonbalanced rows.
Their Markov term is t d times variance4r(1-r), giving (2).
Each balanced row contributes at most -sigma+xi/2. Therefore (3) gives
a strictly negative common-temperature row sum at most -xi |B|/2.
All edge variables use the SAME t; the usual graph Finner compiler
therefore applies without separately estimating mixed cross terms.

For fixed macroscopic total variance, |B| must be macroscopic too.
Indeed -log r>=1-r and -log(1-r)>=r give h(r)>=2r(1-r),
hence4r(1-r)<=2h(r), and (2)--(3) imply

    sum_(i not in B)4r_i(1-r_i)
                         <=2(sigma-xi)|B|/p.

Thus a total normalized row variance at least delta m forces |B|>=c delta m.
The exponential Finner margin dominates the exp(O(m log m)) choices of
integer means, balanced sets, and polarities. For smaller total variance,
the deterministic bulk operator bound gives O(delta N^(3/2)) energy.
Let order grow first and then delta decrease to zero. This supplies
the uniform additive-o formulation of (4). Balanced-column repair has
uniform o(N^(3/2)) cost and exact annihilation of fibre constants.

## An explicit useful subset

Suppose every nonbalanced fibre has r_i<=r0<=1/2 and a fraction alpha of
all fibres is balanced. Condition (3) holds whenever

    (1-alpha)[p h(r0)+4t d r0] <= (sigma-xi) alpha.       (5)

Thus for small r0 a balanced fraction of order r0 log(1/r0) suffices;
the remaining rows may contain arbitrarily many finer minority scales.
This is a partial mixed-profile result derived from one actual common
temperature. It does not handle low balanced mass failing (3), general
intermediate biases, or prove the full variance inequality.

## Safe numerical specialization

Use L32,k31,p31/32,t97/20 and the director's directed certificate
E_t(nu_p)<=-4/5, with the exact stratification fee. Its balanced cap
coefficient is below .497761196. Set b=499/1000 and xi=1/200.
Outward interval arithmetic gives

    sigma > .011827156412724424,

and verifies (5) for alpha=1/400 and r0=1/1000000, with slack greater
than .00000240838643247892. Thus the same actual ensemble satisfies (4)
at coefficient .499 for every word having at least0.25 percent exactly
balanced fibres and all other fibres within minority fraction10^(-6)
of their own constant state. Their majority signs and tiny scales are
unrestricted. This is explicitly a constraint-set cap theorem.

The65-digit outward-rounded calculation and its source certificate are
recorded by `../computations/principle_construct_2026_09_07_mixed_cone_constants.py`
and the corresponding results JSON. The earlier proof line h(r)>=2r
was incorrect with natural logarithms; the director caught it. The
direct inequality h(r)>=2r(1-r) above proves exactly the needed bound.
