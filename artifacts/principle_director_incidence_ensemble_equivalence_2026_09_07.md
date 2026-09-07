# Entropy-universal balanced incidence factorization of actual signings

2026-09-07. New consequence of the independently reconstructed signed
Eulerian count. Complete proof below; separate review requested. This is
an all-order representation and a uniform measure comparison, not a proof
that the represented optimum or its disorder pressure converges.

## 1. The exact ensemble comparison

Let G be a connected loopless graph with positive even degrees 2k_v and e
edges. At each vertex independently choose a uniformly random balanced word
on its incident halfedges. Set S_vw=u_(v,vw)u_(w,vw). Write P_G for the law
of S and U_G^* for uniform edge signs conditioned on

    product_e S_e=(-1)^e.

The signed Eulerian count proves that the support of P_G is EXACTLY this
parity class. More quantitatively, at every point of that support,

    1/(2 product_v k_v) <= dP_G/dU_G^*
                       <= (1/2) product_v (2k_v+1).          (1)

Indeed put D=product_v binom(2k_v,k_v), the number of balanced incidence
arrays. The signed count gives C_G(S)>=D/(2^e product k_v), while trivially
C_G(S)<=2^e (choose one incidence sign per edge). Since P_G(S)=C_G(S)/D
and U_G^*(S)=2^{-(e-1)}, the lower ratio is immediate. For the upper use
binom(2k,k)>=2^(2k)/(2k+1), as the central binomial term is the largest
of 2k+1 terms summing to 2^(2k).

Consequently, for EVERY nonnegative function f of the ENTIRE signing,

    E_U* f/(2 product k_v) <= E_P f
                            <= (product (2k_v+1))/2 E_U* f. (2)

This has no local-statistic, continuity, concentration, or positive-probability
cutoff assumption. In particular it retains extremely rare events down to a
single signing. When e is much larger than sum_v log(2k_v+1), the normalized
log probabilities and partition functions coincide at speed e.

The comparison is logarithmic, not a claim of small total variation. The
row variables are independent before taking products. The resulting edge
variables generally are not independent.

## 2. An all-order full-sign representation

For odd n>=3 take G_n=K_n. For even n>=4 take G_n=K_n minus a fixed perfect
matching. Its degree d_n is even and G_n is connected. Independently fill
the omitted matching edges with fair signs. This produces a law P_n of
actual hollow symmetric full sign matrices at EVERY such order.

Its support consists of all full signings satisfying the single G_n parity
condition. Thus every full signing can be put in the support by changing
at most one G_n edge. Since one edge changes Q by at most 2,

    M_n <= min_(A in support P_n) Q(A) <= M_n+2.             (3)

The support parameterization uses n independent balanced words (and an
independent matching). This is a genuine exact row constraint with an
entropy-rich realization theorem, but it is not a small-information state:
the words still carry order n^2 bits in total.

## 3. Negative-disorder pressure is preserved at the critical scale

Let U_n be uniform on all hollow signings, let Q be the original absolute
Boolean cap, and fix lambda>0. Define

    Z_n(lambda)=E_U_n exp[-lambda sqrt(n) Q(A)],
    Z_n^inc(lambda)=E_P_n exp[-lambda sqrt(n) Q(A)].

Flipping one fixed G_n edge bijects its two parity classes. Since
|Q(A')-Q(A)|<=2, conditioning the uniform law on either class changes
the logarithm of this positive partition function by at most
2lambda sqrt(n). Applying (2) and k_v=d_n/2 then gives, for n>=4,

    |log Z_n^inc(lambda)-log Z_n(lambda)|
       <= n log(n-1)+log 2+2lambda sqrt(n).                 (4)

The deliberately loose bound covers both ratio directions and even orders.
At fixed lambda its error is o(n^2). The comparison holds for lambda
depending on n too, with the displayed extra cost retained.

For m_n=M_n/n^(3/2), the elementary minimum/atom sandwich is

    m_n <= -log Z_n(lambda)/(lambda n^2)
         <= m_n + binom(n,2)log2/(lambda n^2).              (5)

Therefore convergence of -n^(-2)log Z_n^inc(lambda) for every fixed lambda>0
would imply convergence of m_n: use (4), then (5), and let lambda increase
only AFTER the n limsup/liminf comparison. No claimed pressure limit is
being inferred merely from the independent-row representation.

For sublevel counts the direct event comparison (2) similarly preserves
their speed-n^2 rate within the parity class. Passing between parity classes
requires a cap threshold shift 2/n^(3/2); it must not be discarded without
checking the intended continuity or open/closed-set statement.

## 4. What this principle explains, and what it does not

Independent row constraints can look globally restrictive yet represent
every underlying edge signing, up to one parity bit, with only O(n log n)
log-density distortion. Thus replacing the full orthogonal Boolean row code
in a weave by merely balanced types can lose the distinguishing seed
information: it has enlarged the feasible incidence matching set almost
to the entire disorder space. This is a proved scope statement, not a claim
that every representation or every growing hierarchy fails.

Conversely, the theorem offers a concrete row-independent representation
for a rare-event pressure problem without sacrificing atom-scale events.
A useful next step would have to exploit that row structure to prove a
pressure comparison or limit. The remaining nested Boolean maximum is still
present; this representation is NOT yet certified simpler than the original
optimization, and the original convergence question remains open.

Dependencies: `principle_director_signed_eulerian_entropy_2026_09_07.md`,
its independent transition-count audit, and the ordinary Schrijver
Eulerian-orientation inequality cited there. All other steps are finite
counting, a one-edge bijection, and the minimum/atom sandwich.
