# Exponential cost of good principal selectors in actual bounded-cap parents

Date: 2026-09-06. Status: proved. This sharpens a typical-selector
obstruction into a quantitative bound on the number of good selectors.
It does not exclude their existence or settle convergence of M_n/n^(3/2).

Conventions: A is symmetric and hollow, Q(A)=max_x |x^T A x|/2.
All selector laws below are on subsets, not on new independent signings.

## 1. A finite lower-tail inequality for principal caps

Let A be any real symmetric hollow D-by-D matrix, L=||A||op>0, and
let T be a uniform n-subset of [D], 1<=n<D. If

    Pr{Q(A_T)>=b} >= 1/2,

then for every a<b,

```math
\Pr\{Q(A_T)\le a\}
\le 2\exp\left\{-\frac{(b-a)^2}{576L^2n}\right\}.       (1)
```

The imported result is the multislice convex-distance inequality of
[Sambale--Sinulis, arXiv:2010.16289v1, Proposition 1.7, pp. 5--6](https://arxiv.org/pdf/2010.16289v1):
for a subset E of the uniform multislice,

    Pr(E) E exp(d_T(u,E)^2/144) <= 1,

where d_T is the supremum, over Euclidean-unit coordinate weights, of
the minimum weighted coordinate-mismatch distance to E. The ordinary
binary slice is the case with alphabet {0,1} and multiplicities D-n,n.
There is no requirement that E be invariant under coordinate permutations.

Here is the exact mapping, including the convexity issue. For signs x and
sigma in {+1,-1}, let

    G_(sigma,x)=sigma diag(x) A diag(x)+L I,
    f_(sigma,x)(u)=u^T G_(sigma,x)u/2,
    F(u)=max_(sigma,x) f_(sigma,x)(u).

Every G is positive semidefinite and has norm at most 2L. Consequently
F is convex. On the binary n-slice,

    F(1_T)=Q(A_T)+Ln/2.

Fix a high-cap indicator u and an active pair (sigma,x). The vector
g=G_(sigma,x)u is a subgradient of F at u and has norm at most
2L sqrt(n). For every v in the low-cap event E={Q(A_T)<=a},

    b-a <= F(u)-F(v) <= g dot (u-v)
          <= sum_i |g_i| 1_(u_i != v_i).

The middle inequality is the supporting-plane inequality, not a claim
that the Boolean maximum itself is convex in unrestricted weights.
After normalizing g, this proves d_T(u,E)>=(b-a)/(2L sqrt(n)).
The imported exponential expectation is therefore at least
(1/2)exp((b-a)^2/(576L^2n)), proving (1).

No globally Lipschitz bound for a quadratic function is assumed.
Only its active subgradient at a slice point is used. The shift Ln/2
is constant on this slice, which is essential.

## 2. Actual bounded-cap parents: uniformly small fixed retentions

For every fixed C0>0 and 0<c<2/pi there exist rho>0, kappa>0 and n0
depending only on C0,c such that EVERY hollow sign parent with
Q(A_D)<=C0 D^(3/2), and EVERY n0<=n<=rho D, obeys

```math
\boxed{\quad
\Pr\{Q((A_D)_T)\le c n^{3/2}\}
\le 3\exp\{-\kappa n^2/D\}.\quad}                       (2)
```

In particular, for any fixed p in (0,rho], n=floor(pD), the fraction
of such good selectors is exponentially small in n. The statement is
uniform over actual optimizing parents and all sufficiently large
orders; it does not assume those parents have bounded operator norm.

### Uniform positive-probability input

The fixed-degree local-moment proof in
`transfer_seed_sparse_bounded_cap_parent_theorem_2026_09_06.md`
has the following immediate bounded-operator form. For every C and
c1<2/pi there are rho1 in (0,1/2],m0 such that, whenever B_d is a
hollow signing, ||B_d||op<=C sqrt(d), and m0<=m<=rho1 d,

    Pr{Q((B_d)_S)>=c1 m^(3/2)} >= 3/4.                    (3)

To see the quantifiers directly, choose a FIXED Chebyshev degree q with
2q/[pi(q+2)]>c1. The mean-square errors of its normalized trace and
trace-energy, and its expected mean-square diagonal error, are bounded
by constants times rho1+1/m+d^-1/2. The trace of R_q^2/m
has a bounded expectation. First choose a large fixed bound on this last
quantity, then small fixed tolerances for the other errors; Markov and
Cauchy--Schwarz give failure probability below 1/4 by choosing rho1
small and m0 large. These constants depend only on C,c1. There is no
exchange of a growing polynomial degree with the order limit.

The later stronger K/n estimate in
`transfer_adversary_fixed_retention_random_loss_2026_09_06.md`
also implies (3), but is not needed for the exponential argument.

### Spectral core and its sampled size

Choose c<c1<2/pi and fixed epsilon in (0,1/4) with

    g=c1(1-2epsilon)^(3/2)-c > 0.

The proved simultaneous diagonal-majorant deletion gives a core R of
size d>=(1-epsilon)D with

    ||A_R||op <= L0 sqrt(D),
    L0=4 K_G C0/epsilon,  K_G=pi/(2 asinh(1)).

Thus (3) applies to this core with C=L0/sqrt(1-epsilon).
Take rho<=(1-epsilon)rho1 and n0 large enough that
(1-2epsilon)n0>=m0.

For uniform T let m=|T intersect R|. Its expectation is at least
(1-epsilon)n and

    Pr{m<(1-2epsilon)n} <= exp(-2epsilon^2 n).             (4)

For completeness, expose an ordered sample without replacement. The
Doob martingale for the total number of core hits has, at reveal k,
conditional range length (D-n)/(D-k)<=1. The bounded-interval
exponential estimate gives E exp(t(m-E m))<=exp(nt^2/8), hence (4).
This avoids an independence assumption on the sampled indicators.

Conditional on m, T intersect R is uniform in R. For m above the
threshold, (3) and then (1), with b=c1 m^(3/2), a=c n^(3/2), give

    Pr{Q(A_(T intersect R))<=c n^(3/2) | m}
      <=2 exp(-g^2 n^2/(576 L0^2 D)).                     (5)

Indeed b-a>=g n^(3/2), the sampled size in (1) is m<=n,
and the core operator norm squared is at most L0^2 D.
Principal monotonicity gives Q(A_T)>=Q(A_(T intersect R)). Combining
(4)--(5) proves (2) with

    kappa=min(2epsilon^2, g^2/(576 L0^2)).

The constants are not optimized. When n is much smaller than sqrt(D),
(2) is vacuous. The simultaneous sparse theorem still gives probability
tending to zero in that range; the separate K/n estimate is also useful.

## 3. Counting, sampling, and relative-entropy consequences

Let G be the set of n-subsets meeting the cap target and U the uniform
law on all such subsets. Equation (2) says

    |G| <= 3 exp(-kappa n^2/D) binom(D,n).

For any alternative selector law nu giving G probability gamma,
data processing of relative entropy under the map T->1_(T in G) gives

```math
D(\nu\Vert U)
\ge \gamma\kappa n^2/D-\gamma\log3-h(\gamma).             (6)
```

If U(G)=0 and gamma>0, the divergence is infinite. Otherwise use the
exact binary divergence and discard the nonnegative term
-(1-gamma)log(1-U(G)). For fixed positive retention and fixed gamma>0,
successful selector laws therefore require an extensive entropy change.
This is a theorem about selections of ACTUAL good parents, not an
assumption about a surrogate Gibbs law.

A union bound also limits success from K samples whose individual
unconditional laws are uniform by 3K exp(-kappa n^2/D); independence
between those samples is unnecessary. This does not cover an arbitrary
adaptive algorithm whose samples have different conditional distributions.

## What remains possible

An exponentially small fraction can still contain many good selectors.
The theorem supplies no lower count, no algorithm to find one, and no
all-order cap-preserving realization. It therefore sharpens the failure
of ordinary sampling but is not a proof of nonconvergence or of
impossibility of a designed selector construction.

Finite algebra replay:
`computations/transfer_director_selector_convexity_2026_09_06.py`.

Independent adversarial reconstruction checked the primary theorem,
active-subgradient argument, conditional hypergeometric martingale ranges,
core normalization, fixed-parameter choices, and binary KL inequality.
The finite integer replay passed all 139264 tangent-pair tests; it is a
regression check and not a replacement for the general proof above.
