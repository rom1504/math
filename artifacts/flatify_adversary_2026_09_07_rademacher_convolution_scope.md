# Rademacher-convolution kernel tests and actual Walsh scope

2026-09-07. These distinguish the proved noisy-single-Rademacher sector
from a proposed broader kernel lemma. They do not obstruct a whole
rank-two construction.

## Exact failure of the proposed uniform pair bound

At t=2 sqrt(15), lambda=15/2, the tilted four-variable kernel is

    E exp[t sqrt(2)(AC+AD+BC-BD)-lambda(A^2+B^2+C^2+D^2)].

The noisy-single-Rademacher theorem proves this is below 1/15 for its
reference family. It does NOT extend to every normalized Rademacher sum.
For A,B,C,D independently distributed as (R1+R2)/sqrt(2), the event all
four vanish contributes 1/16. The events with one nonzero coordinate on
the left side and both nonzero on the right side, with favorable orientation,
have total probability 1/16 and exponent 8 sqrt(30)-45 > -2.
Therefore the kernel exceeds 1/16+1/144=5/72>1/15, using e<3.
The exact full numerical value is approximately .0817011660.

More strongly, if A,B are single Rademachers and C,D are normalized
two-Rademacher sums, the favorable events where one of C,D vanishes and
the other aligns with the matching Hadamard sum have probability 1/8.
Their exponent is 8 sqrt(15)-30>0, so this mixed kernel exceeds 1/8.
Its numerical value is approximately .3343475569. These are exact local
kernel counterexamples, not global profile counterexamples.

## Preserved bounded searches

`computations/flatify_adversary_2026_09_07_rademacher_convolution_kernel.py`
and its initial and `_extended.json` results preserve the floating tests.
They include unequal two-, three-, and four-term sums, equal sums through
32 terms, and Gaussian smoothing. Full-spin-entropy mixture optimization
gave a largest found threshold .4934790611 with weights on equal sums
of 1, 2, and 4 Rademachers. This mixture crosses Walsh-order parities and
is not an actual common-order obstruction. Restricting to the tested
even parity family (1,4,16 terms) gave .4845062464; the odd parity family
(2,8,32 terms) gave .4884461636. These are floating optimization outputs,
not certified universal upper bounds or certified global optima.

## Exact common-order realizability of the parity-compatible families

Let L=2^l. For all sufficiently large n with n=l modulo 2, put
d=(n+l)/2 and assume d>=L-1. On F2^l, use the L-1 nonempty squarefree
monomials as coordinates of a map T into F2^d, padding with zeros.
The first l coordinates are the original variables. Its image is a
graph over F2^l. The L graph points are affinely independent: adjoining
the constant coordinate gives the invertible evaluation matrix of all
squarefree monomials on the Boolean cube.

Translate this graph in its last d-l coordinates to partition F2^d into
fibres of a map pi:F2^d -> F2^(d-l). For x in F2^(d-l), y in F2^d, and
any Boolean g(y), define f(x,y)=x dot pi(y)+g(y). Its Walsh coefficient
at frequency (u,v) is

    2^(d-l) sum_{y:pi(y)=u} (-1)^(g(y)+v dot y).

After division by sqrt(2^n), this is 1/sqrt(L) times a sum of L signs.
As v varies uniformly, affine independence makes these signs uniform
modulo their common global sign. Thus their ABSOLUTE sum has exactly
the distribution of |(R1+...+RL)/sqrt(L)|. The same holds in every fibre
u, so the entire absolute Walsh profile has exactly that distribution.

Consequently all fixed powers-of-four term counts coexist at sufficiently
large even n, and all twice-powers-of-four counts coexist at sufficiently
large odd n. This gives actual growing-order reference profiles, not
generic surrogate measures. Varying g alone gives only 2^(2^d) rows,
which is subexponential in 2^n for fixed L; no positive physical entropy
is inferred from that count. Independent bit noise can be applied as in
the separate noisy-bent/profile argument, with its explicitly paid noise
entropy and Gaussian smoothing, but this does not fix the failed local
uniform-kernel bound.
