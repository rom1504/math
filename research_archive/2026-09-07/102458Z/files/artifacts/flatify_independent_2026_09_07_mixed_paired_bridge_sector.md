# A pinned-center bridge paying a union of paired-noise exceptions

Later simplification: EVERY output of this pinned ensemble already has
localized center fields, and therefore obeys the elementary all-spin
overlap bound in
`flatify_adversary_2026_09_07_dephased_hadamard_energy_window_audit.md`.
The kernel and entropy calculations below remain valid, but are not
needed for the fixed-center child-energy-window application.

Date: 2026-09-07. Proved actual-sign sector operation. Unlike simply
discarding nongeneric rows, this operation keeps their paired covariance
and pays their physical entropy. It controls a union of G/P type choices
with one fixed ensemble. It does not cover arbitrary translations,
arbitrary profile exceptions, or the whole original parent cube.

The independent researcher derived the pinned-center kernel; the
adversarial researcher independently reconstructed it, identified the
same-ensemble type union, and observed the stronger .72 threshold with
the same rational tilt. Their separate all-paired deterministic theorem
is in `flatify_adversary_2026_09_07_support_adapted_paired_bridge.md`.

## 1. Physical frames, pairings, and a single fixed pin layout

Let k=2^d=2m, n=mk, and fix arbitrary physical centers x0,y0. In each
fibre use the center-adapted Walsh basis diag(x0_i)H_k, or its right-side
analogue. The center's coefficient vector is sqrt(k) at frequency zero.

Choose a nonzero physical Walsh translation a_i in EVERY row in advance.
Choose a frequency b_i with a_i dot b_i=1, and pair frequencies j,j+b_i,
listing the active frequency a_i dot j=0 first. In particular the center
pair is (0,b_i). Different rows may have different fixed translations.

A physical target row can later be labeled G (unrestricted) or P (its
relative noise eta=x0_i*x_i satisfies eta(u)=eta(u+a_i)). For a P row,
every inactive frequency vanishes exactly. The construction does not
depend on these later labels.

Reserve two G rows a,b on the left and two G rows c,d on the right. Pin
the center pair of left row i to neighbor phi(i), and the center pair
of right row j to neighbor psi(j), where

    phi(i)=c for i!=b,   phi(b)=d,
    psi(j)=a for j!=c,   psi(c)=b.                        (1)

Every center pair points to a reserved G row and no edge has center pairs
at both ends. Permute all other frequency pairs uniformly among the
remaining neighbor slots. Give each non-center column an independent
fair sign; fix the sign of the frequency-zero column positive. These are
literal restricted signed column permutations of the Walsh bases.

Use the usual rank-two tile (1/2)F_i,pair H_2 G_j,pair^T. The resulting
C is a full rectangular sign matrix and C/sqrt(n)=UKV^T is orthogonal.
Because each center has support only at its pinned first coordinate,

    x0^T C y0=0                                         (2)

for EVERY output of this ensemble. There are exactly 2m distinct pinned
edges. All choices in (1), including reserved rows and translations, are
fixed before considering the G/P labels or target spins.

## 2. The pinned paired-profile sector

Set alpha=4/5 and s=3/5. In a physical target row, let u0 be its coefficient
at frequency zero, u1 its coefficient at b_i, and let rho_i be the
empirical law of the absolute values of the remaining k/2-1 ordered
coefficient pairs. The two reference laws are

    nu_G=law(|sZ1|,|sZ2|),
    nu_P=law(|sqrt(2)sZ|,0),                             (3)

with independent standard Gaussians. If the row is labeled P, its exact
physical pair-equality constraint is also required.

Define its squared profile error by

    d_i^2=(u0-alpha sqrt(k))^2/k+u1^2/k
              +[(k/2-1)/k] W_2(rho_i,nu_type)^2.         (4)

Require the average d_i^2 on each side to be at most tau^2. In particular
the reserved G rows must satisfy the G profile definition. The parameter
tau=tau_k tends to zero in the asymptotic theorem.

The covariance difference in (3) is essential. P rows have one active
Gaussian coordinate of variance 2s^2 and one identically zero coordinate;
they are not replaced by independent Gaussians of variance s^2.

## 3. Exact-law coupling, including the pinned coordinates

For a reference row, set its pinned first coordinate to alpha sqrt(k).
Its pinned partner is N(0,s^2) for G and zero for P. All other ordered
pairs are independent with the signed versions of (3). Condition their
empirical absolute-pair law to be W_2-close to (3), and condition the
pinned partner to have magnitude at most epsilon_k sqrt(k). For a suitable
epsilon_k->0 these two events have joint probability at least 1/2,
uniformly in the two types. This uses only fixed-dimensional empirical
W_2 convergence and a scalar Gaussian tail.

Optimal pair matching, followed by a uniform common pair permutation and
the allowed independent coordinate signs, gives exactly the physical
construction law with its center pair pinned. The pinned first-coordinate
error and the actual pinned partner's energy are paid explicitly in (4).
The reference partner's additional error is at most epsilon_k sqrt(k).
Thus, with eta=tau+2epsilon_k, the full physical feature vectors u,v and
their comparison vectors g,h can be coupled so that

    ||u-g||,||v-h||<=eta sqrt(n),
    ||u||=||v||=sqrt(n),
    ||g||,||h||<=(1+eta)sqrt(n).

Putting e=eta(2+eta), this implies

    |u^T K v-g^T K h|<=e n,
    | ||u||^2-||g||^2 |, | ||v||^2-||h||^2 |<=e n.      (5)

Dropping all row conditioning costs at most 2^(2m). The same coupling
works for every fixed physical target pair and type assignment; it need
not couple different target pairs simultaneously before the union bound.

## 4. Exact Gaussian kernels and the pinned mean contribution

Take t=10/3, lambda=35/18 and use the exponential

    exp[t g^T K h-lambda(||g||^2+||h||^2-2n)].           (6)

Let a=1+2lambda s^2=12/5, a_P=1+4lambda s^2=19/5. On an ordinary
non-pinned edge, exact Gaussian integration gives

    GG: (a^2-t^2 s^4)^(-1)=(108/25)^(-1),
    GP: [a(a a_P-2t^2 s^4)]^(-1/2)=(1872/125)^(-1/2),
    PP: (a_P^2-2t^2 s^4)^(-1/2)=(289/25)^(-1/2).        (7)

These follow by keeping O=H_2/sqrt(2) and the actual covariance matrices
s^2 I_2 and diag(2s^2,0). In particular the GP and PP effective noise
singular value is sqrt(2), with the appropriate zero direction retained.

Every pinned edge has one deterministic mean mu=alpha sqrt(k) at one
end and a G row at the other. Rotating the G pair by O splits the mean
interaction from the first endpoint's Gaussian partner. The exact mean
factor is

    exp[mu^2(t^2 s^2/(2a)-lambda)],
    t^2 s^2/(2a)-lambda=-10/9.                          (8)

This coefficient is independent of the OWN row's G/P type. Relative to
the ordinary kernel for that edge, the remaining determinant factor is
sqrt(9/5) for a pinned G row, and sqrt(13/5) for a pinned P row. Hence
all 2m replacements cost at most exp[m log(13/5)].

The norm compensation and all mean factors combine exactly to

    2lambda n+2m alpha^2 k(-10/9)=(37/15)n.              (9)

No center-center term is omitted: it is identically zero by (1). The
potentially order-n mean-noise contribution is retained in (8)--(9).

## 5. Physical entropy is paid with the actual pin constraint

Let eps_L,eps_R be the fractions of P rows on the two sides. For a
physical row write alpha_i=u0/sqrt(k), delta_i=(1-alpha_i)/2. At that
overlap, a G row has at most exp[k h(delta_i)] spins, while a P row has
at most exp[(k/2)h(delta_i)] because its independent physical variables
are the k/2 translation pairs.

Equation (4) gives average (alpha_i-alpha)^2<=tau^2. By Cauchy-Schwarz,
average |delta_i-1/10|<=tau/2. The elementary binary-entropy continuity
bound |h(p)-h(q)|<=h(|p-q|), followed by Jensen, shows that the spin-pair
count for a fixed assignment is at most

    exp[n{(2-(eps_L+eps_R)/2)h(1/10)+2h(tau/2)}
          +2m log(k+1)].                               (10)

The last term sums over all possible row overlaps. This avoids the FALSE
claim that a mixed size-one/size-two group system has the same entropy
drop under only a global Hamming constraint. Here the pinned profile
itself forces the needed row-overlap control.

## 6. Uniform margin and a SAME-ensemble union of row types

Put D=108/25, E=1872/125, P=289/25 and

    A=log D-(1/2)log E,
    B=-(1/2)log D+(1/2)log E-(1/4)log P.

From the exact edge counts in (7), the leading log moment per n is

    F(eps_L,eps_R)=37/15-(1/2)log D
                    +(A/2)(eps_L+eps_R)+B eps_L eps_R.  (11)

The exact rational and directed-interval certificate is
`computations/flatify_independent_2026_09_07_mixed_paired_bridge_certificate.py`
with matching results JSON. It proves

    (1/2)log D-2h(1/10)-1/15 >1/75,
    B>0,                 h(1/10)-A-B>1/5.              (12)

Since eps_L eps_R<=(eps_L+eps_R)/2, the Chernoff rate at bridge level
18/25=.72 exceeds the entropy in (10), before its vanishing errors, by

    >1/75+(eps_L+eps_R)/10.                             (13)

Both bridge polarities have the same kernels. The coupling error in (6)
is at most (t+2lambda)e n=(65/9)e n. Its cost tends to zero per n, as does
2h(tau/2). All remaining terms are sublinear: row conditioning costs
2m log2, pinned-edge determinant ratios cost at most m log(13/5), and row
overlap enumeration costs 2m log(k+1).

Most importantly, the actual ensemble (1) is INDEPENDENT of all G/P
labels outside the four reserved rows. There are at most 2^(2m-4) such
assignments; unioning them costs only another (2m-4)log2. The fixed positive
margin (13) therefore proves that, for every tau_k->0 and all sufficiently
large compatible orders, there exists ONE actual sign bridge such that

    x0^T C y0=0,
    |x^T C y|<=(18/25)n^(3/2)                           (14)

simultaneously for every target pair admitting ANY allowed G/P assignment
whose average paired-profile errors (4) are at most tau_k^2.

## 7. Original-value relevance and limitations

This includes typical iid noisy rows and typical translation-paired noisy
rows at noise level .1, with their covariance differences intact. In
particular it treats, rather than discards, the adversary's high-entropy
paired exceptional family when its paired rows avoid the four reserved
G fibres and use the fixed chosen translations. Typicality follows from
the corresponding two-coordinate and four-coordinate Walsh Lindeberg
argument: a pair test uses two coordinates, and its variance compares two
distinct pairs. Orthogonality and bounded coefficients give the same fixed-
dimension replacement estimate. No atypical
cardinality conclusion is required for (14).

Inside the usual actual-child noisy-energy windows, child cap c gives
combined energy at most (32c/25+o(1))n^(3/2). Thus this sector of the actual
parent has value at most [.72+1.28c+o(1)]n^(3/2), below 2sqrt(2)c already
at c=.47. This is a sector inequality, not a global parent bound.

The four reserved rows are a real profile requirement, despite their
vanishing fraction: their noise receives the pinned-center interactions.
The translations and pair partitions are fixed once, not selected after
seeing a target. General mixtures of other constraints, unequal noise
levels, and profiles outside the G/P references remain unbounded by this
theorem. The actual scalable exceptional counterexample remains compatible
with these limitations. Original convergence is still not established.
