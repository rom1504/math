# Random-selector entropy and a strict stability improvement of the actual weave

2026-09-07. **Verified by two independent reconstructions.**
The proposed conclusion is a new unconditional strict improvement of an
existing all-order upper certificate, not convergence or ensemble optimality.

## 1. Sparse spectral rows cannot hide their selector entropy

Let k/m=p_m->p in (0,1). For ANY orthogonal m-by-m matrix O, consider
f in {0,+-1/sqrt(p_m)}^m with exactly k nonzero coordinates and put u=Of.
Every such u has squared norm m. Fix epsilon>0 and V>sqrt(2). Say u is
concentrated if

    sum_{j:|u_j|<=V} u_j^2 <= epsilon m.

For sufficiently small epsilon, its number of possible f is at most

```math
\exp\{m[\xi(\epsilon,V)+o(1)]\},
\qquad
\xi=h(V^{-2})+V^{-2}\log(1+2/\sqrt\epsilon)
            +h(16\epsilon)+16\epsilon\log2.                 (1)
```

The error is uniform over O. Here h is binary entropy, and epsilon<1/32
ensures the Hamming-ball entropy estimate is on its increasing branch.

Proof. The heavy support has size at most m/V^2. Its coordinate projection
lies in a radius-sqrt(m) Euclidean ball and is within sqrt(epsilon m) of u.
For every such support take a sqrt(epsilon m)-net of that ball with at most
(1+2/sqrt(epsilon))^d points. Summing over d<=m/V^2 costs
exp[m h(V^-2)+o(m)]. A selected net point is within2sqrt(epsilon m) of u.
Pull it back by O^T and round each coordinate to the nearest of the three
allowed alphabet values. Any discrepancy from f costs at least1/(4p_m)
in squared distance. Thus at most16p_m epsilon m<=16epsilon m coordinates
can disagree. Each disagreement has at most two alternative symbols.
The Hamming ball gives the last two terms of (1). Enforcing exact support
k can only decrease that bound.

Now choose the retained physical row set T uniformly among all k-subsets,
independently in each fibre and independently of its row basis. Summing
over its 2^k row-spin vectors and averaging T is exactly the count over
these ternary f divided by binom(m,k). Consequently for EVERY row basis,

    E_T sum_{x_T: concentrated} max_j L_t(u without j)
       <= exp[m(-h(p)+xi(epsilon,V)+o(1))],           (2)

because 0<=L_t<=1. This is a direct count, not a Gaussian assumption or
an inference from a posterior channel. Choosing epsilon small and then V
large makes xi arbitrarily small, independently of m or the basis.

## 2. Diffuse rows force many unfrozen incoming fibres

Consider ANY array u_i(j) with sum_j u_i(j)^2=m for every i. Suppose at
least delta m rows are not concentrated, so their light energy exceeds
epsilon m. The total directed light mass is at least delta epsilon m^2.

Choose U,C with

    U^2 >= 4V^2/(delta epsilon),
    C   >= 4V^2/(delta epsilon).

Discard entries whose reverse coordinate has magnitude above U. There
are at most m^2/U^2 such entries, each carrying light weight at most V^2.
Discard columns whose total incoming squared norm exceeds Cm. There are
at most m/C such columns, each carrying light weight at most V^2 m.
Diagonal entries cost at most V^2 m. The remaining off-diagonal light-light
mass is at least delta epsilon m^2/3 for sufficiently large m.

It follows, using the per-column upper bound V^2 m, that at least

    theta m columns, theta=delta epsilon/(8V^2),

have incoming norm squared at most Cm and at least delta epsilon m/8 of
light mass whose reverse coordinates have magnitude at most U. Under the
actual edge tilt, their light variance is therefore at least kappa m, where

    kappa=(delta epsilon/8) sech^2(2t U V)>0.

This assertion is deterministic and holds for EVERY output-column
permutation. It does not assume that row and column variables are independent.

By the independently verified orthant theorem in
`principle_director_stability_orthant_penalty_2026_09_07.md`, each such fibre
has stability probability at most exp(-c0 k), with a fixed c0>0. Finner
therefore supplies an exp(-theta c0 mk/2) penalty for full-spin stability.

## 3. Combine the two sectors before optimizing the cap

Let a be a valid limiting exponent for the existing full row sum, including
its spin count and maximum deleted-coordinate L_t:

    E_{basis,T} Z_i <= exp[(a+o(1))m].

Assume a>-h(p). Choose epsilon,V so that
xi(epsilon,V)<[a+h(p)]/2, and fix delta=1/2. Put
b=-h(p)+xi, so b<a. The full spin sum is split into configurations with
at least delta m diffuse rows and configurations with fewer.

For the diffuse sector the deterministic conclusion of Section2 pulls the
stability factor out before averaging any column permutations. The original
PSD graph-Cauchy--Schwarz/permanent argument then bounds the remaining sum
by the SAME product of the full Z_i. Its logarithmic exponent improves by
at least theta c0 p/2 at scale m^2.

For the other sector at least (1-delta)m rows are concentrated. Specify a
set of such rows; their indicator depends only on the row multiset and is
unchanged by output-column permutations. It therefore passes through the
same graph-Cauchy--Schwarz bound. Independence of selectors and bases gives
a product of the row estimates (2) and the ordinary full-row estimate.
Summing over at most2^m choices changes the logarithm by only O(m).
The exponent in this sector improves by at least (1-delta)(a-b).

For example choose the fixed positive number

    Delta=min(theta c0 p/4, (1-delta)(a-b)/2).

Use the common full-row upper a+zeta, with zeta strictly smaller than the
already halved Delta, and choose a sufficiently large but FIXED depth to
achieve that bound. No attainment of the limiting row exponent is assumed.
The harmless halves absorb this depth error, finite-order errors, k/m-p,
and the union of the two sectors. Counting both energy polarities and only actual locally
stable extrema now gives the all-order bound

```math
\limsup_n\frac{M_n}{n^{3/2}}
 \le \frac{t+a-\Delta}{2t\sqrt p}.                 (3)
```

All physical diagonal contributions are handled through the exact hollow
stability inequalities in the companion proof. In the moment estimate one
drops only nonnegative diagonal defects, exactly as in the original chain.
No macroscopic diagonal block is removed by a triangle inequality.

The required existing row theorem is uniform over fixed selectors, so
averaging independent uniform selectors does not weaken it. Conversely,
the new concentrated-row count is uniform over the entire basis. Therefore
one may first fix epsilon,V,U,C,theta,kappa,c0,Delta; then fix a finite
recursive depth large enough for the old row exponent; then let admissible
orders grow; finally use the existing H2/H12 ratio-density and principal
restriction to obtain every sufficiently large order. No depth is allowed
to grow inside the fixed-depth type calculation.

## 4. Application to the current certified upper endpoint

Take p=24/25, t=97/20, and

    a=p log2-5151/6250.

The old directed-interval certificate and Gaussian-boundary closure give
the required row exponent. Its entropy gap is strictly positive:

    a+h(p)=h(24/25)+(24/25)log2-5151/6250 > 9/1000.

Thus (3) is strictly smaller than the old exact endpoint

    [97/20+(24/25)log2-5151/6250]/[(97/10)sqrt(24/25)].

The tracked outward-rounded checker
`computations/principle_director_2026_09_07_selector_entropy_certificate.py`
verifies the entropy premises with epsilon=1/100000,V=100,delta=1/2,
U=100000,C=10^10: xi<.004<gap/2 and gap>.009. Its JSON records the
exact symbolic positive thermal-variance constant, not a floating zero.
The independent construction-role checker uses an even smaller epsilon.

The elementary constants furnished by this proof may make Delta extremely
small. No changed displayed decimal or numerical optimization is claimed.
The substantive step is a strict actual-sign construction improvement from
joint local stability plus the selector entropy of concentrated rows.

## 5. Stronger theorem: a selector can exclude EVERY concentrated word

**Independent audits: both construction and invention roles PASS.** The
two-sector proof above is valid, but its assumption a>-h(p) is unnecessary.
The following argument improves every fixed p in (0,1), finite t>0 row
certificate to which the existing weave counting theorem applies.

Let N_bad(O,T) count the concentrated words in the fibre selected by T.
The unweighted version of (2) gives, uniformly over EVERY orthogonal O,

    E_T N_bad(O,T)<=exp[-(h(p)-xi+o(1))m].

Choose epsilon,V so xi<h(p). This is possible precisely because p is
strictly between zero and one. Since N_bad is a nonnegative INTEGER,

    P_T(N_bad(O,T)>0)<=exp(-c m)                         (4)

for a fixed c>0 and all sufficiently large m, uniformly in O. Thus the
selector, not a particular spin word, can be chosen to eliminate every
concentrated response simultaneously. All 2^k physical spin words remain
available; none are discarded from the original optimization.

For each independent row basis choose its selector conditional on
N_bad=0. The conditional density is at most (1-exp(-cm))^-1 uniformly
in the basis. Consequently ANY nonnegative old row statistic Z satisfies

    E_new Z <= E_old Z/(1-exp(-cm)).                      (5)

In particular its limiting row exponent a is unchanged. Output signed
column permutations preserve the event N_bad=0: they merely permute/sign
the coordinates of every u. Their conditional uniform distribution is
therefore also unchanged. Independent conditioning in each fibre preserves
fibre independence. These facts are exactly what the graph-CS argument
requires; no distributional symmetry of a generic orthogonal matrix is
assumed beyond the output orbit already present in the construction.

Now Section2 applies with delta=1 to EVERY full spin configuration. Take
theta=epsilon/(8V^2), kappa=(epsilon/8)sech^2(2tUV) and the corresponding
c0>0 from the orthant theorem. The common stability gain is
exp(-theta c0 mk/2). After fixing all these constants, choose a fixed old
recursive depth with row exponent at most a+theta c0 p/8. Let m grow and
absorb the remaining o(m^2) errors. For example the safe conclusion is

```math
\limsup_n M_n/n^{3/2}
 \le \frac{t+a-\Delta_*}{2t\sqrt p},
\qquad \Delta_*=\theta c_0p/8>0.                       (6)
```

The larger numerical slack in (6) leaves room for finite-depth, rounding,
selector-conditioning, and all-order restriction errors. The all-order
argument is the same fixed-depth H2/H12 construction as above. Conditional
sampling is a genuine finite construction law with nonempty support, not
conditioning on the desired parent cap or on an unknown optimum.

Equivalently, keep the original unconditioned selectors. Their probability
of any bad fibre is at most m exp(-cm). On its complement the uniform
stability penalty applies; after extracting that penalty, dropping the good
selector indicator only enlarges the original nonnegative moment. The cap
failure probability is then at most m exp(-cm)+exp(-c' m^2). This version
does not even require altering the sampling law.

At p=24/25 one may take epsilon=1/1000,V=20,U=2000,C=2000000.
The directed-interval checker verifies xi<.122<h(p), with
theta=1/3200000 and kappa=sech(388000)^2/8000. These are explicit
positive analytic constants; evaluating them as floating zero is not a
valid estimate. The improvement remains too small to advertise as a new
displayed decimal.

This is a finite-alphabet uncertainty statement followed by an actual
local-stability count. It explains why random retained-coordinate sets
can remove all sparse spectral resonances without losing the independent
row algebra. It does NOT assume that the old annealed exponent is attained.

## 6. What this does and does not explain

An upper certificate which counts every spin can be slack for two different
reasons: energy can hide in a low-description spectral sector, or diffuse
states can violate many local optimality inequalities. Random selectors
charge the former sector; convex product concentration charges the latter.
They are complementary estimates and are combined before the absolute
maximum, not separately paid energy channels.

This supplies a positive principle in the SUCCESSFUL recursive sign family.
It does not show that that family realizes arbitrary favorable exact seeds,
that its true cap has a limit, or that it is asymptotically optimal among all
signings. Those remain separate obligations. In particular, strictness of
the old certificate is not a proof of the original convergence theorem.
