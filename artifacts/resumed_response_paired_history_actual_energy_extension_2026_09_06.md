# Strict extension of the actual-energy optimum of every finite paired history

Date: 2026-09-06. Status: the director and bound-audit agent independently
reconstructed compactness, the largest-index argument, strict extension,
and deletion of asymptotically redundant queries.
This is a fixed-Haar-query-history theorem, transferred only to the
exact-flat/comparator class already audited. It is NOT a theorem for
arbitrary original sign matrices or a claim about the optimal depth limit.

## 1. The finite moment body optimizes actual retained energy

Let H=(q_0,p_0,...,q_d,p_d) be an orthonormal finite paired scalar history
of the adaptive Haar-query construction. Each p_j=N_j is a fresh
independent standard Gaussian. Each q_j is measurable before N_j is
exposed, and all earlier q_i,p_i are also measurable then. The finite
initial side information is independent of all fresh N_j. The limiting
involution acts on the exposed span by the swap matrix J.

Let f range over bounded jointly odd functions of the complete currently
available row history and seed information, with |f|<=1. Put

    c(f)=E[H f],
    E_frame(f)=c(f)^T J c(f)/2.

The paired-query theorem identifies this expression with the actual
normalized half-energy of that output: the unexposed residual has a
fresh independent Gaussian response, whose self-contribution is zero.
This is not the previous terminal lower functional J(F,H), and no
conditional-Jensen reduction of the energy is being made.

The moment body

    M={c(f): |f|<=1, f odd}

is compact. Indeed the feasible L-infinity ball is weak-* compact and
closed under the oddness constraint, and every coordinate of H is in
L1. Its finite-dimensional moment image is compact. Therefore

    C_frame=max_(c in M) c^T Jc/2                (1)

is attained. Nonconvexity causes no problem for attainment, but this
does not assert efficient global optimization. The moment bound is
||c||^2<=E f^2<=1.

An arbitrary bounded measurable maximizer is an attained POPULATION
moment-body/ordered-L2 optimum. It is not silently treated as a fixed
globally Lipschitz implementation. Actual finite algorithms use fixed
approximants close enough to preserve the strict margin established below.

## 2. No maximizer exhausts the norm in the exposed span

Let c maximize (1), with representative f. First-order optimality in
the convex feasible function domain gives

    f=sign(k) wherever k!=0, k=H^T Jc.          (2)

To see this, the derivative toward any feasible v is E[k(v-f)]<=0;
the jointly odd pointwise maximizer sign(k) is feasible. Null ties need
not be resolved for this argument.

Suppose ||c||=1. Equality in the L2 projection inequality forces

    f=H^T c almost surely, |f|=1 almost surely.

Write c=(a_0,b_0,...,a_d,b_d). Conditional on all side information and
the Gaussians preceding N_d, boundedness of f forces b_d=0: otherwise
f contains the nonzero affine Gaussian term b_d N_d. Now f is measurable
before N_d. But k contains a_d N_d. If a_d were nonzero, conditionally
sign(k) takes both values with strictly positive probability, contrary
to (2) and the already fixed Boolean value of f. Thus a_d=0 as well.

Repeat for d-1,d-2,...,0. All coefficients vanish, contradicting
||c||=1. Hence EVERY maximizer satisfies ||c||<1. Compactness of the
argmax set gives the strict, frame-dependent quantity

    epsilon_frame=1-max_{c maximizes (1)}||c||^2>0. (3)

There is also delta_frame>0 such that every feasible c with energy
at least C_frame-delta_frame satisfies
||c||^2<=1-epsilon_frame/2. Otherwise compactness produces a norm-one
maximizer. No explicit uniform bound over different histories follows.

## 3. Boolean realization and one actual new query

If f is fractional, conditionally round it to a Boolean u using fresh
independent row seeds, so E[u|old history]=f. This preserves c and hence
the actual limiting energy. A pathwise spin-equivariant construction
uses a fresh even uniform gate U (the CDF of |G_new|), and a separate
fresh odd Rademacher seed R:

    u=sign(f) 1{U<|f|}+R 1{U>=|f|}.

Every multiplication input remains jointly odd. The gate is independent
of the entire old history. Bounded measurable rules are passed through
L2 Lipschitz approximation before the matrix limit; the independent
uniform gate makes rounding probabilities stable under L1 errors.

Query the actual vector u, not an unrelated surrogate. The paired
conditional kernel yields

    Bu -> k+sigma N_new,
    k=H^T Jc, sigma^2=1-||c||^2>0.              (4)

The fresh Gaussian is independent of all row history and side seeds.
Thus the actual stability gain

    g=E|k+sigma N_new|-E[u k]

is at least Gamma_sigma(sqrt(1-sigma^2)), with
Gamma_sigma(z)=E|z+sigma N|-z. This follows from the convex-in-z^2
Jensen calculation in the mixed-charge feedback note. In particular,
for a true old-frame maximizer,

    g>=g_frame:=Gamma_sqrt(epsilon_frame)
                       (sqrt(1-epsilon_frame))>0. (5)

Here Gamma_sigma increases with sigma and decreases with its nonnegative
argument, which licenses using the lower floor (3).

The feasible population-closure damped endpoint

    m=(1-eta)u+eta sign(Bu), eta=g_frame/4,

has normalized energy at least C_frame+g_frame^2/8, by the exact
operator-norm-one quadratic inequality. Independent final Boolean
rounding preserves that limiting energy. Select fixed Lipschitz/gate
approximants before the matrix limit, with total energy error less than
g_frame^2/16. This produces a genuine fixed finite algorithm of energy
strictly greater than C_frame+g_frame^2/16. Softsign first, then its
limit, avoids any finite-n zero convention; (4) also supplies the needed
no-atom fact. No exact implementation of a pathological Borel maximizer
is required.

The new algorithm uses only finitely many extra seeds and one new
matrix query to construct m. To EXPOSE its energy in a state-evolution
test one may of course use a further matrix multiplication; this is
an analysis step, not an extra input to the output rule.

## 4. Consequence and the essential scope limit

Every fixed finite paired-query history has an attained population
optimum for actual energy, and that optimum is strictly improved by extending the history
with one genuine query and its feedback. The same holds uniformly for
sufficiently near-optimal old-history policies, using the half-sized
floor from Section 2 and a fixed positive gain depending on that frame.

This is a strict-extension theorem for retained-energy algorithms. It
does not recycle the final terminal functional and therefore is not
subject to the terminal family's .45 bound merely by construction.
It does not imply a common positive increment along successively larger
histories. The constants in (3) may collapse, global moment-body
optimization need not be practical, and the separate fixed-GFOM Haar
upper ceiling remains applicable to every individual extension. No
conclusion for the original min-max convergence problem follows without
an additional comparison or dimension-growing argument.

## 5. The fixed-algorithm supremum cannot have a finite-history maximizer

Let C_fixed be the supremum of limiting normalized positive energy over
all fixed finite equivariant algorithms in the paired-query/GFOM class.
This is a well-defined bounded number. Passing from absolute to positive
energy causes no difficulty: running a fixed rule at -B reverses its
limiting energy measured with B. The rigorous available bracket is

    .4333221116640807 <= C_fixed <= sqrt(15)/8.

The lower endpoint comes from the already certified terminal algorithm
and its exact Hadamard endpoint law; the upper endpoint is the separate
Haar/GFOM union-bound theorem. The numerical feedback evidence for
C_fixed>.45 is not substituted for a rigorous lower endpoint here.

No finite paired-history algorithm attains C_fixed. Otherwise its finite
history would have C_frame at least its value and at most C_fixed, hence
C_frame=C_fixed. Section 3 gives a finite further algorithm with larger
energy, a contradiction.

The same conclusion applies to fixed globally Lipschitz GFOMs after
deleting asymptotically redundant queries. Indeed a query with zero
limiting residual variance is an L2 linear combination of previously
exposed query/response vectors; multiplication by the bounded operator
replaces its response by the swapped combination. A finite Lipschitz
continuation propagates this L2 equivalence. Removing such queries
leaves precisely a finite nondegenerate paired history (or a seed-only
zero-energy rule). This reduction does not require a nonsingular Gram
matrix at every original step, and it does not merely approximate a
maximizer by frames whose strict-extension constants could disappear.

Nonattainment is not a numerical value theorem. It is consistent with
strictly increasing attainable values converging to any limit inside
the displayed bracket, with vanishing fresh innovation. In particular
it neither identifies the original-signing infimum nor proves convergence
of its normalized finite-size sequence.
