# Boolean BH boundedness: closing synthesis

Campaign: 2026-09-18, 04:27:50--07:27:50 UTC. This synthesis is being
checked during the final verification period. The original quadratic-signing
convergence problem remains paused and its reported bounds are unchanged.

## Outcome

For COMPLEX degree-at-most-m Walsh polynomials, neither uniform boundedness
of B_m nor a diverging sequence has been proved. The strongest direct
statements about the unrestricted constants obtained in this campaign are

```math
B_4>2+2^{-18},\qquad
\liminf_{m\to\infty}B_m>2+2^{-207}.
```

These are strict CONSTANT lower bounds, not growth estimates. They do not
resolve the user's dichotomy and are not asserted for the real field.
No claim of external novelty is made. Exact finite arithmetic and complete
analytic limit arguments, not decimal optimization, support both results.

## 1. Direct lower-bound mechanism and dependencies

The [counterexample proof](bh_boundedness_counterexamples_2026_09_18.md),
Sections 5, 7 and 8, supplies an explicit five-bit, degree-two polynomial f
with sixteen flat coefficients, sixth-root values and ratio exactly two.
An exact pair a,b has degrees 1 and 3 and b=f^2 conjugate(a). On two copies,

    F=f tensor f,  h=a tensor b-b tensor a,
    Re(conjugate(F)h)=0,
    |h|^2 in {0,48}.

The tangent h creates 36 new coefficients, with squared mass 45/4, without
first-order loss on the old flat support. For epsilon=2^-10, new critical
coefficient mass has order epsilon^(8/5), while the cap cost is quadratic.
Convexity and rational bounds give the displayed B_4 inequality.

For the persistent gap, fix epsilon=2^-104, R=|h|^2/48, and
lambda=log(1+48epsilon^2)/2. The normalized phase G=(F+epsilon h)exp(-lambda R)
has entropy exceeding 8log2 by a quantified epsilon^2 log(1/epsilon) term.
Joint polynomial normalization of k copies has degree rate 4+16lambda,
not the rate 8 of naive exact normalization, and uniform error tending to
zero. The final exponent q_m and bounded gaps between degree budgets are
handled explicitly; this proves a liminf over EVERY sufficiently large m.

Two separately coded Eisenstein-integer replays verify the seed and tangent.
Eleven rational checks verify the analytic margins used for normalization;
500-digit calculations are supplementary diagnostics. Both other researchers
and the director reconstructed the entire proof, including the order of limits.
The director's final independent physical-cube verifier works instead in
Z[sqrt(-3)], evaluates all 1,024 vertices, and transforms back to coefficients.
It again certifies the degree-eight normalization term and rational B_4 margin;
its source and exact JSON output are tracked in computations/.

## 2. Strongest positive class theorems

The [positive proof](bh_boundedness_positive_2026_09_18.md) establishes:

| Class | Verified dimension-uniform conclusion | Essential hypothesis |
| --- | --- | --- |
| Permutation-symmetric complex polynomials | BH ratio at most sqrt(6e) | Full permutation invariance, not merely transitivity |
| Unit-circle-valued polynomials affine in each of k disjoint blocks, k>=1 | Support at most 4^k/2 and q_k norm at most 2^(1-1/(2k)); both sharp | Exact constant modulus and block-affinity |
| Root-of-unity-valued degree-m polynomials | Support at most 4^m and BH ratio at most 2 | Values are roots of unity; alphabet size may grow |

The third proof is in the [barrier artifact](bh_boundedness_barrier_2026_09_18.md),
Section 7. It uses algebraic integrality, all Galois embeddings and Holder.
The sharp block-affine proof uses odd Boolean postprocessing and an even
integer lattice; it covers arbitrary phases, not only roots of unity.
Exact degree-two examples show that this Boolean postprocessing need not
preserve degree outside the block-affine class.

A quantitative one-block stability theorem also holds: an affine p with
|||p|^2-1||_infinity<=eta<=1/4 has sup distance at most 9sqrt(eta) from an
exact affine unitary, independently of dimension. Compatibility across
several blocks remains unproved. None of these results gives a decomposition
of arbitrary bounded polynomials into the stated positive classes.

## 3. A reusable approximation theorem, and its quadratic application

The late combined result is in the [director proof](bh_boundedness_director_2026_09_18.md),
Sections 7 and 9, and barrier Sections 31--32. For a unitary seed G on n bits,
let V(G) be the largest absolute expectation of M_G^*JM_G-J over density
matrices with uniform physical diagonal, where J=(1/2)sum X_i. If Q is a
degree-D polynomial on k seed copies, C=||Q||_infinity>0, and
delta=||Q-G^(tensor k)||_2, then the EXACT finite inequality is

    D/k >= [V(G)-n delta]/C.

Proof ingredients: operator trigonometric Bernstein; a tensor product of
flat-diagonal states; two Cauchy--Schwarz bounds for the approximation error.
The diagonal constraint, not an assumed spectral identity, makes the error
the ordinary uniform-cube L2 norm.

For every real degree-at-most-two P, at arbitrary time and coefficients,

    H(|widehat(exp(iP))|^2)<=64 V(exp(iP)).

The proof uses simultaneous product-state row witnesses. Averaging their
auxiliary signs leaves the density diagonal exactly uniform, even though
their other coordinates were chosen globally to align all row contributions.
Both independent researchers reconstructed this compatibility step.

Consequently, for any FIXED nonconstant quadratic-phase seed on n bits,
polynomials on EXACTLY nk bits with L2 error tending to zero and cap tending
to one satisfy

    limsup_k R_(actual degree of Q_k)(Q_k)<=exp(32).

For independent linear phases the bound improves to exp(1). These constants
are not optimized. This closes the actual-ratio version of that fixed-seed
approximation construction, not just its Shannon-entropy lower-bound formula.
Extra auxiliary variables and simultaneous growing-seed limits are explicitly
excluded unless their additional error factors are controlled.

## 4. Falsifiers that changed the argument

- The weighted bootstrap's level-q_r norm has an unavoidable linear loss
  at EVERY fixed retained level. The same examples have bounded full q_m
  norm. This is a method obstruction, not a lower bound on B_m.
- Ordinary disjoint tensoring cannot increase the largest critical ratio
  of its factors. Selector gates at their natural degree budget do not
  amplify a constant above two. Actual degree collapse is a distinct case.
- The normalized ten-bit seed loses its nonconstant asymmetric
  degree-transfer pairs at total degree budget eight, certified by exact
  finite-field minors. Its symmetric (4,4) space is the one-dimensional
  span of P and has zero antisymmetric transfer. The trivial constant
  pair also survives; a universal iteration no-go is not inferred.
- Full operator velocity is NOT an L2 approximation lower bound. A fixed
  sixteen-bit rare-AND sign has velocity four, but cap-one L2 approximants
  to its tensor powers have degree rate at most 1/sqrt(2). Exact rare-mark
  enumeration from de Wolf's primary quantum-query construction proves this.
  The replacement flat-diagonal certificate respects that counterexample.
- Fixed-degree unimodular functions are finite juntas, but the proved
  Ramsey bound is huge. It does not supply an exponential-in-degree support
  bound or bounded BH constants. For any fixed error delta<1, there is no
  approximation of ALL supremum-norm-at-most-one degree-m polynomials, uniformly over m and
  ambient dimension, by unitary atoms with
  arbitrary finite dimension-independent degree allowance D(m) and total
  mass C(m). Individual decompositions can succeed; this is not a claim
  that every use of such atoms fails.

## 5. Literature and verification scope

The actual linked [Slote--Volberg paper](https://arxiv.org/abs/2609.07758)
was reconstructed separately from Ivanisvili's 2609.12427. The expert's
later exponent reports were not supplied proofs. The precise linear
auxiliary-norm obstruction above is proved here; we do not claim to have
reconstructed every optimization behind the reported exponent-one barrier.

Entropy/radial normalization has an antecedent in the complex-polydisc
paper 2608.16584; its domain and constants are different. Boolean endpoint
bounds and operator degree/commutator estimates also have primary antecedents.
New combinations and adapted proofs are identified without claiming priority.

In particular, [Aaronson--Ben-David--Kothari--Rao--Tal, arXiv:2010.12629](https://arxiv.org/abs/2010.12629),
Section 4.3, already contains the operator commutator/degree ingredient
and its pointwise Boolean approximation application. The new proof here
adds the explicit flat-diagonal L2 certificate and the quadratic-phase
entropy witnesses; its novelty has not been established by the limited search.
Rotation averaging gives B_m(real)<=B_m(complex)<=c_(q_m)^(-1)B_m(real),
where c_q=(E|cos theta|^q)^(1/q)>=2/pi. Thus boundedness is equivalent
over the two fields, but a strict complex lower bound above two does NOT
imply a real lower bound above two.

The closing replay runner executes thirteen selected certificate/diagnostic
programs, hashes their sources and preserves their full outputs. Every replay
passed. Diagnostic eigenvalues, searches and floating-point checks remain
explicitly separate from exact certificates and analytic proofs. New code
uses explicit seeds/constants or tracked sibling modules, not ignored input data.

## 6. Remaining obstacle and continuation judgment

The successful tangent step is not yet a scalable amplification operation.
Its fixed ambient/effective-degree ratio imposes a constant ceiling; naive
normalization destroys the useful asymmetric transfer pairs. A divergent
construction needs a growing family that retains coefficient entropy while
paying proportionately less degree, with cap control and changing exponents
proved, not assumed. The present work supplies no such family.

On the upper side, uniform control for torsion, symmetric and block-affine
classes, and for the stated fixed quadratic-phase compiler families, does
not cover arbitrary complex bounded polynomials.
In particular H(G)<=constant V(G) for ALL unitary G would include the open
Boolean Fourier-entropy/influence conjecture, since V(b)=I(b). It is not an
available shortcut. Shannon estimates alone also do not bound all critical
Renyi tails of unrestricted polynomials.

Director judgment: do NOT automatically extend the campaign or continue
tuning the fixed constant. A longer run on the tested tensor, address,
normalized-seed or fixed quadratic-phase constructions is not justified by
the present evidence. A new campaign needs a concrete growing-degree
operation or a genuine extension beyond those proved classes. The original
signing problem has not been reopened, solved, or improved by these results.
