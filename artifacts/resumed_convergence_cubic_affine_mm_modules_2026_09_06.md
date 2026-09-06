# Cubic affine-MM modules and complete binary Fourier-sign closure

Date: 2026-09-06. Status: exact finite signed-module identities for
cubic phases and cubic controlled-index gates; a separate tensor
quotient argument gives all binary Fourier signs as limiting tests.
No arbitrary index-permutation or gate-composition theorem follows.

## 1. The cubic derivative identity

All polynomials are binary algebraic-normal-form polynomials. For
a homogeneous squarefree cubic P on F_2^r, put F=grad(P), with
formal partial derivatives. Then

    a dot [F(v+a)+F(v)]=P(a)       for every a,v.           (1)

It suffices to check P(v)=v_i v_j v_k. Every mixed term involving
v occurs twice after contracting its three derivatives with a;
the pure a monomial occurs three times and hence once.

More generally every P of degree at most three with P(0)=0
admits a degree-at-most-two vector map F satisfying (1):

* use the full gradient for every cubic monomial;
* for a quadratic monomial a_i a_j, add v_j only to F_i;
* for a linear monomial a_i, add v_i to F_i.

All these contributions add over F_2.

## 2. An exact finite signed module for every cubic Fourier phase

On F_2^r x F_2^r define

    h(x,y)=(-1)^(x dot y),     Phi(x,y)=x+F(y).

The label Phi is exactly uniform, since x is uniform for each y.
For every a the signed component is a Maiorana--McFarland bent
function, with exponent

    q_a(x,y)=x dot y+a dot [x+F(y)].

In its normalized Walsh transform, summing first over x forces
y=u+a. Thus its dual exponent at (u,v) is

    v dot (u+a)+a dot F(u+a)
       =u dot v+a dot [v+F(u)]+P(a).                      (2)

Swap the two output coordinate blocks and write W for this
symplectic Walsh involution. Equation (2) says exactly

    W[h (-1)^(a dot Phi)]
             =(-1)^P(a) h (-1)^(a dot Phi).

For every real scalar profile f on F_2^r, therefore,

    W[h f(Phi)]=h [T_P f](Phi),
    T_P=H_r diag((-1)^P) H_r.                             (3)

This is an invariant SIGNED module with no limiting argument,
no removed constant channel, and exactly uniform labels. For
every Boolean profile array F_profile and symmetric seed B,

    R(B)>=(1/2) E_t ||[T_P F_profile](t) B||_1.           (4)

The physical Walsh order is 2^(2r), an allowed even order in
the bilinear regularization. A constant term in P only changes
the global sign of T_P and does not affect (4).

The construction proves overlapping cubic Fourier-phase gates
in one exact finite Walsh module. It does not concatenate
operators by assuming their separate realizations compose.

## 3. Sharp scope of the one-shift identity

The condition (1) cannot provide higher-degree phases even if
one allows an arbitrary nonlinear F. Subtract F(0), harmless
in (1), and assume F(0)=0. Setting v=0 gives P(a)=a dot F(a).
Writing (1) for a=u+v then gives

    P(u+v)+P(u)+P(v)=u dot F(v)+v dot F(u).                (5)

The right side has only mixed bidegrees (1,d) or (d,1) in the
separate groups of variables u and v. A monomial of degree
d>=4 in P would produce on the left a mixed split of sizes
2 and d-2, with neither degree one. Its coefficient cannot
cancel with another monomial, so all such monomials vanish.
Hence deg(P)<=3. This is a scoped maximality result for (1),
not a restriction on all affine bent families.

## 4. A genuinely nonlinear controlled-index permutation

Let P,F satisfy (1). On the STATE space with coordinates
(u,v,w) in F_2 x F_2^r x F_2, define the permutation

    pi_(a,c)(u,v,w)
       =(u, v+u a, w+c u+a dot F(v)).                     (6)

It depends AFFINELY on the INDEX parameters (a,c), although
it is generally nonlinear in the state variables. Solving
for the inverse and using (1) gives

    pi_(a,c)^(-1)=pi_(a,c+P(a)).                          (7)

Indeed, after v is replaced by v+u a, the extra last-coordinate
term is u a dot [F(v+a)+F(v)]=u P(a).

Consequently the affine MM family

    q_(a,c)(x,y)=x dot pi_(a,c)(y)

has, after the standard swap of output state blocks, the
nonlinear affine dual-index permutation

    sigma(a,c)=(a,c+P(a)).                                (8)

This is the controlled-cubic index gate itself, not merely a
diagonal Fourier phase. The associated profile operator is
H_(r+1) P_sigma H_(r+1).

To check equidistribution for the signed module, write
q_(a,c)=x dot y+(a,c) dot Phi(x,y). For every nonzero index,
the character expectation of Phi equals the fraction of y
for which pi_(a,c)(y)=y. If a!=0, the branch u=1 has no fixed
point because v shifts by a, so that fraction is at most 1/2.
If a=0,c=1, exactly the branch u=0 is fixed, again fraction
1/2. Thus direct sums equidistribute all labels at a rate
at most 2^(-s) in each nonzero character. The affine signed
module lemma applies to ALL Boolean profiles, including
unbalanced ones, with physical base order 2^(2r+4).

For a cubic P the family is not quadratic in the physical
variables, so the inverse-POLAR-matrix information bound does
not apply to it. The number of cubic choices in a fixed
label dimension is nevertheless only exponential in a
polynomial of that dimension. Thus this alphabet expansion
alone does not remove the independent-cloud matching barrier.

## 5. All binary Fourier signs already follow from reflections

The earlier simplex quotient construction makes the uniform
averaging reflection Q_N=2E_N-I on N=2^r atoms available as
a limiting probability-operator test for R. Identify those
atoms with F_2^r. For b in F_2^r let D_b multiply a profile
by the Boolean character chi_b. The gauged operator

    S_b=-D_b Q_N D_b

has Fourier multiplier -1 at b and +1 at every other
frequency. Independent row and column Boolean gauges are
legitimate in the bilinear R test.

Given ANY prescribed sign function epsilon:F_2^r->{+1,-1},
tensor one S_b for each b with epsilon(b)=-1, on independent
copies of the label group. Restrict to the SUM quotient

    L(x_1,...,x_t)=x_1+...+x_t.

A character chi_a pulled back along L is the product of
chi_a in every factor. The tensor operator therefore has
multiplier product_b S_b_hat(a)=epsilon(a) on this quotient.
Characters span its function algebra, so the quotient is
invariant and its exact limiting operator is

    H_r diag(epsilon) H_r.                                (9)

For fixed r the number of factors is finite. Approximating
each uniform N-point reflection by the earlier weighted
simplex tests and then taking their finite tensor product
preserves every fixed bilinear profile objective. The sum
quotient itself is exactly uniform on the limiting uniform
atom spaces. If there are no negative signs, use the identity
profile test. These observations justify the finite-then-limit
order of operations in (9).

Thus (3) adds an exact finite Walsh implementation for cubic
phases, but does NOT enlarge the already available limiting
binary diagonal alphabet. The controlled-INDEX map (8) is a
different operation and is not supplied by (9).

The diagonal closure is a legitimate simultaneous tensor
quotient. It is not a claim that arbitrary separately
available operators on overlapping coordinates compose.

## 6. Connection and current gap

The arbitrary Fourier-permutation matching theorem computes
the majorant T(B) in an enlarged model. The present modules
still do not realize those data-dependent matching indices.
Arbitrary diagonal signs do not change Fourier feature
directions, and a fixed-dimensional cubic-index family does
not have the exponential index freedom used by that proof.

No R=T claim, general polar realization, improved original
constant, or convergence conclusion is asserted here.

Exact integer replay is in
`computations/resumed_convergence_cubic_mm_verify_2026_09_06.py`.
It passes for all 128 zero-constant polynomials on three bits
and 17 four-bit polynomials, including overlapping cubics.
It checks (1), every state permutation and its inverse index,
all component Walsh identities, the exact uniform phase
labels, and the controlled-index label contraction bound.
