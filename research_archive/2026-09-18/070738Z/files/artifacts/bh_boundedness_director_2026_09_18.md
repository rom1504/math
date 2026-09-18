# Director derivations: Boolean BH boundedness

Status: separate three-hour campaign, 2026-09-18. The signing problem is
paused. Results below concern the actual scalar Boolean supremum, not
a completely bounded or torus norm. External novelty is not asserted.

## 1. Uniform control for open and closed Hadamard chains

**Proved; independently reconstructed by the barrier researcher.**
Let d>=2 and let H_1,...,H_(d-1) be real or complex Hadamard matrices
of common order N. Thus U_j=H_j/sqrt(N) is unitary and every entry of
U_j has modulus N^(-1/2). Give each slot its own N Boolean variables.
The normalized open-chain polynomial is

    T(x_1,...,x_d)
      =N^(-1) x_1^T U_1 D(x_2) U_2 ... D(x_(d-1)) U_(d-1) x_d.

Its N^d distinct coefficients have modulus N^(-(d+1)/2). Consequently
its coefficient q_d norm is exactly one. Its scalar Boolean supremum is
at least 1/sqrt(2) over the real field, and at least sqrt(2)/pi over
the complex field. Thus its BH ratio is at most sqrt(2), respectively
pi/sqrt(2), uniformly in BOTH d and N.

Proof: for d>=3 fix the variables from x_3 onward and write the tail
as v, with ||v||_2=sqrt(N). Randomize x_2. Every coordinate of
U_1 D(x_2)v has squared Rademacher coefficient norm one. Khintchine's
L1 inequality gives expected coordinate modulus at least 1/sqrt(2).
In the real case choose x_1 to align all coordinates. In the complex
case use the elementary inequality

    max_(epsilon_i=+-1) |sum_i epsilon_i z_i|
       >=(2/pi) sum_i |z_i|.

Indeed maximize the real part after a global phase, and average that
phase. This proves the stated bound after division by N. For d=2
randomize x_d directly. The complex Rademacher L1 constant 1/sqrt(2)
follows from its real version: average real projections, then minimize
the concave average sqrt(lambda_1 cos^2 t+lambda_2 sin^2 t) at a
rank-one covariance with fixed trace. No tensor polarization is used.

The same conclusion holds for the cyclic polynomial

    C(x_1,...,x_d)=Tr(D(x_1)U_1 ... D(x_d)U_d),

where all d edges are normalized Hadamards. Its coefficient q_d norm
is sqrt(N). Fix x_3,...,x_d. The remaining bilinear coefficient matrix
has entries (U_1)_(ij) V_(ji), with V unitary. Every row has squared
norm 1/N. Random x_2 and the same sign/phase argument give Boolean cap
at least sqrt(N/2), respectively sqrt(2N)/pi. This includes the closed
cycle; closing a path does not create a degree-growing BH ratio here.

Scope: slots are disjoint, all edges are square of the same size, and
each has flat unitary entries. Identifying variables can create Walsh
collisions and degree collapse and is NOT covered by this proof.

Literature connection: the primary September 14 paper
[Hadamard Rigidity and Sharp Stability in the Completely Bounded
Bohnenblust--Hille Inequality](https://arxiv.org/abs/2609.16329)
identifies Hadamard chains as its exact extremizers. Our scalar estimate
above is elementary and does not depend on that classification proof.
Its fixed-degree stability theorem is in coefficient norm, with
degree-dependent constants. Neither fact supplies uniform scalar BH
control for arbitrary polynomials or even a dimension-free transfer of
coefficient perturbations to the scalar supremum. No such transfer is
claimed.

## 2. A useful exact diagnostic, not a solution

For nonzero f define pi_S=|fhat(S)|^2/||f||_2^2. For every positive
degree budget m the exact identity is

    log R_m(f)=log(||f||_2/||f||_infinity)
                  +H_(m/(m+1))(pi)/(2m).

Hence support at most exp(Cm), or ambient dimension at most Cm,
immediately gives a uniform ratio. A divergent family must escape these
simple support bounds and must compensate for its L2-to-supremum loss.
This identity alone is merely a reformulation, not a lower or upper
bound for the unrestricted constants.

The counterexample track now has a finite degree-four example with
ratio strictly above two. It escapes torsion-valued outputs through a
tangent perturbation and actual inter-block degree transfer. Ordinary
tensor amplification cannot turn this into divergence. The next test
is whether degree transfer can operate in a growing ambient/degree
geometry, such as balanced address lifts, without losing its cap or
spending the degree saved in normalization. See the counterexample
artifact for the exact coefficients and current evidentiary status.

## 3. Entropy transfer without an ambient-dimension error factor

**Proved below; independently reconstructed by both other researchers.** This is a
transfer lemma, not an assertion that its hypotheses produce divergence.
Let G be any fixed complex function on a finite cube with |G|=1. Put
pi_S=|Ghat(S)|^2 and H=-sum pi_S log pi_S. Suppose Q_k is a polynomial
on k disjoint copies of that cube, of degree at most D_k, where
D_k/k -> d>0, and

    ||Q_k-G^(tensor k)||_2 <=delta_k,
    ||Q_k||_infinity <=1+gamma_k.

If limsup delta_k<=delta<1 and limsup gamma_k<=gamma, then

    liminf_k R_(D_k)(Q_k)
       >= [(1-delta)/(1+gamma)] exp[H/(2d)].

In particular if both errors tend to zero, no ambient-dimension factor
is needed in transferring the full Shannon-entropy gain. All-degree
conclusions require a separate bounded-gap/padding argument, as supplied
for the explicit normalization construction in the counterexample file.

Proof: fix eta>0. Let a be the Walsh coefficient vector of G^(tensor k)
and restrict it to the usual spectral typical set T_k on which

    |a_S|^2 <=exp[-k(H-eta)].

By the law of large numbers applied to -log pi, its squared mass
s_k=sum_(S in T_k)|a_S|^2 tends to one. Zero probabilities are omitted.
For b=Q_khat, Parseval and Cauchy--Schwarz give

    Re sum_(T_k) bar(a_S)b_S >=s_k-delta_k sqrt(s_k).

The conjugate exponent to q_D is q_D'=2D/(D-1). For D>1,

    ||a 1_T||_(q_D')
      <=||a 1_T||_infinity^(1/D)||a 1_T||_2^(1-1/D)
      <=exp[-k(H-eta)/(2D)] s_k^((1-1/D)/2).

Holder now supplies a lower bound for ||b||_(q_D). Divide by the
supremum bound, take liminf, and send eta to zero. Since d>0,
D_k>1 eventually. No coefficient count or dimension enters the proof.

The same argument is a finite certificate: a unit L2 target with at
least 1-eta of its squared Fourier mass on coefficients of magnitude
at most exp(-h/2) forces a nearby polynomial of degree D to have a
coefficient q_D norm at least

    [(1-eta)-delta sqrt(1-eta)]
       exp[h/(2D)]/(1-eta)^((1-1/D)/2),

provided delta<sqrt(1-eta). Here the exact retained mass may instead
be used throughout; with only a lower bound on mass, the displayed
bound follows because (sqrt(s)-delta)s^(1/(2D)) increases with s.

This sharpens the error bookkeeping of the explicit tensor normalization
construction. It does not manufacture the approximants. Producing a
family with H/d unbounded and controlled cap/error would prove
unbounded BH constants; no such family is presently available.

## 4. Rejected final-exponent bootstrap proposal

The director proposed treating low levels at the final q_m rather than
their own stronger q_r. This correctly avoids the fixed-level Chebyshev
inflation, but the proposed starting estimate contained a false step:
the polynomial-in-m bound on a homogeneous projection's SUPREMUM was
mistaken for a bound on its coefficient l1 norm. A normalized Hadamard
quadratic already has cap at most one and coefficient l1 norm sqrt(N).
The barrier researcher caught this error before it entered any proof.

There is a second exact limitation. Mixed-norm interpolation at q_m
requires child exponent budgets to ADD to m. At a full-degree split
d+e=m, dimension-independent admissibility forces child budgets exactly
d,e. Keeping q_m unchanged on both children produces q_(2m), not q_m.
Thus the proposed contraction does not close. The detailed budget and
noise-cost accounting is in the barrier artifact; this rejects this
particular proof modification, not uniform boundedness itself.

## 5. Fixed-dimensional sphere-valued polynomials are finite juntas

**Proved; the scalar unit-circle case was discovered by the positive
researcher and independently reconstructed by the director and barrier
researcher.** The argument extends to every fixed real target dimension r.
It is qualitative: its Ramsey-size bound does not imply bounded BH
constants as the degree grows. External novelty is unclaimed.

For each m,r there is a finite J(m,r) such that every Walsh polynomial
F:{+-1}^n ->R^r of degree at most m and ||F(x)||_2=1 at every vertex
depends on at most J(m,r) coordinates. For complex scalar outputs take
r=2. For complex target dimension r use real dimension 2r.

Here is an explicit recursive bound. Set J(0,r)=0 and J(1,r)=r.
Choose a finite cover of the unit sphere in R^r by caps such that any
two vectors assigned the same cap have positive inner product. Let t(r)
be the number of caps, and color zero coefficients with one extra color.
For m>=2 set

    L=max(2m, floor(2^(m-1) J(m-1,r))+1),
    J(m,r)=R_m(L; t(r)+1)-1,

where the finite uniform-hypergraph Ramsey number guarantees a
monochromatic L-set in every (t(r)+1)-coloring of the m-subsets.

Proof. Delete irrelevant coordinates and suppose n>J(m,r). Color each
nonzero top-degree coefficient vector by its normalized direction's cap.
A monochromatic nonzero clique cannot contain 2m vertices U: the Fourier
coefficient of ||F||^2 at U is

    sum_(S subset U, |S|=m) <Fhat(S),Fhat(U\S)>.

There are NO omitted lower-degree or external-coordinate terms: the
symmetric difference has size 2m while both factors have degree<=m.
Every displayed summand is strictly positive in a single acute cap,
contradicting ||F||^2=1. Ramsey therefore supplies an L-set A whose
m-subset coefficients are all zero.

Every restriction of the coordinates outside A now has degree<=m-1
on A. No outside factor could have accompanied a top m-subset, by the
original degree bound. By induction every such restriction depends on
at most J(m-1,r) variables.

For each originally relevant i in A choose a nonzero coefficient vector
with i in its support, and project to one nonzero real component. The
coefficient on its intersection with A, after restriction outside A,
is a nonzero polynomial of degree at most m-1 in the outside variables.
A nonzero multilinear polynomial of degree k is nonzero on at least
2^(-k) of the cube (induction on a variable in a top monomial).
Consequently i survives as a relevant variable with probability at
least 2^(1-m). Expected relevance in A is at least L 2^(1-m), exceeding
J(m-1,r), a contradiction. No independence among these survival events
is needed.

For m=1, write F=a_0+sum_i a_i x_i. Constancy of ||F||^2 forces the
nonzero vectors a_i to be pairwise orthogonal, giving at most r of them.
This establishes the induction.

The fixed target dimension is essential. The map
(x_1,...,x_n)/sqrt(n) has degree one and constant Euclidean norm but
target dimension n and depends on every coordinate.

Consequence: the positive track's locally indistinguishable correlated
law obstruction now excludes dimension-independent, bounded-mass
supremum-norm approximation by ALL complex unimodular atoms, not merely
Boolean or torsion atoms, even allowing any finite prescribed degree
inflation. Use K=J(D,2) in the local-marginal bound. This does not exclude
coefficient-norm approximations or boundedness proofs that do not use
that particular endpoint decomposition.

## 6. Phase entropy with a small control block and an arbitrary affine block

**Proved; independently reconstructed by both the barrier and counterexample
researchers.** This is a restricted
positive theorem about the entropy compiler, not a bound for all BH
polynomials and not an assertion about an optimal approximation compiler.
All logarithms are natural. Let

    P(z,x)=p_0(z)+sum_(i=1)^N x_i p_i(z)

be real, of degree at most an integer d>=1, where z has k coordinates and x has N
coordinates. Put Q=||P||_infinity and G=exp(itP). Then

    H(|Ghat|^2)
      <= k log2 + 2|t|Q + 2t^2(d-1)Q^2.                 (A)

In particular, at bounded |t|Q and k=O(d), arbitrary growth of the affine
block N does not make this phase entropy superlinear in d. No Boolean
range or fixed coefficient alphabet is required for the p_i.

Proof. Fourier-transform only the x variables first. Apart from the
irrelevant phase exp(itp_0(z)), the amplitude at a subset S of data
coordinates is

    a_S(z)=i^|S| product_(i in S) sin(tp_i(z))
                      product_(i notin S) cos(tp_i(z)).

For each z, q_S(z)=|a_S(z)|^2 is a product Bernoulli probability law.
Since P is real, maximizing over all x gives
|p_0(z)|+sum_i |p_i(z)|<=Q. The elementary binary entropy bound
h(sin^2 u)<=2|sin u|<=2|u| gives

    E_z H(q(z)) <=2|t|Q.

For completeness h(p)<=2sqrt(p(1-p)) follows by symmetry and the
derivative comparison 2log v<=v-1/v for v>=1. This also implies
h(p)<=2sqrt(p), the version used above.

The full Fourier probability of G, marginalized to S, is
pi_S=E_z q_S(z), by Parseval in z. Hence

    H(pi)=E_z H(q(z))+I(Z;S),

where Z is uniform and the conditional law of S given Z=z is q(z).
Boolean logarithmic Sobolev, applied to each complex amplitude a_S,
implies

    I(Z;S)=sum_S Ent(|a_S|^2)
        <=2 sum_(j,S) E |D_j a_S|^2,

with D_j a=(a(z)-a(z^j))/2. The amplitudes form a unit vector which is
a tensor product of two-dimensional vectors. Their inner product at z
and z^j is product_i cos(t(p_i(z)-p_i(z^j))). Consequently

    sum_S |D_j a_S|^2
       = [1-product_i cos(t(p_i(z)-p_i(z^j)))]/2
       <= t^2 sum_i |D_j p_i|^2.

Here 1-product c_i<=sum(1-c_i) for real c_i in [-1,1], by induction;
and 1-cos u<=u^2/2. Thus no assumption that the cosine factors are
positive was made. Since deg p_i<=d-1,

    I(Z;S)<=2t^2(d-1) sum_i E p_i^2
           <=2t^2(d-1)Q^2.

After the z transform there are at most 2^k control-frequency labels,
so the full Fourier entropy is at most H(pi)+k log2. This proves (A).
For Q=0 the conclusion is immediate. The phase p_0 is absent from the
conditional probabilities; its possible contribution to final control
frequencies is included in k log2.

The logarithmic Sobolev step is the standard hypercube inequality,
equivalently obtained by differentiating the noise hypercontractive
inequality at exponent 2. It holds for complex amplitudes by applying
the real inequality to their moduli and using
|D_j |a||<=|D_j a|. The new feature of (A) is its explicit separation
of the small control block from arbitrarily many affine data variables.
External novelty has not been assessed.

### Scope and a rejected extension

Summing this argument over many affine blocks loses their mutual Fourier
information and can produce a quadratic-in-d bound. It does not prove
the desired universal estimate H(exp(itP))<=C d |t| ||P||_infinity.
The positive track's iterated Boolean selector gives an exact warning:
the sum of coordinate marginal entropies can exceed the true Fourier
entropy by an unbounded factor. One must not replace the joint Fourier
entropy by independent-coordinate information and call the loss harmless.

The proposed universal phase-entropy inequality is itself only a
candidate. Uniform boundedness of BH would imply it, through the actual
Taylor compiler; a family violating it would give BH divergence. Neither
implication furnishes such a family. The quadratic theorem, disjoint
Boolean-block theorem, and (A) are verified/testable subclasses, not a
solution by reformulation.

### 6.1 A stronger all-time bound by rare-output information

**Proved; independently audited by the barrier researcher.** The same family in
fact obeys a linear-in-time estimate without a degree assumption on the
control functions:

    H(|widehat(exp(itP))|^2)
        <= (6+2k log2) |t| ||P||_infinity.                 (B)

Consequently the proposed linear phase-entropy inequality holds uniformly
for all such families with k=O(deg P), at every time t, not just bounded
time. This still leaves families with much richer overlapping control
information uncovered.

Put a=|t|Q. For a<=1, under the partial-Fourier conditional law above,
the event E that S is nonempty has probability at most

    Pr(E|Z=z)<=sum_i sin^2(tp_i(z))<=a^2.

The label S is constant when E fails. The chain rule therefore gives

    I(Z;S)<=h(Pr E)+Pr(E)H(Z|E)
           <=2a+a^2 k log2.                            (C)

The second inequality uses the WORST-CASE bound Pr(E|Z=z)<=a^2,
not only its average. To verify it, the subprobability vector
u_z=Pr(Z=z,E) is coordinatewise bounded by a^2 Pr(Z=z).
The homogeneous entropy S(u)=sum_z u_z log(sum u/u_z) is
coordinatewise increasing, so S(u)<=a^2 H(Z)=a^2 k log2.
Also h(v)<=2sqrt(v) and Pr E<=a^2. This is the same elementary
rare-output estimate underlying the strengthened Fano check in the
barrier artifact, used here in the opposite information direction.

Let T denote the control-frequency label in the FULL Fourier probability.
By Parseval and partial transformation,

    Pr(T!=empty)=E_x Var_z exp(itP(z,x))
               <=E_(x,z)|exp(itP(z,x))-1|^2<=a^2.

Thus H(T)<=2a+a^2 k log2. Combining this with (C) and
E_z H(S|Z=z)<=2a gives

    H(S,T)<=H(S)+H(T)<=6a+2a^2 k log2
                              <=(6+2k log2)a.

For a>=1 use the simpler bounds I(Z;S)<=H(Z)=k log2 and
H(T)<=k log2. They give H(S,T)<=2a+2k log2, which is also
at most the right-hand side of (B). The case a=0 is a constant
phase and has zero entropy.

The conditional physical label Z and the final Fourier label T are
different objects. The proof compares their separate joint laws only
through the common data-frequency marginal S; it does not assume a
simultaneous physical/Fourier measurement of the control variables.

For uniform k-bit controls the bound Pr(E)H(Z|E)<=a^2 k log2 also
follows directly from H(Z|E)<=k log2 and Pr(E)<=a^2; worst-case
conditioning and homogeneous entropy are not needed for this special
case. Their use above is valid but should not be advertised as essential.

### 6.2 Shared control partitions and Boolean inner blocks

**Proved; director extension combining the rare-output argument with the
counterexample track's Boolean-block entropy estimate, independently
audited by the barrier researcher.** Let C_1,...,C_L be a partition of a Boolean control cube into
nonempty cells. Let K be the cardinality of the union of the Walsh
supports of their indicators. Let b_i be nonconstant Boolean-valued
functions on disjoint data blocks, with degree at most D>=1. Suppose

    P(z,x)=p_0(z)+sum_i p_i(z)b_i(x^(i))

is real and every p_i, including p_0, is constant on every C_j. Then,
for Q=||P||_infinity and every real t,

    H(|widehat(exp(itP))|^2)
       <= [2D+4+log L+log K] |t|Q.                      (P)

The number and dimensions of the data blocks, their biases, and the
number of coordinates in the control cube do not occur in this estimate.

Proof. Set a=|t|Q. All b_i attain both signs independently, so
|p_0(z)|+sum_i|p_i(z)|<=Q. Conditional on a control cell, the squared
data-Fourier amplitudes form a product. For one block with u=tp_i(z),
its Fourier distribution is EXACTLY

    cos^2(u) delta_empty+sin^2(u)|bhat_i|^2.

At the empty coefficient the two contributions are perpendicular in
the complex plane, since b_i is real; bias causes no omitted cross term.
Boolean granularity gives H(|bhat_i|^2)<=2(D-1)log2. Entropy of the
mixture is therefore at most2D|u|, by the calculation in counterexample
Section15. Summing gives H(S|cell)<=2Da. Its nonempty-label probability
is at most sum_i sin^2(tp_i(z))<=a^2, on every cell.

For a<=1, the same event decomposition yields
I(cell;S)<=2a+a^2 log L. The control-frequency label T has support
inside the indicated union of size K: for each fixed data input the
function of z is a linear combination of the cell indicators.
Moreover Pr(T!=empty)<=a^2 by the previous Parseval argument. Thus
H(T)<=2a+a^2 log K and

    H(S,T)<=(2D+4)a+a^2(log L+log K).

For a>=1 use I(cell;S)<=log L and H(T)<=log K instead. The two
ranges prove (P). Again, cell and control-Fourier labels belong to
different experiments sharing only the same data-frequency marginal.

For literal k-bit controls take singleton cells: L=K=2^k. This gives
(2D+4+2k log2)a and recovers (B) at D=1. For a common depth-h
decision tree of controls, L<=2^h. Each leaf indicator has at most2^h
Fourier coefficients, so K<=4^h; hence the bound is
(2D+4+3h log2)a. The real polynomial P then has a representation
degree budget h+D. An ACTUAL-degree cancellation must still be
accounted for: this is a uniform linear-in-degree estimate when
D+log L+log K=O(deg P), not an assertion that every low-degree
polynomial has such a small shared control partition.

## 7. All-time quadratic phase entropy versus rotation speed

**Proved; independently reconstructed by the barrier researcher.** Let

    G(x)=exp(i sum_(i<j) theta_ij x_i x_j),
    J=(1/2)sum_i X_i,   v(G)=||M_G^* J M_G-J||_(2->2),

where M_G is pointwise multiplication, X_i flips coordinate i, and the
Hilbert space uses uniform cube measure. Then, with natural logarithms,

    H(|Ghat|^2) <= (152/3) v(G) < 52 v(G).             (Q)

There is no small-time, dimension, sign, or coefficient-size restriction.
This is an entropy inequality, not a bound on every Renyi entropy or on
the unrestricted BH constants. It includes arbitrarily large quadratic
phases, rather than merely a specified Taylor approximation procedure.

### Proof

Reduce every theta_ij modulo pi/2 to [-pi/4,pi/4]. Multiplication by the
removed edge phase is a constant times x_i x_j, so the original function
is a constant times chi_D times the residual phase. Put

    sigma_i=(-1)^(1_(i in D)),
    c_i=product_j cos(2 theta_ij)>=0,   q_i=(1-c_i)/2<=1/2.

The i-th marginal of the ORIGINAL Fourier probability is p_i=q_i if
sigma_i=1 and p_i=1-q_i otherwise. Therefore entropy subadditivity gives

    H(|Ghat|^2)<=sum_i h(q_i).

Moreover, evaluating M_G^*JM_G-J at the all-plus X eigenvector gives
minus the Fourier mean degree. Consequently

    v(G)>=sum_i p_i>=sum_i q_i.                        (Q1)

For a row with q_i<=1/8 we have c_i>=3/4 and all its cosine factors are
positive. Put u_j=tan(2 theta_ij) and T_i=sum_j u_j^2. Then

    4q_i<=T_i<=c_i^(-2)-1<=7/9<1.                     (Q2)

The upper bound follows from product_j(1+u_j^2)=c_i^(-2), and the
lower bound from c_i>=exp(-T_i/2) and 1-exp(-T_i/2)<=T_i/2.

Choose independent signs z_1,...,z_n. The product qubit state having
Bloch coordinates Z_j=z_j/2 and X_j=0 can freely choose each Y_j as
either sign of sqrt(3)/2. Define

    M_i(z)=c_i Im product_j (1+i u_j z_j/2).

Equivalently, for ALL rows, including ones with a zero cosine, use the
nondivided expression

    M_i(z)=Im product_j [cos(2theta_ij)+i(z_j/2)sin(2theta_ij)].

Conjugation of X_i by the residual quadratic phase produces X_i times
a cosine of its neighbor field and Y_i times its sine (the overall
sign depends on the Pauli convention and is immaterial). The character
adds sigma_i. Choose Y_i with sign sigma_i sign M_i and the appropriate
common convention. These choices are simultaneous: changing Y_i does
not change any Z_j or X_j. Since the unconjugated J expectation is zero,

    v(G)>=(sqrt(3)/4)sum_i |M_i(z)|
         >=(1/4)sum_i |M_i(z)|                       (Q3)

for every z. This is the joint product-state witness; it does not bound
different rows by incompatible maximizing states.

For a small row, the imaginary product equals L+R, where
L=(1/2)sum_j u_j z_j and R contains odd Walsh levels at least three.
The elementary fourth-moment Khintchine estimate gives

    E|L|>=sqrt(T_i)/(2sqrt(3))>=sqrt(T_i)/4.

Walsh orthogonality and elementary symmetric sums give

    ||R||_2^2<=sum_(s>=3) (T_i/4)^s/s! <=(T_i/4)^3.

The last inequality holds for T_i/4<=1/4, for example by bounding the
tail by a^3 exp(a)/6<=a^3. Hence

    E|M_i|>=c_i(sqrt(T_i)/4-T_i^(3/2)/8)
            >=3sqrt(T_i)/32>=3sqrt(q_i)/16.

Averaging (Q3) now gives

    v(G)>=(3/64)sum_(q_i<=1/8) sqrt(q_i).             (Q4)

Finally h(q)<=2sqrt(q) for small q, while h(q)<=log2<=1<=8q
for q>1/8. Combining (Q1) and (Q4),

    H<=2sum_small sqrt(q_i)+8sum_large q_i
       <=(128/3+8)v=(152/3)v.

This also covers zero rows, endpoints theta=+-pi/4, and character
gauges with q_i=0. No continuity argument at a zero cosine is needed.

### Exact scope of the amplification obstruction

For any FIXED such G, a sequence Q_k uniformly approximating G^(tensor k)
with error tending to zero and degree D_k must obey

    liminf_k D_k/k>=v(G).

This follows independently from the rotation-speed lower bound proved
in the barrier artifact; no converse is assumed. Thus any limiting
degree rate d available to such a uniform compiler has

    H(|Ghat|^2)/(2d)<=76/3.

In particular a family of quadratic-phase seeds cannot make the
Shannon-entropy transfer LOWER BOUND diverge using fixed-seed uniform
compilation, even with arbitrary times and optimal such compilers.
This does NOT upper-bound the entire BH ratio of an approximation,
does NOT rule out L2-only approximations with controlled supremum, and
does NOT justify an interchange of a growing-seed and approximation
limit. It is a precise obstruction to the specified entropy mechanism.

The proof uses only Pauli conjugation, Fourier entropy subadditivity,
an explicit simultaneous product-state witness, and a fourth-moment
inequality. External novelty has not been established. Numerical tests
in the companion director quadratic-phase script are diagnostics only;
the argument above is independent of them.

The barrier researcher subsequently extended this to arbitrary real
degree-at-most-two P, including linear terms: H(exp(iP))<=64v(exp(iP)).
The director and counterexample researcher each reconstructed that
extension. Its complete proof is in the barrier artifact Section 31;
the diagnostic script also tests those linear terms. The uniform versus
L2 scope distinction remains unchanged.

## 8. The remaining cap-controlled L2 gap: failed attempted extensions

**Historical failed attempts, superseded for fixed quadratic seeds by
Section 9 below.** The entropy-transfer
lemma only needs L2 approximation with a controlled supremum. Section 7
instead controls uniform approximation via rotation velocity. It was
tempting to identify these two resources; the barrier artifact Section 29
now rigorously refutes that identification even for a fixed Boolean seed.
For b=1-2 AND_16, v(b)=4, whereas cap-one polynomials approximating k
copies in L2 have degree rate at most 1/sqrt(2). The construction uses
exact rare-mark enumeration and the ordinary polynomial method. Its
Fourier entropy is also small, so it is NOT a BH counterexample.

There is a real positive exception: the counterexample artifact Section 17
proves H(G)<=64d_eff for independent linear phases G, even for cap-bounded
L2 approximants. Complex evaluation equals correlation with the product
phase, and a univariate Bernstein--Walsh inequality charges degree for
each dyadic angle band. That proof, including conditional fixing of
unused coordinates, was independently reconstructed by the director.

The director tested three possible extensions to interacting quadratic
phases, without obtaining a theorem:

1. Replace velocity witnesses by their L2 matrix elements. The product
   witness need not have bounded density relative to the uniform cube;
   a vanishing uniform-cube L2 error need not control its expectation.
   The rare-AND example proves that some such loss is necessary in general.
2. Restrict to independent vertex sets, obtaining linear phases. This
   can lose a graph-coloring factor, unbounded on dense graphs. It also
   discards conditional Fourier information; no universal entropy recovery
   inequality from those restrictions was proved.
3. Factor a quadratic phase into independent edge phases. These factors
   share coordinates. Their Fourier labels combine by symmetric difference,
   so the exact independent-phase complex-evaluation identity no longer
   applies. Treating edge labels as independent would change the cap and
   the degree model. No such replacement is licensed.

One discriminating future question is therefore whether all quadratic
phases obey H(G)<=C d_eff for cap-controlled L2 tensor approximation,
or whether a growing family violates it. The known velocity estimate
alone does not answer it. This is a genuine restricted question, not an
equivalence asserted with full BH boundedness. A failure with unbounded
H/d_eff and bounded approximation caps would, by Section 3, give a
diverging BH family; proving it would close that particular amplification
class only. The unrestricted Renyi coefficient problem remains broader.

## 9. A flat-diagonal certificate closes the quadratic L2 gap

**Proved jointly; independently reconstructed by all three research
tracks and the director.** This section changes the assessment in
Section 8. It closes a specific amplification class, not BH boundedness.

For any unit-circle-valued function G on n signs define

    V(G)=sup_rho |Tr rho (M_G^* J M_G-J)|,
    J=(1/2)sum_i X_i,

where rho ranges over positive trace-one matrices whose diagonal in the
physical cube basis is exactly uniform, rho_(xx)=2^(-n). Ordinary matrix
trace is used. This constraint is essential: it makes the norm of a
diagonal multiplication error in the rho-inner product equal to its
ordinary uniform-cube L2 norm. The full operator velocity lacks this
property and has the explicit rare-AND counterexample in Section 8.

### 9.1 Exact finite lower bound on cap-controlled L2 degree

Let Q be a degree-at-most-D polynomial on k disjoint copies of the seed
cube. Put C=||Q||_infinity>0 and delta=||Q-G^(tensor k)||_2. Then

    D/k >= [V(G)-n delta]/C.                         (F)

No limiting interchange, approximation algorithm, or conjecture enters
this finite inequality. Its usefulness for a growing seed depends on
controlling n delta; for each fixed seed this term vanishes with delta.

Proof. Fix an admissible rho and use rho^(tensor k), still flat on its
diagonal. Let U=M_(G^(tensor k)), K=sum of the k copies of J, and E=U-M_Q.
Tensor factorization gives

    Tr rho^(tensor k) U^*[K,U]
       =k Tr rho (M_G^* J M_G-J).

Expand [K,U]=[K,M_Q]+K E-E K. The first term has absolute trace at most
||[K,M_Q]||. Simultaneous spin rotation and the Banach-valued trigonometric
Bernstein inequality give ||[K,M_Q]||<=D C: the rotated multiplier is a
trigonometric polynomial of degree at most D and constant operator norm C.

For the first error term, Cauchy--Schwarz in the rho^(tensor k) inner
product gives

    |Tr rho^(tensor k) U^* K E|<=||K|| delta.

For the second use the same inequality with A=E^*U and B=K. Because
E and U are commuting diagonal operators and U is unitary, A^*A=|E|^2.
Flatness of the density diagonal therefore gives the SAME error norm
delta, and

    |Tr rho^(tensor k) U^* E K|<=||K|| delta.

Since ||K||=nk/2, these inequalities prove
k|Tr rho (M_G^*JM_G-J)|<=DC+nk delta. Take the supremum over rho to get (F).

In particular, if delta_k->0 and limsup C_k<=C0<infinity, then

    liminf_k D_k/k >=V(G)/C0.                        (F-limit)

The important special case C_k->1 has no cap factor. The more elaborate
centered-energy proof first found by the barrier researcher is retained
in that artifact's exploration notes; the direct commutator proof above
is sufficient and avoids its unnecessary second-moment bookkeeping.

### 9.2 The quadratic entropy witnesses have the required diagonal

All density matrices actually used in Section 7 and its linear-term
extension are flat after the stated averaging:

* The all-plus X eigenstate has uniform physical diagonal and proves
  V(G)>=sum_i p_i>=sum_i q_i.
* For the product-state witness, Z_i=z_i/2 and X_i=0 while Y_i has the
  optimally chosen sign, which may depend on the ENTIRE sign vector z.
  Its diagonal at a physical vertex x is
  product_i (1+x_i z_i/2)/2. The Y signs do not appear. Averaging over
  independent uniform z makes this diagonal exactly 2^(-n).

Its averaged energy is the same positive sum of absolute row responses
used in Section 7. Thus the full entropy proof works with V, not merely
the unrestricted operator norm v. We have, for every real degree-at-most-two P,

    H(|widehat(exp(iP))|^2)<=64 V(exp(iP)).           (QF)

The hollow case retains constant 152/3. Combining (QF) with (F-limit),
every fixed quadratic-phase seed and cap-one-asymptotic L2 compiler obey

    H(G)<=64 liminf_k D_k/k.

Therefore the Shannon-entropy transfer lower bound exp[H(G)/(2d)] is
at most exp(32) throughout this class, even for L2-only compilers.
This does not upper-bound arbitrary additional Fourier tails of Q_k.
It does not show boundedness for all polynomials. It does close the
previously live quadratic-phase entropy-amplification loophole.

### 9.3 Validation and limitations of the certificate

For a real Boolean-valued b, V(b)=I(b) exactly. The operator M_b^*JM_b-J
has entry -1 on each directed sensitive cube edge and zero elsewhere.
Positivity and the flat diagonal force |rho_(xy)|<=2^(-n), so its trace
has magnitude at most I(b); the all-plus state attains that bound.
Thus this certificate does not contradict the rare-AND example, and a
universal inequality H(G)<=constant V(G) would INCLUDE the unresolved
Boolean Fourier-entropy/influence conjecture. No such general theorem
is asserted.

For independent phases G=product_j exp(i theta_j x_j), one has the exact
formula V(G)=sum_j |sin theta_j|. The operator is a sum of disjoint 2x2
blocks with eigenvalues +- |sin theta_j|, each with a uniform-diagonal
extremal eigenvector; their product is an admissible witness. Hence
H(G)=sum_j h(sin^2 theta_j)<=2V(G). This yields a sharper constant than
the separate dyadic complex-evaluation proof when delta_k->0 and C_k->1.
That separate proof still has its distinct strength: it allows error
bounded strictly below one and subexponential, rather than bounded, caps.

The certificate is additive under disjoint tensor products. Partial
traces of a flat-diagonal density remain flat, giving the upper bound.
For the reverse bound choose product optimizing densities with positive
energy. A sign reversal of the energy is always available by conjugating
with the full parity multiplier, which anticommutes with J and commutes
with M_G while preserving the diagonal. Thus no incompatible sign
choice is hidden in the additivity statement.

This is a reusable theorem connecting approximation, a concrete state
constraint, and degree. It remains a class-specific obstruction to one
counterexample mechanism, not a claimed resolution of the campaign target.

### 9.4 Stronger conclusion when the ambient tensor space is fixed

**Proved; independently audited by the barrier researcher.** For a FIXED
nonconstant quadratic-phase seed G on n bits, let Q_k live on EXACTLY the
nk tensor coordinates, with L2 error tending to zero and cap tending to
one. Take D_k to be its actual degree. Then

    limsup_k R_(D_k)(Q_k)<=exp(32).                   (AR)

This bounds the actual ratios, not merely the value delivered by the
entropy-transfer lower bound. Indeed V(G)>=I(G)>0, so (F) implies
D_k>=(V(G)/2)k eventually. The full counting-measure comparison gives

    ||widehat(Q_k-G^(tensor k))||_(q_Dk)
       <=2^(nk/(2D_k))delta_k<=2^(n/V(G))delta_k->0.

Every subsequence has a further subsequence with D_k/k->d in [V(G),n].
Since the seed Fourier distribution is fixed and finite, its Renyi
entropy at D_k/(D_k+1) tends to its Shannon entropy. Consequently the
actual ratio along this subsequence tends to exp[H(G)/(2d)], which is
at most exp(32). This proves (AR).

For an independent linear phase seed the same argument gives the
sharper upper bound exp(1), using H<=2V. The constants are not optimized.
No extra variables may be silently introduced: a tiny L2 error on a
much larger auxiliary cube can have a large critical coefficient norm.
Likewise, the error factor 2^(n/V(G)) is fixed-seed, not uniform in a
simultaneously growing seed. The earlier conservative statements that
did not control extra Fourier tails are refined by (AR) precisely under
these additional, explicit ambient-space and nonconstancy hypotheses.
