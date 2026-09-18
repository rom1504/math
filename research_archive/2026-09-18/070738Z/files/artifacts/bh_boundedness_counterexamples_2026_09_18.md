# Boolean BH boundedness: counterexample track

2026-09-18, separate campaign 04:27:50--07:27:50UTC. The old signing
problem is paused. This file distinguishes reconstructed literature,
proved failure classes, exploratory calculations, and unresolved routes.
No divergence of the Boolean BH constants is claimed.

Write q_m=2m/(m+1), and

    R_m(f)=||fhat||_(q_m)/||f||_infinity,  deg(f)<=m.

All Fourier coefficients are normalized Walsh coefficients on the
probability cube. Degrees below are the ACTUAL reduced Walsh degrees,
unless a degree budget is explicitly named.

Main proved outcomes (complex scalars):

* Section7 gives the explicit finite certificate B_4>2+2^(-18).
* Section8 gives actual polynomials at every sufficiently large degree
  and proves liminf_m B_m>2+2^(-207). This is a persistent constant
  gap, NOT divergence. Both independent research tracks audited it.
* Section11 supplies the general entropy/effective-degree transfer;
  no family with an unbounded entropy/rate ratio is produced.
* Sections10,12,14 rule out the tested addressed tangent, orthogonal
  coupling, and uncompressed inheritance mechanisms. Sections13 and15
  give rigorous phase-Hamiltonian failure classes, including a
  moving-seed optimal-uniform-compiler bound for small quadratic phases.
* Section17 independently rules out fixed independent-coordinate phase
  targets even under cap-controlled L2 compilation. The cross-audited
  flat-diagonal state theorem in the barrier artifact extends this
  conclusion to every fixed real quadratic phase at arbitrary time.
* Section18 gives an exact rare-marker spectral identity and entropy
  bound. Its compiler comparison is deliberately limited to the
  displayed paid budget, not an unproved optimal-degree lower bound.

The unrestricted boundedness question remains unresolved. Statements
about real scalar constants or external novelty are not inferred.

Final replay checkpoint, 2026-09-18 06:47UTC: all eight deterministic
scripts in this track passed again (`checks`, `tangent`, `normalization`,
`address_tangents`, `skew_coupling`, `inheritance`, `hamiltonian`, and
`phase_l2_compiler`). Numerical searches are excluded from this proof
checkpoint. Exact lower-bound certificates and analytic estimates,
not a floating-point near-unitarity test, support the stated theorems.

## 1. Primary-source baseline and Boolean-valued obstruction

The three September16 audit/synthesis files were read before starting.
The polynomial upper theorem does not settle boundedness. Its weighted
proof losses must not be treated as lower bounds for R_m.

The exact Boolean-valued result is already known: Arunachalam et al.,
[*A cb-Bohnenblust--Hille inequality with constant one and its applications
in learning theory*](https://ir.cwi.nl/pub/35177/35177.pdf), Section4.1,
Lemma4.1 and Proposition4.2. The complete decisive proof was read directly.
If f takes values in {-1,1} and has degree<=m, then

    R_m(f)<=2^((m-1)/m).

An address function attains equality. Reconstruction: after writing
g(z)=(1-f(1-2z))/2 on {0,1}^n, its multilinear monomial coefficients
are integers by Mobius inversion. Its degree is<=m, so every Walsh
coefficient of f is a multiple of 2^(1-m). Parseval implies at most
4^(m-1) nonzero coefficients. Counting-measure Holder yields the result.
This rules out ALL Boolean-valued bent/addressing/recursive constructions,
not merely the usual examples. This statement is imported, not new here.

## 2. Exact entropy bookkeeping for tensors

Let f be nonzero, d=deg(f)>=1, pi_S=|fhat(S)|^2/||f||_2^2, and
alpha_m=m/(m+1). Direct algebra gives, for every degree budget m,

    R_m(f)=(||f||_2/||f||_infinity)
             exp[H_(alpha_m)(pi)/(2m)].                 (1)

For disjoint tensor powers f^[tensor k], degree kd and coefficient
probability law pi^[tensor k], this becomes EXACTLY

    R_(kd)(f^[tensor k])
      =(||f||_2/||f||_infinity)^k
          exp[H_(alpha_(kd))(pi)/(2d)].                 (2)

Thus any fixed non-unimodular seed has ratio tending to ZERO. For a
unimodular seed it tends to exp[H(pi)/(2d)], a finite constant; indeed
the displayed Renyi factor decreases with k. A fixed-degree advantage
cannot be exponentiated while retaining its old q_d exponent.

### 2.1 Arbitrary injective binary recoding still cannot amplify a seed

This extends(2), and is proved here. Translate the spectral support to
have full affine span F_2^r, where r>=1. Define

    epsilon=min_(nonzero ell in F_2^r)
               min{P_pi(ell.S=0),P_pi(ell.S=1)}>0.

Full affine span makes epsilon positive. Under pi^[tensor k], every
nonconstant affine parity has both probabilities at least epsilon:
its bias is, up to sign, a product of nontrivial seed-parity biases,
each of absolute value at most 1-2epsilon.

Now apply ANY injective affine map to all spectral indices, into an
arbitrarily larger physical coordinate space. The resulting coefficient
list is only relabeled, and the cube supremum norm is unchanged (the
dual linear map is surjective). The physical coordinate rows have rank
kr, so at least kr coordinate bits are nonconstant. Consequently their
expected Hamming weight, and therefore the maximum supported weight D_k,
is at least epsilon kr. Formula(1) then gives

    R_(D_k)(recoded tensor)
      <=(||f||_2/||f||_infinity)^k
          exp[H_(alpha_(D_k))(pi)/(2epsilon r)].         (3)

The right side is bounded uniformly in k, and tends to zero for a
non-unimodular seed. Translation by a character is included. An affine
span of rank zero is a single coefficient and trivially has ratio one.

IMPORTANT exclusion: a noninjective quotient may merge Walsh indices.
Then both coefficient probabilities and the supremum norm change, and
(3) is not asserted. Such collisions require a fresh exact calculation.

### 2.2 Even changing seeds cannot help by ordinary tensoring

There is a stronger no-amplification statement before recoding. Let the
nonzero factors f_j live on disjoint cubes, with positive degree budgets
d_j, and D=sum d_j. Complex interpolation of counting-measure norms gives

    ||fhat_j||_(q_D)
       <=||fhat_j||_(q_(d_j))^(d_j/D)
           ||fhat_j||_2^(1-d_j/D),

because 1/q_D=(d_j/D)/q_(d_j)+(1-d_j/D)/2. Tensor coefficients and
supremum norms multiply, while ||f_j||_2<=||f_j||_infinity. Hence

    R_D(tensor_j f_j)
       <=product_j R_(d_j)(f_j)^(d_j/D)
       <=max_j R_(d_j)(f_j).                           (8)

The factors may vary with D and have arbitrarily increasing dimension.
This is not just the fixed-seed limit in(2). Constants can be deleted
from the tensor and do not affect the argument.

### 2.3 Selector recursion has a uniform closure bound

Let z be a fresh bit and put

    f=(1+z)g/2+(1-z)h/2,   D=1+max(deg g,deg h).

The two branches may share all their variables. Their Fourier
coefficients combine exactly as (g_A+h_A)/2 and (g_A-h_A)/2. For
1<=q<=2, concavity of t^(q/2) and the parallelogram identity give

    |a+b|^q+|a-b|^q<=2(|a|^q+|b|^q).

Thus ||fhat||_q<=2^(2/q-1)max(||ghat||_q,||hhat||_q).
Normalize max(||g||_infinity,||h||_infinity)=||f||_infinity=1.
If both branch BH ratios at their degree budgets are at most C>=1,
the same interpolation as above gives

    R_D(f)<=2^(1/D) C^((D-1)/D)<=max(2,C).            (9)

Constants as branches have coefficient norm<=1. Equations(8)--(9)
therefore show that an ARBITRARY tensor/selector formula cannot produce
an unbounded ratio from primitives with bounded ratios, evaluated at
its recursively assigned additive/max-plus-one degree budget.
All branch coefficient collisions in(9) have been retained.

The important qualification is not cosmetic: if the actual degree
falls substantially BELOW the formula budget through cancellation,
q_actual is smaller and(9) does not control its BH ratio. Such a degree
collapse is a necessary new ingredient for this amplification route.

## 3. Growing phases: three scalable failure classes

These elementary propositions apply to genuinely complex outputs; they
do not follow merely by calling them Boolean-valued functions.

### 3.1 Arbitrary complex selector leaves, with depth paid

A binary decision tree of depth at most d, with arbitrary complex leaf
values of modulus<=1, represents a polynomial of degree<=d. A path of
length l contributes a product of l literals, with at most 2^l Walsh
terms. There are at most 2^d leaves, hence at most 4^d Fourier terms.
Parseval and Holder therefore give

    ||fhat||_(q_d)<=2||f||_infinity.                    (4)

The same conclusion holds for any explicitly constructed polynomial
whose support count is at most 4^d and whose degree budget is d. Growing
or even continuous phase alphabets at the leaves do not alter(4).
If the ACTUAL polynomial degree is much smaller than the decision-tree
depth, this argument cannot replace d by that smaller degree.

### 3.2 Full affine Fourier supports

Suppose the nonzero Fourier support is an affine subspace a+V of rank r,
and its largest Hamming weight is m. Choose r independent physical
coordinate functionals on V. Their joint projection maps a+V ONTO
F_2^r, so some supported word is one on all r coordinates. Thus r<=m and

    R_m(f)<=|a+V|^(1/(2m))<=sqrt(2).                  (5)

No flatness of the coefficients is needed. In particular, bent and
quadratic phase spectra, whenever their full affine support description
applies, cannot escape through an injective linear code or a character
translation. Mere CONTAINMENT in an affine subspace is insufficient.

### 3.3 Generic coded products of phases

Fix arbitrary characters chi_(v_1),...,chi_(v_t), and set

    f_theta(x)=product_j[cos(theta_j)+i sin(theta_j)chi_(v_j)(x)].

This is exactly unimodular. Away from cos(theta_j)=0, the coefficient at
v is the common cosine product times

    sum_(J: xor_(j in J)v_j=v) i^|J| product_(j in J)tan(theta_j).

For every v in span{v_j}, this is a NONZERO polynomial in the independent
tangent variables: distinct J give distinct monomials, even when the
Walsh character indices collide. Thus outside a finite union of proper
algebraic zero sets, the Fourier support is the ENTIRE linear span.
By(5), its actual-degree BH ratio is at mostsqrt(2). This is a generic-parameter
theorem for every circuit size and code, not a probabilistic numerical
observation. Fine-tuned angles lying on the exceptional zero sets remain
unexcluded, unless another result (e.g. fixed root alphabet) applies.

## 4. Precision and fixed algebraic alphabet bounds

If f is bounded by1, has degree<=m, and all its values belong to
delta(Z+iZ), then Mobius inversion on {0,1}^n gives

    fhat(S) in delta 2^(-m)(Z+iZ).

Every nonzero coefficient has magnitude at least delta2^(-m), so
Parseval and Holder imply

    R_m(f)<=2 delta^(-1/m).                            (6)

The scalar delta is part of the hypothesis. A diverging sequence on a
uniform value grid must therefore have log(1/delta)/m tending to infinity
along a subsequence. Four-phase-valued functions have delta=1, giving2.
This does not control arbitrary irrational circle-valued outputs.

The following initial pointwise algebraic estimate is valid but is
strictly superseded by the Galois--Holder theorem below. Let all values
be Kth roots of unity. Put D=phi(K).
For K>=3, every nonzero coefficient a has 2^m a in Z[zeta_K]. Every
Galois conjugate is the corresponding Fourier coefficient of another
unit-modulus degree<=m function, hence has modulus<=1. The algebraic
norm is a nonzero integer. Pairing conjugate complex embeddings yields

    |a|>=2^(-mD/2),   and consequently R_m(f)<=2^(D/2). (7)

Indeed the norm of 2^m a is at most 2^(mD)|a|^2 because all other
normalized conjugates have modulus<=1. Thus each FIXED root alphabet
is excluded from a divergence construction. The bound deteriorates with
K. For K=2, Section1 is stronger.

The barrier track then supplied a MUCH stronger uniform argument,
independently reconstructed here. Every Galois conjugate f_sigma has
the SAME Fourier support and Parseval norm1. Nonzero a_S satisfies

    product_sigma |a_(S,sigma)|^2 >= 2^(-2mD).

Generalized Holder, now SUMMED over coefficients before discarding
conjugate information, yields

    #supp(fhat) 2^(-2m)
      <=sum_S product_sigma |a_(S,sigma)|^(2/D)
      <=product_sigma (sum_S |a_(S,sigma)|^2)^(1/D)=1.

Thus root-of-unity-valued degree<=m functions have Fourier support at
most4^m and R_m<=2, UNIFORMLY over every growing root alphabet.
See `bh_boundedness_barrier_2026_09_18.md` for the canonical proof and
its limitations. Arbitrary circle-valued functions cannot simply be
approximated by torsion-valued functions while preserving degree; the
director/barrier track gives an explicit failure of that density claim.

## 5. Exact sixth-root extremizer: ratio TWO already at degree TWO

This finite construction was found in this track and then certified
exactly, independently replayed by the barrier track. No external
novelty claim is made. It is an actual lower example, NOT divergence.

Let zeta=exp(i pi/3), with zeta^2=zeta-1. On five Boolean variables,
identify subsets with binary masks, lowest bit corresponding to x_1.
The polynomial is

    f(x)=1/4 sum_(s in M) zeta^(k_s) x^s,

with the complete table

    M   =[0,1,2,3,4,5,6,8,9,10,12,16,17,18,20,24],
    k_s =[0,2,0,5,1,4,2,4,3, 1, 4, 1, 2, 2, 3, 2].

M is EXACTLY the Hamming ball of radius2 in five coordinates, with
16 entries. Every coefficient has modulus1/4 and degree is exactly2.
Here is a finite proof of its supremum, with no floating-point premise.
In the integer basis (1,zeta), the sixth roots are

    (1,0),(0,1),(-1,1),(-1,0),(0,-1),(1,-1).

The32 Walsh sums 4f(x), in binary input order0,...,31, are

    (-4,4),(0,4),(-4,0),(4,-4),(0,4),(4,0),(-4,4),(4,-4),
    (0,4),(0,4),(0,4),(4,0),(4,0),(4,-4),(0,4),(4,-4),
    (4,-4),(4,0),(0,-4),(4,-4),(4,-4),(4,-4),(-4,0),(0,-4),
    (4,0),(0,4),(0,4),(0,4),(4,-4),(0,-4),(-4,4),(-4,0).

Every displayed pair is4 times a sixth root. Thus |f(x)|=1 exactly,
and

    R_2(f)=16^(3/4)/4=2.                             (10)

This is genuinely beyond the Boolean-valued bound sqrt(2) at degree2.
Its Fourier support is NOT affine, so it does not conflict with(5).
It also saturates the uniform torsion theorem, showing its bound2
cannot be improved even for the fixed alphabet of sixth roots.

Taking k disjoint copies has degree2k, support16^k, coefficients of
modulus4^(-k), supremum1, and R_(2k)=2 EXACTLY. Thus all even-degree
constants are at least2. This tensor construction deliberately does
not claim an increasing ratio. Adding one fresh character gives the
odd-degree lower bound2^((m-1)/m), tending to2.

Any injective recoding of the tensor's frequencies still has16^k
nonzero coefficients and remains sixth-root-valued, so the torsion
support bound forces its degree to be at least2k. For THIS exact
extremizing sequence, not even a constant-factor injective degree
compression is possible. Noninjective substitutions remain sixth-root-
valued as well and have ratio<=2 after all collisions are reduced.

## 6. Current boundary

The still-open counterexample directions after these filters are not
ordinary fixed-seed tensoring, generic phase products, Boolean outputs,
or depth-paid selectors. The all-torsion theorem also excludes ANY
root-valued construction, including growing phase alphabets and
noninjective coded substitutions. Possible escape requires genuinely
non-torsion values or non-unimodular functions, with degree collapse or
another mechanism beyond the closure theorems above. These are research
targets, not asserted mechanisms or evidence of divergence.

## 7. A genuine nontorsion escape: B_4 > 2 + 2^(-18)

The constant2 from the root-of-unity theorem does NOT extend to all
bounded complex polynomials. This track proves the following explicit
finite lower bound, independently audited by the barrier track:

    B_4 >= R_4(P) > 2 + 2^(-18).                     (11)

Here P is a degree4 polynomial on10 Boolean variables, given below.
This is not a divergence theorem, nor an external novelty claim. Its
values have two different moduli; this is exactly where it leaves the
unimodular/torsion class. The construction uses a genuine degree
transfer between two blocks, rather than ordinary tensor amplification.

### 7.1 An exact local degree-transfer pair

Use f and zeta from Section5, and put on its five-variable cube

    a=conjugate(zeta)+zeta x_3+x_4,
    b=f^2 conjugate(a).                               (12)

Although the displayed product for b has nominal degree5, its reduced
Walsh degree is EXACTLY3. The complete nonzero coefficient table for
4b, written as integer pairs u+v zeta, is

    mask : (u,v)
       0 : ( 1,-1)       3 : (-1,-1)
       4 : (-1, 0)       7 : ( 1,-2)
       8 : ( 0,-1)      11 : ( 0,-3)
      17 : (-1,-1)      18 : (-3, 3)
      21 : ( 3, 0)      22 : ( 1,-2)
      25 : ( 2,-1)      26 : (-2, 1).

This table follows by XOR convolution of the finite coefficient list
in Section5, using zeta^2=zeta-1; it is also checked independently by
the barrier replay. Define the antilinear involution J_f u=f^2 bar(u).
Because |f|=1, J_f a=b and J_f b=a exactly.

For two disjoint five-variable blocks x,y, set

    F(x,y)=f(x)f(y),
    h(x,y)=a(x)b(y)-b(x)a(y),
    P=F+2^(-10)h.                                    (13)

Both F and h have degree<=4. Indeed the two h terms have degree1+3,
and no unreduced product degree is substituted for that fact.
Since J_F h=-h, one has Re(bar(F)h)=0 pointwise.

There is also an exact cap calculation. The quantity u=a/f takes
values in {0} union {2 zeta^j:j=0,...,5}. This follows directly from
the Section5 value table and a's four possible values. Therefore

    h/F = u(x)bar(u(y))-bar(u(x))u(y)
          belongs to {0,4i sqrt(3),-4i sqrt(3)}.

The three values occur respectively640,192,192 times among the1024
inputs. In particular, for every real epsilon,

    ||F+epsilon h||_infinity=sqrt(1+48epsilon^2).      (14)

### 7.2 The new coefficient mass is exactly paid

The Fourier support S of F has256 terms, all of modulus1/16.
The polynomial h has60 nonzero coefficients. Exactly36 are outside S:
24 have modulus sqrt(3)/4 and12 have modulus3/4. Thus

    sum_(outside S) |hhat|^2 = 45/4,
    sum_(outside S) |hhat|^(8/5) >= 45/4.             (15)

The latter follows since all these magnitudes are at most1.
For example, the coefficient at block masks(0,21) is3 bar(zeta)/4,
so the existence of genuinely new Fourier terms is not a counting
artifact or numerical support threshold.

Let q=8/5 and A=||Fhat||_q^q=2^q. Convexity on the OLD support gives

    sum_(s in S)|Fhat(s)+epsilon hhat(s)|^q
       >= A + q(1/16)^(q-2)epsilon
                    Re sum_(s in S)bar(Fhat(s))hhat(s)
       = A.

The linear term is EXACTLY zero by Parseval and Re(bar(F)h)=0.
This retains every old coefficient and does not discard a harmful
first-order term. Adding the new support and applying(14)--(15),

    R_4(F+epsilon h)^q
      >= [A+(45/4)|epsilon|^q]/(1+48epsilon^2)^(4/5). (16)

The positive new-support exponent q<2 beats the quadratic cap cost.
At the explicit rational epsilon=2^(-10), epsilon^q=2^(-16), and
concavity of t^(4/5) gives the entirely rational denominator bound

    (1+48epsilon^2)^(4/5)<=1+3/81920.

Consequently

    R_4(P)^q-A
       >= [45/262144-A(3/81920)]/(1+3/81920)
       >= 33/(16*81923) > 2^(-16),                   (17)

where only A<=4 was used in the second inequality. Since
A+2^(-16)<4 and the derivative of t^(5/8) is greater than1/4 on
(0,4], taking the5/8 power in(17) proves(11).

The decimal R_4(P)=2.0000467013... is a diagnostic only; the rational
certificate proves the weaker explicit improvement in(11). A numerical
one-parameter optimization gives about2.00015338, but is NOT used.

### 7.3 What this does and does not settle

This disproves the conjectural extension `all complex B_m<=2`, if
one proposes that particular bound. It does NOT disprove boundedness
by a larger absolute constant. Ordinary tensoring of P cannot amplify
its ratio by Section2.2, and fixed-seed tensor powers eventually lose
through ||P||_2/||P||_infinity<1.

An attempted sum of many such tangent interactions across many seed
blocks introduces a separate cap-versus-variance cost. No estimate in
this file turns that many-block construction into a diverging sequence.
Thus the full boundedness question remains explicitly unresolved.

### 7.4 Exact normalization and recoding do not amplify this example

The positive track suggested normalizing P pointwise to the unit circle.
This can be tested without approximating a reciprocal square root.
Write R=-h^2/F^2, whose values are exactly0 and48, and set

    alpha=((1+48epsilon^2)^(-1/2)-1)/48.

Then the exactly unimodular normalization is the finite polynomial

    (F+epsilon h)(1+alpha R)
       = F+epsilon(1+48alpha)h+alpha C,
    C=FR=-h^2 bar(F).

The integer replay proves that C has751 nonzero Fourier coefficients,
degree EXACTLY8, and coefficient3/4 at the10-bit mask503, which has
weight8. Since F,h have degree4, this coefficient cannot cancel when
epsilon is nonzero. Thus this PARTICULAR exact radial normalization
doubles the degree; it is not a degree4 unimodular counterexample.

There are two further rigorous scalability guards. First, a construction
on k disjoint five-bit seed blocks with degree2k always satisfies

    R_(2k) <= (2^(5k))^(1/(4k)) = 2^(5/4)

by Parseval and the total ambient Fourier support bound. Therefore
arbitrarily many balanced degree-transfer terms on this fixed block
geometry cannot themselves produce divergence.

Second, P contains every coefficient of the original256-term root-valued
F. Indeed on that support |Fhat|=1/16 and
epsilon |hhat|<=epsilon||h||_2<1/16 at epsilon=2^(-10).
For k tensor copies, an injective parity recoding therefore contains the
recoded support of F^tensor k, of size256^k. That latter function remains
root-valued, so the Galois support theorem forces its degree to be at
least4k. Consequently no injective parity recoding can compress P^tensor k
below degree4k. Deliberate coefficient cancellations or noninjective
substitutions are not covered by this last argument.

## 8. Degree-efficient normalization gives a persistent gap above2

There is nevertheless a different, scalable use of the same example.
The following is stronger than an isolated finite-degree lower bound:

    liminf_(m->infinity) B_m > 2+2^(-207).           (18)

This remains a CONSTANT lower bound, not divergence. It uses nonunitary
polynomial approximants and pays all added degrees. Both other campaign
tracks independently checked the normalization mechanism; the elementary
constants below make the seed choice explicit.

### 8.1 Explicit seed and its entropy gain

Retain F,h from Section7, and fix epsilon=2^(-104). Put

    R=|h|^2/48 in {0,1},
    lambda=(1/2)log(1+48epsilon^2),
    P=F+epsilon h,
    G=P exp(-lambda R).

Here deg(P)<=4, deg(R)<=8, and |G|=1 EXACTLY. Let pi_s=|Ghat(s)|^2
and H(G)=-sum pi_s log(pi_s), with natural logarithms. Then

    H(G) >= 8log2 + epsilon^2[(45/4)log(1/epsilon)-256]. (19)

Here is a finite proof, not an asymptotic assertion with an unspecified
smallness constant. Write u=(1+48epsilon^2)^(-1/2), C=48FR, and

    G=F+epsilon h+w,
    w=(u-1)epsilon h+(u-1)C/48.

From the exact value histogram, ||h||_2^2=18 and ||C||_2^2=864.
Since 1-u<=24epsilon^2, for every0<epsilon<=2^(-10),

    ||w||_2<=16epsilon^2,
    ||G-F||_2<= (9/2)epsilon.                        (20)

Each coefficient obeys the same bound as its coefficient l2 norm.
On the old256-term support S, |Fhat|=1/16 and the probabilities
p_s=|Ghat(s)|^2 are at least(1/2)(1/256). Taylor's inequality for
g(p)=-p log p, whose second derivative is -1/p, gives

    sum_S g(p_s)
      >= log256 +(log256-1)sum_S(p_s-1/256)
                     -256 sum_S(p_s-1/256)^2.

The first sum is minus the outside mass, of magnitude at most
(81/4)epsilon^2. Also

    [sum_S(p_s-1/256)^2]^(1/2)
       <= (1/8+||G-F||_2)||G-F||_2
       <= (153/256)epsilon.

Using log256-1<=7, the total negative correction is less than
256epsilon^2. On each of the36 originally new coefficients,
1/4<=|hhat|<=3/4. Equation(20) therefore gives

    |Ghat|^2 >= (1-128epsilon)epsilon^2 |hhat|^2,
    |Ghat|<=epsilon.

Their entropy contribution is at least
2(1-128epsilon)(45/4)epsilon^2 log(1/epsilon), which is at least
(45/4)epsilon^2 log(1/epsilon). All remaining entropy summands are
nonnegative. This proves(19).

### 8.2 Polynomial radial normalization with the degree counted

For k disjoint copies of the10-variable block, set J=sum_(i=1)^k R_i
and P_k=product_(i=1)^k P_i. Define

    L_k=ceil(4lambda k),
    T_k(J)=sum_(j=0)^(L_k-1) (-lambda J)^j/j!,
    Q_k=P_k T_k(J),
    M_k=4k+4(L_k-1).

The smaller added degree is essential and is exact: R_i^2=R_i, so
J^j reduces to products of at most j distinct R_i. Moreover h_iR_i=h_i,
and P_iR_i=C_i/48+epsilon h_i has degree at most8. Each selected block
therefore costs only4 degrees beyond the unselected block's4. Thus
Q_k is an explicit polynomial on10k signs of degree at most M_k,
and M_k/k -> d:=4+16lambda. Taylor's remainder on0<=J<=k, together
with |P_k|<=exp(lambda k), proves

    ||Q_k-G^tensor k||_infinity
       <= exp(lambda k) (e lambda k/L_k)^(L_k)
       <= exp(-lambda k/3)=:rho_k.                   (21)

Indeed L_k>=4lambda k, the exponent is decreasing for L_k>lambda k,
and 5-4log4<=-1/3 follows from log2>=2/3. In particular
||Q_k||_infinity<=1+rho_k. The exponential decay can be extremely
slow for the explicit tiny epsilon, but this has no bearing on the
asymptotic theorem.

Crucially, the uniform approximation also controls the FINAL coefficient
norm, although G^tensor k itself may have degree8k. For every function
E on10k signs and every integer m,

    ||Ehat||_(q_m) <= 2^(5k/m)||E||_2
                  <= 2^(5k/m)||E||_infinity.        (22)

The factor in(22) stays bounded for m>=M_k. No dimension-dependent
error has been silently discarded.

### 8.3 Final changing exponent and every-degree conclusion

For each sufficiently large m choose the largest k with M_k<=m.
The bounded gaps of M_k imply m/k->d. Since |G|=1, the exact Renyi
entropy identity from Section2 gives

    ||(G^tensor k)hat||_(q_m)
       = exp[k H_(m/(m+1))(pi)/(2m)]
       -> exp[H(G)/(2d)].

Equations(21)--(22), divided by1+rho_k, show

    liminf_m B_m >= exp[H(G)/(2d)].                  (23)

Finally lambda<=24epsilon^2, so d<=4+384epsilon^2<5. Equation(19),
epsilon=2^(-104), and log2>=2/3 give

    H(G)-2d log2
       >= [(1170-768)log2-256]epsilon^2
       >=12epsilon^2.

Consequently H(G)/(2d)>log2+epsilon^2, and(23) is strictly greater
than2+2epsilon^2=2+2^(-207), proving(18).

This construction does not contradict tensor non-amplification: Q_k
is NOT the plain tensor P^tensor k. It includes a joint polynomial
in the block modulus indicators, with additional degree16lambda k
asymptotically. Nor does it contradict the torsion theorem: Q_k is
not root-of-unity-valued. The fixed ambient ratio10/d still imposes
a constant ceiling, so this mechanism alone does not prove divergence.

## 9. Arbitrary Boolean block substitution of a fixed seed cannot diverge

The following includes strongly biased inner functions and all Fourier
collisions that those biases create. It is not just the balanced-address
case. Let f be any fixed nonzero complex polynomial on r Boolean inputs,
and let g_i be arbitrary nonconstant Boolean-valued functions on disjoint
blocks, with exact degrees d_i>=1. Write F=f(g_1,...,g_r). Then its actual
degree is

    D=max_(fhat(S)!=0) sum_(i in S) d_i,

and, whenever f is nonconstant,

    ||Fhat||_(q_D)/||F||_infinity <= 2*2^(r/(2D)).   (24)

In particular, for a fixed outer seed, every such family with D tending
to infinity has limsup ratio at most2. Increasing the block dimensions,
using different Boolean constructions in different blocks, and allowing
their biases to approach1 do not evade this statement.

Proof. Put mu_i=E g_i and sigma_i=sqrt(1-mu_i^2)>0. Expand the outer
function in the orthonormal basis for the product measure with these
biases:

    f(y)=sum_S b_S product_(i in S)(y_i-mu_i)/sigma_i.

Thus sum|b_S|^2=E_mu |f|^2<=||f||_infinity^2. Centering each inner function
gives phi_i=(g_i-mu_i)/sigma_i, which has l2 norm1, has no constant Fourier
coefficient, and has at most4^(d_i-1) nonzero coefficients by the Boolean
granularity theorem. Therefore

    ||phihat_i||_(q_D)<=2^((d_i-1)/D).

The nonconstant-block sets of the products product_(i in S)phi_i are
different for different S, so their Fourier supports are DISJOINT.
This is the step that resolves, rather than ignores, all collisions in
the original uncentered substitution. The exact coefficient norm is

    ||Fhat||_(q_D)^(q_D)
       =sum_S |b_S|^(q_D) product_(i in S)
                               ||phihat_i||_(q_D)^(q_D).

Every b_S!=0 has sum_(i in S)d_i<=D, because the biased basis change is
triangular by inclusion. Hence each product factor is at most2, and

    ||Fhat||_(q_D)<=2||b||_(q_D)
       <=2*(2^r)^(1/(2D))||b||_2
       <=2*2^(r/(2D))||f||_infinity.

Each nonconstant Boolean g_i attains both signs, and their input blocks
are independent, so ||F||_infinity=||f||_infinity. Finally, a maximal
positive-weight outer support S has no strict supported superset;
its coefficient in the biased basis is fhat(S)product sigma_i!=0.
The corresponding disjoint product has actual degree sum d_i=D and
cannot cancel against another nonconstant-block set. This proves the
claimed ACTUAL degree and(24). For constant f the ratio is1.

This does not rule out outer seeds whose number of inputs grows fast
relative to their weighted degree; the fixed-seed qualification is
substantive. It also does not cover non-Boolean inner maps, whose ranges
need not preserve the cube supremum norm of f.

## 10. Addressed independent tangent directions: exact variance, l1 cap

This directly tests the proposed infinite-small-perturbation route:
summable squared amplitudes but divergent squared-amplitude entropy.
There really are many independent low-degree tangent directions in an
exponentially large addressed geometry. However, the natural unselected-
branch family pays the l1 norm of its amplitudes in the supremum norm,
not their l2 norm.

### 10.1 Correct addressed seed and its spectrum

The unmodified seed f of Section5 is not balanced: E f=1/4. Addressing
N=2^k independent copies therefore creates a constant-coefficient collision.
Its support is1+15N^2, not16N^2, and its critical ratio at degree k+2
tends sqrt(61)/4, not2.

Use instead the balanced six-variable seed g=s f(x). For each branch j,
let g_j=s_j f_j and A_j=s_j a_j, with a from(12). On k fresh address
bits t, define the exact address indicators I_j(t) and

    F(t,blocks)=sum_(j=1)^N I_j(t)g_j.

Then |F|=1, deg(F)=d=k+3, and its Fourier support consists of16N^2
coefficients, all of modulus1/(4N). Thus R_d(F)=2^((d-1)/d), tending
to2. The extra leaf sign is needed for this exact flat spectrum.

### 10.2 Many exact orthogonal tangent directions

Define, on the SAME addressed cube,

    H_j=A_j-F^2 bar(A_j),
    H=sum_j epsilon_j H_j,             epsilon_j real. (25)

Since |F|=1, J_F H_j=-H_j. Hence Re(bar(F)H)=0 and

    ||F+H||_infinity^2=1+||H||_infinity^2.            (26)

For N>=2, each nonzero H_j has actual degree k+6=d+3. Its top term
comes from an unselected branch l!=j: the full k-address monomial,
a degree4 term of f_l^2, and a degree2 term of bar(A_j). No other
branch produces the same physical monomial. Distinct H_j supports
are disjoint: every term of H_j contains the distinguished sign s_j
and no other distinguished leaf sign, whereas F^2 contains none.

The local exact identities are

    E f^2=-1/8,   E a^2=0,   <a,b>=0,
    ||a||_2^2=||b||_2^2=3,
    ||projection_(|S|<=2) b||_2^2=9/8.

They follow directly from the Section7 coefficient table. Consequently

    ||H||_2^2=6 sum_j epsilon_j^2.                  (27)

Even the genuinely NEW coefficient mass is exactly orthogonal. If
S_F is the Fourier support of F, then

    sum_(s outside S_F)|Hhat(s)|^2
       = [(189-69/N)/64] sum_j epsilon_j^2.          (28)

For completeness, conditional on the address selecting j, the old-
support projection of H_j is A_j minus the old part of s_jb_j, with
squared norm3+9/8=33/8. Conditional on any other selected branch, it
is A_j+(1/8)bar(A_j), with squared norm3(1+1/64)=195/64. Thus its old
mass is(195+69/N)/64; subtract from6 to obtain(28).

### 10.3 The exact alignment obstruction

Each A_j takes values in {0} union {2zeta^l:0<=l<6}. Once the selected
branch l and its physical inputs are fixed, F is a sixth root. The
unselected A_j remain independently selectable, and for every such F,

    max Im(A_j/F)=sqrt(3),   min Im(A_j/F)=-sqrt(3).

Choose each unselected branch to align with the sign of epsilon_j.
The selected branch contributes at worst the opposite amount. Taking
l with minimum |epsilon_l| proves

    2sqrt(3)(1-2/N)||epsilon||_1
       <=||H||_infinity<=2sqrt(3)||epsilon||_1.      (29)

The lower bound is useful for N>=4. Thus(27)--(28) do NOT justify a
quadratic supremum cost: the different branch variables align at a
vertex. In particular epsilon_j^2 proportional to1/[j(log j)^2] has
bounded variance but an unbounded supremum cost.

If ||epsilon||_1<=L, its direction-label entropy is bounded independently
of N:

    sum_j epsilon_j^2 log(L^2/epsilon_j^2)<=2L^2/e.

This follows by writing u_j=|epsilon_j|/L and using
u log(1/u)<=1/e. So the proposed divergent label-entropy gain cannot
occur in this natural independent addressed tangent family at bounded
cap cost. More correlated tangent families are not excluded.

There is also a separate coarse support ceiling. F^2 has at most16N^2
coefficients, A=sum epsilon_jA_j has at most3N, and H=A-F^2bar(A).
Thus F+H has at most67N^3 coefficients. Its actual degree for any
nonzero epsilon and N>=2 is k+6, so Parseval gives

    R_(k+6)(F+H)<2^(3/2).

This support estimate alone rules out divergence within(25), while
the stronger alignment calculation identifies the missing resource in
the suggested summable-variance/divergent-entropy mechanism.

Even a hypothetical improvement from l1 to l2 cap control would not by
itself give divergence through the N direction labels. Their normalized
Shannon entropy is at most log N=k log2, already proportional to the
baseline address degree. In particular probabilities proportional to
1/[j(log j)^2] have truncated entropy only of order log log N, which is
o(k). A successful version would need additional within-direction Fourier
entropy or a genuinely cheaper way to access the labels.

## 11. General entropy / effective-degree transfer

This section isolates the mechanism behind Section8. The typical-set
dual argument was supplied by the director and independently reconstructed
here. It removes any dependence of the approximation error on ambient
dimension. No literature-novelty claim is made.

Primary comparison: the complete proof of Theorem8.6 in
[Pellegrino--Teixeira, arXiv:2608.16584v1, Sections8.1--8.2](https://arxiv.org/html/2608.16584v1#S8)
was read directly. It establishes an entropy/radial-index lower principle
for analytic-polydisc inner functions by tensoring, controlling a Taylor
tail through analytic continuation, retaining coefficient entropy, and
homogenizing. Thus the general entropy-plus-normalization strategy has
clear prior precedent. The mechanism proved here is specifically Boolean:
it corrects the physical log-modulus using a finite polynomial, accounts
for Walsh degree after reduction, and needs no analytic radial Cauchy
estimate or homogenization. Equation(30) is a robust l2 transfer formulation,
not a claim that the underlying entropy strategy is new.

### 11.1 A finite dimension-free coefficient witness

Let G be any unit-modulus function on a finite Boolean cube, a=Ghat,
and let T be a set of Fourier indices satisfying

    M=sum_(s in T)|a_s|^2>0,
    |a_s|^2<=exp(-L) for every s in T.

Let Q be a nonzero polynomial of degree at most D>=1, with
||Q-G||_2<=delta. Then

    R_D(Q) >= M^(1/(2D)) (sqrt(M)-delta)_+
                 exp[L/(2D)] / ||Q||_infinity.     (30)

Proof. Truncate a to a_T. Its pairing with Qhat has real part at least
M-delta sqrt(M), by Parseval and Cauchy--Schwarz. For the conjugate
exponent q_D'=2D/(D-1), with infinity understood when D=1,

    ||a_T||_(q_D')
       <=||a_T||_infinity^(1/D)||a_T||_2^(1-1/D)
       <=exp[-L/(2D)] M^((D-1)/(2D)).

Holder gives(30). Importantly, a_T can include frequencies ABOVE D;
the pairing automatically ignores them on the Q side. No assertion
that the approximated target has small degree is being used.

For a fixed unitary seed G with probability law pi=|Ghat|^2 and entropy
H(pi), take k independent tensor copies and the set where
-log(pi_s1 ... pi_sk)>=k(H(pi)-eta). Its mass tends to1 by the weak law
of large numbers. Therefore any Q_k satisfying

    deg Q_k<=D_k,   D_k/k->d>0,
    ||Q_k-G^tensor k||_2->0,   ||Q_k||_infinity<=1+o(1)

obeys

    liminf R_(D_k)(Q_k)>=exp[H(pi)/(2d)].           (31)

Only l2 approximation is needed in(31), together with the separately
stated cap control. For moving seeds, equation(30) remains exact, but
one must verify the stated information-tail mass rather than silently
use a uniform law of large numbers. This permits, but does not supply,
a divergence route with entropy/effective-degree ratio tending to infinity.

### 11.2 A general radial-normalization criterion

Let P be a fixed nonvanishing complex polynomial of actual degree d>=1
on n signs. Put

    a=min_x |P(x)|>0,
    lambda=log(max_x|P(x)|/a),
    G=P/|P|.

For lambda>0 define R=log(|P|/a)/lambda, a real function with0<=R<=1.
Assume an integer c>=0 satisfies

    deg(P R^j)<=d+c for every integer j>=1.         (32)

There is always the elementary choice c=n-d, because these are functions
on n signs; a special radial algebra can give a much smaller c.
Condition(32) is finitely checkable. If R takes K distinct values,
including0, its value-annihilating polynomial has degree K and zero
constant term. Therefore every positive power R^j is a linear combination
of R,...,R^(K-1), so it suffices to check1<=j<=K-1 in(32). Then

    liminf_m B_m >= exp[H(|Ghat|^2)/(2(d+4c lambda))]. (33)

If lambda=0, the same statement has denominator2d and follows from the
unitary tensor itself. For lambda>0, let J=sum_(i=1)^k R_i and
L_k=ceil(4lambda k). The explicit polynomial is

    Q_k=product_(i=1)^k (P_i/a)
                    sum_(j=0)^(L_k-1) (-lambda J)^j/j!.

Every monomial of J^j involves at most j distinct blocks. In a block
with a positive exponent, (32) limits its degree to d+c, independent
of that exponent. Hence deg Q_k<=D_k:=dk+c(L_k-1), and D_k/k tends
to d+4c lambda. Taylor's remainder gives exactly the estimate(21),
since0<=J<=k and |product(P_i/a)|<=exp(lambda k). Thus Q_k approximates
G^tensor k uniformly with error at most exp(-lambda k/3), so(31) applies.
The bounded gaps of D_k yield the every-degree liminf in(33).

For Section8, R is an indicator and P R=C/48+epsilon h has degree8;
therefore c=4, and(33) has denominator2(4+16lambda), agreeing with the
fully explicit certificate. Ordinary fixed-seed tensoring omits this
radial correction, and therefore does not contradict(33).

The criterion is quantitative but conditional on a favorable entropy
versus radial-degree budget. For every fixed seed its conclusion is
only a constant. Neither this criterion nor the explicit seed establishes
an unbounded sequence of Boolean BH constants.

## 12. Skew-orthogonal coupling does not supply a cheap addressed tangent

The director proposed correlating unselected branches through a skew-
Hadamard or orthogonal matrix instead of adding independent directions.
The following explicit calculation tests that mechanism without assuming
anything about a dense signing optimization problem.

Let u_i=a_i/f_i on independent five-bit seed blocks. Their exact moments
and available values are

    u_i in {0} union {2zeta^j:0<=j<6},
    E u_i=1/4,   E|u_i|^2=3,   E u_i^2=0.           (34)

Every one of the six nonzero values occurs. For a real skew-symmetric
N-by-N matrix K, put

    Q_K=sum_(i<j) K_ij(u_i bar(u_j)-bar(u_i)u_j).

This is purely imaginary. If F is the balanced addressed unitary of
Section10, then H=F Q_K is a tangent to F regardless of their dependence.
The local Walsh degree of u is3, so deg H<=k+9 when N=2^k.

### 12.1 Exact variance and a self-contained supremum lower bound

Write u_i=R_i+iI_i. Then Q_K=2i I^T K R. Equation(34) gives
E R_i=1/4, E I_i=0, E R_i^2=E I_i^2=3/2, and E R_iI_i=0.
Expanding the four indices, using K_ii=0 and independence, yields

    ||Q_K||_2^2=(69/8)||K||_F^2+(3/8)||K 1||_2^2. (35)

There is also the dimension-uniform deterministic inequality

    ||Q_K||_infinity >= sum_i ||K_(i,*)||_2.        (36)

Proof. For any cut I,J of the index set, select u_i=2x_i for i in I
and u_j=2zeta y_j for j in J. Both signs x_i,y_j can be chosen
independently because all six phases occur. Then

    Q_K=-4i sqrt(3) sum_(i in I,j in J)K_ij x_i y_j.

For any real Rademacher sum Z with variance v, E Z^4<=3v^2 and
Holder imply E|Z|>=sqrt(v/3). Optimizing the x_i after randomizing
the y_j thus gives cap at least4 sum_(i in I)(sum_(j in J)K_ij^2)^(1/2).
For a random independent fair cut, conditional on i in I, the expectation
of the row square root is at least half its full row norm: use
sqrt(S)>=S/sqrt(v) for0<=S<=v. Paying the probability i in I gives
one quarter of each full row norm. Some cut realizes this average,
proving(36). The sharp L1 Khintchine constant would improve the factor,
but is not needed or imported here.

If K is also orthogonal, ||K||_F^2=||K1||_2^2=N and each row norm is1.
Therefore

    ||Q_K||_2=3sqrt(N),   ||Q_K||_infinity>=N,
    ||Q_K||_infinity/||Q_K||_2>=sqrt(N)/3.           (37)

This includes normalized skew-Hadamard matrices whenever available.
The orthogonal coupling does not create a bounded cap-to-variance ratio.

### 12.2 Bounded-cap normalization cannot hide unbounded entropy here

Fix c<infinity and choose a real epsilon with
||epsilon Q_K||_infinity<=c. For orthogonal K, (37) forces
|epsilon|<=c/N. The exact unitary phase normalization

    G_N=F(1+epsilon Q_K)/sqrt(1+epsilon^2|Q_K|^2)

satisfies

    ||G_N-F||_2<=|epsilon|||Q_K||_2<=3c/sqrt(N).    (38)

Moreover its Fourier entropy is O_c(log N). Here is a complete useful
support argument that does not rely on an entropy continuity estimate
in exponentially large ambient dimension. The seed u has24 Fourier
coefficients, so Q_K has at most576N^2, and F has16N^2. Put
t=-i epsilon Q_K, so |t|<=c. The function

    phi(t)=(1+it)/sqrt(1+t^2)

has polynomial approximants of degree2j+1 with uniform error at most
C_c rho_c^j, where rho_c=c^2/(1+c^2)<1. Explicitly expand

    (1+t^2)^(-1/2)
      =(1+c^2)^(-1/2) sum_(l>=0)
        [binom(2l,l)/4^l] [(c^2-t^2)/(1+c^2)]^l.

The coefficients in brackets are at most1, so the geometric remainder
bound is elementary. Let S_j be the union of the Fourier supports of
F Q_K^l for0<=l<=2j+1. Then these sets are nested, their cardinalities
are at most exp[C(j+1)log(eN)], and the Fourier probability of G_N
outside S_j is at most C_c^2 rho_c^(2j). Assign each Fourier index
the first level j where it appears. This integer label has bounded
mean and entropy depending only on c. Conditional on level j, its
Fourier-index entropy is at most log|S_j|. The entropy chain rule proves

    H(|Ghat_N|^2)<=C_c log(eN).                     (39)

The case c=0 is simply G_N=F. Uniform convergence of the series shows
there is no missing Fourier support outside the union of the S_j.

Finally, F's address-frequency degree is Bin(k,1/2): its coefficients
are flat and all address masks occur. Its Fourier mass above degree k/3
tends to1. Equation(38) implies the same for G_N. Hence any l2-o(1)
polynomial approximation to r_N tensor copies of G_N, for ANY sequence
r_N>=1, must have degree at least(1/3-o(1))k r_N. To see the uniform
quantifier, the expected fraction of copies with degree below k/3 tends
to0; Markov bounds this fraction even when r_N varies. The Fourier
mass below the stated total degree therefore tends to0.

Thus the entropy/effective-degree ratio of these targets stays bounded
by a constant depending only on c, even after subsequent joint radial
normalization or another l2-accurate preparation. This rules out the
proposed bounded-cap skew-orthogonal coupling as a divergence mechanism
through the Section11 entropy transfer. It does not rule out arbitrary
nonlocal tangent constructions or phase amplitudes growing without bound.

### 12.3 Exact parity restriction prevents compressing the address bits

There is also an exact, stronger degree obstruction for these targets.
The five-bit truth word11 in the binary indexing of the certificate
has a=0 and f=1. Fix every leaf's five-bit input to that word. The
balanced leaf g_j=s_j f_j becomes s_j, and every tangent generator
A_j=s_j a_j vanishes. Thus all the Section10 tangents, the skew
coupling Q_K, and their normalized phase feedbacks vanish or become
their base multiplier1 on this subcube.

The surviving target is the ordinary address function
sum_j I_j(z)s_j. Set the signs s_j to the parity of the binary address
label j. The target then becomes the parity of all k address bits.
For r disjoint copies it becomes a parity character of degree kr.
Coordinate fixing never increases polynomial degree or uniform error,
and a polynomial of degree below kr is orthogonal to that character.
Therefore any uniform approximation with error strictly below1 has
degree at least kr, for EVERY k and r, even if both vary.

This applies more generally to any addressed feedback construction
which is a fixed unit scalar times F whenever all leaves have a=0,
f=1. It does not require small feedback amplitude or proximity to F.
It improves the preceding k/3 spectral argument for the actual tested
targets and rules out compression of their address-bit cost. It does
not bound the entropy of an arbitrary feedback multiplier or exclude
a family with genuinely superlinear entropy in k.

## 13. Real quadratic Hamiltonian phases have linear entropy cost

Final cross-track upgrade: the independently audited flat-diagonal
state argument in the barrier artifact, Sections31--32, proves
`H(exp(itP))<=64 v_flat(exp(itP))`. Every fixed-seed tensor compiler
with cap tending to1 and L2 error tending to0 has degree rate at least
`v_flat`. Thus even arbitrary-time quadratic phases cannot give a
divergent entropy/effective-degree criterion through such L2 compilers.
The proofs below remain independent, more elementary restricted
estimates. The new fixed-seed conclusion must not be promoted to a
moving-seed theorem without paying its finite `n*L2_error` term.

This is a new proved failure class for the phase-Hamiltonian attack, not
a bound for all complex Boolean BH constants. Let

    P(x)=sum_(i<j) a_ij x_i x_j,  a_ij real,
    Q=||P||_infinity,  G_t(x)=exp(i t P(x)).

The hollow symmetric matrix A has entries a_ij. For every real t,

    H(|Ghat_t|^2)<=8 sqrt(3) |t| Q.                 (40)

There is no restriction on the number of variables, rank, coefficient
magnitudes, or spectral norm. Entropies use natural logarithms.

### 13.1 Complete proof

Let p_i be the probability that Fourier index i belongs to a random
index distributed as |Ghat_t|^2. Parseval and the coordinate difference
identity give the exact formula

    p_i=E sin^2(t sum_(j!=i) a_ij x_j)
       <=t^2 sum_(j!=i) a_ij^2.

Shannon entropy is at most the sum of its coordinate marginal entropies.
For 0<=p<=1, the Shannon/Renyi comparison and log(1+u)<=u give

    h(p)<=2 log(sqrt(p)+sqrt(1-p))
         =log(1+2sqrt(p(1-p)))<=2sqrt(p).

Consequently H(|Ghat_t|^2)<=2|t| sum_i ||A_i||_2.
For a real Rademacher sum Z, E Z^4<=3(E Z^2)^2, and Holder yields
E|Z|>=sqrt(E Z^2/3). Hence

    sum_i ||A_i||_2<=sqrt(3) E_y ||Ay||_1
                   <=sqrt(3) max_(x,y Boolean) |x^T A y|.

For any Boolean x,y put u=(x+y)/2 and v=(x-y)/2. Then
x^T A y=2P(u)-2P(v). The multiaffine polynomial P has the same
supremum on [-1,1]^n as on its vertices, so the last maximum is at
most 4Q. This proves (40), including Q=0 or t=0.

### 13.2 What the direct Taylor compiler does and does not rule out

Fix P and t with tau=|t|Q>0. On k disjoint copies the target phase is
exp(i t sum_j P(x^(j))). Truncate its exponential series before
L_k=ceil(4 tau k). The real-axis integral remainder is bounded by
(tau k)^L_k/L_k! <=(e tau k/L_k)^L_k, which tends exponentially to zero.
The approximating polynomial has degree at most 2(L_k-1), so its
asymptotic degree rate is at most 8tau. Its cap tends to one.
Equation (40) shows that inserting THIS declared rate into the Section11
entropy transfer gives a lower-bound output at most exp(sqrt(3)).

This is deliberately a statement about the direct Taylor compiler.
It does not lower-bound the optimal uniform approximation degree, nor
exclude a sharper compiler from producing a better entropy/rate ratio.
The entropy bound alone does not prove bounded BH constants. A separate
degree lower bound would be needed to close that stronger claim.

The elementary proof and finite Fourier identities are replayed by
`computations/bh_boundedness_counterexamples_2026_09_18_hamiltonian.py`.

The barrier agent independently audited the entire proof and332 finite
cases. Its homogenization observation extends this to every real
degree-at-most-two polynomial P, including linear and constant terms:
replace it by P_tilde(x_0,x)=P(x_0 x)-P(0). This is a hollow homogeneous
quadratic, with cap at most2||P||. The phase Fourier coefficients are
unchanged in modulus, with index S relabeled by adjoining x_0 exactly
when |S| is odd. Thus the general quadratic bound is
H(exp(itP))<=16sqrt(3)|t|||P||.

### 13.3 Small quadratic phases cannot escape through an optimal compiler

There is a stronger conclusion when P is homogeneous quadratic and
0<a=|t|Q<=1/2. Every uniform-o(1) tensor approximation degree rate
d_eff for G=exp(itP) satisfies

    d_eff>=a(1-a)>=a/2,
    H(|Ghat|^2)/(2d_eff)<=8sqrt(3).                (40a)

Thus in this small-amplitude class even an optimal compiler cannot
make the Section11 entropy/degree output unbounded. This is stronger
than the Taylor-only statement above and uses no converse theorem for
rotation velocity.

For completeness, let X_i,Y_i,Z_i be the Pauli operators, with Z_i
diagonal in physical cube coordinates, and J=(1/2)sum_i X_i. Put

    C=i[J,M_P]=sum_(i<j) a_ij(Y_i Z_j+Z_i Y_j).

Given a Boolean maximizer x of |P|, the product pure qubit state with
Bloch expectations <X_i>=0 and <Y_i>=<Z_i>=x_i/sqrt(2) has
expectation <C>=P(x). Hence ||C||>=Q. The rotation-velocity operator is

    A(t)=M_G^*J M_G-J
        =integral_0^t exp(-isM_P) C exp(isM_P) ds.

Since ||[M_P,C]||<=2Q||C||, the norm difference between A(t) and
tC is at most t^2 Q||C||. Therefore

    ||A(t)||>=|t|(1-|t|Q)||C||>=a(1-a).

The barrier artifact Section21.3 proves the valid lower bound that
every uniform tensor degree rate is at least this rotation velocity.
Combining it with (40) proves (40a). The barrier researcher independently
reconstructed this entire argument: PASS. Its separate Section23 shows
why one must NOT assert that the velocity is itself an attainable rate.
The amplitude restriction in (40a) is essential to this proof; no
optimal-compiler conclusion for unrestricted |t|Q is being claimed.
As in the velocity theorem, the target P is fixed while the number of
tensor copies tends to infinity. The resulting numerical bound is
uniform over such fixed targets. A simultaneous moving-seed use needs
the explicit finite velocity error n*delta to vanish; mere uniform
errors delta->0 with uncontrolled growing n are not silently covered.

### 13.4 A dimension-free finite estimate covers moving quadratic seeds

The ambient error caveat just stated can be removed for this particular
analytic target class by using finite differences instead of directly
comparing commutators. For P homogeneous quadratic, a=|t|Q>0, and
any polynomial Q_k of degree D uniformly within delta>0 of
G^(tensor k), one has the exact finite lower estimate

    D/k >= [v(G)-4sqrt(delta(a^2+a/k))]/(1+delta). (40b)

There is no ambient-dimension term. In particular, for arbitrary
moving targets, if 0<a<=1/2, D>=1, and delta<=1/2048, then

    k H(|Ghat|^2)/(2D)<=16sqrt(3).                 (40c)

Thus small-amplitude quadratic phases cannot produce divergence via
the entropy/effective-degree criterion even by a moving-seed optimal
uniform compiler. This is not an upper bound on the BH ratio of an
arbitrary nearby polynomial, nor does it cover merely L2 approximation.

Proof of (40b). Write P(s)=exp(isJ)M_P exp(-isJ). It is an
operator-valued trigonometric polynomial of degree2, with constant
operator norm Q. Bernstein gives ||P'||<=2Q and ||P''||<=4Q.
Differentiating exp(itP(s)) by the unitary Duhamel formula gives

    ||F'||<=2a,  ||F''||<=4a+4a^2,
    F(s)=exp(isJ)M_G exp(-isJ).

For k tensor copies, the second derivative is bounded by
4ak+4a^2 k^2. The first derivative at0 has norm k v(G). Taylor's
integral remainder therefore gives, for every real phi>0,

    ||F(phi)^(tensor k)-F(0)^(tensor k)||
       >=phi k v(G)-2phi^2(ak+a^2k^2).

The rotated approximating polynomial has derivative norm at most
D(1+delta), again by Bernstein. Its two endpoint errors cost at most
2delta. Combining these estimates and dividing by k phi gives

    v(G)<=D(1+delta)/k+2delta/(k phi)
                           +2phi(a+a^2k).

Choose phi=sqrt(delta)/sqrt(ak+a^2k^2) to prove (40b).
For (40c), if ak<=1 then (40) already gives kH<=8sqrt(3), and D>=1
suffices. If ak>=1, (40b) and v(G)>=a/2 give

    D/k >= a[1/2-4sqrt(2delta)]/(1+delta)>=a/4

at the stated delta threshold. Equation (40) now proves (40c).
The barrier researcher independently reconstructed all derivative
constants, the optimized step, and the two-regime argument: PASS.

## 14. The normalized seed does not inherit nonconstant offset pairs

This is a focused exact failure of the proposed recursive finite-feature
mechanism, not a classification of all feature algebras. Let G be the
exact normalized ten-bit unitary from Section8, at either
epsilon=2^(-10) or epsilon=2^(-104), and put P=F+epsilon h. Define
the anti-linear involution J_G A=G^2 conjugate(A). Then:

* There is no nonzero complex polynomial A with either
  (deg A,deg J_G A)<=(1,7), (2,6), or (3,5).
* The space satisfying (deg A,deg J_G A)<=(4,4) is exactly the
  one-complex-dimensional space spanned by P.
* J_G P=P, so this last pair has no antisymmetric transfer:
  P(x)J_G P(y)-J_G P(x)P(y)=0.

Here a complex scalar multiple of P is sent to its conjugate scalar
multiple; the dimension statement is about the underlying complex
linear solution space after conjugating the input variable. Thus the
nonconstant degree-saving/spending pairs of total degree at most8
have disappeared. By involutivity the same assertions cover the
reversed degree pairs.

### 14.1 Exact rank certificate

Write q=h/F, R=|q|^2/48 in {0,1}, epsilon=1/E. Pointwise,

    G^2=F^2 [E^2+48-96R+2E q]/[E^2+48].          (41)

All entries lie in Q(zeta), zeta^2-zeta+1=0. Compute its exact
Walsh transform using integer Eisenstein pairs. To require
deg(G^2 B)<=8-s for an unknown B of degree at most s, form the matrix

    M_s(T,S)=widehat(G^2)(T symmetric_difference S),
    |S|<=s, |T|>8-s.

Conjugating A to B imposes no restriction on its possible complex
coefficients. Reduce the common-denominator matrix modulo1009 and
substitute zeta=375; indeed375^2-375+1=0 modulo1009. The factors
1024(E^2+48) are nonzero modulo this prime for both chosen E.
Exact Gaussian elimination gives:

    s             0     1     2      3      4
    columns       1    11    56    176    386
    rows         11    56   176    386    638
    rank mod p    0    11    56    176    385.

A nonzero minor modulo p remains a nonzero algebraic minor over
Q(zeta), and hence over C. This proves the asserted upper bounds on
kernel dimensions. The exact identity J_G P=P supplies the one
dimension for s=4; P has degree4. The replay also checks over the
integers, before reduction, that G^2 has degree exactly8.

The self-contained certificate is
`computations/bh_boundedness_counterexamples_2026_09_18_inheritance.py`.
It uses no numerical rank threshold or approximate equality.

The barrier researcher independently reconstructed the formula and
finite-field implication and reran the exact certificate: PASS. In
particular the claimed complex rank is not inferred from a floating
near-nullspace calculation.

### 14.2 Remaining trivial pair and scope of the recursion barrier

The constant pair (1,G^2), of degrees(0,8), DOES survive. Therefore
the rank certificate must not be reported as excluding every possible
next tangent. Nor does it exclude deliberately increasing the total
degree budget or changing the target.

There is a separate elementary ceiling on an uncompressed binary
recursion: if each generation replaces two disjoint copies by a new
unitary on their union and charges at least twice the previous
effective degree, then ambient dimension/effective degree never
increases. Starting with this ten-bit target and its degree budget
at least4, one always has H/(2d)<=5log(2)/4. Such a recursion cannot
make the Section11 entropy/effective-degree output diverge, even if
a feature pair could be inherited perfectly. A successful recursive
attack therefore needs genuine degree compression or a seed family
whose ambient/effective-degree ratio already grows.

## 15. Growing-degree disjoint Boolean Hamiltonians also obey linear cost

Let b_j be arbitrary nonconstant real Boolean-valued polynomials on
disjoint variable blocks, with degrees d_j>=1. Put

    P=c+sum_j a_j b_j,  a_j real, d=max_(a_j!=0) d_j.

Then, uniformly in the number and sizes of the blocks, their phase
alphabets, the degrees, and the real weights,

    H(|widehat(exp(itP))|^2)<=2|t|d||P||_infinity. (42)

This includes growing addressed Boolean blocks, not just fixed finite
seeds. Ignore terms with a_j=0 and the constant phase exp(itc). For one
Boolean b and u real,

    exp(iu b)=cos(u)+i sin(u)b,
    |widehat(exp(iu b))|^2
       =cos^2(u) delta_empty+sin^2(u)|bhat|^2.

The probability identity at the empty index uses that b has real
coefficients. Entropy of a mixture is at most the entropy of its label
plus the conditional entropy. Boolean granularity from Section1 gives
H(|bhat|^2)<=2(d_b-1)log2. Since
h(sin^2u)<=2|sin u|<=2|u| and sin^2u<=|u|,

    H(|widehat(exp(iu b))|^2)
      <=2[1+(d_b-1)log2]|u|<=2d_b|u|.

The disjoint phase factors have product Fourier probabilities, so
their entropies add. Finally each nonconstant b_j takes both signs
independently on its own block, so ||P||=|c|+sum_j|a_j|. This proves
(42). No claimed bound for overlapping Boolean Hamiltonians follows.

### 15.1 Why the unrestricted phase conjecture is a substantive target

The parent asked whether every real degree-d P might satisfy
H(exp(itP))<=C|t|d||P|| with an absolute C. A counterexample with an
unbounded ratio would already prove unbounded Boolean BH constants:
the k-copy Taylor compiler has degree rate4|t|d||P|| and uniform error
tending to zero, so Section11 gives

    liminf_m B_m >= exp[H(exp(itP))/(8|t|d||P||)]  (43)

for each fixed nonconstant P and nonzero t. Conversely, a hypothetical
uniform BH bound B would force the displayed phase entropy to be at
most8|t|d||P||log B. These are rigorous implications, not a solution of
the unrestricted conjecture. The quadratic proof and (42) are genuine
restricted classes where the required dimension-free linear estimate
can be proved directly.

## 16. Independent cross-audit notes for the final synthesis

* Director Section6, affine data with a small control block: PASS.
  The exact partial-Fourier conditional product distribution, complex
  logarithmic Sobolev step, tensor amplitude inner product, and final
  control-frequency entropy charge were independently reconstructed.
  The estimate is klog2+2|t|Q+2t^2(d-1)Q^2 with d>=1; a constant
  phase should be handled separately if using degree budget0.
* Barrier Section23, failure of the rotation-velocity converse: PASS.
  For f=(1-x1x3+i x2(x1+x3))/2 the exact velocity is sqrt3, certified
  by A^3=3A and tr(A^2)=12. Identifying x3=x1 gives i x1x2, so k
  copies require degree at least2k for uniform error below1. The
  exact symbolic replay passed. This prevents promoting the velocity
  lower bound to an unproved optimal compiler theorem.
* Barrier Section26, correlated-state quantum scalarization: PASS.
  The fixed-suffix anticommuting uncertainty bound, input-dependent
  Parseval budget, stabilizer one-core restriction, and exact real
  norm-saturation argument were independently reconstructed. The
  actual source construction in [Slote, Section4](https://arxiv.org/html/2608.01424#S4)
  was read directly: its terms are an anticommuting core tensored with
  commuting diagonal suffixes, exactly as needed by the audit. Coherent
  Fourier collisions and renormalization by a smaller actual scalar
  cap remain explicit loopholes, not silently excluded.
* Barrier Section28, controlled disjoint Boolean blocks and a shared
  control partition/tree: PASS. Conditional block spectra are exact
  mixtures, including biased Boolean blocks. The nonempty-data event
  bounds control/data information; the rare final control-frequency
  event bounds its separate entropy. Shared-partition support counts
  are paid explicitly.
* Barrier Section29, the rare-AND16 L2/compiler counterexample: PASS.
  The actual primary de Wolf0802.1816v2 AppendixA exact-search loop
  and query-to-degree statement were read directly. Verification of
  candidates maintains the decreasing-count loop invariant even when
  a search has the wrong guess. The query bound, composition degree,
  and binomial L2 tail all check. This defeats the unrestricted
  rotation-velocity lower bound for cap-bounded L2 compilers; it does
  not give an unbounded BH ratio.
* Barrier Sections31--32, quadratic entropy and flat-diagonal velocity:
  PASS after full reconstruction. For a flat-diagonal state, two
  separate state-Cauchy--Schwarz inequalities bound the commutator
  approximation error by `nk*delta`; the main term is at most
  `degree(Q)*cap(Q)`. Hence the exact finite bound is
  `degree(Q)/k >= [v_flat(G)-n*delta]/cap(Q)`. Averaging the quadratic
  product-state witnesses over their auxiliary signs makes their
  physical diagonal exactly uniform; adaptive Y signs do not affect
  that diagonal. This gives `H<=64v_flat`, and closes the fixed-seed
  cap-1+o(1), L2-convergent quadratic compiler class. The Boolean
  identity `v_flat=I` and independent-phase identity
  `v_flat=sum |sin(theta_j)|` check separately. No moving-seed L2
  conclusion is inferred without control of `n*delta`.
  The director's further Section9.4 actual-ratio corollary also passes:
  if the nonconstant seed is fixed and Q_k lives on EXACTLY its nk
  tensor coordinates, the positive degree rate makes
  `2^(nk/(2D_k))*delta_k=o(1)`. Thus its actual critical coefficient
  ratios, not merely the entropy-transfer value, have limsup at most
  `exp(32)` for quadratic phases. Extra auxiliary coordinates are not
  covered by this last comparison.

Field guard: the source paper's statement about no known scalar lower
bound asymptotically above2 is made for real-valued functions. The
strict lower bounds in this artifact are complex. No realification
preserving those strict constants is proved or claimed here.

These audit results do not add an unrestricted upper or divergent
lower bound for the Boolean BH constants.

## 17. Independent phase products cannot exploit an L2-only compiler

This closes a genuine gap left by the uniform-approximation velocity
argument. Fix any finite independent-phase target

    G(x)=exp(i sum_(j=1)^n theta_j x_j),  theta_j real.

Suppose Q_k has degree at most D_k, uniformly bounded supremum norms,
and limsup_k ||Q_k-G^(tensor k)||_2<1. Set
d_eff=liminf_k D_k/k. Then

    H(|Ghat|^2)<=64 d_eff.                         (44)

In particular the Section11 entropy/effective-degree output is at most
exp(32) in this class, even with L2 rather than uniform approximation.
The theorem is uniform over all fixed finite seeds and angles. A
simultaneous moving-seed argument requires the finite estimates below;
no unstated uniformity of the limiting restriction argument is claimed.
This does not bound the BH ratio of every polynomial near the target.

### 17.1 Exact complex evaluation lower bound

Changing angles by integer multiples of pi changes only a global phase,
and input sign flips replace theta by -theta. Thus assume
0<=theta_j<=pi/2. Consider a subset B of coordinates whose angles lie
in (u/2,u], with 0<u<=pi/4. Write n_B=|B|. In k tensor copies, fix
every input outside B in every copy to a value whose conditional squared
L2 error is at most the original squared error delta_k^2. Such a value
exists by averaging. Remove the target's resulting global phase from
the restricted polynomial q_k. It still has degree at most D_k and
cap at most C_k=||Q_k||_infinity, and it approximates the independent
product phase on B^k in L2 to error at most delta_k.

The target is unitary, so its inner product with q_k has real part at
least1-delta_k. Expanding its exact product Fourier coefficients gives

    1-delta_k <= product_(j in B)(cos theta_j)^k
                   |q_k((-i tan theta_j)_(j in B^k))|.

Define the single-variable polynomial

    p(z)=q_k(((tan theta_j)/(tan u)) z)_(j in B^k).

For real z in [-1,1], every argument belongs to [-1,1]. Multiaffinity
and the triangle inequality extend the cube cap to this full box, so
|p(z)|<=C_k there. The elementary Bernstein--Walsh estimate yields

    |p(-i tan u)|<=C_k exp[D_k asinh(tan u)].

For a self-contained proof, put z=(w+w^(-1))/2. The polynomial
w^D p((w+w^(-1))/2) has degree at most2D and boundary modulus at most
C_k on |w|=1. Apply the maximum principle to its reciprocal polynomial.
At z=-iy, choose the root with |w|=y+sqrt(1+y^2), proving the estimate.

Consequently the following finite inequality is exact whenever
delta_k<1:

    (D_k/k) asinh(tan u)
       >= sum_(j in B) log(sec theta_j)
            +(1/k) log[(1-delta_k)/C_k].           (45)

For the stated fixed-seed assumptions the last term tends to zero.
Since log(sec theta)>=theta^2/2 and
asinh(tan u)=integral_0^u sec(s) ds<=sqrt(2)u, we obtain

    d_eff>= n_B u/(8sqrt(2)).                      (46)

This proof needs neither a uniform approximation error nor a velocity
converse. In fact the cap assumption can be weakened to
log C_k=o(k), while keeping limsup delta_k<1.

### 17.2 Dyadic angle summation controls the whole entropy

Put u_l=(pi/4)2^(-l), l>=0, and let n_l count angles in
(u_l/2,u_l]. The independent spectral probabilities are Bernoulli
with success probabilities sin^2 theta_j. For theta<=u<=pi/4,

    h(sin^2 theta)<=u^2 log(e/u^2).

Indeed h(p)<=p log(e/p), sin^2 theta<=u^2, and the last function
increases for0<=p<=1. Equation (46) bounds every n_l u_l by
8sqrt(2)d_eff. Therefore the entropy from angles at most pi/4 is
at most8sqrt(2)d_eff times

    sum_(l>=0) u_l log(e/u_l^2)
      =(pi/2)[1+2log(8/pi)]<5.

The last elementary bound follows from pi>3, pi<16/5, and e>8/3.
For the n_large remaining angles above pi/4, their entropy is at most
n_large log2. The ordinary Fourier-degree weak law gives
d_eff>=I(G)=sum_j sin^2 theta_j>=n_large/2: below this mean, the
k-copy Fourier mass tends to zero, and the best possible L2 error
tends to1. Hence their contribution is at most2log2 d_eff.
Adding the two parts gives

    H(G)<[40sqrt(2)+2log2]d_eff<64d_eff,

which proves (44). Zero angles contribute nothing. If d_eff=0 the
same bounds force zero entropy, so there is no divergent criterion.

The barrier researcher independently audited the conditional fixing,
complex evaluation identity, Bernstein--Walsh proof, dyadic constants,
and the large-angle degree argument: PASS.
The boundedness researcher independently reconstructed the entire
proof as well: PASS. The deterministic replay is
`computations/bh_boundedness_counterexamples_2026_09_18_phase_l2_compiler.py`;
400 complex evaluation/Bernstein cases and the dyadic sum pass.

The later flat-diagonal state theorem in barrier Section32 gives the
sharper `H<=2d_eff` under cap tending to1 and L2 error tending to0.
The complex-analytic theorem here remains independent and allows a
bounded nonvanishing error strictly below1 and any subexponential cap.

## 18. Rare-marker insertion does not create a free entropy gain

The rare-event quantum compiler in the barrier artifact motivates a
different scalable attempt: mark only a small fraction of blocks, put
an arbitrary complicated phase V inside them, and search for the marked
blocks before paying for V. The following exact identity records the
entropy cost; no query-complexity converse is assumed.

Let E be a nonzero {0,1}-valued polynomial of degree r>=1 on one cube,
with p=E E<=1/4. Let V be any unitary-valued function on a disjoint
cube. Put s=E|V-1|^2<=4 and

    G=1+E(V-1).

Then G is unitary, and

    H(|Ghat|^2)<=p H(|Vhat|^2)+16pr.              (47)

If s=0 this is immediate. Otherwise define the probability vectors
mu=|Ehat|^2/p and nu=|widehat(V-1)|^2/s. Since ps<=1, the ENTIRE
Fourier probability law, including its empty coefficient, is

    |Ghat|^2=(1-ps)delta_empty+ps(mu tensor nu).    (48)

Every nonempty coefficient has this form by the disjoint product
expansion. At the empty index write c=E V. The right side is
1-ps+p^2|c-1|^2, while the left side is |1+p(c-1)|^2; they agree
because s=2-2Re(c). This proves the identity, rather than just a bound
obtained by discarding interference.

Boolean granularity for1-2E implies that all nonzero coefficients of
E have magnitude at least2^(-r), and that p>=2^(-r). Parseval gives
|supp Ehat|<=p4^r, hence H(mu)<=2r log2. Mixture entropy yields

    H(G)<=h(ps)+ps H(mu)+p s H(nu).

For0<s<=4,

    h(ps)<=ps log(e/(ps))<=4p log(1/p)+p.

The second inequality uses s(1-log s)<=1. If c=E V, comparison of
the sole changed coefficient at the empty index gives

    s H(nu)-H(V)
      =|c|^2 log|c|^2-|c-1|^2 log|c-1|^2+s log s
      <=1/e+4log4<6.

Combining these estimates with log(1/p)<=r log2 gives
H(G)<=pH(V)+12pr log2+7p<=pH(V)+16pr, proving(47).

Here is the precise compiler limitation, separated from the entropy
theorem. Suppose a particular rare-marker construction charges a
declared effective degree budget

    B>=c_0 r sqrt(p)+p d_V,   c_0>0, d_V>0,

where d_V is its charged cost for the inner phase. Then

    H(G)/B<=max{H(V)/d_V,16sqrt(p)/c_0}.           (49)

Both numerator and denominator are sums with the displayed weights.
Thus a search-and-evaluate compiler paying its marker-search and
active-inner costs cannot amplify an existing entropy/cost ratio by
rare insertion alone. In the literal quantum implementation the inner
cost must be an actual phase-evaluation query cost; a small abstract
Walsh degree must NOT be substituted for query cost without a compiler.
Equation(49) is not an optimal-degree lower bound and does not exclude
an additional degree collapse or a different implementation. The
proved finite contribution is the exact identity(48) and bound(47).

The barrier researcher independently read and reconstructed the full
finite argument, including the empty-coefficient interference and
declared-budget limitation: PASS.

## 19. Bounded primary-literature check of the strict complex lower bound

Checked on 2026-09-18, during the final audit. This is a short targeted
search, not a systematic review or a novelty certificate. I did not
identify, among the sources checked, a previously stated strict
`B_m>2` or `liminf B_m>2` for COMPLEX-VALUED Walsh polynomials with the
supremum taken on the Boolean cube. That limited negative search result
does not establish that no such result exists.

The important source distinctions are as follows.

* [Defant--Mastylo--Perez, arXiv:1706.03670, Section1.1](https://arxiv.org/html/1706.03670#S1.SS1)
  defines its cube constants for real-valued functions. Its comparison
  with complex polynomials concerns the polytorus, a different domain.
* [Volberg--Zhang, Theorem1.1](https://research-explorer.ista.ac.at/download/13318/17299/2024_MathAnnalen_Volberg.pdf)
  explicitly extends the cube inequality to complex-valued functions,
  noting that the earlier real-valued proof also works there. This is
  an upper-bound statement, not a sharp complex lower-bound statement.
* [Arunachalam--Dutt--Escudero Gutierrez--Palazuelos, Proposition1.5](https://link.springer.com/article/10.1007/s00208-025-03142-5)
  proves the exact bound `2^((d-1)/d)` for Boolean-VALUED functions,
  with equality for an address function. Its adjacent discussion of
  the best known multilinear lower bound is explicitly real-scalar.
  It is not a classification of all complex-valued cube polynomials.
* [Slote, arXiv:2608.01424, Section2](https://arxiv.org/html/2608.01424#S2)
  states the hypercube theorem for `f:{+-1}^n -> R` immediately before
  its report that no asymptotic lower bound above2 is known to the
  author. Therefore our COMPLEX result must not be announced as
  resolving or improving that real-scalar claim.
* [Pellegrino--Teixeira, arXiv:2608.16584v2, definitions(1.1) and Theorem9.4](https://arxiv.org/html/2608.16584v2)
  supplies a directly relevant earlier mechanism: tensorized inner
  functions, coefficient entropy, and radial approximation yield a
  persistent lower bound `liminf D_m>1.27`. But `D_m` is defined for
  homogeneous holomorphic polynomials with their cap on the complex
  polydisc. This is not a strict-above2 Boolean-Walsh result. Its
  entropy--radial idea is prior literature, so generic use of that
  idea is not claimed new by this campaign.

The introductions, relevant theorem/definition passages, and lower-bound
discussion were inspected directly. The newer cyclic-group upper papers
arXiv:2609.07758v2 and arXiv:2609.12427 and the support-sensitive paper
arXiv:2607.05594 were also checked for the relevant distinction; no
strict complex Boolean lower bound was located in those checked passages.
Targeted searches included complex/Boolean/hypercube/Walsh lower bounds,
small-degree constants, and circle/root-of-unity-valued extremizers.
Search results for complex BH frequently refer instead to polydisc caps;
results using the word unimodular frequently mean unit-modulus
COEFFICIENTS, not unit-modulus VALUES.

Conclusion: retain the exact proofs and independent replays of Sections7--8,
label them complex scalar, and make no external priority claim. Any
publication-level novelty claim would require a fuller bibliography and
expert checking of this precise field/domain convention.

For clarity, boundedness itself is nevertheless equivalent over the two
scalar fields. If `c_q=(E_theta |cos(theta)|^q)^(1/q)`, rotation averaging
of the real inequality gives

    B_m(real)<=B_m(complex)<=B_m(real)/c_(q_m).

Indeed the averaged q-th power of the coefficient norm of
`Re(exp(i theta)f)` equals `c_q^q ||fhat||_q^q`, whereas every rotated
real part has cap at most `||f||_infinity`. Since `c_q>=2/pi` for
`1<=q<=2`, this compares boundedness by an absolute factor. But
`c_(q_m)->1/sqrt2`, so a complex lower bound just above2 does not yield
a real lower bound above2. This elementary comparison is not asserted
as a new result.
