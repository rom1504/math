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
