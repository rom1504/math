# Exact center iteration and the fresh-spin obstruction

Date: 2026-09-06. Status: exact finite-dimensional dynamics and exact
Hadamard examples proved below. The twin-root limit uses the already
independently audited fixed-polynomial old/star CLT and own-root removal
module. This does not assert an iterative Gaussian state evolution.

The one-step center theorem is recorded in
`resumed_bound_audit_restricted_channel_center_update_2026_09_06.md`.
Its restricted tested identities do **not** identify the entire
conditional law of BF. In particular, throughout Section 4 the center
spin is the explicit tested spin sign(S K+s Z), not sign(BF).

## 1. A nonperturbative iteration that is always feasible

Let B be any real symmetric hollow matrix, fix a deterministic or random
mask h in [0,1]^n, and let

    |f_i| <= 1-h_i,                 |c_i| <= h_i.

Then both f+c and -f+c belong to the cube. Put E(f,c)=f^T Bc and
q_B(x)=x^T Bx/2. Exactly,

    [q_B(f+c)-q_B(-f+c)]/2 = E(f,c).

Consequently the absolute same-spin maximum is at least |E(f,c)|. This
uses only symmetry and cube feasibility, not any stochastic independence.

Starting with any feasible f^0,c^0, alternate

    c^{t+1}_i = h_i sign((Bf^t)_i),
    f^{t+1}_i = (1-h_i) sign((Bc^{t+1})_i).

Choose sign(0) by retaining an existing endpoint sign whenever one is
available. Both updates are exact maximizations of the bilinear
objective over their respective boxes. The gains are exactly

    E(f^t,c^{t+1})-E(f^t,c^t)
      = sum_i [h_i |(Bf^t)_i|-c^t_i(Bf^t)_i] >= 0,

    E(f^{t+1},c^{t+1})-E(f^t,c^{t+1})
      = sum_i [(1-h_i)|(Bc^{t+1})_i|-f^t_i(Bc^{t+1})_i] >= 0.

These formulas remain valid at h_i=0 or 1 and at zero fields: there is
no division by a radius. After the first full round, each vector takes
only endpoint values, so there are finitely many possible states. With
the stated tie rule every changed coordinate has strictly positive
gain. Hence the iteration terminates at a coordinatewise best-response
pair in finitely many rounds, although no useful uniform bound on the
number of rounds follows.

This provides a valid iterative feasible algorithm, but no positive
lower bound on a later gain. Such a bound is false in general, including
at bilinear values greater than .43, as Section 3 shows.

## 2. A fixed unbalanced binary mask has a genuine spectral ceiling

Write R_4=J_4-2I_4. It is a symmetric Hadamard matrix, R_4^2=4I_4,
with row sum 2 and diagonal -1. For k>=1 put

    n=4^k,  R=R_4^{tensor k},  d=(-1)^k,  A=R-dI_n.

A is a genuine symmetric hollow signing, and

    ||A||op=sqrt(n)+1,
    max_{x in {+1,-1}^n} |x^T A x|/2 = n(sqrt(n)+1)/2.

For odd k the all-one vector attains the positive spectral endpoint.
For even k take a balanced sign eigenvector of one R_4 factor and the
all-one vectors of the other factors; the resulting sign eigenvector
attains the negative endpoint. This proves the displayed equality,
not just its upper bound.

If h is binary and has pn nonzero coordinates, every pair with
|f|<=1-h and |c|<=h satisfies

    |f^T A c|/n^(3/2)
       <= (1+n^(-1/2)) sqrt[p(1-p)].

Thus no fixed-mask bilinear iteration can approach 1/2 for a fixed
unbalanced limiting density p, even on a signing sequence whose actual
absolute same-spin maximum approaches 1/2. For example the banked
mask density .52970552247711725 gives a ceiling strictly below 1/2.
This is a restriction on that bilinear certificate and fixed mask,
not a limitation on the actual available signings or their energies.

## 3. An exact stationary pair at 7/16

Use the lexicographic Kronecker ordering for R_16=R_4 tensor R_4, and set

    h=c=(1,1,1,0, 0,0,0,0, 1,0,1,1, 0,0,0,0),
    f=  (0,0,0,-1, 1,1,1,1, 0,-1,0,0, 1,1,1,1).

Direct integer multiplication gives

    R_16 f=(4,6,4,2, -2,0,-2,0, 4,2,4,6, -2,0,-2,0),
    R_16 c=(0,2,0,-2, 2,4,2,4, 0,-2,0,2, 2,4,2,4).

Therefore all support-relevant fields have strictly correct signs and

    c=h sign(R_16 f),  f=(1-h)sign(R_16 c),  f^T R_16 c=28.

Lift by tensoring f,c,h with the all-one vector of length m=4^ell,
and use A=R_16 tensor R_4^{tensor ell}-(-1)^ell I. Since the second
factor has row sum sqrt(m), the lifted pair remains strictly
stationary. The diagonal deletion does not affect either support-
relevant field. Its bilinear value, divided by (16m)^(3/2), is

    28/64 = 7/16 = .4375 > .43.

Thus fixed-mask alternating center/response updates cannot have a
universal strictly positive increment above .43. This example is NOT
a counterexample to improvement of the actual best signing energy:
here f^T R_16 f=-12, c^T R_16 c=4, and
(-f+c)^T R_16(-f+c)=-64, so that particular endpoint signing already
attains the asymptotic value 1/2. The distinction is essential.

### 3.1 A stronger stationary example whose endpoints also stay below 1/2

For R_64=R_4^{tensor 3}, use the following lexicographic vectors, with
spaces only separating groups of eight entries:

    f: 0000+000 --+00++- -0--0++0 -+0+-0-+ +--+000+ -00+0+-+ 00+0-0-0 -+0+-0--
    c: +-+-0-+- 000-+000 0-00+00+ 00+00+00 0000+--0 0--0+000 ++0+0-0- 00+00+00

Here + means 1 and - means -1; set h=|c|. The exact standard-library
integer replay `computations/resumed_response_iteration_certificate_2026_09_06.py`
checks R_64^2=64I, its row sums, all entries of both response fields,
and the following identities and inequalities:

    sum h=27,  |f|+|c|=1,
    min_{c_i!=0} c_i(R_64 f)_i=3,
    min_{f_i!=0} f_i(R_64 c)_i=1,
    f^T R_64 c=223,  f^T R_64 f=35,  c^T R_64 c=-5.

The same tensor lift, now to N=64m, is therefore strictly stationary
under both fixed-mask updates. Its normalized bilinear value is

    223/512=.435546875 > .43.

The two endpoint vectors f+c and -f+c have full-Hadamard quadratic
forms 30+446=476 and 30-446=-416 at the base level. Under tensor lift,
diagonal deletion has normalized effect O(N^-1/2). Consequently the
better absolute same-spin endpoint energy converges to

    max(476,416)/(2*512)=119/256=.46484375 < 1/2.

Nevertheless the actual absolute same-spin maximum for that very same
signing sequence converges to 1/2, by Section 2. Thus this elementary
alternation can genuinely stop with both its endpoint signings below
the available optimum. No claim is made that this particular stationary
state is reachable from the specific randomized Gaussian-mask initial
state of the one-step theorem. What it rules out is a universal positive
subsequent gain based only on cube feasibility, high bilinear value,
bounded operator norm, and coordinatewise best-response updates.

### 3.2 Even variable-mask coordinate ascent has strict suboptimal plateaus

Allow h to change, equivalently optimize over the full coordinatewise
diamonds |f_i|+|c_i|<=1. Because B is hollow, holding all other
coordinates fixed makes the objective affine in (f_i,c_i). The exact
best choice maximizes

    f_i(Bc)_i+c_i(Bf)_i,

whose maximum is max(|(Bc)_i|,|(Bf)_i|). In endpoint coordinates
x=f+c and y=-f+c, this is exactly independent coordinate ascent on
q_B(x) and coordinate descent on q_B(y):

    E(f,c)=[q_B(x)-q_B(y)]/2,
    |u+v|/2+|u-v|/2=max(|u|,|v|).

This more permissive algorithm also has strict suboptimal fixed points.
For R_64, take the following two sign vectors:

    x: -+++-+++--+---+-----+++++-+-+-+-++-+----+++++----+-+-+-+----++++
    y: -++-+--++-+--+-+---------+----+-+++++++++++++++++--+-++---++++--

The same exact integer certificate verifies all 64 local inequalities:

    x^T R_64 x=448,       min_i x_i(R_64 x)_i=4,
    y^T R_64 y=-464,      min_i -y_i(R_64 y)_i=4.

Thus x is a strict single-spin local maximum of its quadratic form,
and y is a strict single-spin local minimum. With
f=(x-y)/2, c=(x+y)/2 and h=|c|, the pair is a strict local maximum
against every individual coordinate diamond update, including mask
switches. Tensor lift by an all-one vector of length m=4^ell preserves
strictness even after diagonal deletion: the signed margins are at
least 4sqrt(m)-1>0. The resulting normalized limits are

    bilinear value:       (448+464)/(4*512)=57/128=.4453125,
    best endpoint energy: max(448,464)/(2*512)=29/64=.453125,
    actual optimum:       1/2.

The normalized operator norm tends to 1. This falsifies a universal
strict positive drift based only on high value, cube feasibility,
operator norm near 1, and coordinatewise optimization, even if the
mask can change. It does not rule out nonlocal updates or an escape
theorem for the particular distribution of states reached from the
Gaussian-mask construction.

## 4. Marginally fresh tested centers need not be product spins

### 4.1 The one-root coupling

In the audited one-root Gaussian limit, S is a fair sign independent
of the old Gaussian fields X and a star channel Z=sqrt(v)N. Let K=K(X)
be the even first-chaos adjoint, s>0, and v>0. Define the explicit
tested center spin

    S'=sign(S K+s Z).

Conditionally on X, S K+s Z is symmetric. Hence S' is a fair sign
independent of the entire old one-root field X. It is nevertheless
correlated with both the old spin and the new channel:

    E[S S'|X] = 2 Phi(K/(s sqrt(v)))-1,
    E[Z S'|X] = 2 sqrt(v) phi(K/(s sqrt(v))) > 0.

Equivalently, away from the zero-probability threshold,

    S'=S sign(K)     if |s Z|<|K|,
    S'=sign(Z)       if |s Z|>|K|.

Single-root freshness relative to X is therefore insufficient to
reuse S' as a vector of independent input signs. The next construction
shows actual macroscopic pair dependence in genuine sign matrices.

### 4.2 A vertex-transitive signing sequence with twin rows

Let R=R_4^{tensor k}, n=4^k, N=2n, d=(-1)^k, and

    A=R tensor J_2-d I_N,  B=A/sqrt(N-1),  m=N-1.

This is a symmetric hollow signing with ||B||op tending to sqrt(2).
Each pair (a,1),(a,2) consists of twins: the two rows agree outside
their own two positions. The matrix is invariant under swapping any
twin pair, and the product of the within-factor S_4 permutations acts
transitively on its vertices. Consequently all rootwise distributions
and root-removal errors below are the same at every vertex.

With P=I_n tensor J_2, one has exactly

    Q=B^2=[2n P-2d(R tensor J_2)+I_N]/m.

Its diagonal is 1, its within-twin off-diagonal entry is 1-1/m, and
all entries between distinct twins have absolute value 2/m. For each
fixed odd r>=3,

    ||Q^(circ r)-P||op
      <= r/m+(N-2)(2/m)^r = o(1).

For a fixed normalized finite odd mixture h=sum_{r>=3} b_r h_r,
sum b_r^2=1, let Z=B h(BS) and
R_h=sum b_r^2 Q^(circ r). Then R_h=P+o_op(1), uniformly for this fixed
mixture. Direct multiplication gives

    (B P B)_ii=(4n-3)/(2n-1) -> 2,
    v_i=(B R_h B)_ii -> 2.

### 4.3 Derivation of the two-root limit; no assumed multiroot law

Let i,j be twins and X_i a fixed finite old tree family. Under the
permutation swapping i and j,

    X_j(S)=X_i(S with S_i,S_j interchanged).

The exact marked tree fields at root i exclude S_i. Their fixed-degree
nonroot influences are O(N^-1/2) in every fixed Lp. Thus
X_i-X_j=o_Lp(1). This also follows directly by deleting the terms
containing the other twin label; the remaining injective tree sums
coincide under the twin automorphism.

Put G=BS. Exactly,

    G_i-G_j=B_ij(S_j-S_i),
    Z_i-Z_j=B_ij[h(G_j)-h(G_i)].

For fixed polynomial h, hypercontractive moment bounds and the first
identity imply Z_i-Z_j=O_Lp(N^-1). The audited old/star theorem gives
the one-root limit (X_i,Z_i) -> (X,sqrt(2)N_0), with independence of
X and N_0, together with vanishing own-output-spin removal error for
Z_i. If that removal result is formulated as an average over roots,
vertex transitivity makes it rootwise here. Twin coincidence makes
the S_j-removal error vanish as well. Old X_i is own-i-free and its
S_j-removal error vanishes by the preceding influence bound.

Conditioning both X_i,Z_i on all spins except S_i,S_j therefore changes
them by o_L2(1). The conditioned pair is exactly independent of the
two fair signs. Single-root convergence and twin coincidence now give

    (S_i,S_j,X_i,X_j,Z_i,Z_j)
      -> (S_1,S_2,X,X,sqrt(2)N_0,sqrt(2)N_0),

where S_1,S_2,N_0,X are mutually independent. This is the required
two-root test limit, derived from the proved one-root theorem,
automorphism, and influence/removal estimates; a general multiroot
CLT is neither assumed nor needed.

For a fixed finite continuous K(X), the continuous mapping theorem
applies to the tested new spins: conditional on X, the positive-
variance Gaussian term makes either sign argument zero with probability
zero. Gaussian-a.e.-continuous K follows by the same approximation
argument. Their limiting pair is

    S'_a=sign(S_a K(X)+s sqrt(2)N_0),  a=1,2.

Conditionally on X,N_0, average independently over S_1,S_2 to obtain

    E[S'_1 S'_2|X,N_0]
      = 1_{|s sqrt(2)N_0|>|K(X)|}.

Hence

    lim E[S'_i S'_j]
      = P(|s sqrt(2)N_0|>|K(X)|) > 0.

If ||K||_2<=1, the right side is at least Phi(-1/s): restrict to
|K|<=sqrt(2), an event of probability at least 1/2, and use independence.
Each S'_a is fair and independent of X, yet the pair is not independent.
This is an actual signing-sequence obstruction to blindly treating
the updated center as fresh product input.

### 4.4 What is not inferred about actual best response

For any bounded response vector F, twin rows give the deterministic
estimate

    |(BF)_i-(BF)_j| <= 2/sqrt(N-1).

This alone does not imply sign((BF)_i)=sign((BF)_j) with high
probability: an anti-concentration estimate near zero is additionally
needed. The restricted one-root identities do not identify the full
conditional law of BF because coherent first-chaos Q-spin terms may
remain. Accordingly the pair-correlation conclusion above is asserted
only for the explicit tested center sign(S K+sZ).

## 5. Consequences and remaining freedom

There are two distinct successful statements: exact finite alternation
is feasible and monotone; the Gaussian tested center gives a uniform
one-step improvement over the entire original marked-mask class. They
do not combine automatically into an iterative improvement theorem.

The twin example identifies the missing datum: the joint dependence
of the new center spins across roots, not their one-root fairness or
independence from old one-root Gaussian fields. Enlarging the one-root
Gaussian dictionary without tracking this dependence is insufficient.
Conversely, this counterexample has asymptotic operator norm sqrt(2),
not 1; it does not rule out a more specialized near-optimal regime
with an extra incoherence hypothesis.

Resampling independent signs independent of the entire previous
history cannot preserve a gained center correlation: conditional on
that history, the expected new bilinear center contribution is zero.
A successful iterative scheme would instead require a proved joint
coupling that retains the useful correlation, or a different feasible
comparison whose new independent inputs also modify the odd response.
Neither construction is established here. Variable-mask updates and
schemes tracking coherent Q-spin blocks remain genuinely open as
global or distribution-sensitive mechanisms; purely coordinatewise
variable-mask ascent has the rigorous obstruction in Section 3.2.
