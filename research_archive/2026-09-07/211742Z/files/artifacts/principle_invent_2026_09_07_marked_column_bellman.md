# A marked-column Hadamard recursion for biased Boolean slices

2026-09-07. **Proved one-row exponent theorem; independently reconstructed
PASS by the construction agent. No mixed-profile cap certificate is asserted
here.** This is an actual ensemble operation preserving one DC
column, not an application of the ordinary signed-input law to an asymmetric
word. It removes the precise dephasing dependence left in
`principle_invent_2026_09_07_balanced_transform_compiler.md`.

The original information envelope throughout this note is

    E_t(nu)=sup_L {g_t(E Var(X|L))-I(X;L)}.

The Gaussian reward is applied AFTER averaging the conditional variance.
No alternate envelope with E[g_t(Var(X|L))] is substituted.

## 1. Exact marked ensemble and disappearance of input phases

Use the recursive normalized Hadamard ensemble of
`continued_convergence_recursive_orbit_bound_2026_09_06.md`:

    U_s=diag(U_q^1,U_q^2) H_2 g,  s=2q,
    H_2=2^(-1/2)[I I; I -I],

with a fresh independent uniform signed input permutation g. Mark a fixed
output row, always in the first child, and multiply each physical input
coordinate by that row's sign. Denote the resulting matrix by Ubar_s. Its
marked row is exactly 1/sqrt(s) in every coordinate. Thus sqrt(s) Ubar_s^T
is an actual sign Hadamard with a marked all-ones column.

Write r_1 for the marked row signs of U_q^1. Before g the marked row signs
are (r_1,r_1). If Pi is the underlying UNSIGNED permutation of g, direct
coordinate evaluation gives

    g D_(marked row signs after g)=D_(r_1,r_1) Pi.

Consequently the exact dephased matrix is

    Ubar_s=diag(Ubar_q^1, U_q^2 D_(r_1)) H_2 Pi.       (1)

Conditional on the first child, U_q^2 D_(r_1) has exactly the ordinary
unmarked child law: its own fresh signed input permutation absorbs this
deterministic diagonal sign matrix. This conditional law does not depend
on the first child, so the two displayed child matrices are independent.
Pi is independent and uniform. Hence (1) is a recursion with ONE marked
child, ONE ordinary child, and an UNSIGNED input permutation. At terminal
nodes the marked basis is a fixed dephased Hadamard followed by a uniform
unsigned input permutation. All identities are exact at finite order.

## 2. Mixed orbital estimate at a marked terminal node

Let P_m denote the unsigned permutation group and G_m the signed group.
Write L_G(v)^2=E_(g in G_m) exp(-t||v-gv||^2), and define L_P analogously.
The Gaussian-Fock proof of the ordinary orbital estimate gives, uniformly
over orthogonal V and ||v||^2<=Cm,

    E_(Pi in P_m) L_G(V Pi v)
        <=exp(O_(t,C)(sqrt(m))) L_P(v).             (2)

Indeed the unsigned orbit covariance has operator norm L_P(v)^2: its
finite orbit Gram matrix is positive entrywise with constant row sum.
The OUTPUT projection remains onto G_m-invariant Fock vectors, whose
degree-truncated rank is exp(O(sqrt(m))). The same Poisson truncation and
the lower bound L_P(v)^2>=exp(-4t||v||^2) finish the proof. No unsigned
output invariance or column-phase independence is being presumed.

For a finite raw empirical law nu,

    m^(-1) log L_P(v_m) -> Phi_t(nu)=-F_t(nu)/2,
    F_t(nu)=inf_(coupling nu,nu) [D+t E(X-Y)^2].     (3)

This is the ordinary raw Gaussian-kernel permanent type calculation.

## 3. Exact finite-depth marked type recursion

Let b_r=B^r Phi be the ordinary unmarked exponent on symmetric laws,
with the Bellman operator from the recursive-orbit source. Define marked
exponents a_r on arbitrary finite raw laws by a_0=Phi and

    a_r(nu)=sup_pi { [a_(r-1)(law((A+B)/sqrt(2)))
                   +b_(r-1)(sym law((A-B)/sqrt(2)))]/2
                    -D(pi||nu tensor nu)/2 },      (4)

where the average FULL marginal of pi is nu. Unlike the ordinary
recursion, this is not merely an average absolute-marginal constraint.

For an unsigned input permutation, the number of source words is
s!/prod n_a!, while the number producing an ordered pair table is
(s/2)!/prod n_ab!. Stirling's formula gives exactly the exponent
-D(pi||nu tensor nu)/2 under the average-full-marginal constraint.
The marked and ordinary child matrices in (1) are independent. The
concatenation bound L_G(v_+,v_-)<=L_G(v_+) L_G(v_-) and (2) at marked
terminals therefore prove, at EACH FIXED depth r and finite alphabet,

    limsup_m m^(-1) log E L_G(Ubar_m v_m)<=a_r(nu).  (5)

All type counts have polynomial size at fixed r and terminal Fock losses
sum to O_r(sqrt(m)). This does not assert a growing-depth type estimate.

An important exact simplification is allowed: average pi under input
swap (A,B)->(B,A). This preserves the raw PLUS law and the symmetrized
MINUS law, while decreasing relative entropy. After averaging, BOTH
input marginals equal nu and the minus law is already symmetric. Thus
(4) can be optimized over swap-invariant couplings with identical
marginals, and its entropy charge is I(A;B)/2. If nu has mean zero,
both child laws also have mean zero. Their average second moment equals
that of nu, so every descendant at depth ell has moment at most
2^ell m_2(nu).

## 4. The original envelope is a mixed supersolution

The precision-Schur theorem
`decisive_director_precision_schur_supersolution_2026_09_06.md`, equation
(2), was proved for ARBITRARY finite-second-moment input laws:

    E_t(U)+E_t(V)-I(A;B)<=E_t(A)+E_t(B),
    U=(A+B)/sqrt(2), V=(A-B)/sqrt(2).               (6)

No symmetry is used in its conditional-label/precision factorization.
For the swap-invariant coupling above, both right-side laws are nu and
V is symmetric. Hence E is a supersolution of the mixed recursion (4)
when E is used for BOTH children. This uses the original information
envelope, not a newly enlarged one.

## 5. Extinction of the marked spine: no asymmetric boundary theorem

The already-proved ordinary theorem gives b_d(nu)->E_t(nu), uniformly
over symmetric laws with second moment at most each fixed B. For clarity,
this uniformity follows directly from the stopped-tree estimate in
`continued_convergence_terminal_gap_reduction_2026_09_06.md`: with
K(s)=-g_t(s), for any C>B,

    delta_d(B):=sup_(symmetric m_2<=B) [b_d-E_t]
       <=epsilon+B sup_(s>C) K(s)/s
                    +K(C)K(B)/(d kappa(epsilon,C)).      (7)

First C increases, then epsilon decreases, then d increases. The precise
Schur theorem supplies the temperature-alignment hypothesis of that
stopped-tree estimate. Its compactness/rigidity constants do not depend
on the chosen source inside the moment ball.

For every centered nu with v=m_2(nu), every r>L>=1, expansion of just the
first L marked nodes and the mixed supersolution (6) give

    a_r(nu)<=E_t(nu)
       +sum_(ell=1)^L 2^(-ell) delta_(r-ell)(2^ell v)
       +2^(-L) K(2^L v).                          (8)

To verify the last term, every Bellman value a_j is nonpositive, since
its terminal value and negative entropy charges are nonpositive; while
E_t(mu)>=g_t(Var(mu))>=-K(m_2(mu)). The remaining marked descendant has
moment at most 2^L v and total tree weight 2^(-L). Every sibling at depth
ell is ordinary, so (7) bounds its excess. Prefix envelope terms telescope
by (6). The argument is uniform over all root policies, so it also bounds
the supremum defining a_r.

For fixed L the sum in (8) tends to zero as r tends to infinity. Then
2^(-L) K(2^L v) tends to zero, since K(s)=O_t(log(1+s)). We obtain

    limsup_(r->infinity) a_r(nu)<=E_t(nu),          (9)

uniformly over each fixed centered bounded-second-moment source class.
Only one branch remains marked; there is no need to prove an asymmetric
Gaussian rigidity or a new controlled CLT. Equality in (9) is not needed
and is not asserted.

## 6. Actual biased-slice input and its row exponent

Let q/m->p in (0,1), choose q physical rows T, and take a Boolean word x
of fibre mean a=q^(-1)sum x_i. Its centered transformed coordinates are

    Ubar_m v,  v=(1_T(x-a))/sqrt(p).

The raw source law is therefore EXACTLY

    nu_(p,a)=(1-p) delta_0
       +p(1+a)/2 delta_((1-a)/sqrt(p))
       +p(1-a)/2 delta_((-1-a)/sqrt(p)).            (10)

It has mean zero and variance 1-a^2. The marked spectral coordinate is
identically zero, because sum v=0. The number of words in this physical
slice is binom(q,q(1+a)/2), with exponent p h((1+a)/2). Thus (5), (9),
and fixed-depth-then-order limits give the one-row upper exponent

    p h((1+a)/2)+E_t(nu_(p,a)).                   (11)

More precisely every positive exponent tolerance is obtained at a
sufficiently large FIXED recursion depth, common to all a in [-1,1].
The bounded moment uniformity above supplies the common depth. At that
fixed depth alphabet size is bounded uniformly in a, multinomial/Stirling
errors are O_r(log(m)/m) uniformly even at zero masses, and the integer
table maximum is bounded by the continuous Bellman supremum. The envelope
is uniformly W_2 continuous on the resulting bounded source supports.
These observations justify the stated uniformity for order-varying biases.
There are at most q+1 physical slice types, so summing over those types
has no leading entropy cost.

The DC column may be placed in the fibre's self-port and omitted. The
other output columns retain independent signed permutations. The usual
orbital deletion bound loses at most sqrt(2m), hence no leading row
exponent. Exact balanced-column repair then transfers the centered
transform in operator norm as proved in the balanced-transform artifact.

Equation (11) is a proved and usable biased-row certificate. It is NOT
yet a strict coefficient times (1-a^2): the entropy, scalar envelope,
single-temperature choice, heterogeneous fibre profiles, and full
physical diagonal must still be combined in the actual weave pressure.
In particular an assertion that the old symmetric scalar certificate
automatically scales with the variance would not follow from this note.

## 7. Exact profile normalization and a scalar falsifier

Let v_i=1-a_i^2 and omit each zero DC port. For the unrepaired centered
weave its transformed row has squared norm m v_i. If H_bulk is the
half-Hamiltonian, the squared edge defect satisfies

    D_sigma/(2q)=m sum_i v_i-2 sigma H_bulk/q.      (12)

Therefore a sufficient single-temperature row condition for the uniform
variance-weighted bound

    |H_bulk(x)|<=b sqrt(N) sum_i ||P_i x_i||^2+o(N^(3/2))

is, up to suitable uniform strictness away from zero total variance,

    p h((1+a)/2)+E_t(nu_(p,a))
           +t(1-2b sqrt(p))(1-a^2)<=0.             (13)

The prefactor in (13) includes BOTH the entropy of the physical biased
slice and the full variance-dependent energy budget. It uses a common t
over edges, not independently chosen endpoint temperatures.

At fixed p in (0,1) and t>0, (13) has a rigorous rare-minority obstruction.
Write a=1-2r and reveal only whether X is the minority atom of (10), which
has probability p r. The resulting conditional-variance average is

    V_r=4(1-p) r^2(1-r)/(1-p r).

The label has information h(p r), so the ORIGINAL envelope obeys

    E_t(nu_(p,1-2r))>=g_t(V_r)-h(p r).

Since g_t(V)=-tV+O_t(V^2) at zero,

    p h(r)+E_t(nu_(p,1-2r))
           >=p r log p+O_(p,t)(r^2).              (14)

The variance is 4r(1-r). Hence a necessary condition for the particular
scalar certificate (13), uniformly over positive fixed minority fractions,
is

    b >= [1+p log p/(4t)]/[2 sqrt(p)].             (15)

In particular this floor exceeds 1/2 whenever

    t> -p log p/[4(1-sqrt(p))].                    (16)

The right side approaches 1/2 as p increases to one. The previously
successful p=.96,t=4.85 strict symmetric certificate lies well inside
this obstructed regime. Thus the new marked exponent DOES NOT by itself
certify a strict variance-weighted bulk bound at that temperature. The
obstruction holds at arbitrarily small but FIXED minority fractions, so
it is not removed merely by allowing an additive o(N^(3/2)) remainder.

This is a falsifier of the scalar information/row-relaxation route, not
an actual mixed-word lower witness for the signing ensemble. It leaves
possible stronger geometric row estimates, multitype edge kernels, or a
nonlinear response sufficient for a specified constant-mode seed.

There is a stronger exact statement allowing the temperature to depend
on the bias. Put v=4r(1-r), d=V_r/v=(1-p)r/(1-pr),
H=p h(r), J=h(pr)-H>=0, and z=t v. The constant label and the minority
label, respectively, imply

    A_t:=p h(r)+E_t(nu_(p,1-2r))
       >=max {H-z, -d z-J},

because g_t(s)>=-t s. Minimizing the maximum of the two resulting affine
functions of 1/z gives

    inf_(t>0) [1+A_t/(t v)]/(2 sqrt(p))
       >=(1-d) p h(r)/(2 sqrt(p) h(pr)).            (17)

The two branches meet at z=(H+J)/(1-d); the displayed minimum can also
be checked directly on either side. Its right side tends to
1/(2 sqrt(p)) as r decreases to zero. Conversely E_t<=0, so taking
t to infinity bounds the same infimum above by 1/(2 sqrt(p)). Hence

    lim_(r->0) inf_(t>0) [1+A_t/(t v)]/(2 sqrt(p))
                 =1/(2 sqrt(p))>1/2.             (18)

This is an analytic all-temperature falsifier, independent of floating
optimization. It concerns this scalar certificate, not the actual cap.
