# Support-free critical-scale realization and anchored extension

2026-09-17. Joint localization/Bernoulli-track deduction. The nuclear
normalization and the rank-one matrix-MGF proof were independently
reconstructed by both tracks. External novelty is not asserted. This
artifact preserves exact physical signs and the exact old quadratic block.

The applicable positive conclusion is a **center-cover theorem**, not a
finite-entropy-slope theorem for the whole near-level set. The latter has
been proved inapplicable; see Section 5.

## 1. A single empirical list controlling every response matrix

Let nu be any probability law on {+-1}^k with E hh^T=I_k. It need not
be centered. Define phi(T)=E_nu ||Th||_2 for any real matrix with k columns,
and write ||T||_* for its nuclear norm. Suppose

    q >= 3k log(4k),
    delta = sqrt(3k log(4k)/q) <= 1,
    C0 = 16 sqrt(e*pi/2).

There exists one deterministic list h_1,...,h_q from the support of nu
such that

    sum_j h_j h_j^T <= q(1+delta) I_k,                 (1)

and, simultaneously for EVERY real matrix T with k columns,

    (1/q) sum_j ||T h_j||
       <= phi(T) + C0 sqrt(k/q) ||T||_*
       <= (1+C0 k/sqrt(q)) phi(T).                    (2)

There is no restriction on the support size, rank of T, number of query
matrices, or hierarchy depth. In particular, (2) is not obtained by
paying separately for different response levels.

### 1.1 Nuclear norm is controlled by the true response

Put A=(T^T T)^(1/2). Isotropy, Cauchy--Schwarz pointwise in h, and
||h||=sqrt(k) give

    ||T||_* = tr A = E h^T A h
             <= sqrt(k) E ||A h|| = sqrt(k) phi(T).   (3)

Every response and nuclear norm depends only on T^T T. Therefore it is
enough to control square k by k matrices; rectangular output dimensions
do not add complexity.

### 1.2 Sharp aggregate covariance from a scalar matrix MGF

Sample H_1,...,H_q independently from nu and write G=sum_j H_j H_j^T.
Since (hh^T)^2=k hh^T,

    E exp(theta hh^T) = [1+(exp(theta k)-1)/k] I_k.

The Golden--Thompson inequality, conditioning on the final summand and
iterating, consequently gives

    E tr exp(theta G)
       <= k [1+(exp(theta k)-1)/k]^q
       <= k exp[(q/k)(exp(theta k)-1)].               (4)

This is a specialized matrix Chernoff proof, with no noncommuting MGF
product asserted. The scalar identity after conditioning is essential.
The source matrix trace inequality and iterative method are recorded in
Tropp, [User-friendly tail bounds for sums of random matrices,
Sections 2.4 and 3.7](https://tropp.caltech.edu/papers/Tro11-User-Friendly-preprint.pdf).
The same paper's Corollary 5.2 gives the general Chernoff statement.

Markov with theta=log(1+d)/k yields, for 0<=d<=1,

    Pr{lambda_max(G) >= (1+d)q}
       <= k exp[-(q/k)((1+d)log(1+d)-d)]
       <= k exp[-q d^2/(3k)].                        (5)

Thus (1) fails with probability at most 1/4. Jensen applied to
exp(theta lambda_max(G))<=tr exp(theta G), with theta=1/k, also gives

    E lambda_max(G) <= k log k+(e-1)q <= e q.         (6)

The last inequality only requires q>=k log k.

### 1.3 One nuclear-unit-ball empirical event

Let

    D = sup_(||T||_*<=1) [(1/q)sum_j ||T H_j||-phi(T)].

The set is compact and contains zero, so D is measurable and nonnegative.
Symmetrization, and Jensen for the independent absolute values of standard
Gaussians, give

    E D <= (2/q) E sup_T sum_j epsilon_j ||T H_j||
         <= (2/q)sqrt(pi/2) E sup_T sum_j gamma_j ||T H_j||.

Condition on the H_j. The Gaussian process on the final line has squared
increment at most sum_j ||(T-U)H_j||^2. Gaussian comparison therefore
bounds its expected supremum by

    E sup_(||T||_*<=1) sum_(j,a) gamma_(j,a)(T H_j)_a
       = E ||sum_j gamma_j^(vector) H_j^T||_op.

Conditional on the sample, the matrix in this operator norm has the law
Z G^(1/2), where Z is a standard k by k Gaussian matrix. Its expected
operator norm is at most 2sqrt(k) sqrt(lambda_max(G)). To see the first
factor directly, compare the process u^T Z v on two unit spheres with
u^T g+v^T g': the increment inequality follows from
(1-u^T u')(1-v^T v')>=0. The comparison supremum has expectation at most
2sqrt(k).

Combining this with (6) and Jensen proves

    E D <= 4 sqrt(e*pi/2) sqrt(k/q).                  (7)

Markov shows D<=C0 sqrt(k/q) with probability at least 3/4. Its
intersection with (1) has probability at least 1/2. Select one list in
that intersection, use homogeneity, and then (3). This proves (2).

## 2. Exact signs and sharp concentration at the critical aspect ratio

Write n=kp+ell with 0<=ell<k. Fix labels from Section 1, and independently
sample fair scalar signs g_(j,t), 1<=t<=p, and fair signs on the ell
leftover old coordinates. Define the jth bridge column

    C_j=(h_j tensor g_j, xi_j).                       (8)

Every entry is a sign. For w in R^n, identify its core with a k by p
matrix W. All physical signs within one k-coordinate bundle share the
same scalar driver; these dependencies are retained. For every new word
y, the linear response w^T C y has subGaussian proxy

    V(w)=sum_j ||W^T h_j||^2+q||w_left||^2
          <=q(1+delta)||w||^2.                       (9)

For the exact new-spin maximum B_w=sum_j |w^T C_j|, scalar-sign bounded
differences give the same centered MGF proxy V(w), not twice that proxy.
Indeed changing a scalar driver changes B_w by at most twice the
absolute value of its linear coefficient. Thus

    E exp[t(B_w-E B_w)] <= exp[t^2 V(w)/2],
    E B_w <= sqrt(q V(w)).                            (10)

The nuclear response estimate additionally gives

    E B_w <= q(1+C0 k/sqrt(q)) E_nu||W^T h||
                +q||w_left||.                       (11)

Suppose now k=k_n=O(sqrt(n)) and q=floor(epsilon n) for any fixed
epsilon>0. Then q/(k log(4k)) tends to infinity, delta tends to zero,
and 1+C0 k/sqrt(q) is bounded (with a constant depending on epsilon).
Consequently every center family F_n with

    log|F_n|=o(n),
    sup_(f in F_n) E_nu||F_core^T h||/sqrt(kp) -> 0   (12)

satisfies, with probability tending to one,

    max_(f in F_n,y) |f^T C y| = o(n^(3/2)).         (13)

Leftovers cost at most q sqrt(ell)=o(n^(3/2)). The multiplicative
empirical-response error is paid at the center scale, before taking any
near-level hierarchy, and tends to zero at each fixed epsilon.

## 3. Entire low-response center codes are permitted

For any mu_n tending to zero, take ALL core Boolean matrices X satisfying

    Phi_nu(X)=E_nu||X^T h||/sqrt(kp) <= mu_n,

and all leftover words. The scalar information theorem in
[the discrepancy track](paper_discrepancy_scalar_response_entropy_2026_09_17.md)
gives log|F_n|=o(n), with no aspect-ratio or support restriction.
Specifically, the physical isotropic sign law h tensor g satisfies
E|<h tensor g,X>|<=sqrt(kp) Phi_nu(X); the theorem bounds the normalized
logarithmic cardinality by O(mu_n^2 log^2(1/mu_n)). The ell leftovers
cost ell log 2=o(n). Hence this ENTIRE low-response code satisfies (12).

The code can be empty for an unsuitable law/cutoff; a cover hypothesis
below necessarily rules that out. No existence of favorable nu_n for
an actual optimizer is presumed.

## 4. Anchored full-sign extension, including k of order sqrt(n)

For actual full-sign A_n, write Q_n=max_x |H_(A_n)(x)| and let E_n(eta)
be its full absolute eta n^(3/2) near-level set. Choose arbitrary
isotropic laws nu_n in dimension k_n=O(sqrt(n)) and center families as
in (12), including the whole low-response code of Section 3 if desired.
Assume that for small eta,

    rho(eta)=limsup_n max_(x in E_n(eta))
                    min_(f in F_n) d_H(x,f)/n < 1/2,
    K_H=limsup_(eta down to0) rho(eta)h(rho(eta))/eta < infinity.

Then for every tau>2K_H and all sufficiently small fixed epsilon>0,
there exist actual full-sign extensions W retaining A_n exactly, of
order n+floor(epsilon n), such that

    limsup_n [Q(W)-Q(A_n)]/n^(3/2)
                         <=tau epsilon+epsilon^(3/2). (14)

Here is the residual payment, to make the change from the impossible
whole-nearcode entropy criterion explicit. For x within r n flips of
a center f, put w=x-f. Its squared norm is at most 4rn. The number
of center/deviation pairs is at most exp[o(n)+n h(r)+o(n)]. Equations
(9)--(10), followed by a union bound only over these deviations, give

    B_(x-f)/n^(3/2)
       <=2epsilon sqrt(r)+sqrt(8epsilon r h(r))+o_n(1). (15)

The exact new-spin maximization is already inside B_w; there is no
extra 2^q entropy. Equation (13) is paid once for the entire center code.

For completeness choose K_H<kappa and a finite geometric hierarchy
eta_j=eta_0 R^(-j), R>1, ending at eta_J of order epsilon^2. Strict
radius majorants obey r_j h(r_j)<=kappa eta_j. For a nonterminal shell,
its old deficit is at least eta_j/R, and

    sup_(eta>=0) [sqrt(8epsilon kappa eta)-eta/R]
                                            =2kappa R epsilon.

The first term of (15) has coefficient at most 2sqrt(r_0), which tends
to zero with eta_0. The terminal shell costs O(epsilon^(3/2)), absorbed
by the strict tau margin. The exterior old deficit eta_0 pays the
global O(sqrt(epsilon)) bridge bound for sufficiently small epsilon,
using (9) and a union over both spin words. A new principal signing
with cap at most q^(3/2) costs the final epsilon^(3/2). Fix epsilon and
its finite hierarchy before taking n to infinity. This proves (14).

This is a genuine improvement in the deployment range for arbitrary
isotropic nu: previous absolute empirical-response control required
k=o(sqrt(n)); relative nuclear control includes k=O(sqrt(n)). The
orthogonal, deterministically balanced special case already allowed
k=o(n), and is not claimed to be improved here. Nor does (14) assert
that its center-cover hypothesis holds for actual minimizers.

More generally the same proof only needs k log(4k)=o(n), together with
maximal center response mu_n tending to zero and k mu_n=o(sqrt(n)).
These conditions make the empirical correction in (11) vanish at each
fixed epsilon and make delta tend to zero. Thus the theorem can also
reach k larger than sqrt(n) when the protected centers have sufficiently
small response; no such stronger optimizer geometry is asserted.

## 5. Why the direct whole-near-level shortcut was rejected

Let x maximize sigma H_A(x)=Q(A). Flip a uniformly chosen set of r
coordinates to obtain Y. Exactly

    E[Q(A)-sigma H_A(Y)]
            =4r(n-r)Q(A)/[n(n-1)].                  (16)

The deficit is nonnegative. If Q(A_n)/n^(3/2) tends to c>0 and
4c rho<eta, Markov retains a positive fraction of all binom(n,r)
words with r/n tending to rho inside E_n(eta). Consequently

    liminf_n log|E_n(eta)|/n >= h(rho).

Thus the full near-level entropy divided by eta diverges as eta tends
to zero for every bounded positive-cap family. An extension criterion
requiring a finite such slope is valid as an implication but has no
application to the original matrices. It cannot justify (14).

The stronger response-rigidity argument of the director/discrepancy
tracks additionally shows that vanishing-response whole-nearcode
profiles cannot be o(sqrt(eta)). This does not invalidate the anchored
construction: centers are cheap, but the unavoidable neighborhood
deviations are explicitly retained and paid through (15).

## 6. Audit status and reproducibility

The finite inequalities and critical-scale center-cover argument were
self-audited before circulation. The director independently read all of
Sections 1--4 and returned PASS, explicitly checking the scalar matrix-MGF
iteration, all-matrix Gaussian comparison, one-event selection, exact-sign
construction, bounded-difference constant, and order of limits. The replay
`computations/paper_localization_2026_09_17_critical_realization.py`
passed 600 nuclear-response checks, 600 Gaussian increment checks,
128 exact centered-absolute MGF diagnostics, and 45 exact fixed-radius
flip identities, using isotropic laws at k=2,3,4,8. The largest numerical
error in the scalar matrix-MGF identity was 4.45e-16. Output is preserved
at `tmp/paper_portfolio_2026_09_17/localization/critical_realization_audit.json`.
These diagnostics do not substitute for uniform proof. The constants are
intentionally not optimized; the essential feature is one relative-response
certificate for all matrices and all levels, with covariance 1+o(1).
