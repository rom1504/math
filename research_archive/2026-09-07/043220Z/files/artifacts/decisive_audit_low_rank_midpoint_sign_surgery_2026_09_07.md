# Actual sparse-sign midpoint repair from a quadratic ramp certificate

Status: exact finite conditional construction, with a nonvacuous actual
sign example. The existence of its ramp certificate for width minimizers
is OPEN. A stronger all-high-state capture hypothesis is shown below to
be asymptotically circular and is not the main theorem.

Let A be a hollow symmetric sign matrix of order n. Write

    H_A(x)=x^T A x/2,       P=max H_A,       R=max(-H_A),
    w=(P+R)/2,             I=(P-R)/2.

Complement A if needed so I>=0. The case I=0 already has Q(A)=w.
All statements below assume I>0. Let Pi be an orthogonal projector of
rank r, and suppose the following RAMP inequality holds for every
Boolean x:

    H_A(x)<=w+(I/n) x^T Pi x+e,       e>=0.                 (1)

Unlike requiring near-complete projection of every state with H_A>w,
(1) only requires projection in proportion to that state's excess.

## 1. General low-rank theorem, with no incoherence hypothesis

Choose mu>0 such that t*mu<=1, where t=2I/n. Delete from the projector
the coordinates with Pi_ii>mu: let Z be the diagonal indicator of the
remaining coordinates, k=n-Tr Z, and Pi_0=Z Pi Z. Then

    k<=r/mu,     0<=Pi_0<=Id,     Tr Pi_0<=r.

There exists a full hollow sign matrix A' on the SAME n vertices with
at most 4Ir changed undirected edges and

    Q(A') <= w+e+2I sqrt(k/n)+Ir/n
               + (4 K_G I mu w)/n
               + sqrt(8Ir a)+(4/3)a,                      (2)

where a=(n+2)log 2 and K_G is the real bilinear Grothendieck constant.
One may replace k in (2) by min(n,r/mu). This construction changes old
coefficients; it is not an auxiliary-only extension or an unpriced block
completion.

### 1.1 Masking the projector

For every Boolean x,

    ||Pi Zx|| >= ||Pi x||-sqrt(k).

Since both projection norms are at most sqrt(n), this implies

    x^T Pi_0 x >= x^T Pi x-2 sqrt(nk).

Thus (1) gives

    H_A(x)-(I/n)x^T Pi_0 x <= w+e+2I sqrt(k/n).              (3)

Also H_A>=I-w and 0<=x^T Pi_0 x<=n, so the left side of (3) is at
least -w. Removing the forbidden diagonal from Pi_0 adds the harmless
constant I Tr(Pi_0)/n. Therefore

    Q(A-t offdiag(Pi_0))
       <= w+e+2I sqrt(k/n)+Ir/n.                           (4)

Here and below Q applies equally to a hollow real coefficient matrix.

### 1.2 Implementing the penalty by valid sign flips

Factor Pi=sum_(l=1)^r v_l v_l^T with orthonormal v_l. Set z_l=Zv_l and

    K=sum_l |z_l| |z_l|^T.

Then K is PSD, K_ii=(Pi_0)_ii<=mu, and |(Pi_0)_ij|<=K_ij<=mu.
Independently flip the existing edge sign A_ij with probability

    p_ij=(t/2)[K_ij+A_ij(Pi_0)_ij].                         (5)

These probabilities lie in [0,t mu] subset [0,1]. Their exact mean
matrix, on off-diagonal entries, is

    Abar=A-t offdiag(Pi_0)-t(A circ K).                     (6)

The extra Schur term in (6) is real and must be paid for; it is not
silently discarded as unbiased rounding.

Hollow polarization gives beta(A)<=4w. The simultaneous diagonal
majorant D>=A,-A has Tr D<=K_G beta(A). Schur multiplication by K
preserves these PSD inequalities. Hence, for Boolean x,

    |H_(A circ K)(x)|
       <=(1/2)sum_i D_ii K_ii
       <=2 K_G mu w.                                     (7)

Combining (4), (6), and (7) accounts for the first five terms of (2).
No operator norm assumption on A was used.

### 1.3 The full-cube rounding error is subleading when r=o(sqrt(n))

Since sum_l ||z_l||_1^2<=n sum_l ||z_l||_2^2<=nr,

    sum_(i<j) p_ij <= t sum_(i<j)K_ij <= Ir.                 (8)

The centered random edge coefficients are independent, bounded in
absolute value by 2, and have total variance at most 4Ir. For any fixed
spin, Bernstein's inequality bounds its two energy tails by

    2 exp[-u^2/(2(4Ir+2u/3))].

Take u=sqrt(8Ir a)+(4/3)a, a=(n+2)log2. Then
u^2>=2a(4Ir+2u/3), so the union bound over the entire cube gives

    Pr(Q(A'-Abar)>u)<=1/2.

Markov's inequality and (8) give probability at most 1/4 of more than
4Ir edge flips. Both desired events therefore hold together with
probability at least 1/4. This proves (2) and the stated actual edit count.

## 2. A simpler entrywise-incoherent variant with larger allowed rank

Let kappa=max_(i!=j)|Pi_ij|, assume t*kappa<=1, and do not mask Pi.
Independently flip edge ij with probability

    p_ij=(t/2)[kappa+A_ij Pi_ij].                           (9)

Again these are valid probabilities. Now the exact mean is

    Abar=(1-t kappa)A-t offdiag(Pi).

Using (1), the forbidden-diagonal correction, and triangle inequality,

    Q(Abar)<=w+e+Ir/n+t kappa Q(A).

Writing d=n(n-1)/2, the total flip expectation is at most t kappa d.
The same Bernstein argument proves existence with

    Q(A')<=w+e+Ir/n+t kappa Q(A)
                   +sqrt(8t kappa d a)+(4/3)a.             (10)

At most 4t kappa d edges need change. This variant uses no Grothendieck
inequality and permits any rank r=o(n), provided
kappa=o(n^(-1/2)) and the ramp error e=o(n^(3/2)).

### 2.1 Sign-compatible ramps have no mean surcharge

There is a different useful sufficient hypothesis which does not demand
entrywise incoherence. Suppose

    A_ij Pi_ij>=0 for every i<j,       t max_(i<j)|Pi_ij|<=2.

Flip edge ij independently with probability p_ij=t A_ij Pi_ij/2.
The mean is EXACTLY A-t offdiag(Pi), with no Schur or scalar bias.
Let

    S=sum_(i<j)p_ij=t Tr(A Pi)/4>=0.

The same proof provides a signing with at most 4S changed edges and

    Q(A')<=w+e+Ir/n+sqrt(8Sa)+(4/3)a.                     (10a)

Also S<=t ||A||op r/4. Therefore for ||A||op=O(sqrt(n)), bounded
normalized w,I, e=o(n^(3/2)), and r=o(n), the rounding cost in (10a) is
O(n sqrt(r))+O(n)=o(n^(3/2)). A projector supported on positive clique
blocks can satisfy this condition even when its off-diagonal coherence
is not o(n^(-1/2)). The sign-compatibility requirement is essential to
the exact mean formula; it is not automatic for a general projector.

In fact bounded operator norm is not necessary for this asymptotic
conclusion. Sign compatibility and Cauchy--Schwarz give

    Tr(A Pi)=sum_(i!=j)|Pi_ij|
       <=sqrt(n(n-1) Tr(Pi^2))<=n sqrt(r),
    S<=I sqrt(r)/2.                                       (10b)

If normalized I is bounded and r=o(n), (10b) makes S=o(n^2), so the
normalized rounding error in (10a) is O((r/n)^(1/4))+O(n^(-1/2)).
Together with e=o(n^(3/2)) and t*kappa<=2, this proves midpoint repair
for arbitrary bounded-width parents satisfying the sign-compatible ramp,
without a spectral hypothesis. The operator estimate above is a sharper
rate when a bounded operator norm is additionally available.

Explicitly, writing i=I/n^(3/2), the normalized bound is

    Q(A')/n^(3/2) <= w/n^(3/2)+e/n^(3/2)+i*r/n
       +2sqrt(i(1+2/n)log2)*(r/n)^(1/4)
       +(4/3)(n+2)log2/n^(3/2).                           (10c)

The edit count is at most 2I sqrt(r)=2i sqrt(r/n)*n^2. The entrywise
feasibility condition is still explicitly I*kappa<=n; the trace bound
does not imply it and does not remove it.

### 2.2 Positive contractions and a trace budget also suffice

Sections 1--2.1 remain valid with 0<=Pi<=Id in place of an orthogonal
projector, and r=Tr Pi in place of its rank. In the factorization use
Pi=sum_l z_l z_l^T from its spectral square root, whose total squared
column norm is Tr Pi. In the masking argument use Pi^(1/2) rather than
Pi in the two norm expressions. Every probability, diagonal, Schur,
and total-variance estimate is otherwise identical. This makes the
ramp a finite semidefinite certificate with a small trace budget; no
claim that such a certificate exists for width minimizers is added.

## 3. Quantitative asymptotic implications and the unproved premise

Suppose w<=C n^(3/2), e=o(n^(3/2)), and r=o(sqrt(n)). Put

    alpha=r/sqrt(n),       mu=alpha^(1/3)/sqrt(n).

Then k/n<=alpha^(2/3), t mu<=2C alpha^(1/3)<=1 eventually, and (2)
gives

    Q(A')/n^(3/2) <= w/n^(3/2)+e/n^(3/2)
      +(2C+4K_G C^2)alpha^(1/3)+C r/n
      +sqrt(8C alpha(1+2/n)log2)+O(n^(-1/2)).               (11)

Thus Q(A')<=w+o(n^(3/2)), by an actual same-order sign construction.
For a sequence of width minimizers, existence of these ramp projectors
would imply M_n-W_n=o(n^(3/2)). The same conclusion follows from (10)
with r=o(n) and kappa=o(n^(-1/2)).

Neither global width minimality nor the exact extremal-cut identity has
been shown to provide such projectors. Equation (1) is a concrete
low-complexity quadratic certificate that suffices for recovery; it is
not asserted for arbitrary width minimizers. In particular this theorem
does not establish convergence of the original optimum.

## 4. Why the earlier stronger high-set capture hypothesis is circular

Suppose T=w+h<P, and every spin with H_A(x)>T satisfies

    ||Pi x||^2 >=(1-epsilon)n.                              (12)

Let x* attain P and put v=(P-T)/(2P), b=sqrt(1-v). Choose independent
spins y_i of means b x*_i. Then

    E H_A(y)=(1-v)P=(P+T)/2.

Since H_A<=P, this implies Pr(H_A(y)>T)>=1/2. On the other hand,

    E[y^T Pi y]=b^2 x*^T Pi x*+v r<=n-v(n-r),
    Var(y^T Pi y)<=2n.                                    (13)

For completeness, write y=b x*+z with independent centered coordinates.
Since y_i^2=1, the centered expression has exactly two orthogonal chaos
levels: the linear part is 2b[(Pi-diag Pi)x*] dot z, and the quadratic
part is 2 sum_(i<j) Pi_ij z_i z_j. Since -Id<=Pi-diag Pi<=Id and
Var(z_i)=v, their total variance is at most
4b^2 v n+2v^2 r<=(4v-2v^2)n<=2n. This proves (13) directly.

Chebyshev and (12)--(13) imply the exact necessary bound

    v(1-r/n)<=epsilon+2/sqrt(n).                           (14)

Equivalently,

    I <= h+2P[epsilon+2/sqrt(n)]/(1-r/n),                   (15)

when r<n. Consequently r=o(n), epsilon=o(1), h=o(n^(3/2)), and
P=O(n^(3/2)) already force I=o(n^(3/2)) without doing any surgery.
That stronger premise is therefore not a noncircular route to a leading
midpoint repair. More quantitatively, when T=aP for fixed a in [0,1),
equation (14) forces

    r>=n-2epsilon*n/(1-a)-4sqrt(n)/(1-a).

Thus vanishing-error capture of a whole fixed-fraction positive energy
layer requires almost full rank. The ramp in (1) avoids this obstruction: after small
independent noise around x*, the required fraction decreases at rate
P/I, which is larger than the loss of projection fraction.

## 5. A nonvacuous actual-sign ramp example

Let n=4^s and H=(J_4-2Id_4)^(tensor s). It is a symmetric regular
Hadamard matrix, with H1=sqrt(n)1 and operator norm sqrt(n). There is a
flat Boolean vector g orthogonal to 1 with Hg=-sqrt(n)g. Hollow H.
Fix theta>1 and independently flip its negative off-diagonal entries
with probability p=(theta-1)/sqrt(n). This is valid for large n.

The mean matrix is

    (1-p)H+pJ-[ (1-p)(-1)^s+p ]Id.

Its centered random perturbation has independent upper-triangle entries,
uniform bound 2, and row variance O_theta(sqrt(n)). Matrix Bernstein
therefore gives operator norm O_theta(n^(1/4)sqrt(log n)+log n)=o(sqrt n)
with probability tending to one. In particular there are actual full
signings A with

    P(A)/n^(3/2)->theta/2,
    R(A)/n^(3/2)->1/2,
    w/n^(3/2)->(theta+1)/4,
    I/n^(3/2)->(theta-1)/4.                                (16)

The positive and negative lower witnesses in (16) are 1 and g; the
matching upper bounds are spectral. For Pi=J/n put
q=x^T Pi x/n. The spectral bound on the orthogonal complement of 1 gives
uniformly over all spins

    H_A(x)/n^(3/2)<=1/2+((theta-1)/2)q+o(1)
      <=(theta+1)/4+((theta-1)/4)q+o(1),                    (17)

because q<=1. Thus (1) holds with rank one, kappa=1/n, and error
o(n^(3/2)), despite a fixed positive normalized midpoint. The stronger
all-high-state hypothesis (12) cannot hold, by (15).

The entrywise-incoherent surgery gives cap at most (theta+1)/4+o(1),
strictly below the original theta/2. These matrices are not claimed to
minimize width or absolute cap. Iterating this example's bias reduction
only approaches the half constant; it is not a strict-subhalf or a
convergence construction. Its role is to demonstrate that the ramp
criterion is genuinely realizable rather than implicitly assuming an
already vanishing midpoint.

The matrix Bernstein input here is the independent self-adjoint sum
bound in [Tropp, *User-Friendly Tail Bounds for Sums of Random Matrices*](https://tropp.caltech.edu/papers/Tro11-User-Friendly-preprint.pdf),
Theorem 1.4, printed page 4, with variance parameter max row variance;
the theorem was checked directly in the primary author PDF. The scalar
Bernstein step in Sections 1--2 follows directly from independent
bounded centered variables. The finite ramp/probability algebra is
checked separately by the accompanying exact-rational script.

The reproducible run
`computations/decisive_audit_midpoint_ramp_sign_surgery_2026_09_07.py`
passed 5120 masked constructions, 81920 exact mean-energy identities,
and 5120 exact product-law projector-variance checks. The independent
full proof audit by the bridge agent passed the finite construction,
the asymptotic constants, and the example; the subsequent variance
sharpening in (13) is the exact two-chaos calculation displayed above.
The checker also passed 300 sign-compatible exact-mean constructions,
including the trace/Frobenius bound S^2<=I^2*r/4.

## 6. Precise remaining actual-minimizer obligation

For each sufficiently large order n, one would need to SELECT a signing
A_n with W(A_n)<=W_n+o(n^(3/2)), orient its midpoint nonnegatively, and
prove one of the ramp certificate hypotheses above with vanishing
normalized error. For example, a sufficient package is a positive
contraction Pi_n of trace o(n), the ramp (1) with e=o(n^(3/2)), and

    A_ij(Pi_n)_ij>=0,       I*max_(i<j)|(Pi_n)_ij|<=n.

Equations (10a)--(10c) would then produce actual same-order signings
with Q<=W_n+o(n^(3/2)), proving M_n-W_n=o(n^(3/2)). The low-trace
or incoherent alternatives in Sections 1--3 give other sufficient
packages. If I=o(n^(3/2)) already, no certificate is needed.

No such package has been proved for selected width minimizers. Arbitrary
bounded-cap signings need not satisfy the first two regimes, as the
independently audited clique example shows; the sign-compatible regime
remains a hypothesis as well. The all-high-set capture substitute in
Section 4 is circular and must not be used. Finally, even equality of
the two normalized optimization problems up to o(1) would still require
a separate all-order limit theorem for width to establish the original
convergence question.
