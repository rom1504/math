# Soft flatness: an exact uniform variational bridge

Status: elementary proofs independently audited PASS by transfer_seeds.
This is an original-value reduction, NOT a thermodynamic-limit theorem.
The standard ingredient is the strong-variance smooth-max/type-2 bound.
The earlier archive `ar_matrix_rounding_literature_toolkit.md`, section 4,
already records variance-sensitive biased rounding with a Bernstein O(n)
remainder. The proof below removes that remainder and gives a dimension-
uniform penalty approximation. No novelty is claimed for the underlying
smooth-max inequality.

Let d=n(n-1)/2, n>=2, and let B be a hollow symmetric real matrix with
|B_e|<=1. Put q_n(B)=Q(B)/n^(3/2) and

    delta(B)=d^(-1) sum_e (1-B_e^2).

## 1. Exact strong-variance rounding

For independent centered real random variables Z_e of finite variance,
and a finite set V of vectors whose coordinates have absolute value at
most one,

    E max_(v in V) v.Z <= sqrt(2 log|V| sum_e Var Z_e).       (1)

Proof: let F_lambda(z)=lambda^(-1) log sum_v exp(lambda v.z).
Its second derivative in coordinate e is lambda times the softmax
variance of v_e, hence at most lambda. Starting from z=0 and revealing
one independent coordinate at a time, second-order Taylor with the
global curvature bound and mean zero gives

    E F_lambda(Z) <= log|V|/lambda + (lambda/2) sum Var Z_e.

The linear term has zero conditional expectation. Optimize lambda.
No boundedness or subgaussian assumption on the Z_e is needed: this
uses the SUM of coordinate variances, a strong variance parameter.
For zero variance the assertion follows by continuity or directly.

Take V={s(x_i x_j)_(i<j): x_1=1, s=+-1}. It has at most 2^n elements
(duplicates do not matter). Independently round B_e to signs A_e with
E A_e=B_e. Then Var(A_e-B_e)=1-B_e^2, and (1) yields

    E Q(A-B) <= sqrt(2 n log2 sum_e(1-B_e^2)).              (2)

In particular some rounded signing satisfies this bound. Consequently

    m_n <= q_n(B) + sqrt(a_n delta(B)),
    a_n=((n-1)/n) log2,  m_n=M_n/n^(3/2).                 (3)

This is exact at every n; there is no additive O(n^(-1/2)) term.

## 2. Uniform soft-penalty approximation

For tau>0 define

    e_n(tau)=min_(B in [-1,1]^d) [q_n(B)+tau delta(B)].

Compactness gives attainment. Flat signings are feasible, so e_n<=m_n.
Equation (3) and sqrt(a delta)-tau delta<=a/(4tau) give

    0 <= m_n-e_n(tau) <= a_n/(4tau) <= log2/(4tau).        (4)

Thus a proof that e_n(tau) has an all-order limit for every fixed tau
would prove the original convergence, by a uniform approximation
argument. This implication does not assume the original limit.

At any global soft minimizer B_tau, comparison with a flat optimum and
(3) gives tau delta <= sqrt(a_n delta), and hence

    delta(B_tau) <= a_n/tau^2.                            (5)

The inequality is automatic when the right side exceeds one.

For example, the optimum restricted to amplitudes in [1-epsilon,1]
differs from m_n by at most sqrt(a_n(2epsilon-epsilon^2)). The soft
penalty is stronger because it needs only average variance deficit.

## 3. Finite-temperature exactness and its limitation

Define the two-sided, projectively reduced pressure

    f_(beta,n)(B)=(beta n)^(-1)
      log sum_(x_1=1,s=+-1) exp[(beta/sqrt(n)) s q_B(x)].

Here q_B(x)=sum_(i<j) B_ij x_i x_j.

Then q_n(B)<=f_beta(B)<=q_n(B)+log2/beta and

    partial_e^2 f_beta <= beta/n^2.

Therefore f_beta(B)+tau delta(B) is concave in every individual
coordinate whenever

    tau >= beta(n-1)/(4n).                               (6)

Successively replacing each coordinate by an endpoint that does not
increase the function proves the exact identity

    min_(B in [-1,1]^d) [f_beta(B)+tau delta(B)]
       = min_(A flat) f_beta(A)                          (7)

under (6). This gives a second proof of (4), choosing beta=4tau n/(n-1).
It also explains why this regularization is not a free convexification:
the penalty is concave, and precisely when coordinate curvature is
dominated, the original discrete optimization is recovered exactly.
No min-max interchange is justified by (7).

## 4. A continuous integral with uniform zero-temperature control

One can also remove the outer minimum without a dimension-dependent
log n entropy error. Let U_d be uniform probability on [-1,1]^d and

    L_n(alpha,tau)=-(alpha d)^(-1) log
      integral exp[-alpha d (q_n(B)+tau delta(B))] dU_d(B).

Clearly L_n>=e_n. For 0<epsilon<=1, let B_* minimize the soft objective,
let U have independent uniform[-1,1] entries, and set
B_epsilon=(1-epsilon)B_*+epsilon U. By (1),
E q_n(U)<=sqrt(log2/3), so with probability at least one half,
q_n(U)<=K:=2sqrt(log2/3). On that event

    q_n(B_epsilon)<=q_n(B_*)+epsilon K,
    delta(B_epsilon)<=delta(B_*)+4epsilon.

The affine image is a box inside the cube, with U_d-volume epsilon^d.
Hence its good half has volume at least epsilon^d/2, and

    0 <= L_n-e_n <= epsilon(K+4tau)
                    +log(1/epsilon)/alpha+log2/(alpha d). (8)

For alpha(K+4tau)>=1, choosing epsilon=1/[alpha(K+4tau)] gives

    L_n-e_n <= [1+log(alpha(K+4tau))]/alpha
                  +log2/(alpha d).                      (9)

Thus all-order limits of L_n(alpha,tau) at fixed alpha,tau would also
prove original convergence: first alpha tends to infinity at fixed tau,
then tau tends to infinity. The limit theorem for this continuous
integral has NOT been proved or imported.

## 5. Why normalized copying does not supply that theorem

For a hollow flat A of order n and integer k>=1, the matrix

    C=(A tensor J_k)/sqrt(k)

has order kn, lies in the coefficient cube, and satisfies q_(kn)(C)=m(A).
The quadratic equality holds because the block magnetizations lie in
[-k,k] and the original hollow quadratic maximum is attained at the
cube vertices. But its variance occupancy is

    1-delta(C)=(n-1)/(kn-1),

so delta(C) tends to one as k grows. Its soft penalty approaches tau,
not zero. This defeats this particular all-order copying operation;
it does not prove that every possible local-penalty realization fails.
The positive variance-density / flat-sign recovery problem survives.

## 6. Genuine remaining analytic obstacle

The continuous model has a convex homogeneous norm q_n, but the term
-tau d^(-1)||B||_2^2 is concave. Ordinary log-concave matrix-potential
arguments do not apply. The q_n observable is also a Boolean-spin
supremum, not a normalized spectral trace; standard spectral-measure
large-deviation principles do not identify it. The Gaussian radial-
collapse defect of the unpenalized outer model is excluded by (5),
but this does not establish interpolation, subadditivity, or recovery.

The strongest honest output here is the exact uniform bridge (4), the
quantitative occupancy (5), and the continuous-integral estimate (9).
They do not themselves prove a new inequality comparing different orders.

## 7. Exact oversaturation / variance-relocation barrier

Suppose N=m+l and B is ANY block-diagonal real symmetric matrix on that
split. For any flat hollow signing A, let D=A-B and let C=A_(S,T) be its
m by l cross block. The exact half-width identity gives

    Q(D) >= W(D) >= ||C||_(infinity->1).                  (10)

Indeed q_D(x)-q_D(y), with y obtained by flipping all spins in T, is
twice x_S^T C x_T. Averaging one of the two independent spin vectors
and optimizing the other yields the finite bound

    Q(A-B) >= max{m E|S_l|, l E|S_m|},                   (11)

where S_k is a sum of k independent fair signs. This is independent of
the internal blocks of B. For a fixed proportional split it is

    (sqrt(2/pi)+o(1)) max{m sqrt(l),l sqrt(m)}.

For balanced blocks this is (1/(2sqrt(pi))+o(1))N^(3/2).
Thus no extension of the difference-norm rounding estimate (2) can
round a zero-bridge endpoint at subleading cost, even when there are
only TWO blocks, their magnitudes are bounded by sqrt(2), and they are
scaled actual optimizers. Any successful comparison must exploit
cancellation in Q(A) itself, not prove Q(A-B)=o(N^(3/2)).

The desired normalized children B=diag(sqrt(N/m)A_m,
sqrt(N/l)A_l) have total squared entry mass approximately that of a
flat parent, but lie outside the coefficient cube. Unbiased flat
rounding is then impossible coordinatewise. The signed variance
deficit 1-average(B_e^2) can be nearly zero only because the missing
cross-edge variance is canceled by internal oversaturation.

Nor does any nonnegative local penalty p(|b|), with p(1)=0 and
p(0)>0, vanish on this endpoint: the zero cross block alone contributes
asymptotic average penalty 2(m/N)(l/N)p(0). This rules out direct use
of such a local penalty to make the zero-bridge endpoint cheaply
feasible. It is not a no-go theorem for every nonlocal amplification
or for a TOTAL-cap comparison of two optimized children.
