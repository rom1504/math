# Actual optimizer noise response: exact identities and a scope limit

This note addresses the child-specific route, not another universal kernel.
Let A be a hollow symmetric full-sign matrix of order n, and switch and,
if necessary, negate it so that H_A(1)=Q(A)=q. Define d_i=sum_j A_ij.

## 1. All-cut positivity and exact response

For every vertex subset S,

 cut_A(S)>=0, H_A(1-2 1_S)=q-2 cut_A(S),
 0<=d_i<=n-1, and sum_i d_i=2q.

Let z_i be independent signs with E z_i=rho, and put xi_i=z_i-rho,
v=1-rho². The exact orthogonal-chaos expansion is

 H_A(z)=rho²q+rho sum_i d_i xi_i+sum_{i<j}A_ij xi_i xi_j.

The three displayed terms are pairwise orthogonal after centering. Hence

 E H_A(z)=rho²q,
 Var H_A(z)=rho² v sum_i d_i²+v² binom(n,2).           (1)

In particular, q<=C n^(3/2) implies sum d_i²<=2(n-1)q and variance
O_C(n^(5/2)). For any fixed delta<1/4, Chebyshev gives

 P(|H_A(z)-rho²q|>n^(3/2-delta))=O_C(n^(-1/2+2delta)).

Thus the deterministic old-child slack (1-rho²)q is recovered with a
power-saving error for one random noisy ground state, uniformly over
actual optimal children. This is NOT uniform over an exponentially large
cloud. Turning this fact into an all-spin parent cap cannot discard the
probability/entropy distinction.

For uniform noise of exactly s flipped vertices, the mean has the equally
explicit form

 E H_A(z)=[1-4s(n-s)/(n(n-1))]q.

No uniform positive lower cut bound at a fixed Hamming distance follows
from the stated all-cut constraints alone.

## 2. Near-optimality cannot give a stronger variance scale

Take an ACTUAL optimal sign matrix A_s, switch a positive ground state to
the all-ones vector, and adjoin r vertices with every newly introduced
edge positive. Call the resulting order n=s+r matrix C. Then

 Q(C)=M_s+rs+binom(r,2).                              (2)

The triangle inequality gives the upper bound and the all-ones vector
attains it, so (2) is exact. Also M_n>=M_s: restriction to any s vertices
is cap-contractive by averaging all outside spins. Consequently

 0<=Q(C)-M_n<=rs+binom(r,2).

Taking r=o(sqrt n) makes C asymptotically near-optimal with o(n^(3/2))
excess. Nevertheless every adjoined vertex has ground degree n-1, so

 sum_i d_i(C)²>=r(n-1)².

For any epsilon<1/2, choose n^epsilon<<r<<sqrt n. At any fixed
0<|rho|<1, (1) then has variance much larger than n^(2+epsilon).
Therefore asymptotic near-optimality, all-cut positivity, and the correct
leading cap scale do NOT imply an O(n^(2+epsilon)) noisy-ground variance
bound. The best general n^(5/2) exponent can be approached up to an
arbitrarily slowly vanishing factor within this near-optimal class.

This is not an actual-minimizer counterexample. It shows precisely why a
new variance regularity theorem would have to use exact optimality, not
merely an o(n^(3/2)) optimization gap. An error n^(3/2)/log n is also too
large to be silently treated as the requested fixed power saving.

## 3. A bridge cannot cancel the fluctuation in both child orientations

This section keeps the two internal child matrices fixed. It constrains
bridge-only completion, not the authorized broader flatification operation
which may also change old internal edges globally.

Switch selected positive ground states of children A,D to1. For a bridge
C and independent biased noises u,v with common mean rho, let

 X=H_A(u)+H_D(v), Y=u^T C v.

The two independently reversible child orientations give parent energies
X+Y and X-Y. Exactly,

 [Var(X+Y)+Var(X-Y)]/2=Var(X)+Var(Y)>=Var(X).          (3)

In expanded form, with v0=1-rho², row degrees d_A,d_D and bridge row/column
sums c=C1,e=C^T1,

 Var(X+eta Y)=rho² v0[||d_A+eta c||²+||d_D+eta e||²]
              +v0²[binom(m,2)+binom(n,2)+mn].

Thus bridge row sums can cancel local-field fluctuations for one selected
orientation, but the opposite child reversal restores or worsens them.
The average over the two orientations has only positive extra terms.
This is a direct algebraic obstruction to paying for noisy extremizers
by a one-orientation variance cancellation.

There is also an exact center constraint. If z=1^T C1, then the two ground
centers alone force

 Q(parent)>=Q(A)+Q(D)+|z|.

For equal children, choosing bridge row sums -d_A and column sums -d_D
would give z=-2Q(A), hence a parent cap at least4Q(A), already exceeding
the desired leading weighted-child level2sqrt2 Q(A). This is the unpaid
bridge-cancellation problem in a concrete local-field design.

## 4. What remains unproved

The exact mean response supplies real, optimizer-specific slack. What is
missing is an exponential-scale, joint cloud bound that retains this
slack and controls both reversals. Polynomial probability concentration
for a fixed ground state, a positive-cut condition alone, and cancellation
of one orientation's linear response do not provide it. No construction
satisfying the favorable flatification target is claimed in this note.
The bridge-only constraints above are not no-go statements for global
old-edge surgery, which remains within the user's authorization.

## 5. An actual minimizer can have an indefinite ground-state Laplacian

One possible strengthening would turn all-cut nonnegativity into a
quadratic lower bound for nonuniform conditional means. It is false even
for an actual minimizer. Take

 A=[[0,1,1,1,1],
    [1,0,1,1,-1],
    [1,1,0,-1,-1],
    [1,1,-1,0,1],
    [1,-1,-1,1,0]].

Its all-ones energy is4 and degrees are(4,2,0,2,0). Every cut has weight
between0 and4: by complements one checks only singletons and pairs,
using cut({i,j})=d_i+d_j-2A_ij. Hence Q(A)=4. Any order5 sign matrix has
E H²=10, while H is even-valued, so Q<4 would force Q<=2 and contradict
the second moment. Thus this is an ACTUAL optimal signing, not a planted
near-minimizer.

For L=diag(d)-A and v=(-2,-2,-2,-1,0), one has v^T L v=-2.
In general, nonuniform independent means m_i satisfy the exact identity

 q-H_A(m)=(1/2)sum_i d_i(1-m_i²)+(1/2)m^T L m.

Taking m=v/2 makes the second term -1/4. The actual slack is1/2, whereas
the putative nonnegative-Laplacian lower bound would give3/4. Hence
cut positivity cannot justify dropping this signed quadratic term, even
under exact optimality. This finite example does not refute a suitably
weakened asymptotic statement; it excludes the direct PSD shortcut.

## 6. Projection balancing: a real small-rank correction, no cap inequality

We also tested whether a child could be made row-balanced before using it
as a balanced port seed. For a Boolean vector v let P_v=I-vv^T/n.
The difference A-P_v A P_v has rank at most2. Under a uniform random v,
E||Av||²=||A||F²=n(n-1), so its nuclear norm is typically O(sqrt n),
using

 A-P_v A P_v=(Av)v^T/n+v(Av)^T/n-(v^TAv)vv^T/n².

The previously proved target-contraction surgery could round such a
correction with an O(n^(5/4)) cap-error scale. The missing point is the
cap of the target P_v A P_v itself: P_v is not a cube contraction.

For exact checking, its full quadratic energy at Boolean x is

 H_A(x)-(x^Tv)(x^TAv)/n+(v^TAv)(x^Tv)²/(2n²).

We exhaustively evaluated this finite expression for the single stored
optimal representatives of orders3 through10. Minimum projected caps
were respectively1.333333,2.5,4,5,6.367347,8,11.160494,13; original caps
were3,4,4,5,9,10,12,13. These numbers are a diagnostic, not a theorem or
an exhaustive census of optimal signings. In particular every ground
direction of the stored A5 gives projected cap4.64>4, and every ground
direction of stored A8 gives10.875>10. Thus projecting away a selected
ground state is definitely not automatically cap-contractive, even for
actual minimizers.

Choosing some other good v might conceivably help balancing, but no
uniform existence inequality was proved. Nor does a small-rank balancing
operation itself fill the leading number of zero bridge edges in the
original favorable flatification task. We did not promote this missing
projection inequality to a new convergence claim.
