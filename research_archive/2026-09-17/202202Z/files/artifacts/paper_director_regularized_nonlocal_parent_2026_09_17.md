# A favorable restricted parent from actual signs and nonlocal columns

2026-09-17. Director combination. **Independently reconstructed by the
discrepancy researcher, including both preparation events and exponents.**
This combines the same-order information regularizer,
nonlocal Boolean eigenlaws, Gaussian--Poisson response comparison, and
the scheduled exact-sign covariance repair. It solves a restricted
parent VALUE problem on a concrete family; it does not solve the full
parent problem or claim that this family minimizes the original cap.

## 1. Precise result

Let p=2^a, n=p^2, and let A=F-I be the quadratic Walsh signing in the
[nonlocal law theorem](paper_localization_nonlocal_sign_response_2026_09_17.md).
For all sufficiently large n there is a genuine same-order full signing W
such that, with natural logarithms and universal constants,

    Q(W)<=Q(A)+C n^(5/4)(log n)^(1/2),
    T_n=n^(3/4)/(log n)^(1/2),
    log|E_W(T_n)|<=C n^(3/4)(log n)^(1/2).          (1)

Here E_W(T)={x:Q(W)-|H_W(x)|<=T}; changing the fixed constant multiplying
T_n only changes C. For EVERY fixed epsilon>0, every integer sequence
q/n->epsilon, and EVERY full-sign child D of order q, one can choose an
actual sign bridge B so that

    max_(x in E_W(T_n), y Boolean)
       |H_W(x)+H_D(y)+x^T B y|
    <= Q(W)+(1/sqrt(pi)) q sqrt(n)+Q(D)
              +C_epsilon n^(11/8)(log n)^(1/4).    (2)

All new spins, both absolute polarities, and complete sign support are
retained. The restriction on old spins is explicit and essential. The
coefficient 1/sqrt(pi)<3/4 is obtained from actual nonlocal columns;
it is neither an independently Gaussian bridge nor a covariance-only
Jensen bound. The order sequence here is n=4^a, not every n.

## 2. Prepare one actual signing without losing its old energy geometry

Set k=floor(sqrt(n) log n). The read-two construction gives W by
rewriting a k-vertex star, with its random rewritten-edge cap at most
C n sqrt(k), and its conditional fresh-vertex increment at most
C n/sqrt(k). The two events intersect in that proof.

Here we need the extra original-energy comparison; it is NOT automatic
for an arbitrary signing. Since ||A||op=p+1, the removed original cross
block has cap at most (p+1)sqrt(n k), and its removed internal block
has cap at most (p+1)k/2. Therefore this SAME W satisfies

    Q(W-A)<=C n sqrt(k).                            (3)

The cap Q is also used for the weighted difference matrix in (3).
Triangle inequality gives |Q(W)-Q(A)|<=Q(W-A). Thus every x in
E_W(T_n) belongs to the original A-nearcode at relative deficit

    eta_n=[T_n+2Q(W-A)]/n^(3/2)
          <=C n^(-1/4)(log n)^(1/2).                (4)

The regularizer's width bound b(E_W(T))<=T+C n/sqrt(k), followed by
its Sauer counting bound, gives (1). No count of all Hadamard
ground states or all fixed-relative-window nearstates is assumed.

## 3. Physical response and simultaneous fluctuations

Use the fair fixed-point-free pole mixture P_0 for q-r columns and
the identity-law repair for r=round(q/p) columns, with the deterministic
schedule from the nonlocal theorem. At p|q the aggregate covariance
is exactly qI; otherwise its operator discrepancy is O(p). No averaging
over fractional columns is substituted for the physical bridge.

The [Gaussian--Poisson theorem](paper_director_matching_gaussian_poisson_2026_09_17.md),
with its weighted-Palm refinement, gives uniformly on the code in (1)

    E_(P_0)|h dot x|/sqrt(n)
       <=1/sqrt(pi)+4sqrt(eta_n)+O(n^(-1/8)).        (5)

The harmless O(1/p) pole-weight difference is absorbed in that error.
The all-query exponential-moment bound for P_0 is uniform in p. Hence,
putting H=log|E_W(T_n)|, a union bound over this ENTIRE code gives

    max_x sum_(j<=q-r)|h_j dot x|
      <=(q-r)sqrt(n)[1/sqrt(pi)+4sqrt(eta_n)+O(n^(-1/8))]
         +C sqrt(n)[sqrt(q(H+1))+H+1]              (6)

with a fixed positive probability. The separate identity columns have
their ENTIRE bipartite cap at most O_epsilon(n^(5/4)), also with a
fixed high probability. These events intersect; no small-code bound
is needed for the repair part.

Using (1),(4), the leading response error and the square-root term in
(6) are both O_epsilon(n^(11/8)(log n)^(1/4)). The linear H term,
the identity repair, and the preparation cost have smaller exponents.
Pay the exact Q(D), maximize every new spin, and bound old |H_W| by
Q(W). This proves (2).

## 4. What the combination removes, and what it does not

On this prepared actual family, there is now a favorable leading mean,
a quantitative power-saving response error, physical sign support,
all-new-spin control, and a subexponential full microscopic nearcode.
These did not follow from covariance rounding, entropy regularization,
or scalar Gaussian replacement alone.

There remains an independent ESCAPE obligation: (2) does not bound
old words outside E_W(T_n). Its error is larger than T_n, so an
outside word cannot be excluded just by comparing these two scales.
All-order realization of arbitrary near-minimizing seeds is also absent.
No new upper bound on M_n, convergence theorem, or nonconvergence theorem
is inferred. The published full-parent certificate keeps the missing
all-energy counts visible instead of silently dropping those words.
