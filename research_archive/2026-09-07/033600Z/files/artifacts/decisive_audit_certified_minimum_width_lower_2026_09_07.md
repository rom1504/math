# The certified .4333221116640807 lower bound also bounds minimum width

Status: proved by the freshly reconstructed marked-response theorem.
This is a consequence for the original half-range optimization, not a
claim that its optimum equals the minimum absolute cap.

For a hollow symmetric sign matrix of order n, write

    H_A(x)=x^T A x/2,
    P(A)=max_x H_A(x),       R(A)=max_x -H_A(x),
    Q(A)=max(P(A),R(A)),     W(A)=(P(A)+R(A))/2,
    M_n=min_A Q(A),          W_n=min_A W(A).

All maxima here are over Boolean spins. The same endpoint extrema apply
on the full cube by independent coordinate rounding and hollowness.

## 1. The result and its remaining gap

Let

    c_* = .433322111664080753415812928897579346558634648033693413106996.

The finite rational interval certificate proves

    liminf_n W_n/n^(3/2) >= c_* > .4333221116640807.             (1)

In particular, together with the all-order original upper bound,

    0 <= limsup_n (M_n-W_n)/n^(3/2)
      <= .494515125-.4333221116640807
       = .0611930133359193.                                 (2)

Equation (2) is only a numerical bound on the possible gap. It does not
prove that gap vanishes. It also does not prove convergence of either
normalized optimum.

The proof below identifies exactly why every step of the cap lower
bound continues to hold for width. The analytic and interval details
were freshly reconstructed in
`decisive_audit_fresh_full_lower_chain_2026_09_07.md`; no earlier audit
verdict or optimality hypothesis is used as a replacement for them.

## 2. Spectral deletion is directly controlled by width

The mean of H_A on the cube is zero. Thus P,R>=0 and

    W(A)<=Q(A)<=2W(A).

There is a slightly sharper direct route to the relevant majorant.
For Boolean x,y put u=(x+y)/2 and v=(x-y)/2. Both lie in the cube and

    x^T A y = 2[H_A(u)-H_A(v)].

Consequently the bilinear norm beta(A) obeys

    beta(A)<=4W(A).                                         (3)

The simultaneous Grothendieck diagonal majorant therefore gives

    D>=A,-A,       Tr D<=K_G beta(A)<=4K_G W(A).              (4)

For a sequence with W(A)<=C n^(3/2), delete the coordinates with
D_ii>4K_G C sqrt(n)/epsilon. At most epsilon*n are deleted. The retained
principal signing A' of order n'>=(1-epsilon)n has

    ||A'/sqrt(n'-1)||op<=L(C,epsilon),                        (5)

with a fixed finite L when epsilon is fixed. It remains a hollow sign
matrix with exact row squared norm n'-1.

Principal width monotonicity is exact. For any retained spin x, average
the omitted coordinates as independent fair signs. The conditional
energy mean equals H_A'(x). Hence

    P(A)>=P(A'),       R(A)>=R(A'),       W(A)>=W(A').         (6)

There is no cancellation or midpoint assumption in (6).

## 3. The marked response certifies a range, not just an absolute cap

Fix L and consider any sequence satisfying (5). The full nonlinear
marked-response theorem, reconstructed in Sections 2--5 of the fresh
lower audit, says the following. For each fixed admissible Gaussian
pair F,H, with F odd, H even and nonnegative, and |F|+H<=1, the literal
finite signing constructions satisfy

    liminf_n (1/n) sum_i E[H_i |(BF)_i|]
      >= J(F,H):=E H E_N |K_F+tau N|,                        (7)

where B=A/sqrt(n-1), K_F is the inverse first-chaos projection, and
tau^2=||F-P_1F||_2^2. Ordered finite approximations are used before the
matrix limit. In particular, (7) is a theorem for every fixed-operator-
bounded signing sequence, not merely original cap minimizers.

The two literal feasible mean vectors

    mu_+=F+H sign(BF),       mu_-=-F+H sign(BF)

lie in the cube pointwise. Their energy difference is exactly

    H_A(mu_+)-H_A(mu_-)
      =2 sqrt(n-1) sum_i H_i |(BF)_i|.                       (8)

Both energies lie in [-R(A),P(A)], so the left side of (8) is at most
P(A)+R(A)=2W(A). Taking expectations and then the limit gives

    liminf_n W(A_n)/n^(3/2)>=J(F,H).                        (9)

No demand that the two energies have opposite signs has been made.
Their difference alone proves (9). This is exactly the factor of two
used in the original cap proof, now with the sharper available interval.

## 4. The same frozen rational policy and the same limits

The certified policy is the fixed 21-anchor/degree-200 rich core and
its second two-Gaussian ternary rectangle response at theta=1. Its
conditional inverse correlation, weighted hole-mask Jensen bound,
degree-20 remainder, and rational outward interval arithmetic are all
unchanged by replacing Q by W. They supply J(F,H)>=c_*.

The exact replay used for the fresh audit was

    .venv/bin/python computations/resumed_response_rich_core_birth_certificate_2026_09_06.py --policy computations/results/resumed_response_rich_core_optimized_rectangle_policy_2026_09_06.json --theta 1 --target 4333/10000 --output /home/math/quadra/tmp/decisive_audit_rich_lower_replay_2026_09_07.json

Its output is byte-identical to the recorded certificate, with SHA256

    c6ab5199c366952bbfc831d26d1527149f2adaeb7120c4c5ec9b1d3ee7f86dcb.

To remove (5), first fix epsilon, and hence L, and finish all finite
construction / matrix-order / polynomial and smoothing approximation
limits from the lower theorem. Apply (6) and (9) to the retained cores:

    liminf_n W(A_n)/n^(3/2)>=(1-epsilon)^(3/2)c_*.

Only then send epsilon to zero. Sequences of unbounded normalized width
do not threaten a lower bound; on any subsequence with bounded width,
the preceding argument applies. Equivalently, one may choose actual
width minimizers, which are bounded because W_n<=M_n, without importing
any stationarity or near-minimality property of those matrices.

Thus (1) is a universal signing lower bound. Width minimality is used
only to name its optimized consequence, not in the proof mechanism.
