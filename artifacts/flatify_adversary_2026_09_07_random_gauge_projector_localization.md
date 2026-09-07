# Random-gauge projector localization through almost-linear rank

2026-09-07. New actual full-sign bridge operation. The strongest version is Section 6: relative-leverage contraction improves the error exponent to 1/4, with no masking. Sections 1--5 preserve the original random-gauge/masking construction and the director's 1/6 optimization as a separately checked alternate proof. This is not a payment of the internal child energies.

Independent and construction agents separately reconstructed the argument and passed the relative-leverage refinement. The construction audit is preserved in `flatify_construct_2026_09_07_nearlinear_projector_audit.md`.

## Theorem at Hadamard orders

Let H be any order-n full sign Hadamard. Let P,Q be arbitrary orthogonal projections with total rank r>=1, and set

    a=r log(24n²)/n.

When a tends to zero, there exists a FULL sign matrix C such that, simultaneously for every Boolean x,y,

    |x^T C y|<=sqrt(n)||(I-P)x||||(I-Q)y||
                 +O(n^(3/2)a^(1/6)+n).                (1)

Constants are absolute. In particular every r=o(n/log n) is allowed; r<=n^(1-delta) gives a power saving, up to the displayed logarithm. The projections are fixed before the random gauges are selected but may otherwise be arbitrary and depend on the actual children.

## 1. Low-leverage rows and columns

Choose 0<mu<=1. Let I={i:P_ii<=mu} and J={j:Q_jj<=mu}. The total number b of omitted rows and columns satisfies

    b<=r/mu.

Let L=log(24n²), choose independent uniform sign vectors s,t, and put H'=D_s H D_t. This remains a full sign Hadamard. Write

    T=(I-P)H'(I-Q).

It has operator norm at most sqrt(n), kills both designated subspaces, and has exact Frobenius deficit

    n²-||T||F²=n r-||PH'Q||F²<=nr.                  (2)

## 2. Three simultaneous entrywise estimates

For i in I and every j, (PH')_ij is, up to the factor t_j, a Rademacher sum with squared coefficients sum_u P_iu²=P_ii<=mu. Hoeffding and a union bound give

    max_(i in I,j) |(PH')_ij|<=sqrt(2mu L)             (3)

except with probability at most 2n² exp(-L). Likewise, except with the same failure bound,

    max_(u,j in J) |(H'Q)_uj|<=sqrt(2mu L).            (4)

For the double term, first expose t and define V_uj=sum_v H_uv t_v Q_vj. The event (4) is exactly the bound |V_uj|<=sqrt(2mu L), independent of s. Conditional on that event, for i in I,j in J,

    (PH'Q)_ij=sum_u P_iu s_u V_uj

is a Rademacher sum with squared coefficient norm at most 2mu²L. Another Hoeffding union bound, now conditional on t, gives

    max_(i in I,j in J)|(PH'Q)_ij|<=2mu L             (5)

except with conditional probability at most 2n² exp(-L). The event (3) need not be independent of (5); a union bound suffices. All three estimates therefore hold with probability at least 1-6n²exp(-L)>=3/4.

Fix such gauges. On I x J,

    |T_ij-H'_ij|<=epsilon,
    epsilon=2sqrt(2mu L)+2mu L.                       (6)

This is the only place randomness of the Hadamard phases is needed.

## 3. Masking, scaling, and actual sign recovery

Let R_I,R_J be coordinate projections onto I,J, put U=R_I T R_J, and set B=U/(1+epsilon). Every entry of B lies in [-1,1], including its zero entries on the masked rows or columns.

Let b_L,b_R be the two mask sizes. Since ||T||op<=sqrt(n),

    beta(T-U)<=n(sqrt(b_L)+sqrt(b_R))<=n sqrt(2b).

Moreover beta(U-B)<=epsilon n^(3/2), because ||U||op<=sqrt(n). Each T row and column has squared norm at most n, so (2) gives

    n²-||U||F²<=n(r+b).

Independently round EVERY B_ij to a sign with mean B_ij. The variance budget, including the masked zero means, is

    V=sum_(i,j)(1-B_ij²)
      <=2epsilon n²+n(r+b).                          (7)

For u=(2n+2)log 2, Bernstein and a union bound over all 2^(2n) Boolean bilinear pairs give, with positive probability,

    beta(C-B)<=sqrt(2Vu)+(4/3)u.                     (8)

The factor two for the two tails is included in this choice of u. All entries of C are actual signs. Combining (6)--(8),

    beta(C-T)/n^(3/2)
      <=sqrt(2b/n)+epsilon
          +O(sqrt(epsilon+(r+b)/n)+n^(-1/2)).        (9)

Finally x^T T y=((I-P)x)^T H'((I-Q)y), which proves the desired residual bound once (9) is paid. No norm bound on P x in the coordinate cube is assumed.

## 4. Optimizing the threshold

Set a=rL/n and choose

    q=mu L=a^(2/3),  mu=a^(2/3)/L.

For sufficiently small a this is in (0,1]. Then

    b/n<=a/q=a^(1/3),
    epsilon=O(sqrt(q)+q)=O(a^(1/3)).

Every leading term in (9) is O(a^(1/6)), proving (1). This improves the initial unoptimized a^(1/8) choice. For rank zero simply take the original Hadamard bridge.

## 5. All-order extension and exact scope

Embed the two subspaces by zeros into any Hadamard order h>=n with h/n tending to one. Apply the construction at h and restrict C to the original n rows and columns. Its estimate holds on the real cube by separate linearity, so the zero-padded vectors are allowed and their original projection norms are unchanged. The leading factor becomes sqrt(h); its excess costs at most n(sqrt(h)-sqrt(n)).

The already audited all-order Paley/Baker--Harman--Pintz construction supplies h-n=O(n^(21/40)), hence that excess is O(n^(41/40)), smaller than the new error for r>=1. Thus at EVERY sufficiently large order the theorem holds with total error

    O(n^(3/2)(r log n/n)^(1/6)+n^(41/40)).

One may use the older nuclear/projector rounding bound O(n^(5/4)sqrt(r)+n^(41/40)) when it is smaller. No computation, certificate, or favorable child-energy assumption is required to select the bridge in this theorem.

The independent Rademacher-width theorem, strengthened by the fourth spectral moment bound, has an obstruction for every r=o(n) to paying child energies and this residual norm separately. Accordingly (1) is a stronger actual bridge construction, NOT a closure of the original recurrence, and not a claim that a scalar residual envelope now succeeds.

## 6. Stronger theorem: relative leverage contraction, without masking

The SAME exposure argument proves, simultaneously over ALL i,j, with probability at least 3/4,

    |(PH')_ij|<=sqrt(2L P_ii),
    |(H'Q)_ij|<=sqrt(2L Q_jj),
    |(PH'Q)_ij|<=2L sqrt(P_ii Q_jj).                 (10)

Indeed use the individual row/column variance instead of mu in (3)--(5). After exposing t, every V_uj is bounded by sqrt(2L Q_jj); the conditional variance of the final sum is at most 2L P_ii Q_jj. A zero leverage makes the relevant row or column identically zero and requires no tail estimate. The same three union bounds apply.

Define nonnegative leverage amplitudes

    a_i=sqrt(2L P_ii), b_j=sqrt(2L Q_jj).

Equation (10) and the full-sign entries of H' give the EXACT factorized majorant

    |T_ij|<=1+a_i+b_j+a_i b_j=(1+a_i)(1+b_j).

Consequently D_ii=(1+a_i)^(-1), E_jj=(1+b_j)^(-1), and B=DTE are feasible sign means, without masking any row or column.

For every Boolean x,

    ||(I-D)x||²=sum_i[a_i/(1+a_i)]²<=2L r_L,

and similarly on the right. Since ||T||op<=sqrt(n),

    beta(T-B)<=n(sqrt(2Lr_L)+sqrt(2Lr_R))
                 <=2n sqrt(Lr).                     (11)

The Frobenius deficit of T is still at most nr. Each of its rows and columns has squared norm at most n, so

    V=n²-||DTE||F²
      <=nr+n sum_i(1-D_ii²)+n sum_j(1-E_jj²).

Using 1-(1+a)^(-2)<=2a and Cauchy--Schwarz on the leverage sums yields

    V<=nr+4n sqrt(Lnr).                              (12)

Independent mean-preserving sign rounding and (8) now give

    beta(C-T)<=2n sqrt(Lr)
                 +sqrt(2u[nr+4n sqrt(Lnr)])+(4/3)u.

For a=rL/n tending to zero this is

    O(n^(3/2)a^(1/4)+n),                             (13)

which is stronger than the earlier a^(1/6) estimate. In particular the actual full-sign bridge satisfies

    |x^T C y|<=sqrt(n)||(I-P)x||||(I-Q)y||
                  +O(n^(3/2)(r log n/n)^(1/4)+n)

uniformly over all Boolean x,y. Every r=o(n/log n) is still permitted, and r<=n^(1-delta) now gives normalized error O(n^(-delta/4)(log n)^(1/4)).

The all-order zero-padding argument of Section 5 gives the same statement at arbitrary large n with the additional O(n^(41/40)) term. This is a direct amplitude-cover contraction of a randomly gauged projected Hadamard target. It does not invoke an unproved incoherence property of the children or their chosen subspaces.
