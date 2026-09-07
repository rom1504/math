# Independent audit: universal near-stability entropy deficit

2026-09-07. **PASS.** Reconstruction of the director's uniform-field core lemma and its total-violation strengthening. Constants below are deliberately nonoptimal. All large-order qualifications are uniform in the external field.

## 1. Operator-bounded core with an arbitrary external field

Let B be symmetric hollow of order k, with off-diagonal entries of magnitude 1/sqrt(n), k>=n/2, and ||B||op<=L, where L>=1 is fixed. For a fixed arbitrary b in R^k and a uniform Boolean x, define

    I(x)=sum_i [-x_i(Bx+b)_i]_+.

Then, for some delta(L)>0 and all sufficiently large n,

    Pr{I(x)<=k/512}<=exp(-delta(L) k).                 (1)

Set T=128L and H={i:|b_i|>T}.

If h=|H|>=k/4, every mismatch x_i != sign(b_i) either has |(Bx)_i|>T/2, or contributes at least T/2 to I. The first set has size at most 4L^2 k/T^2=k/4096. On the event in (1), the second has size at most 2I/T<=k/(256T). Their sum is at most 9k/32768<k/1024. Thus the number of possible x on that event is at most

    2^(k-h) sum_(j<=k/1024) binom(h,j).

Since h>=k/4, the relative Hamming radius is at most 1/256, giving a fixed positive entropy deficit proportional to k. This estimate is deterministic; no independence between the field Bx and x is assumed.

If h<k/4, condition on all spins in H. Put R=H^c, k'=|R|>3k/4, and b'=b_R+B_(R,H)x_H. Then

    ||b'||_2 <=129L sqrt(k),  ||B_R||op<=L.

The remaining spins are still independent fair signs. Define the convex function

    F(x_R)=||B_R x_R+b'||_1.

Its Euclidean Lipschitz constant is at most L sqrt(k'). Each coordinate of B_R x_R is symmetric, so shifting its mean cannot decrease its expected absolute value. The elementary Khintchine lower bound gives

    E F >= k' sqrt((k'-1)/(2n)) >=k/8

for all sufficiently large n. Apply the previously audited Talagrand convex zero-event estimate to (F-k/16)_+. Its mean is at least k/16 and its Lipschitz constant is at most L sqrt(k), whence

    Pr{F<=k/16}<=exp[-k/(4096L^2)].                    (2)

The identity

    F-[x_R^T B_R x_R+b' dot x_R]=2 I_R

is exact. On I<=k/512 and F>k/16, its bracket is at least 15k/256>k/32. Hence either x_R^T B_R x_R>k/64 or b' dot x_R>k/64.

The quadratic form is centered because B_R is hollow. Its operator norm is at most L and squared Frobenius norm at most L^2 k. Hanson--Wright therefore gives an exp[-c(L)k] upper tail at k/64. The linear form has squared coefficient norm at most 129^2 L^2 k, so the same conclusion follows from the elementary Hoeffding MGF. Adding these tails and (2), then reducing the exponent constant to absorb the finite prefactor, proves (1). Every estimate is uniform in the conditioned heavy spins and in b.

I independently checked the primary [Rudelson--Vershynin theorem](https://arxiv.org/pdf/1306.2872), Theorem 1.1: independent centered subgaussian coordinates give a quadratic tail controlled by the minimum of the squared-threshold/Frobenius-squared and threshold/operator-norm terms. Fair signs meet the hypotheses, and hollowness makes the mean zero. No absolute-entry operator norm is substituted.

## 2. Every low-cap actual matrix has such a core

Let A be hollow full signs with Q(A)<=C n^(3/2), and put B=A/sqrt(n). The simultaneous diagonal Grothendieck majorant gives D>=B,-B with Tr D<=8Cn. Thus the set J={i:D_ii<=32C} has size at least 3n/4, and ||B_J||op<=32C. Set L=max(1,32C).

Condition on all spins outside J. Any prescribed full external field, after normalization by sqrt(n), merely adds to the arbitrary deterministic external field in the core lemma. If the total physical violation is at most n^(3/2)/1024, the normalized core violation is at most n/1024<=|J|/512. Equation (1) therefore proves

    # {x:sum_i[-x_i(Ax+b)_i]_+ <=n^(3/2)/1024}
       <=2^n exp[-delta(C)n].                         (3)

Here b is ANY prescribed physical external field. Both energy polarities follow by applying the same argument to A and -A separately. Exact one-spin-stable states are included as the zero-violation case.

## 3. A metric consequence for cleanup procedures

Let S be the set in (3). Fix rho in (0,1/2) such that h(rho)<delta(C)/2. A union of Boolean Hamming balls gives

    Pr_(X uniform){d_H(X,S)<=rho n}
      <=exp[-delta(C)n/2+o(n)].                       (4)

Thus almost every uniform starting word is a linear Hamming distance from EVERY low-total-violation word. Any cleanup map producing such a word must change at least rho n spins on all but an exponentially small fraction of uniform starts. This uses no assumptions about a particular dynamics, its update order, or its basin structure.

Equation (3) is an actual-matrix statement, stronger than a polynomial count reduction valid for arbitrary quadratic functions. It does not by itself give an annealed exponential-moment improvement: energy tilting changes the relative weights of stable words. In particular (4) is not asserted for arbitrary Gibbs starting laws.

## 4. Sharp universal threshold: independent audit PASS

The director's stronger threshold is valid. For every fixed C and every fixed 0<=v<1/sqrt(2 pi), there is delta(C,v)>0 such that, uniformly over all actual A with Q(A)<=Cn^(3/2) and all prescribed external fields b,

    Pr{sum_i[-x_i(Ax+b)_i]_+ <=v n^(3/2)}
       <=exp[-delta(C,v)n]                                  (5)

for all sufficiently large n. The field may depend on A, but not on the sampled inside spins.

Indeed choose small fixed epsilon,theta>0 so rho=(1-epsilon)(1-theta) satisfies sqrt(2/pi)rho^(3/2)>2v. The diagonal majorant gives a core of size k>=(1-epsilon)n and operator norm L=max(1,8C/epsilon). Normalize A and its field by sqrt(n). If at least theta k core coordinates have field magnitude greater than T, the mismatch count on the event in (5) is at most

    4L^2 k/T^2 + 2v n/T.

A sufficiently large fixed T makes this at most one quarter of the number of heavy coordinates. The associated Hamming count has a uniform exponential deficit.

Otherwise condition on the fewer than theta k heavy spins, leaving k'>=rho n independent signs. The remaining deterministic field has norm at most (T+L)sqrt(n). For each remaining coordinate the unshifted row sum has EXACT law S_(k'-1)/sqrt(n), up to signs, where S_j is the sum of j independent fair signs. Therefore

    EF >= k' E|S_(k'-1)|/sqrt(n)
       >=[sqrt(2/pi)rho^(3/2)+o(1)]n >2vn+a n

for some fixed a>0. Convex concentration makes F<EF-a n/2 exponentially unlikely. On the complementary event, F-(quadratic+linear)=2I_R and I_R<=vn force the centered quadratic-plus-linear form above a n/2. The same Hanson--Wright and linear Hoeffding bounds finish the proof. All prefactors and CLT errors are uniform because k'>=rho n, and the external field affects only a uniformly bounded coefficient norm after conditioning.

The threshold cannot be increased. For b=0 and EVERY hollow full-sign A, not requiring a cap bound,

    E I_A,0 = n E|S_(n-1)|/2
            =[1/sqrt(2pi)+o(1)]n^(3/2).

Consequently for any fixed v>1/sqrt(2pi), Markov's inequality gives Pr(I_A,0<=v n^(3/2)) bounded below by a positive constant for all sufficiently large n. No assertion at exact equality is made.
