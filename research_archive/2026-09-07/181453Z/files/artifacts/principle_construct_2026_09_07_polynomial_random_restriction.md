# Polynomial-scale random restrictions of every bounded-cap signing

2026-09-07. Status: proved below; root independently read and reconstructed the complete proof, audit PASS. The exact rational checker also passes. This is a quantitative obstruction to RANDOM downsampling, including from selectable exact minimizers. It is not an obstruction to rare selectable restrictions, arbitrary parent construction, or convergence of the naked sequence M_n.

**Archive attribution update.** Subsequent tensor-route collision checking found `transfer_seed_intermediate_random_tensor_ranges_2026_09_06.md`, which already uses this degree-28 greedy-moment architecture for tensor powers. Its last-column bilinear estimate in fact proves the correlations needed here directly, without Schatten-four. Section 8 gives that simpler proof, including the weighted extension. The universal actual-input and weighted statements remain valid, but the moment/minorant architecture is not new; the initial spectral derivation below is retained as a valid alternate proof, not claimed as necessary progress.

## 1. Statement

Let A_n be any deterministic hollow symmetric full signing satisfying Q(A_n)<=C n^(3/2), with C fixed. Let I_n be a uniformly random m_n-subset of vertices.

**Polynomial-size theorem.** If m_n tends to infinity and m_n=o(n^(1/28)), then

    Pr{ Q(A_n[I_n])/m_n^(3/2) >= .5027 } -> 1.       (1)

In particular m_n=floor(n^(1/30)) is allowed. A slightly stronger constant supplied by the exact certificate is .50276269428122185772... minus an arbitrary fixed positive error.

**Subpolynomial-size theorem.** If m_n tends to infinity and log(m_n)/log(n) tends to zero, then, for every epsilon>0,

    Pr{ Q(A_n[I_n])/m_n^(3/2)
           >= (2/3)sqrt(2/pi)-epsilon } -> 1.         (2)

The coefficient in (2) is .5319230405352436.... Both statements are uniform over the choice of A_n subject to the fixed cap bound. They therefore apply to ACTUAL exact minimizers at order n. Since the proved all-order upper bound is below .493609, (1) already gives a fixed normalized loss on those inputs.

This improves the scope of the earlier full-total-variation argument in `transfer_adversary_random_restriction_2026_09_06.md`, whose allowed child order was only of square-root-logarithmic size, by extending the archived tensor greedy-moment method to all bounded-cap inputs. Uniform moment control survives adaptive greedy signs. The initial fourth-moment derivation below is valid but unnecessary, as Section 8 explains.

## 2. Spectral input and random column correlations

Set F=A+I, making every entry a full sign. For the greedy principal-submatrix process the diagonal never enters the actual energy. The verified Schatten-four inequality gives

    Tr(F^4) <= K_G^2 beta(F)^2
             <= K_G^2(4Q(A)+n)^2 <= D_C n^3.         (3)

The proof and primary Grothendieck normalization are in `flatify_construct_2026_09_07_spectral_fourth_moment.md`. Only the bound (3) is needed below; alternatively it can be made a hypothesis of the theorem.

For an ordered distinct k-tuple J of columns define

    c_J = (1/n) sum_i product_(j in J) F_ij.

For every fixed k>=2, uniformly in n,

    E_J c_J^2 <= [D_C+k(k-1)]/n.                    (4)

Indeed, for independent column sampling with replacement, the expectation of the row-pair product is rho_uv^k, where rho_uv=(F F^T)_uv/n and |rho_uv|<=1. Thus its absolute average over u,v is at most

    (1/n^2)sum_uv rho_uv^2 = Tr(F^4)/n^4 <= D_C/n.

Replacing sampling with replacement by distinct sampling changes an expectation of a [-1,1]-valued product by at most twice the collision probability, at most k(k-1)/n. Expanding c_J^2 proves (4). Consequently E|c_J|<=sqrt((D_C+k(k-1))/n).

The case k=1 is intentionally NOT asserted. Only even k>=2 will occur below, so the weaker available row-sum estimate causes no gap.

## 3. Uniform even moments for every signing of a random column set

Let J=(J_1,...,J_j) be a uniformly sampled ordered distinct column set. For v in {+1,-1}^j put

    X_(J,v)(i)=j^(-1/2)sum_(a=1)^j F_(i,J_a) v_a,

where i is uniform over all n rows. Let Z_j=j^(-1/2)sum_(a=1)^j epsilon_a for independent fair signs. For each fixed integer r>=1 define

    Delta_(2r)(J)
      =sup_v | E_i X_(J,v)(i)^(2r) - E Z_j^(2r) |.

Then

    E_J Delta_(2r)(J) <= C_(C,r) j^r/sqrt(n).        (5)

To prove this, expand the 2r-th moment into ordered index tuples. A tuple in which all multiplicities are even contributes exactly the same value to both moments, since every F entry and v entry is a sign. Every other tuple has a nonempty parity-support of EVEN cardinality k, where 2<=k<=2r. Its row average is c_T times a product of selected v coordinates. Taking absolute values BEFORE maximizing v removes all dependence on v. There are at most j^(2r) tuples, the normalization is j^r, and (4) bounds the expected absolute value of each surviving correlation by C_(C,r)/sqrt(n). This proves (5).

This simultaneous bound is crucial: the signs v may be arbitrary functions of all entries previously exposed among the selected vertices. No independence between v and those entries is invoked.

## 4. A global polynomial minorant and its exact certificate

There exists an even polynomial P of degree 28 with P(x)<=|x| for all real x and

    mu_P := E P(G) = .75414404142183278658...,
    (2/3)mu_P > 5027/10000,                          (6)

where G is standard Gaussian. Here is an exact rational specification, so (6) does not rely on a floating-point optimization.

Take seven positive rational numbers

    x_i in {799,1607,2432,3289,4196,5190,6364}/1000,
    y_i=x_i^2.

Let h be the degree-at-most-13 Hermite interpolant to f(y)=y^(-1/2) at these seven nodes, matching f(y_i)=1/x_i and f'(y_i)=-1/(2x_i^3). Set P(x)=x^2 h(x^2).

For every y>0 distinct from the nodes the Hermite remainder is

    f(y)-h(y)=f^(14)(xi)/14! product_i(y-y_i)^2 >=0

for some xi in the positive interval spanned by y and the nodes. The derivative is positive. At the nodes equality holds. Multiplication by y and continuity at zero prove P(x)<=|x| globally, including beyond the largest interpolation node.

The checker `computations/principle_construct_2026_09_07_polynomial_restriction_minorant.py` constructs h and P over exact rationals, checks all interpolation equations, and computes E P(G) using E G^(2r)=(2r-1)!!. It verifies both rational inequalities in (6). The displayed decimals are nonessential.

More generally, for every epsilon>0 there is an even finite-degree polynomial P_epsilon<=|x| with E P_epsilon(G)>=sqrt(2/pi)-epsilon. One elementary construction uses

    |x|=(2/pi) integral_0^infinity (1-cos(tx))/t^2 dt.

For each integer k, the Taylor polynomial of cos through degree 4k is a GLOBAL upper bound on cos. This follows by twice integrating the preceding global lower bound, starting with cos(z)<=1 and cos(z)>=1-z^2/2. Replace cos by that upper polynomial and integrate only over [0,T]. The resulting even polynomial is a global lower bound for |x|. For fixed T its Gaussian expectation tends to the exact truncated integral, because the Gaussian cosine series is exp(-t^2/2) and converges uniformly on the compact t interval after division by t^2. Letting T grow recovers sqrt(2/pi). Thus arbitrary precision needs only some fixed degree, not a degree growing with n.

## 5. Adaptive greedy construction on the random principal restriction

Reveal vertices J_1,...,J_m in uniformly random order without replacement. Assign x_(J_1)=1. Having assigned the first j spins, choose x_(J_(j+1)) to make

    x_(J_(j+1)) sum_(a<=j) A_(J_(j+1),J_a) x_(J_a)

nonnegative. The final induced energy is H=sum_(j=1)^(m-1) Z_(j+1), where each Z_(j+1) is the absolute value of that fresh field.

Write P(x)=sum_(r=0)^R p_r x^(2r) and define the nonnegative error

    D(J)=sum_(r=1)^R |p_r| Delta_(2r)(J).

The constant coefficient, if present, contributes no moment error. By (5), E D(J)<=C_(C,P) j^R/sqrt(n). By standard finite moment expansion of a normalized Rademacher sum,

    E P(Z_j)=mu_P+O_P(1/j).

Consequently the average absolute field over ALL n vertices, uniformly over the adaptively selected signs, is at least

    sqrt(j)[mu_P-C_P/j-D(J)].                        (7)

Excluding the already selected j vertices changes this lower bound by at most j^2/(n-j): each excluded field has absolute value at most j. Thus (7), less that deterministic error, is a lower bound on the conditional expectation of Z_(j+1) given the entire exposed history.

Summing the conditional expectations yields a compensator at least

    (2/3)mu_P m^(3/2)
      -O_P(sqrt(m))
      -sum_(j<m) sqrt(j) D(J_1,...,J_j)
      -O(m^3/(n-m)).                                (8)

No simultaneous-good-prefix event is required. Each prefix remains a uniformly random ordered subset, although its assigned spins are adaptive. Taking expectations of the NONNEGATIVE random error in (8) gives at most C_(C,P) m^(R+3/2)/sqrt(n). Markov therefore makes that error o(m^(3/2)) in probability whenever m^R/sqrt(n)->0. The other errors are also negligible in the ranges claimed.

Finally, the conditional second moment of Z_(j+1) is bounded above by

    [n/(n-j)] j[1+Delta_2(J_1,...,J_j)].

Using (5) at r=1 gives a total expected conditional variance O_C(m^2) whenever m=o(sqrt(n)). The martingale difference between H and its compensator therefore has variance O_C(m^2), hence is o(m^(3/2)) in probability. This proves

    Q(A[I])/m^(3/2) >= (2/3)mu_P-o_P(1)

whenever m tends to infinity and m^R/sqrt(n)->0. Taking R=14 and the exact polynomial proves (1). For m=n^(o(1)), every fixed R is admissible; then use the arbitrary-precision minorants and obtain (2).

## 6. Scope for the original convergence problem

This is an actual-input theorem, not a proxy calculation: no optimality hypothesis is required, so no choice among exact minimizers evades RANDOM restriction at the stated scales. It strengthens the range over which a natural all-order recovery operation is known to lose a leading normalized amount.

It does not assert that every m-subset is bad. The proof neither rules out exceptional selectable near-optimal induced subgraphs nor constructs a full parent from selected children. It supplies no separated subsequences for M_n and no failure of convergence. Coding/free-energy reinterpretations alone were not used to infer an all-order realization theorem.

## 7. Bounded-amplitude row-square-regular extension

The same conclusions (1) and (2) hold for hollow symmetric W_n satisfying

    |W_ij|<=K,   sum_j W_ij^2=n-1 for EVERY i,
    Q(W_n)<=C n^(3/2),

where C,K are fixed. The output in this assertion is the weighted principal restriction itself; it is not asserted to be a full signing or exactly row-regular after restriction.

Set F=W+I and replace K by max(K,1). Every row and column of F has square sum n. The bounded-entry Schatten-four theorem again gives Tr(F^4)<=D_(C,K)n^3. For a random j-column set and arbitrary signs v, let X_(J,v) be as in Section 3. Then, directly relative to Gaussian moments,

    E_J sup_v |E_i X_(J,v)^(2r)-(2r-1)!!|
       <=C_(C,K,r)[j^(-1)+n^(-1/4)+j^r/sqrt(n)].    (9)

Here are the details needed because powers of weighted entries can no longer simply be reduced modulo two. Classify each ordered 2r-tuple in the moment expansion by its d distinct columns, their positive multiplicities a_1,...,a_d, and the number s of multiplicities equal to ONE.

* If s>=2, square the corresponding row-averaged product and average its distinct column labels. Under independent column sampling, each singleton contributes rho_uv=(F F^T)_uv/n. All other factors are bounded by constants depending only on K,r. Since |rho_uv|<=1 and its mean square over row pairs is at most D_(C,K)/n, the expected squared row-average is O_(C,K,r)(1/n). The distinct-sampling collision correction has the same order. There are at most O_r(j^(2r)) such tuples; after normalization this contributes at most C j^r/sqrt(n), uniformly over v.

* If s=1, parity of the total 2r forces at least one other multiplicity to be odd and at least three. Therefore d<=r. The same squared-product calculation now has only one rho factor. Its mean absolute value over row pairs is O(n^(-1/2)) by Cauchy-Schwarz and the fourth moment bound. The row-average consequently has expected absolute value O(n^(-1/4)). Since d<=r, the total normalized contribution remains O(n^(-1/4)). No first-power row mean is assumed small uniformly.

* If s=0 and the multiplicities are not all two, then d<=r-1. The entries are bounded, so the total normalized contribution is O(1/j) without a spectral estimate.

* For the pure paired tuples, put V_i=j^(-1)sum_(a<=j)F_(i,J_a)^2. Their contribution is

      (2r-1)!! (1/n)sum_i V_i^r +O_(K,r)(1/j).

  Each V_i is in [0,K^2], E_J V_i=1, and E_J(V_i-1)^2=O_K(1/j). Moreover (1/n)sum_i V_i=1 EXACTLY for every selected column set, since all columns have square sum n. Taylor expansion of z^r on [0,K^2] therefore cancels the averaged linear term and gives

      E_J |(1/n)sum_i V_i^r-1|=O_(K,r)(1/j).

  This contribution is independent of the signs v.

These four cases prove (9). They also explain exactly why the unweighted parity argument cannot be invoked unchanged.

Apply the same fixed polynomial minorant to (9). Its additional compensator error is O(m^(3/2)n^(-1/4)), which is negligible. Excluding selected rows costs at most K j^2/(n-j), still negligible. For the martingale variance the r=1 calculation is sharper: the pure paired term is exactly one, and only the two-singleton terms remain, with expected error O(j/sqrt(n)). Thus the total variance is again O_(C,K)(m^2). Sections 4 and 5 now prove the weighted extension with exactly the same allowed orders and constants.

## 8. Simpler direct bilinear proof and precise archive collision

For a full signing F and an ordered distinct column tuple J of length q, condition on its first q-1 columns. Their rowwise product is a Boolean vector f. Average the absolute remaining-column correlation over the n-q+1 available columns. Then

    E(|c_J| | first q-1 columns)
       <= ||F^T f||_1/[n(n-q+1)]
       <= beta(F)/[n(n-q+1)].                        (10)

Thus beta(F)=O(n^(3/2)) yields the same O(n^(-1/2)) correlation estimate for EVERY fixed q>=1. This is stronger than needed in Section 3 and avoids both a collision comparison and spectral estimates. It is the direct analogue of Section 1.1 of the archived intermediate-tensor artifact, where product factorization supplied an additional tensor-power decay. Here the global bounded cap already supplies the decay at the ambient order.

For the weighted extension, consider any multiplicity pattern having a singleton column. Condition on all the other columns and choose that singleton as the last sampled column. The rowwise product of the other powers has infinity norm at most K^(2r-1). The bilinear norm inequality, valid on the whole cube, gives the same bound (10) multiplied by K^(2r-1). Summing all these terms proves directly

    E_J sup_v |E_i X_(J,v)^(2r)-(2r-1)!!|
       <= C_(C,K,r)[j^(-1)+j^r/sqrt(n)].             (11)

The patterns with no singleton are handled exactly as in Section 7: nonpaired ones have at most r-1 labels, and pure pairs use the exact row/column square sums. Equation (11) improves (9) by removing its unnecessary n^(-1/4) term and proves the weighted theorem using only bounded entries, row-square regularity, and beta(F)=O(n^(3/2)).

The exact rational minorant checker in this note is an independent certificate. The mechanism of a degree-28 global minorant and adaptive greedy sums was already present in the cited tensor artifact. No literature or archive novelty is claimed for that mechanism.
