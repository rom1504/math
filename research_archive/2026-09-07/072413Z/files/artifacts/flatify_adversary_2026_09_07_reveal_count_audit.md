# Independent audit: posterior-count cleanup for absolute pressure

Date: 2026-09-07. Verdict: the theorem in
`flatify_independent_2026_09_07_reveal_count_concentration.md` is correct,
including uniformity for all beta>0 and the N^(2/3) pressure exponent.
This audit reconstructs the estimates; it does not infer an endpoint drift.

## 1. The potentially dangerous comparison is valid

Use F(v)=min_A log sum_x cosh(beta H_A(v)/sqrt(N-1)), with the ambient
denominator fixed even after deletion. Under the joint orientation/spin
law let Y_e=sigma A_e x_i x_j and lambda_e=beta sqrt(v_e)/sqrt(N-1).
Flipping edge e multiplies the partition sum by
cosh(2lambda_e)-E(Y_e)sinh(2lambda_e). Global sign optimality says this
ratio is at least one, hence E(Y_e)<=tanh(lambda_e). Therefore every
active smooth branch has partial derivative at most beta^2/[2(N-1)].
The finite lower envelope inherits the integrated one-sided bound along
coordinatewise increasing paths. Zero coordinates are handled by positive
approximation and continuity; no lower derivative bound is required.

Separately, F(cv)<=F(v) for 0<=c<=1 because each fixed-sign sum of coshes
is radially nondecreasing. If w>=(1-epsilon)v and the two edge masses are
T, integration from (1-epsilon)v to w costs at most
beta^2 epsilon T/[2(N-1)]. This is precisely the needed upper comparison.
Reverse scaled domination supplies its reverse. Coordinatewise monotonicity
is neither true in general nor used here.

## 2. Deletion and profile arithmetic

Summing one inserted spin gives a factor 2cosh(row field), proving an
increment of at least log2. Fixing an optimal remainder and averaging the
new independent edge signs gives exactly the old partition sum times
2 product cosh(lambda_e). Thus the increment minus log2 is in
[0,beta^2/2] whenever the row variance is at most N-1. This argument works
for the absolute partition sum, including its shared orientation variable.

Consecutive posterior counts differ at one known vertex. After deleting
that vertex, the known labels agree, zero patterns agree, and both edge
masses are (N-2)(N-1)/2. The remaining unknown population is u in both
profiles. Known/unknown weights change by relative at most 2/(delta u).
For unknown/unknown weights, direct subtraction gives

    c(R-1)-c(R)=2[b(u-R)-a(R-1)]/[u(u-1)].

Since a,b<=2/delta and a,b>=1, its absolute value is at most
4/(delta u), while c(R)>=1/3 for u>=4. Thus bidirectional relative
domination holds with epsilon=12/(delta u). The mass comparison costs
at most 3 beta^2 N/(delta u). Restoring the changed vertex costs at most
beta^2/2, not beta^2: both insertion increments lie in the same interval.

## 3. Uniform expectation and conversion

For the hypergeometric revealed count K, Var(K)=p(1-p)ku/(N-1)<=u/2.
Nearest feasible rounding r0 of pk differs by at most 1/2. Hence
E|K-r0|<=2sqrt(u). On |K-r0|<=delta u/4 the consecutive-count bound
telescopes entirely within the central band. The complement has probability
at most 32/(delta^2 u) by Chebyshev. All row-regular pressures lie in
[Nlog2,Nlog2+beta^2 N/4], bounding its contribution. These give
O_delta(beta^2 N/sqrt(u)), uniformly in reveal time.

In the final u-vertex window, compare each posterior profile to any
consistent completed partition, then delete the u unknown vertices from
both. Their common known submatrix gives distance at most u beta^2/2
from the same optimized terminal value. Two count profiles consequently
differ by at most u beta^2. Splitting at u=N^(2/3) proves the claimed
O_delta(beta^2 N^(2/3)) bound. All constants are beta-independent.

The elementary softmax inequalities then give normalized cap error
O_delta(beta N^(-1/3)+1/beta); beta=N^(1/6) is admissible and gives
O_delta(N^(-1/6)). There is no fixed-temperature limit hidden in this step.

## 4. Endpoint scope, including a stronger target hidden by optimization

The terminal inequality T_(N,m)<=Psi_m+Psi_(N-m) is correct: choose the
better of the two child relative global polarities, then average their
diagonal and crossed orientation products. This is an absolute-pressure
inequality and needs no width-to-cap theorem.

Nevertheless, minimizing at the terminal profile makes the proposed
endpoint comparison stronger than a selected-child cap recurrence. For
equal children N=2m there is an exact identity

    T_(2m,m)(beta)=2 W_m(beta),
    W_m(beta)=min_A [log Z_A^+ + log Z_A^-]/2.

Indeed write each child's paired pressure as W_A and its imbalance as
I_A=(log Z_A^+-log Z_A^-)/2. Parent absolute pressure is
W_A+W_B+log cosh(I_A+I_B), at least 2W_m. Taking identical paired-pressure
minimizers with opposite global polarities makes the last term zero.
Thus an absolute-envelope endpoint submartingale estimate directly
compares parent cap-pressure to optimized child WIDTH-pressure. This is
not a defect in the cleanup theorem or the sufficient inequality: it is
an extra structural demand on any proof of that sufficient inequality.

In particular the concentration bound must not be summed over reveal
times. Initial and terminal profiles have no count randomness at all;
the bound alone says nothing about their signed difference.

## 5. Stronger vertex-scaling cleanup: independent verification

After the preceding audit, the independent agent found a stronger local
scaling argument. It is valid and improves the pressure error to
O_delta(beta^2 sqrt(N)), and normalized cap error to O_delta(N^(-1/4)).

Fix a signing and scale all coefficients incident to vertex i by t>=0.
For fixed remaining spins write the energy as H_rest+t x_i h. Summing over
x_i in the absolute partition sum gives

    cosh(H_rest+t h)+cosh(H_rest-t h)
      =2 cosh(H_rest) cosh(t h),

which is nondecreasing in t. Optimizing signs preserves this inequality.
Consequently for a diagonal matrix D with entries in [0,1],
F(DvD)<=F(v): apply successive coefficient-star factors sqrt(D_ii).
This does NOT assert monotonicity of individual edges.

In the consecutive-count remainders above, only edges incident to the u
unknown vertices change. Put d_i=1-epsilon on those vertices and d_i=1
on known vertices, with epsilon=12/(delta u)<=1/2. The bidirectional
relative comparison already proved implies DvD<=w: known-known entries
agree, known-unknown entries lose epsilon, and unknown-unknown entries
lose at least epsilon. Since the two original edge masses agree,

    sum_e [w_e-(DvD)_e] = sum_e [v_e-(DvD)_e]
      <=epsilon sum_(i unknown) sum_(j != i) v_ij
      <=epsilon u(N-1).

The double counting of unknown-unknown edges is intentional and valid:
1-(1-epsilon)^2<=2epsilon. Integrating the one-sided edge derivative
from DvD to w, and using F(DvD)<=F(v), now costs only
6 beta^2/delta. Reverse domination and reinsertion give

    |F_(k,r+1)-F_(k,r)| <= (6/delta+1/2) beta^2.

The same hypergeometric good event therefore costs O_delta(beta^2 sqrt(u)).
The previous Chebyshev bad-event estimate costs O_delta(beta^2 N/u).
For u>=sqrt(N), both are O_delta(beta^2 sqrt(N)); for u<sqrt(N), the
common-terminal deletion bound costs beta^2 u. All estimates remain
uniform over beta>0. Softmax conversion yields

    sup_k E|Q_(k,K)-Q_(k,r0)|/N^(3/2)
       <= C_delta beta/sqrt(N)+2log(2)/beta.

Take beta=N^(1/4). The new error is O_delta(N^(-1/4)). An endpoint drift
bound with the same pressure error would imply an E_N recurrence defect
O_delta(N^(3/4)); this is conditional, not proved by count cleanup.
