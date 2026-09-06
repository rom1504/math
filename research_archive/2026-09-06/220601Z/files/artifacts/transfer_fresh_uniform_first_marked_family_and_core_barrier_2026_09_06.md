# Uniform compact response families and the quantitative core-loss barrier

2026-09-06. Quantitative follow-up to the proved strict first-marked gain.
Two results are established below:

1. A fixed-complexity compact rectangle family admits a uniform finite
   original-degree probe catalog and a fixed positive gain at every fixed L.
2. The PARTICULAR conservative gain certified by the current proof never
   pays its principal-core deletion loss for baselines and parent cap bounds
   at least 1/4. This is not an upper bound on the actual possible gain.

No terminal functional is numerically optimized in this note.

## 1. Response family, topology, and the useful uniform statement

Let F be a nonempty family of functions `f:R^2->{-1,0,1}` which are odd
under joint reversal. Assume a fixed total of at most J>=1 vertical and
horizontal threshold lines suffices to make each f constant on the grid
cells. Fixed unions of rectangles are included after passing to their
coordinate grid. Labels come from the fixed finite set {-1,0,1}.
Thresholds may range in a compact interval; alternatively use their
Gaussian CDF coordinates in [0,1], allowing infinite thresholds in the
closure. Impose the closed constraints

```math
 v_0\le\mathbb E_\gamma f^2\le1-\mu_0,
 \qquad 0<v_0\le1-\mu_0<1.                         (1)
```

One may impose any further closed near-optimality constraint from a
continuous, independently justified benchmark. No value of that benchmark
is asserted here. In particular the current unrestricted decimal lower
bound is not automatically a benchmark for this particular two-coordinate
first-marked family.

Use a fixed odd convention on cell boundaries, for instance value zero
on the symmetric grid boundary. Boundary conventions have zero Gaussian
mass and vanishing averaged actual mass in the old marginal limit.
Coalescing thresholds and zero-area cells are allowed. The parameter map
to Gaussian L2 is continuous: changing thresholds changes f only inside
the corresponding vertical/horizontal strips. With labels fixed,

```math
 \|f_\theta-f_{\theta'}\|_{L^2(\gamma_2)}^2
 \le4\sum_{\text{all cuts }j}
       |\Phi(t_j)-\Phi(t'_j)|.                     (2)
```

There are only finitely many label patterns. Their compact parameter
images, followed by the closed conditions (1) and oddness, form an
L2-compact ambient family. Replace F by its closure in this ambient
family when necessary; all the assertions for the closure imply those
for F. Rational thresholds alone are not a closed compact
parameter set; the uniform result applies to their compact real closure
and hence to the rational subfamily.

For the ACTUAL first marked fields associated with B, define

```math
 F_f=f(G,Y),\quad H_f=1-F_f^2,\quad
 j_n(B,f)=\frac1n\sum_i\mathbb E H_{f,i}|(BF_f)_i|,
 \Lambda_n(B)=\max_x\frac{|x^TBx|}{2n},
 \quad B=A/\sqrt{n-1},\quad\|B\|_{op}\le L.
```

For every fixed L>=1 there is a number Delta_F(L)>0 such that, uniformly
over all such signing sequences and over arbitrary choices `f_n in F`,

```math
 \liminf_n[\Lambda_n(B_n)-j_n(B_n,f_n)]\ge\Delta_F(L).
                                                               (3)
```

Equivalently, one may put `sup_(f in F) j_n(B_n,f)` in (3). Approximate
maximizers suffice, so attainment of that supremum is unnecessary.
Sections 2--5 prove this with a fully specified conservative constant.

## 2. A uniform amount of original degree greater than three

Write the orthonormal Gaussian Hermite expansion in h_a(g)h_b(y), with
ORIGINAL degree `P=a+3b`. For an odd f, the subspace of original degree
at most three is exactly

```math
 S_3=\operatorname{span}\{g,y,h_3(g)\}.
```

There is a positive eta=eta(v0), independent of the rectangle thresholds,
such that every odd |f|<=1 with `E f^2>=v0` satisfies

```math
 \operatorname{dist}_{L^2}(f,S_3)^2\ge\eta^2>0.       (4)
```

Here is a finite-dimensional quantitative specification. Put
`r0=sqrt(v0/2)` and

```math
 \kappa(v_0)=\min_{p\in S_3,\ \|p\|_2=r_0}
                 \mathbb E(|p|-1)_+^2,
 \qquad\eta^2=\min\{v_0/2,\kappa(v_0)\}.            (5)
```

The sphere is compact, the displayed integral is continuous there, and
every nonzero member is an unbounded polynomial on a Gaussian full-support
space. Thus kappa>0. If the projection p of f to S3 has squared norm at
most v0/2, orthogonality proves (4). Otherwise its norm is at least r0,
and radial monotonicity of `E(|p|-1)_+^2`, together with |f|<=1, proves
the same bound. This does not require compactness of the whole response
class; boundedness and the variance floor suffice for (4).

If an entirely explicit lower constant is preferred to the three-dimensional
minimum in (5), set

```math
 T_0=28/r_0,\quad h_0=1/[8(T_0+1)^2],\quad
 \eta^2=\min\{v_0/2,\ h_0^2 e^{-(T_0+1)^2}/(2\pi)\}. (6)
```

To verify (6), write `p=a g+b y+c h3(g)=A g+C g^3+b y`, with coefficient
norm r0. If |p|<=2 at (0,T0),(T0,0),(T0/2,0), then
`|b|<=2/T0`, `|C|<=8/T0^3`, `|A|<=6/T0`, hence
`|a|<=12/T0`, `|c|<=2sqrt(6)/T0` and
`r0<=sqrt(172)/T0<r0`, a contradiction. At one point |p|>2.
On its square of coordinate radius h0 the derivative bound
`|partial_g p|+|partial_y p|<=4(T0+1)^2` keeps |p|>=3/2.
That square has Gaussian density at least
`e^{-(T0+1)^2}/(2pi)` and area 4h0^2. Its contribution to the positive
part integral is at least the second quantity in (6).

## 3. Fixed rectangle complexity supplies a finite Hermite catalog

Let X,Y be standard two-dimensional Gaussians with coordinatewise
correlation rho in [0,1). A disagreement of grid cells crosses at least
one of the J threshold lines. Each threshold-crossing probability is at
most `arccos(rho)/pi`, maximized by threshold zero. Thus

```math
 \mathbb E[f(X)-f(Y)]^2\le4J\arccos(\rho)/\pi.        (7)
```

Writing c_ab for the Hermite coefficients, Mehler's identity gives
`E[f(X)-f(Y)]^2=2 sum_(a,b)(1-rho^(a+b))c_ab^2`.
For integer d>=2 choose rho=1-1/d. Since
`arccos(1-1/d)<=pi/sqrt(2d)` and
`1-(1-1/d)^(d+1)>1/2`,

```math
 \sum_{a+b>d}c_{ab}^2\le3J/\sqrt d,
 \qquad
 \sum_{a+3b>3d}c_{ab}^2\le3J/\sqrt d.                (8)
```

Fix once and for all

```math
 d=\left\lceil\max\{2,36J^2/\eta^4\}\right\rceil,
 D=3d,\qquad\tau_0^2=\eta^2/(6d).                   (9)
```

Equations (4),(8) imply that every f in the family has SOME odd
`P in {5,7,...,D}` with

```math
 \tau_P^2=\sum_{a+3b=P}c_{ab}^2\ge\tau_0^2.          (10)
```

The choice P may depend on f, but the catalog and lower bound do not.
Subtracting the first projection b0*g+b1*y leaves exactly these P>3
coefficients unchanged. Thus these are precisely the residual blocks
required by the actual high-degree feedback theorem. No increasing
original degree is exchanged with the signing limit.

## 4. Fully specified uniform downstream constants

The source proof `transfer_seed_high_degree_actual_feedback_probe_2026_09_06.md`
uses holes mu, a selected tau_P, and fixed L. Since mu>=mu0 and
tau0<=tau_P<=1, all its scalar choices can be replaced by the following
uniform, conservative ones. Write phi for the standard Gaussian density.

```math
 c_*={\mu_0\tau_0^2\over4L}
              \phi\!\left({4L\over\tau_0\sqrt{\mu_0}}\right),
 \delta_g=\mathbb E(c_*|N|-L)_+^2,
 K=\left\lceil\max\{2,36L^4/\delta_g^2\}\right\rceil,
 a_0^2=\delta_g/(2K),
 \epsilon=\min\{1/4,(a_0/2)^{2/3},a_0^4/(64L^8)\},
 a_*={a_0^2\over8L^2},\qquad
 M_0=L^{2K+2}\epsilon^{-K}.                          (11)
```

All are finite, positive and independent of f and n. The finite probe
catalog consists of odd P<=D and odd k<=K. The source first-coefficient
bound is at least c_*, its response amplitude at most L, and its two
small-variance removal inequalities follow from the choice of epsilon.
The actual FULL return correlation is therefore at least a_* n+o(n),
and the actual probe second moment is at most L^2 n+o(n). Its ideal row
variance is at most M0. The normalized covariance/variance errors and
the joint literal-query comparison are taken at FIXED P,k,epsilon.

Insert (11) into the proved endpoint implication as follows:

```math
 a=a_*/2,\quad M=2\max(1,M_0),\quad
 \delta=(a/(8L))^2,\quad R=128ML^4/a^2,
 \rho=a/(4\sqrt{MR}),\quad s=a/(4\sqrt M),\quad
 h=\mathbb E(sN-2\sqrt R)_+,\quad
 d_0=\min\{L,\rho h/4\},\qquad
 \Delta_F(L)=d_0^2/(2L)>0.                           (12)
```

This is the uniform constant promised in (3). The return and old-field
raw second moments have average at most L^2 each. The endpoint proof
retains a fraction at least rho of roots where probe variance is in
[delta,M], return-plus-old-field second moment is at most R, and the
return covariance is at least a/4. Regression coefficients are bounded
by sqrt(R/delta). On those SAME roots the literal W has a fixed
L-dependent moment cap and `E Z_i^2<=2R+2[L^2+(729/2)L^4]`.
Thus the uniformity statement includes the actual root caps, not merely
a Gaussian row-variance calculation.

## 5. Why this is uniform for ACTUAL varying responses

Gaussian L2 compactness alone would not justify a uniform actual-response
theorem: functions can be changed on an n-dependent finite set without
changing their Gaussian equivalence class, while that set can contain
the whole finite actual marginal support. Fixed rectangle complexity
prevents this issue.

Specifically, along any sequence of such actual matrices with fixed L,
the old marginal `(G_i,Y_i)` converges uniformly in the averaged sense
to its continuous Gaussian marginal. Its one-dimensional CDF convergence
is uniform in threshold: if
`F_n^G(t)=n^-1 sum_i P(G_i<=t)` and analogously for Y, then
`zeta_n=max(sup_t|F_n^G(t)-Phi(t)|,sup_t|F_n^Y(t)-Phi(t)|)->0`.
This follows from weak convergence to a continuous CDF by a fixed finite
CDF grid, and does not assert a new quantitative rate. The maximal
averaged atom mass is at most `2 zeta_n`, so moving boundary conventions
also have vanishing cost. Formula (2), with Gaussian strip probabilities
replaced by actual ones plus o(1), therefore implies

```math
 \theta_n\to\theta\quad\Longrightarrow\quad
 {1\over n}\mathbb E\|F_{f_{\theta_n}}-F_{f_\theta}\|^2\to0.
                                                               (13)
```

For ternary vectors, `|H_f-H_g|<=|F_f-F_g|`. Without making ANY Gaussian
replacement of BF, Cauchy--Schwarz and bounded-op transport give

```math
 |j_n(B,f)-j_n(B,g)|
 \le2L\left({1\over n}\mathbb E\|F_f-F_g\|^2\right)^{1/2}.
                                                               (14)
```

Now suppose (3) failed. Choose a violating signing/response sequence.
Pass to a subsequence with one fixed cell-label pattern and convergent
threshold parameters. Its limit f belongs to the compact family. The
FIXED-f actual theorem applies with the uniform constants (9)--(12).
Equations (13)--(14) transfer its gain to the varying sequence, a
contradiction. This avoids asserting a simultaneous increasing-family
Boolean diagram closure. The hard sign in C disappears from j_n, so
the uniform transfer itself needs no new hard-threshold comparison.

The order is: choose the family parameters J,v0,mu0 and fixed L; choose
eta,d,D,tau0; choose the finite probe catalog and all constants in
(11)--(12); for the fixed limiting response apply the already audited
ordered approximations with n first; then use (13)--(14). No threshold
complexity, primitive degree, variance floor, or operator cap grows with n.

## 6. Exact principal-core propagation and its finite slack

Suppose an actual N-vertex hollow signing has
`Q_abs(A)<=C N^(3/2)`. The banked simultaneous diagonal-majorant deletion
uses `beta(A)<=4Q_abs(A)` and yields a diagonal `D_A>=+/-A` with
`tr D_A<=4K_G Q_abs(A)`. Deleting diagonal entries above
`4K_G C sqrt(N)/epsilon` leaves a principal core R with

```math
 m=|R|\ge(1-\epsilon)N,
 \|A_R\|_{op}\le4K_G C\sqrt N/\epsilon,
 \left\|{A_R\over\sqrt{m-1}}\right\|_{op}
 \le {4K_G C\over\epsilon}\sqrt{N\over m-1}.         (15)
```

This keeps the child exactly an original hollow signing. Principal cap
monotonicity follows by averaging the omitted independent spins.
For p=1-epsilon, its asymptotic normalized operator cap is at most

```math
 L_\epsilon={4K_G C\over\epsilon\sqrt p}.            (16)
```

At finite order (15), not (16), is the exact formula. To apply a theorem
requiring a fixed literal cap, use any L>L_epsilon, with epsilon fixed
before N tends to infinity. Alternatively, if C is STRICTLY above the
limsup parent normalized cap, the slack absorbs the `(m-1)` correction
and L=L_epsilon itself is valid eventually. The latter is available for
minimum-realizing sequences with C=1/2 because the audited all-order
upper constant is strictly below 1/2. No continuity through the integer
cutoff K in (11) is silently assumed.

If a uniform family benchmark supplies `liminf j_m>=j_F`, the actual
endpoint theorem and exact principal monotonicity give

```math
 \liminf_N{Q_{\rm abs}(A_N)\over N^{3/2}}
 \ge p^{3/2}\,[j_F+\Delta_F(L)]                     (17)
```

with L as just specified. The exact prefactor before the limit is
`m sqrt(m-1)/N^(3/2)`, not simply m/N. If more than pN vertices survive,
this only increases its limiting lower bound. A negative benchmark
would require additional care with that monotonicity; only positive
benchmarks are considered here.

## 7. An explicit barrier for the CURRENT conservative gain

The certified number in (12), not the true gain, admits a simple upper
bound. First the source constants obey

```math
 c_*\le{1\over4L\sqrt{2\pi}},\quad
 \delta_g\le c_*^2,\quad K\ge2,
 \quad a={a_0^2\over16L^2}
 \le{1\over2048\pi L^4}<1.                          (18)
```

Since `h<=s/sqrt(2pi)`, the constants in (12) give

```math
 d_0\le{a^2\over64M\sqrt R\sqrt{2\pi}},\qquad
 \Delta_F(L)\le{a^6\over2^{21}\pi M^3L^5}
 \le{1\over2^{24}\pi L^5}.                         (19)
```

Here M>=2. Keeping (18) would give the even smaller upper bound
`2^(-90) pi^(-7) L^(-29)`, but that is unnecessary for the comparison.

Take any meaningful benchmark `j_F>=1/4` and parent cap bound `C>=1/4`.
The real Grothendieck constant satisfies K_G>=1, so `D=4K_G C>=1`.
For every epsilon in (0,1) and any allowable `L>=D/(epsilon sqrt(p))`,
equation (19) yields

```math
 p^{3/2}\Delta_F(L)
 \le{\epsilon^5p^4\over2^{24}\pi D^5}
 <j_F\epsilon
 \le j_F(1-p^{3/2}).                               (20)
```

Consequently the lower-certificate expression supplied by this route
always satisfies

```math
 p^{3/2}[j_F+\Delta_F(L)]<j_F.                       (21)
```

This holds for EVERY deletion fraction; there is no numerical epsilon
optimization to perform. Taking epsilon to zero only recovers j_F from
below. Compact-family uniformity fixes a real quantifier issue, but the
current conservative gain still cannot improve its own positive baseline
after this spectral-core passage.

Scope is essential: (19)--(21) bound the explicit SMALL POSITIVE NUMBER
which the current proof guarantees. They do NOT upper-bound the actual
spin-ascent gain, show that stronger estimates cannot pay the deletion
loss, or refute a different method of removing the operator cap.

## 8. Sources and proof-obligation separation

The exact first-marked mixed comparison and literal regression premise
are audited in `transfer_adversary_high_degree_actual_feedback_audit_2026_09_06.md`.
The endpoint constants, both energy orientations, finite hard-tie means,
and row cutoffs are proved in
`transfer_fresh_conditional_actual_feedback_endpoint_gain_2026_09_06.md`.
The diagonal-majorant core and its normalization were rechecked against
`continued_audit_first_return_and_rounding_2026_09_06.md`, Section 5, and
`transfer_director_exponential_selector_cost_2026_09_06.md`, Section 2.
The threshold noise-sensitivity calculation is the elementary one also
used in `transfer_director_inhomogeneous_gaussian_return_rigidity_2026_09_06.md`,
Section 3; only its Gaussian threshold-crossing identity is needed here.

New conclusions of this note: a uniform fixed-complexity response-family
gain with a finite high-original-degree catalog, and an analytic barrier
to using its present certified constant through the stated core deletion.
No unrestricted decimal improvement or convergence theorem is claimed.

Director full-read check reported PASS on 2026-09-06; the explicit
CDF clarification above changes no constant or conclusion.
Independent adversarial audit also reported PASS for the moving-threshold
actual-L2 transfer, explicit eta separation, finite Hermite catalog, and
the all-deletion-fractions certified-gain barrier.
