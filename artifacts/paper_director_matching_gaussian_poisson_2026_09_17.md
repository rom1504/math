# Signed matchings: a quantitative Gaussian--Poisson comparison

2026-09-17. Director proof, combining the localization track's matching
construction with the scalar Stein mechanism and a reconstructed Palm
coupling. **Independently reconstructed by all three researchers.** This is an absolute-response
comparison, not a Gaussian central limit theorem. It allows persistent
large jumps. External novelty has not been established.

## 1. The theorem

Let p>=6 be even. Attach real coefficients c_e to edges of K_p, with

    max_v sum_(e incident v) c_e^2 <= D.

Let M be a uniform perfect matching, let epsilon_e be independent fair
signs, independent of M, and define

    T=sum_(e in M) epsilon_e c_e,
    v=E T^2=sum_e c_e^2/(p-1),
    kappa=sqrt(2/pi).

For every fixed D there is a finite C_D such that

    E|T| <= kappa sqrt(v)+C_D p^(-1/5).              (1)

More generally, after any threshold 0<tau<=1, T is within Wasserstein-1
distance

    C_D[tau+p^(-1/2)+1/(p tau^4)]                   (2)

of a symmetric compound-Poisson variable plus an independent centered
Gaussian, whose TOTAL variance is exactly v. Its characteristic function
is pointwise at least exp(-v t^2/2). This last order, not Gaussianity,
is what proves the absolute-value upper bound. Taking tau=p^(-1/5)
proves (1). Zero coefficients and zero light variance cause no exception.

## 2. Light coefficients: offset-preserving scalar comparison

Split the edges into L={|c_e|<=tau} and B={|c_e|>tau}. Conditional on
M and all heavy signs, the light signed sum has variance

    V_M=sum_(e in M intersect L)c_e^2.

The scalar Wasserstein Stein bound, reconstructed earlier in this
campaign, replaces it by sqrt(V_M)G with error at most

    sum_(e in M intersect L)|c_e|^3 / V_M <=tau.

The ratio is defined as zero if V_M=0. This bound holds for EVERY
Lipschitz-one test, including any shift by the heavy sum.

Write v_L=E V_M. Matching indicators have mean 1/(p-1), zero joint
probability for distinct intersecting edges, and joint probability
1/[(p-1)(p-3)] for disjoint edges. Dropping the negative intersecting
covariances therefore gives

    Var(V_M) <= tau^2 v_L+2v_L^2/(p-3).             (3)

For v_L>0,

    E|sqrt(V_M)-sqrt(v_L)|
      <=sqrt(Var(V_M))/sqrt(v_L)
      <=tau+sqrt(2v_L/(p-3)).                       (4)

Use the SAME independent G on the two sides. The heavy sum can be
arbitrarily dependent on M: the coupling bound (4) remains valid. Since
v_L<=pD/[2(p-1)], this costs O_D(tau+p^(-1/2)). If v_L=0 all light
coefficients vanish, so the replacement costs zero.

We are left with the heavy signed matching sum plus an independent
Gaussian of deterministic variance v_L.

## 3. A finite Palm coupling bound, reconstructed

Let I be a finite vector of zero-one indicators, p_e=E I_e, and let
I^e have the conditional law of I given I_e=1, coupled with I in any
way. Let Z_e be independent Poisson(p_e) variables. For any function
F on nonnegative integer vectors with Lipschitz constant L in the
counting l1 metric,

    |E F(I)-E F(Z)|
      <= L sum_e p_e E||I-(I^e-delta_e)||_1.        (5)

Here is the proof rather than an appeal to a Poisson approximation
theorem with a different metric. Run independent immigration--death
processes with immigration rates p_e and per-particle death rate one.
Their stationary law is Z. The Stein solution

    g(x)=-integral_0^infinity [E_x F(X_t)-E F(Z)]dt

solves Ag=F-EF(Z). This integral converges for a Lipschitz F: couple
all immigrants and let unmatched initial/stationary particles die.
For two added particles at e,f, the mixed difference of F is zero
unless both survive, and is at most 2L otherwise. Thus

    |Delta_e Delta_f g|<=integral_0^infinity 2L e^(-2t)dt=L.

This includes e=f, treating the two particles separately. Telescoping
along additions/deletions proves Delta_e g is L-Lipschitz in l1.
Taking the generator expectation and conditioning the death term on
I_e=1 gives

    E Ag(I)=sum_e p_e E[Delta_e g(I)-Delta_e g(I^e-delta_e)],

which proves (5). Finite truncation followed by the integrable coupling
bound justifies the same argument for unbounded Lipschitz F.

The classical generator/Palm mechanism is in Chen--Xia,
[*Poisson process approximation: From Palm theory to Stein's method*](https://arxiv.org/abs/math/0702820),
Sections2--3, read directly. The unnormalized finite counting metric
and its bound (5) are reconstructed here; their normalized point-process
metric and its special constants are not silently substituted.

## 4. Heavy edges of a uniform matching

The heavy graph has maximum degree d_B<=D/tau^2 and at most
pD/(2tau^2) edges. For a heavy edge e={a,b}, condition a uniform
matching on containing e as follows. If it already contains e, do
nothing. Otherwise its two incident edges are {a,u},{b,w}; replace
them by {a,b},{u,w}. The result is uniform conditional on e: each
target matching containing e has exactly p-1 preimages (itself, and
two preimages for every other edge).

Let I and I^e record ONLY the heavy edges before and after this switch.
The expected discrepancy from the reduced conditional vector satisfies

    E||I-(I^e-delta_e)||_1
      <=1/(p-1)+3 d_B/(p-3).                        (6)

The first term is the event e was already present. On its complement
there are only two removed edges and one added edge. Each removed
partner is uniform outside the distinguished endpoints, and the added
pair is uniform among distinct remaining vertices; their probabilities
of being heavy are bounded by d_B/(p-3). Looser denominators are
intentional and valid for p>=6.

For a Lipschitz-one scalar test f, define

    F(z)=E f(sum_e sum_(j<=z_e)c_e epsilon_(e,j)+sqrt(v_L)G).

All signs and G in this definition are independent. Adding one point
changes its expectation by at most |c_e|<=sqrt(D), so (5)--(6) yield

    |E F(I)-E F(Z)| <= C_D/(p tau^4).               (7)

The original heavy signs have precisely the conditional law used in F.
Equations (3)--(7) prove (2) for every Lipschitz-one test.

## 5. Why large jumps help absolute response rather than obstruct it

The comparison variable Y has characteristic function

    E exp(itY)=exp[-v_L t^2/2
                  +sum_(e in B)(cos(t c_e)-1)/(p-1)].

It is real, nonnegative, and at least exp(-v t^2/2), since
cos u-1>=-u^2/2. For any integrable symmetric real X,

    E|X|=(2/pi) integral_0^infinity [1-E cos(tX)]dt/t^2.

Therefore E|Y|<=kappa sqrt(v). Combining this with (2) proves (1).
No claim of convex order or of comparison for arbitrary convex functions
is made. A nonzero jump part can persist. The scalar absolute-value
functional is exactly the one needed for the physical bridge response.

This infinitely-divisible moment comparison is NOT claimed new. A direct
primary match is Klebanov--Kakosyan,
[*Inequalities for m-Divisible Distributions*](https://www.qeios.com/read/Z7FB80/pdf)
(2024), Theorem1.5's alternate Levy--Khintchine proof and Theorem3.1.
The director read the complete seven-page source and independently checked
these two arguments; the r=1 case is exactly the calculation above.
Other claims in that source are not imported. In particular its displayed
absolute-moment asymptotic (1.6) fails for Gaussian convolution roots, and
Theorem1.6's equality-at-one-nonzero-point assertion fails for a symmetric
lattice compound-Poisson law. Neither issue affects the alternate proof
of Theorem1.5 or the fractional-moment comparison used here.

## 6. Consequence for the actual nonlocal sign law

In the [nonlocal Hadamard construction](paper_localization_nonlocal_sign_response_2026_09_17.md),
the exact-eigensector matching arrays have D<=4. The heavy-row deletion
on an eta-nearlevel set leaves arrays with D<=9. Thus its full nearcode
response bound improves to

    E|h dot x|/sqrt(n)
      <=1/sqrt(pi)+4sqrt(eta)+O(n^(-1/10)),          (8)

uniformly over that entire nearcode, where n=p^2. The separately proved
uniform exponential moments and deterministic identity-column repair
remain unchanged. This is an actual sign-law construction escaping the
Gaussian-angle covariance loss; it is specifically a Hadamard family,
not a theorem about arbitrary exact minimizing children.

Full parent control still needs the nearlevel geometry and all lower
energy words in the finite parent certificate. Equation (8) alone is
neither a convergence recurrence nor a new original cap bound.

## 7. Weighted Palm sensitivity improves the rate to p^(-1/4)

Localization-track refinement, independently suggested by the director.
The localization researcher read and reconstructed Sections 1--6 in
full and returned PASS. The refinement below keeps the coefficient
weights in the same finite coupling; it needs no new primary theorem.

For the function F in Section 4 and an arbitrary scalar Lipschitz-one
test f, adding one point of type e changes F by at most |c_e|. In the
immigration--death coupling the mixed difference for added particles
e and a is therefore at most 2 min(|c_e|,|c_a|) when both survive.
Integrating their survival probability gives the sharper Stein factor

    |Delta_e Delta_a g|<=min(|c_e|,|c_a|)<=|c_a|.

It follows by weighted telescoping that the right side of (5), for
this particular F, can be replaced by

    sum_e p_e E sum_a |c_a| |I_a-(I^e-delta_e)_a|.   (9)

All indices here are heavy edges. Their row and total weighted sums
obey

    max_v sum_(a incident v, a heavy)|c_a|<=D/tau,
    sum_(a heavy)|c_a|<=pD/(2tau),
    number of heavy edges<=pD/(2tau^2).             (10)

For the exact two-switch Palm coupling, the already-present event has
weighted cost |c_e|/(p-1). Otherwise at most two removed edges and
one added edge change. Conditional partner uniformity and the first
bound in (10) show

 E sum_a |c_a| |I_a-(I^e-delta_e)_a|
       <= |c_e|/(p-1)+3D/[tau(p-3)].               (11)

For the added edge, condition first on one endpoint: the other is
uniform among at least p-3 permitted vertices. The same denominator
also safely bounds both removed-edge contributions. This proves
(11) without treating the signs or coefficients as random.

Substituting (10)--(11) into (9) gives the explicit finite Poisson error

    pD/[2tau(p-1)^2]
       +3pD^2/[2tau^3(p-1)(p-3)].                  (12)

Thus the Wasserstein-one comparison to the SAME Gaussian--Poisson
variable in Section 5 has error at most

 (1+kappa)tau + kappa sqrt[pD/((p-1)(p-3))]
       +pD/[2tau(p-1)^2]
       +3pD^2/[2tau^3(p-1)(p-3)].                  (13)

Here the first tau is the light signed-sum replacement, and the
remaining light terms are (4) multiplied by E|G|=kappa. In particular,

    W_1(T,Y)<=C_D[tau+p^(-1/2)+1/(p tau^3)].

Taking tau=p^(-1/4) gives the sharpened theorem

    E|T|<=kappa sqrt(v)+C_D p^(-1/4).               (14)

The matched total variance and characteristic-function comparison of
Y are unchanged. Therefore the complete Hadamard nearcode bound (8)
improves, with n=p^2, to

    E|h dot x|/sqrt(n)
       <=1/sqrt(pi)+4sqrt(eta)+O(n^(-1/8)).         (15)

This is again a Gaussian absolute-response upper bound, not an assertion
of convergence in distribution to a Gaussian: the heavy compound-
Poisson component is retained throughout.

Independent audit: the Bernoulli researcher read and reconstructed all
Sections 1--7, including the weighted second Stein difference, exact
reduced-Palm switch, light variance replacement, and finite constants
in (12)--(13), and returned PASS. The localization researcher separately
reconstructed the director's Sections 1--6 before adding Section 7.

Finite replay: `computations/paper_localization_2026_09_17_matching_palm.py`
exhausts the matching state spaces at p=6,8,10 (15,105,945 matchings).
It verifies all 88 edge-Palm maps have exactly p-1 preimages per target,
3,846 weighted switch bounds, and 480 exact light-variance inequalities.
Output: `tmp/paper_portfolio_2026_09_17/localization/matching_palm_audit.json`.
