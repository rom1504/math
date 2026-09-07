# Actual seed augmentation at the quadratic-logarithmic barrier

2026-09-07. **Proved derivation; independent audit requested.** This strengthens the robust-center construction in `principle_invent_2026_09_07_balanced_valley_augmentation.md`. Every original edge is retained. The conclusion is an actual extension theorem and an optimizer-specific soft-excitation consequence, not convergence.

## 1. Setting and tuned full-sign bridge

Let A_n be symmetric hollow full signs. Let G_n be a nonempty family of Boolean centers, represented as G_n=L_n V_n with column norms of V_n at most one and row norms of L_n at most gamma_n. Assume gamma_n^2=o(n). Write

    r_n(x)=d_H(x,G_n)/n,
    Delta_n(x)=[Q(A_n)-|H_(A_n)(x)|]/n^(3/2).

The primary Gram--Schmidt Walk theorem and conventions were independently checked in `principle_construct_2026_09_07_balanced_augmentation_audit.md`. It produces a full sign coloring with variance proxy 40 times the squared norm of each test direction.

Put g_n=gamma_n^2/n, a_n^2=sqrt(g_n), and b_n^2=1-sqrt(g_n), for all sufficiently large n. Apply the theorem to

    w_i=(a_n (V_n)_:i,b_n e_i).

These vectors have norm at most one. For a nearest center g to x, test against (L_g/a_n,(x-g)/b_n). Thus an actual random full sign column v satisfies

    E exp(lambda v dot x)
       <= exp{lambda^2 [40 n sqrt(g_n)
                           +160 n r_n(x)/(1-sqrt(g_n))]/2}.       (1)

This tuning makes the center term o(n) while removing the factor-two loss of the equally weighted feature sum. No growing-dimensional Gaussian approximation is used.

Sample q=floor(epsilon n) independent columns. Sauer's lemma and VC(G_n)<=gamma_n^2 give log|G_n|=o(n). Union over Hamming shells, all centers, both tails, and all new spins gives one actual full-sign bridge C such that, uniformly in x,y,

    |x^T C y|/n^(3/2)
       <= F_epsilon(r_n(x))+o_n(1),

    F_epsilon(r)=sqrt(320 epsilon r [h(r)+epsilon log 2]).         (2)

Here epsilon is fixed before n tends to infinity. The uniform error follows by the same shell count as the audited robust-center theorem; the normalized variance in (1) converges uniformly on r in [0,1]. Choose a full-sign new child with cap at most q^(3/2). For ANY lower bound Delta_n(x)>=f(r_n(x)), with f(0)=0, the actual parent obeys

    limsup_n [Q(A_n^+)-Q(A_n)]/n^(3/2)
      <= D(epsilon):=epsilon^(3/2)
              +sup_(0<=r<=1) [F_epsilon(r)-f(r)].                (3)

The old principal block is exactly A_n. Equation (3) pays the absolute bridge in full.

## 2. The critical rate

Suppose f is uniformly positive away from zero:

    inf_(delta<=r<=1) f(r)>0 for every delta>0,

and define

    kappa=liminf_(r downarrow 0) f(r)/[r^2 log(e/r)] in (0,infinity].

Then

    limsup_(epsilon downarrow 0) D(epsilon)/epsilon <=80/kappa,  (4)

with 80/infinity=0.

Proof. Since h(r)<=r log(e/r), splitting the square root yields

    F_epsilon(r)
      <=sqrt(320 epsilon) r sqrt(log(e/r))
           +sqrt(320 log 2) epsilon sqrt(r).                    (5)

Fix any kappa'<kappa and choose delta small enough that f(r)>=kappa' r^2 log(e/r) for 0<r<=delta. For w=r sqrt(log(e/r)),

    sqrt(320 epsilon) w-kappa' w^2 <=80 epsilon/kappa'.

The second term of (5) is at most sqrt(320 log 2) epsilon sqrt(delta). On r>=delta, the positive term tends uniformly to zero while f stays uniformly positive, so that region eventually contributes no positive maximum. Divide (3) by epsilon, take epsilon to zero, then delta to zero and kappa' to kappa. The child term epsilon^(3/2) vanishes after division. For kappa=infinity, apply this with arbitrarily large fixed kappa'.

In particular, a barrier growing faster than r^2 log(e/r) near its centers permits actual augmentation at incremental cost o(epsilon)n^(3/2). This includes all fixed power barriers kappa_0 r^alpha with alpha<2, but also barriers such as r^2 log(e/r) log log(e^e/r), which the power theorem does not cover. A sufficiently large finite quadratic-log coefficient already gives profitable dilution.

## 3. Consequence for a liminf-realizing seed family

Let c_inf=liminf_n M_n/n^(3/2)>0, and suppose Q(A_n)/n^(3/2) tends to c_inf along the selected orders. Under the hypotheses of Section 2, necessarily

    kappa <=160/(3 c_inf).                                     (6)

Otherwise (4) gives a small fixed epsilon with D(epsilon)<c_inf[(1+epsilon)^(3/2)-1], contradicting the definition of c_inf using the ACTUAL parents from (3).

The same conclusion has a useful version without a single fixed modulus f. Assume the centers satisfy gamma_n^2=o(n), and are macroscopically isolated in the following precise sense:

    for each delta>0,
    liminf_n inf_{x:r_n(x)>=delta} Delta_n(x)>0.                 (7)

Use inf(empty)=infinity. Define the small-distance stiffness

    kappa_* = lim_(delta downarrow 0) liminf_n
        inf_{x:0<r_n(x)<=delta}
              Delta_n(x)/[r_n(x)^2 log(e/r_n(x))].              (8)

Then kappa_*<=160/(3c_inf). Indeed, if it were larger, choose a smaller fixed coefficient still above this threshold. For each sufficiently small fixed delta the lower bound holds eventually in n throughout 0<r<=delta. Condition (7) pays r>=delta. The same proof of (4), with n taken to infinity before epsilon and delta decrease, gives the contradiction. No uniform exchange of limits is required.

Thus one of two things must occur for a low-factorization-complexity center family along liminf minimizers: there are macroscopically separated states with vanishing normalized deficit, or there are small-relative-distance collective excitations whose deficit is at most a fixed multiple of r^2 log(e/r). This is an optimizer-specific assertion about centers, not the generic linear complexity of a thick energy level.

If G_n contains every exact absolute ground state, the latter excitations have positive energy deficit. Full-sign energies lie in one parity class and have spacing at least two. Any bounded-ratio excitation from (8) consequently satisfies

    d_H(x,G_n)^2 log(e n/d_H(x,G_n)) >=c sqrt(n),

and therefore d_H(x,G_n)>=c' n^(1/4)/sqrt(log(e n)). The forced excitations are not merely one-vertex degeneracies. This lower size statement is conditional on taking the soft-excitation branch and on including all exact grounds among the centers.

## 4. Scope and collision checks

This theorem does not assert that actual minimizers have low-complexity ground centers or a favorable barrier; it proves that sufficiently rigid low-complexity landscapes could be improved by a paid, actual full-sign extension. It does not establish comparable-size flatification.

Random-edge entropy perturbation was considered first and discarded as a new direction because it overlaps `decisive_independent_gaussian_stability_2026_09_06.md` and its cited earlier stratified-noise work. The present proof uses the newly audited balanced bridge, not that perturbation. Targeted searches found no earlier critical quadratic-log barrier theorem. The finite-coefficient and mesoscopic-excitation conclusions above are the additional content beyond the original alpha<2 augmentation.
