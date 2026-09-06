# A global ceiling for the terminal response functional, and what it does not bound

Date: 2026-09-06. Status: structural proof independently reconstructed by
the director and bound-audit agent. The agent's separate exact rational
certificate proves a uniform bound strictly below .45.

## 1. Statement for the whole measurable Gaussian family

Let F be any bounded odd function and H any nonnegative even function
on the countable Gaussian creation space, with |F|+H<=1. Set

    K=U*P1F, r=||P1F||_2=||K||_2,
    J(F,H)=E[H Psi(K,||F||_2^2-r^2)],
    Psi(k,t)=E|k+sqrt(t)N|.

This is the entire terminal full-response functional used in the current
lower-bound campaign, not merely a finite bank or a chosen mask family.
Then

    J(F,H)<=sup_(0<p<1) B(p),
    B(p)=(1-p)Psi(m(p)/sqrt(1-p),p-m(p)^2),
    m(p)=2varphi(Phi^(-1)(1-p/2)).                  (1)

The scalar envelope has floating maximum about .4496994558913. The
independent rational certificate gives the rigorous, slightly looser
uniform bound

    sup_(F,H)J(F,H)
      <=.449998927491931736970926031607336182570891990067547942316383
      <.45.                                      (2)

No actual-purification or finite-frame hypothesis is needed for this
UPPER bound. Put p=E|F| and h=1-p. Replacing H by 1-|F| and replacing
the variance ||F||^2-r^2 by p-r^2 both increase J. These are simply
positivity and Gaussian convexity. If p is zero or one, J is zero.

## 2. The two inequalities responsible for the ceiling

First, boundedness controls the first Gaussian chaos. For r>0 the
Gaussian G=P1F/r is standard and E[FG]=r. Hence

    r<=E[|F||G|]<=m(p).                            (3)

The last inequality is the upper-tail rearrangement of |G| against a
weight in [0,1] of mean p. Also m(p)^2<p, so all following variances
are positive for 0<p<1.

Second, although Psi is convex in k, it is CONCAVE in k^2:

    x -> Psi(sqrt(x),t) is concave for t>0.         (4)

Indeed, with k=sqrt(x) and s(k)=2Phi(k/sqrt(t))-1, its derivative is
s(k)/(2k). This is decreasing because

    s(k)=2 integral_0^(k/sqrt(t))varphi(v)dv
          >=2(k/sqrt(t))varphi(k/sqrt(t)).

The continuous k=0 extension causes no problem. Weighted Jensen under
(1-|F|)/h, followed by E[(1-|F|)K^2]<=r^2, gives

    J(F,H)<=h Psi(r/sqrt(h),p-r^2).                (5)

For fixed p let u=r^2 and k=sqrt(u/h). Differentiation gives

    d/du [h Psi(sqrt(u/h),p-u)]
       =s(k)/(2k)-h varphi(k/sqrt(p-u))/sqrt(p-u)
       >=p varphi(k/sqrt(p-u))/sqrt(p-u)>0.         (6)

Thus insert r=m(p) into (5), proving (1). These inequalities depend
only on Gaussian first-chaos geometry and pointwise feasibility; no
limitation of a particular certificate library is being imported.

The exact upper certificate and its independent proof are in
`artifacts/resumed_bound_audit_full_response_uniform_upper_2026_09_06.md`,
`computations/resumed_bound_audit_full_response_upper_certificate_2026_09_06.py`,
and the correspondingly named results JSON. It uses 109 rational
alpha intervals and monotone outward bounds for all alpha>=0, not a
floating maximizer or a finite grid of unsupported point evaluations.

## 3. Consequence for repeated policy or frame enrichment

Every fixed or infinite measurable terminal pair reduced to this SAME
J obeys (2). Therefore none of the following, by itself, can make this
terminal certificate reach .45:

* enlarging the Hermite bank or the actual Gaussian frame;
* solving the fixed-frame policy optimization globally;
* further causal births, purification gates, or feasible line searches;
* passing to an infinite variational limit of such terminal pairs.

The assertion concerns the value finally CERTIFIED BY J. It does not
say that those algorithms' actual output energies are bounded by J.
It also does not bound a different correlated response formula, a
multi-output energy certificate, or an actual iterative matrix process
whose energy is retained through the updates instead of being reduced
back to a single J. The fixed-GFOM/Haar ceiling from the convergence
track has a different, broader algorithmic scope and must not be
identified with this functional ceiling.

## 4. Exact endpoint energy shows the missing quantity explicitly

For any symmetric matrix B write Q_B(x)=x^T Bx/2. For a feasible pair
F,C with |F_i|+|C_i|<=1, put

    D=Q_B(F)+Q_B(C), E=F^T B C.

Then exactly

    Q_B(C+F)=D+E, Q_B(C-F)=D-E,
    max(|Q_B(C+F)|,|Q_B(C-F)|)=|D|+|E|.           (7)

The terminal response theorem lower-bounds the cross term for
C=H sign(BF). It discards the nonnegative contribution |D|. Therefore
the ceiling (2) does not even cap the energies of its own two actual
endpoint vectors.

The already exact R_16 example in
`resumed_response_iteration_and_twin_obstruction_2026_09_06.md`
illustrates the distinction. Its f,c satisfy

    f^T R_16 c=28, f^T R_16 f=-12, c^T R_16 c=4.

Thus E/16^(3/2)=7/16=.4375, whereas |D|/16^(3/2)=1/16=.0625,
and the better endpoint energy is .5. The tensor lifts are actual
signings after negligible diagonal deletion. This example does not
give a universal lower bound on D, but rules out interpreting (2) as
an upper bound on actual endpoint energies.

To get beyond the terminal barrier, a proof must retain genuinely new
energy information: for example control |Q_B(F)+Q_B(C)|, or exploit
coherent first-chaos/response correlations omitted from the restricted
marked-plus-nonlinear channel, or use another coupled-output mechanism.
No positive universal floor for any such missing term is established
here. In particular the new ceiling does not prove convergence,
nonconvergence, or an upper bound below 1/2 for the original M_n.

A separate, more restrictive discriminator is developed in
`resumed_response_flat_involution_terminal_endpoint_2026_09_06.md`:
for fixed equivariant odd-tree terminal rules on normalized symmetric
Hadamards, matrix-sign charge plus finite-rule universality forces the
common self-energy to vanish. That additional theorem concerns this
randomized architecture on a special matrix class; it is not a
consequence of the scalar inequality (1) alone.
