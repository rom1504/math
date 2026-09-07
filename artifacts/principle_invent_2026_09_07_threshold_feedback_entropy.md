# Quantitative feedback into relaxed stability, and an exact-cleanup obstruction

2026-09-07. **Proved finite statements.** The feedback kernel reaches
threshold-stable states and has an explicit entropy cost. Exact local
maxima are a different target; the final section gives a concrete
delocalized row-regular obstruction to entropy-neutral exact cleanup.

## 1. A deterministic energy-improving kernel with controlled information loss

Let B be any symmetric hollow N by N matrix and H(x)=x^T Bx/2. Write
M=max_x H(x). For a>0 call x a-stable if

    x_i(Bx)_i>=-a for every i.                            (1)

Define T_a by repeatedly flipping the least-index coordinate violating
(1), with any fixed deterministic tie convention. Every actual flip
increases H by more than 2a. Hence the procedure terminates, its output
is a-stable, and, pointwise,

    H(T_a x)>=H(x),
    d_H(x,T_a x)<=number of flips<=(M-H(x))/(2a).          (2)

Repeated flips of one coordinate do not invalidate the second inequality:
final Hamming distance is at most the total number of flips.

Let X have ANY law rho on the cube and Y=T_aX. Put

    r=(M-E H(X))/(2aN),
    hbar(r)=h2(min(r,1/2)),   r>=0.

Then

    H_Shannon(Y)>=H_Shannon(X)-N hbar(r).                 (3)

Indeed let E=X XOR Y be the coordinate error word. Conditional on Y,
E determines X, so H(X|Y)<=H(E)<=sum_i h2(P(E_i=1)). If the mean error
fraction is at most 1/2, entropy concavity bounds this by N times the
binary entropy of that mean, and (2) bounds the mean by r. If r>1/2,
the universal N log2 bound gives exactly hbar(r). As T_a is deterministic,
H(Y)=H(X)-H(X|Y). This proves (3) without an assumed bound on basin size.

The bound depends on the INPUT'S expected energy deficit, not merely on
the full energy range. In particular, for R=max H-min H and R/(2aN)<1/2,
the uniform entropy loss is at most N h2(R/(2aN)), tending to o_a(1)N as
a tends to infinity. The threshold a is fixed before an N-limit.

## 2. A rigorous relaxed-stable partition lower bound

Let A_a(B) be the set in (1), and

    Z_beta,a(B)=sum_(y in A_a(B)) exp(beta H(y)).

The Gibbs variational inequality applied to the pushed-forward law Y
gives, for every rho,

    log Z_beta,a(B)
       >=beta E_rho H+H_Shannon(rho)
           -N hbar((M-E_rho H)/(2aN)).                   (4)

This is a lower bound for the SAME relaxed-stable partition counted by
the actual upper construction. It is not merely the name of a missing
conditional entropy: (4) gives an explicit dimension-free cost.

For the two energy orientations let M_+=max H, M_-=-min H and
w=(M_++M_-)/(2N), the normalized width. Apply T_a separately to a pair
of input laws rho_+,rho_- and put

    j=[E_(rho_+) H-E_(rho_-) H]/(2N),
    s=[H_Shannon(rho_+)+H_Shannon(rho_-)]/(2N).

Concavity of hbar yields

    [log Z_beta,a(B)+log Z_beta,a(-B)]/(2N)
       >= beta j+s-hbar((w-j)/(2a)).                     (5)

For the successful marked product means, the old limit has j>=J(F,H)
and s->s(F,H), with the notation in
`principle_invent_2026_09_07_gaussian_saddle_weak_duality.md`.
If w<=u+o(1), monotonicity in j gives the explicit limit lower bound

    beta J(F,H)+s(F,H)-hbar((u-J(F,H))/(2a)),              (6)

provided u>=J(F,H), as already forced by the zero-temperature lower
theorem. Ordered finite-frame approximations remain unchanged.

## 3. The actual stability-plus-selector upper extends to this exact target

For the weave B=C/sqrt(N), its normalized physical field in the orthant
lemma is exactly sigma x_i(Bx)_i. Thus a-stability replaces each inequality
Z_i>=0 by Z_i+a>=0. This adds a deterministic vector of norm a sqrt(k).
The existing diagonal allowance D=2 becomes D=2+a. Nothing else changes:
the independent edge tilt, light variance condition, convex Lipschitz
constant, Finner step, selector count, and diffuse/concentrated split all
remain valid. For each FIXED a<infinity, the Gaussian hinge lower bound
g_a is still strictly positive, so the full signing construction gives
some strict pressure saving Delta_a>0:

    limsup (1/N)log E Z_beta,a(sigma B)
       <= A_cert(p,t)-Delta_a/p.                        (7)

Here beta=2t/sqrt(p), A_cert=[t+p log2+e_bar]/p, and the finite depth is
chosen AFTER the positive saving is reserved. This is not an estimate
uniform as a grows with N. Every true extremum belongs to A_a, so (7)
also remains an upper-cap certificate.

Combining (5)--(7) gives a genuine quantitative weak duality for the new
filtered object. The cost is explicit and generally nonzero. It must not
be omitted when interpreting the stable partition as a control saddle.
No claim is made that the present constants produce a new useful lower
number: the proved upper saving can be extremely small, whereas the
entropy cost in (6) need not be.

## 4. Exact cleanup cannot be entropy-neutral in the broader row-regular class

Let N=2m and let O be a real orthogonal m by m matrix such that
(Oy)_i and (O^T x)_j are nonzero for every Boolean x,y and every i,j.
Define

    B=[ 0  O ; O^T  0 ].                                (8)

Then B is symmetric hollow, B^2=I, every row has squared norm one,
and ||B||op=1. For a spin pair (x,y), exact one-spin stability requires

    x=sign(Oy),      y=sign(O^T x).                       (9)

In particular y uniquely determines x, so there are at most 2^m exact
stable configurations. Every probability law supported on them has
entropy at most m log2. Starting with the uniform law on the full cube,
ANY kernel whose output is always an exact stable point therefore loses
at least (N/2)log2 entropy, regardless of whether it is deterministic,
randomized, or energy-improving.

This is also a delocalized sequence example. Haar orthogonal O satisfies
the nonzero condition almost surely for each m: each individual forbidden
linear equation has probability zero, and there are only finitely many.
Moreover max_ij |O_ij|->0 in probability, by the elementary spherical
coordinate tail and a union bound over m^2 entries. Thus one may choose
a deterministic sequence with both properties. Its cap obeys Q(B)<=N/2.

The obstruction is for a weighted bipartite involution, not an original
full-sign matrix. It shows precisely that row regularity, bounded norm,
delocalization and a cap at most N/2 do NOT justify negligible-entropy
cleanup to exact local maxima. Any stronger exact-cleanup statement for
the actual full-sign weave would need a further property and proof.

The positive threshold kernel in Sections 1--3 avoids this false claim:
it quantifies its entropy cost and matches a relaxed stability filter
whose actual upper theorem is still strict for every fixed threshold.
