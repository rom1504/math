# Vertex-swap stationarity as an exact overlap likelihood problem

This keeps the full joint cavity in the proposed vertex-swap argument. It
does not assert that pair swaps imply a linear rearrangement inequality.

## 1. Exact two-row insertion comparison

Fix balanced blocks S,T, within-block magnitude a and cross magnitude b,
and raw common scale λ. Let A minimize the two-sided partition over all
signings of this magnitude profile. Delete u∈S and v∈T. For the remaining
spins set

    U_S=Σ_{j∈S\{u}} A_uj xj,   U_T=Σ_{j∈T\{v}} A_uj xj,
    V_S=Σ_{j∈S\{u}} A_vj xj,   V_T=Σ_{j∈T\{v}} A_vj xj,
    f_i=λ(a I_S+b I_T),       g_i=λ(b I_S+a I_T),

where I is U or V. Omitting the direct uv edge, the original insertion
factor is cosh f_u cosh g_v, and swapping the two rows' sign patterns gives
cosh g_u cosh f_v. The direct edge changes either logarithm by at most λb.

Let π be the cavity law tilted by the original insertion factor, and put
L_i=log cosh f_i−log cosh g_i. Global optimality yields

    log E_π exp(L_v−L_u) ≥−2λb.                       (1)

For a quenched global branch field, the exact version averages the outer
logarithm in (1) over g. It does not move that logarithm outside the g
average. All quantities still use the two-row cavity law.

Even as λ=O(n^−1/2), the random row fields and L_v−L_u can be order one.
Consequently (1) does not imply E L_v≥E L_u up to a vanishing error. Its
log-mgf correction can contribute order one per row move, hence order n
over a macroscopic realignment. This differs from a single-edge flip,
whose field is uniformly O(n^−1/2).

## 2. Exact global formulation by hypercube convolution

Let D be the diagonal matrix of a balanced sign vector, defining the two
blocks. Set c=λ(a+b)/2 and d=λ(a−b)/2. The weighted interaction is

    J_D=cA+dDAD.

For raw-temperature one-sided partition functions Z_A^s(t), define
μ_(s,t)(x)=exp(stH_A(x))/Z_A^s(t). If independent X,Y have laws μ_(s,c)
and μ_(s,d), respectively, let p_s(D)=Pr(X⊙Y=D). Direct summation gives

    Z_(J_D)^s = Z_A^s(c) Z_A^s(d) p_s(D).            (2)

No Gaussian approximation or asymptotic limit is used. In particular,

    2^−n Σ_D Z_(J_D)^s=2^−n Z_A^s(c)Z_A^s(d).

For the absolute partition define C=Z_A^+(c)Z_A^+(d)+Z_A^−(c)Z_A^−(d),
choose a common branch s with probabilities proportional to those two
products, and then draw X,Y independently in that branch. Its overlap
law p satisfies

    Z_(J_D)^abs=C p(D).                              (3)

Thus an actual globally optimized anisotropic signing, with its fixed
balanced partition D, minimizes this exact overlap likelihood among all
balanced gauges of the same A. Indeed changing balanced D is equivalent
to a vertex permutation of A in the fixed magnitude profile, which is an
allowed competitor. The local vertex swaps are two-coordinate exchanges
of this balanced overlap vector.

For width pressure, the corresponding exact objective is the geometric
mean sqrt(p_+(D)p_−(D)), multiplied by a D-independent factor. One must
not replace the shared-branch law in (3) by independent branch choices.
With a quenched field hg, the two branch weights are additionally multiplied
by exp(±hg); the optimized overlap objective is E_g log p_g(D), not the
logarithm of an averaged overlap law.

The hoped-for rearrangement theorem therefore requires control of a
minimum balanced overlap point probability for actual optimizing Gibbs
laws. Ordering single-row mean fields is not the exact statement.

## 3. A genuine but insufficient Hubbard lower estimate

For one branch, represent μ_(s,c) and μ_(s,d) independently as product
mixtures with latent means M_i and M_i'. Conditional on the two latent
fields,

    p_s(D)=2^−n E ∏_i(1+D_i M_i M_i').

Each one-coordinate latent mean has a symmetric distribution, because
the branch Ising law has no external spin field; the two copies are
independent. Jensen therefore proves the uniform bound

    log p_s(D) ≥−nlog2 + (1/2)Σ_i E log(1−M_i² M_i'^2). (4)

This is valid for every D, without conditional-PSD assumptions. It is not
sharp enough to pay the reheating loss: it discards all dependence among
latent coordinates and can charge for the auxiliary diagonal noise in
the product-mixture representation. The exact expression (2), not (4),
must be retained in any claimed normalization-matching comparison.

For perspective, balanced D forces the *uniform-spin* covariance

    E_uniform H_A(X)H_A(DX)=−n/2,

whereas E H_A(X)²=n(n−1)/2. Thus these two energies are asymptotically
uncorrelated under uniform spins, but this second-moment fact does not
control their exponentially tilted overlap likelihoods in (2). No
finite-moment-to-pressure inference is made here.

## 4. A finite hot-child comparison that really follows

Write L_A(beta) for the average of its two branch log partitions with
the order-normalized temperature beta, and Psi_N=min_A L_A. For even N,
the convolution implies the unconditional inequality

    2 Psi_(N/2)(beta)
      <= 2 Psi_N(beta/sqrt(2)) - log binom(N,N/2).     (5)

To prove it, choose A attaining the right-hand optimized pressure. In
(2) set c=d=beta/sqrt(2N), and average over balanced D. Each p_s is a
probability distribution on all gauges, so

    average_balanced log p_s(D) <= -log binom(N,N/2).

There is therefore a common balanced D for which the average of the
two branch log likelihoods obeys the same bound. Its interaction is
zero across the two blocks and has within-block raw temperature
2c=beta/sqrt(N/2). The branch-averaged logarithm factors over those
blocks, and is at least 2 Psi_(N/2)(beta), proving (5).

The finite correction is explicit. With psi_N=Psi_N/N this gives

    psi_N(beta) >= (1/2) psi_(N/2)(sqrt(2) beta)
                    +(1/(2N)) log binom(N,N/2).

It improves the elementary restriction inequality by retaining an
entropy term, but still contracts the coefficient of the hot-child
pressure by one half. It is not a same-temperature almost-additive
comparison, and does not supply the temperature payment sought in
the main campaign. No claim of novelty over entropy restriction
inequalities is intended.
