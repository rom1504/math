# Balanced mode bridges: an all-order positive construction

2026-09-17. New deductions within this campaign; external novelty not
established. The deterministic balancing mechanism was derived by the
localization track; the mixed-norm enlargement, Gaussian-information
entropy proof, and direct entropy-profile extension were developed jointly
with the director. The earlier fixed-mode tensor construction and its
near-level hierarchy are reconstructed in
[the Bernoulli track, Section 15](paper_bernoulli_2026_09_17.md).

The theorem below is an actual full-sign extension UNDER AN EXPLICIT
near-extreme geometry hypothesis. It neither establishes that hypothesis
for minimizing matrices nor proves convergence or a new cap constant.

## 1. Exact finite bridge and its balanced covariance certificate

Write n=kp+ell, 0<=ell<k, and let H_k be a real sign Hadamard matrix
with rows h_1,...,h_k. The first kp old coordinates are viewed as a
k by p array. Let q new columns be assigned deterministic mode labels
a(j), with each count q_a either floor(q/k) or ceil(q/k). Independently
for every j and t choose fair signs g_(j,t), and also independent fair
signs on the ell leftover coordinates. Set

    C_j=(h_(a(j)) tensor g_j, xi_j).                     (1)

Every physical entry is a sign. The random variables are independent
k-coordinate bundles, not independent physical edges. For w in R^n,
write w_0 for its array part and w_* for its leftovers, and put

    v_a(w)_t = sum_b h_(a,b) w_(b,t).

Hadamard orthogonality gives sum_a||v_a(w)||^2=k||w_0||^2.
For EVERY y in {+-1}^q,

    E exp[t w^T C y]
      <= exp{t^2[sum_a q_a||v_a(w)||^2+q||w_*||^2]/2}
      <= exp[t^2(q+k)||w||^2/2].                        (2)

The labels are fixed before sampling; the bound holds simultaneously
as a deterministic certificate for every w and y. This removes the
factor k paid when each column's mode is sampled independently.

There is an equally useful exact new-spin statement. Put

    B_w=max_y|w^T C y|=sum_j |w^T C_j|.

Changing one underlying sign changes B_w by at most twice its linear
coefficient's absolute value. The bounded-difference MGF inequality
therefore gives the SHARP proxy from (2), without the extra factor two
of generic absolute-value symmetrization:

    E exp[t(B_w-E B_w)] <= exp[t^2 V(w)/2],
    V(w)=sum_a q_a||v_a(w)||^2+q||w_*||^2
                                           <=(q+k)||w||^2. (3)

For completeness, expose the independent signs one at a time. The
Doob martingale difference has conditional range at most 2|coefficient|;
conditional Hoeffding bounds its MGF by exp(t^2 coefficient^2/2).
Multiply the bounds. No Gaussian concentration theorem is imported.

Cauchy--Schwarz also gives

    E B_w <=sqrt(q V(w))<=sqrt(q(q+k))||w||.              (4)

## 2. Intrinsic low-effective-mode code

Let P_a be the orthogonal projection onto the p-dimensional mode
span(h_a) tensor R^p in the kp-coordinate space. For a Boolean word x
define its effective mode count by

    R_eff(x) = [sum_a||P_a x_0||]^2/(kp)
             = [sum_a||v_a(x)||]^2/(k*kp).              (5)

It lies between 1 and k. For a cutoff u>=1, let F_n(u) consist of ALL
Boolean words with R_eff(x)<=u, with arbitrary leftover coordinates.
No repeated-column or exact mode-support representation is required.

The mean bridge response obeys, uniformly over this whole code,

    E B_x <=(q+k)sqrt(kp*u/k)+q sqrt(ell).              (6)

Indeed sum_a q_a E|v_a(x) dot g| is at most
(q/k+1)sum_a||v_a(x)||; leftovers contribute at most q sqrt(ell).
This uses only E|sum c_i g_i|<=sqrt(sum c_i^2), not a CLT.

The code includes every pure-mode word h_a tensor v, v in {+-1}^p,
and arbitrary leftovers. Uniform measure on that SUBSET has covariance
I_n. Consequently the Euclidean factorization complexity of F_n(u)
has squared value n, by the covariance lower bound and the identity
factorization upper bound. Thus these are not low-factorization centers.

### 2.1 Quantitative cardinality bound by a Gaussian observation

Put n0=kp. For any integer 1<=s<=k set e=u/s+s/k. Whenever e<1/2,

    log |F_n(u)|
      <= log binom(k,s)
         +(sp/2)log(1+k/s)+n0 h(e)+ell log2.            (7)

Here h is binary entropy in nats. This is a direct combination with
the Gaussian observation mechanism reconstructed in the localization
paper track, not an assumption about the code.

Proof: choose X uniformly from the main-coordinate code and let S be
the s modes of largest projection norm, with deterministic tie-breaking.
The elementary decreasing-sequence inequality gives

    ||P_(S complement) X||^2 <= (u/s)n0.

Conditional on S, observe Y=P_S X+G_S, where G_S is a standard Gaussian
on the sp-dimensional subspace. The Gaussian capacity bound yields

    I(X;Y|S) <= (sp/2)log(1+n0/(sp)).

The predictor sign(Y) has at most E||Y-X||^2 erroneous coordinates in
expectation: a sign error costs at least one in squared distance.
Its average error fraction is at most u/s+s/k. Binary conditional
entropy subadditivity and concavity of h give H(X|Y,S)<=n0 h(e).
Finally H(S)<=log binom(k,s) and the leftover coordinates cost ell log2.
The conditional information decomposition proves (7).

An independent net proof checks the same conclusion. Top-s projection
leaves squared distance at most (u/s)n0 to one of binom(k,s) subspaces
of dimension sp. A delta sqrt(n0) net, delta^2=u/s, followed by Boolean
cluster representatives gives an entropy bound

    log binom(k,s)+sp log(1+2/delta)+n0 h(4delta^2)+o(n0),

whenever 4delta^2<1/2. Count leftovers separately by 2^ell; adding them
to the Euclidean net dimension would impose an unnecessary logarithm.

### 2.2 Growing modes make the code cheap in scalar response

Choose Sylvester orders k=k_n with k_n->infinity and k_n=o(n), and
choose 1<=u_n=o(k_n). Set s=ceil(sqrt(u_n k_n)). Formula (7) gives

    log|F_n(u_n)|=o(n).                                (8)

Indeed s/k and u/s tend to zero, and
(s/k)log(1+k/s) tends to zero. The mode-set entropy is at most k log2
and the leftovers cost at most k log2, both o(n).

For every fixed epsilon>0 and q=floor(epsilon n), (6) gives

    max_(f in F_n) E B_f/n^(3/2) ->0.

Using (3) and (8), a union bound then gives with probability tending
to one

    max_(f in F_n,y) |f^T C y|=o(n^(3/2)).             (9)

The construction is all-order: only H_(k_n) is needed, and k_n can
always be chosen as a slowly growing power of two below n.

## 3. Actual full-sign extension under a full near-level cover

Let A_n be actual full-sign matrices. Write
E_n(eta)={x: |H_(A_n)(x)|>=Q(A_n)-eta n^(3/2)} and define

    rho(eta)=limsup_n max_(x in E_n(eta))
                       min_(f in F_n(u_n)) d_H(x,f)/n.

Assume rho(eta)<1/2 for sufficiently small eta and

    K_H=limsup_(eta down to0) rho(eta)h(rho(eta))/eta <infinity.

Then for every tau>2K_H and every sufficiently small FIXED epsilon>0,
there are actual full-sign extensions W_(n+floor(epsilon n)) retaining
A_n as their exact old principal block such that

    limsup_n [Q(W)-Q(A_n)]/n^(3/2)
                                  <=tau epsilon+epsilon^(3/2). (10)

Proof: for words x within Hamming radius r n of a center f, put w=x-f.
Then ||w||^2<=4rn. The number of (f,x) pairs is at most
|F_n| sum_(j<=rn)binom(n,j). Equations (3)--(4) show simultaneously
over those pairs, for every fixed r<1/2,

    B_(x-f)/n^(3/2)
         <=2epsilon sqrt(r)+sqrt(8epsilon r h(r))+o_n(1). (11)

The maximization over new words is already exact in B; no 2^q union
is hidden in (11). The center event (9) is paid once.

Choose K_H<kappa, a small outer level eta0, and a finite geometric
hierarchy eta_j=eta0 R^(-j), R>1, down to eta_J of order epsilon^2.
Strict radius majorants r_j satisfy r_j h(r_j)<=kappa eta_j.
In nonterminal shell j the old energy deficit is at least eta_j/R.
The main residual is paid using

    sup_(eta>=0)[sqrt(8epsilon kappa eta)-eta/R]
                                               =2kappa R epsilon.

The remaining 2epsilon sqrt(r_j) is at most 2epsilon sqrt(r_0),
which has arbitrarily small coefficient as eta0 tends to zero.
The terminal residual is O(epsilon^(3/2)) and is absorbed by the
strict margin in tau. Outside the outer near-level set, (2) and a
union over old and new words give a global bridge cap O(sqrt(epsilon))
n^(3/2), paid by eta0 for small epsilon.

Fill the new principal block with a full signing of cap at most
q^(3/2), available by the elementary Hoeffding union argument. Every
bridge and child entry is a sign, both polarities are retained, and all
events have simultaneous positive probability. Fix the hierarchy after
epsilon, then send n to infinity. This proves (10).

For a liminf-realizing family the necessary consequence is
K_H>=3c_*/4. This is a geometry obstruction, not an assertion that the
mixed-norm cover holds for any minimizing family.

## 4. All-offset transfer to the mode-conditioned Gaussian model

Keep the deterministic labels and replace each scalar g_(j,t) and
leftover xi by independent standard Gaussians. Denote the resulting
bridge C^G. The same k-coordinate bundle has matching covariance
h_a h_a^T, though it is not an isotropic physical-edge Gaussian block.

For ANY deterministic old and new children, the independent-block
theorem gives

    |E Q([A,C;C^T,D])-E Q([A,C^G;(C^G)^T,D])|
                                     <=C_epsilon n^(4/3) k^(2/3). (12)

To check the varying-k dependence, each bundle has dimension k and
Euclidean subGaussian constant sqrt(k). The tilted-moment replacement
bound costs C tau^2 n^2 k^2 exp(4tau^2 k^2), while softmax costs Cn/tau.
Choose tau=n^(-1/3)k^(-2/3). Thus k=o(n^(1/4)) makes (12)
o(n^(3/2)); for k<=n^(1/4-delta) it is a power saving.

Arbitrary child offsets are retained. This is a POSITIVE exact-sign
realization theorem for the singular, mode-conditioned Gaussian model,
not a replacement by iid edge Gaussians. The rank-one bundle covariance
and its deterministic orientation are indispensable data.

The exact-sign cover theorem in Section 3 only needs k=o(n); the
stronger k=o(n^(1/4)) restriction is solely for the separate generic
Gaussian replacement estimate (12).

**Stronger audited replacement:** the subsequent common-Gibbs
[symmetric frame theorem](paper_symmetric_frame_universality_2026_09_17.md)
improves (12) to O_epsilon(n^(5/4)k^(1/2)), valid with vanishing
normalized error for k=o(sqrt(n)). The third-order derivation above
is retained as an independently checked fallback, not the final range.

The archived fixed-frame and bounded-feature Gaussian/sign comparisons
already contain the fixed-k Lindeberg mechanism. The growing-mode,
deterministically balanced geometry and its exact covariance accounting
are the new ingredients here; no novelty is claimed for Lindeberg
replacement itself. The director's direct energy-profile criterion uses
(3), (6), and (7) without a Hamming cover and is maintained separately.

## 5. Audit and replay

Both peer researchers independently audited Sections 1--3 and the
original Section 4 and returned PASS. The replay
`computations/paper_localization_2026_09_17_balanced_modes.py` passed
30 finite covariance/count cases and 540 exact-law centered-absolute
MGF diagnostics, including nondivisible q and leftover coordinates.
Its output is `tmp/paper_portfolio_2026_09_17/localization/balanced_modes_audit.json`,
seed 2026091704. The cardinality examples are numerical evaluations of
the proved bound, not numerical evidence for optimizer applicability.

The discrepancy track's explicit near-half Hadamard family supplies a
strong stress test: its full near-level code is not close to these
low-effective-mode families in the relevant growing-mode range; see
[Sections 17--18 of that track](paper_discrepancy_2026_09_17.md).
Those matrices are not asserted to be actual asymptotic minimizers.
The positive construction remains conditional on the full-code geometry.
