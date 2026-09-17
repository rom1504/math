# Spectral entropy and balanced absolute-response augmentation

2026-09-17. **Proved below; the complete application, entropy argument,
bounded-difference step and limit order were independently reconstructed
by all three researchers and passed.**
This is a combination of a Gaussian information bound, an exact full-sign
construction, and scalar absolute-response concentration. No external
novelty is claimed. It does not establish a hypothesis for minimizing
sequences or prove convergence of the original normalized optimum.

**Applicability correction,18:45UTC:** the finite entropy results remain
valid, but the direct full-nearlevel hypothesis K_s<infinity in Section4
is IMPOSSIBLE for every bounded positive-cap quadratic signing sequence.
Random flips of a ground word force s(eta)>=h(rho) for rho<eta/(4c).
Thus that conditional application is vacuous, not a live optimizer target.
See paper_director_nearlevel_response_rigidity_2026_09_17.md. Anchored
center-cover theorems with a vanishing residual variance are distinct.

## 1. Entropy from concentration in a small union of linear spaces

Let P_1,...,P_k be mutually orthogonal projections of rank p resolving
the identity on R^n, n=kp. For x in the Boolean cube define

```math
R(x)=\frac{(\sum_{a=1}^k\|P_ax\|_2)^2}{n},\qquad 1\le R(x)\le k.
```

For 1<=r<=k let F_r={x:R(x)<=r}. For any integer 1<=s<=k such that
e=r/s+s/k<=1/2,

```math
\boxed{\log|F_r|\le\log\binom{k}{s}
 +\frac{sp}{2}\log(1+k/s)+n h(r/s+s/k).}                 \tag{1}
```

Here h is binary entropy with natural logarithms. The assertion is
vacuous, and need not be used, when F_r is empty.

**Proof.** For each x retain the s largest values ||P_a x||, breaking
ties deterministically, and call their index set S(x). For nonnegative
decreasing v_a, sum_(a>s) v_a^2 <= (sum_a v_a)^2/s. Hence

```math
\|(I-P_{S(x)})x\|_2^2\le rn/s.
```

Take X uniform on F_r. Conditional on S, observe
Y=P_S X+G_S, where G_S is independent standard Gaussian on the
sp-dimensional image of P_S. Its conditional mutual information is
at most (sp/2)log(1+n/(sp)): Gaussian maximum entropy and the
arithmetic-geometric mean bound on the covariance determinant give
this directly, since tr Cov(P_S X|S)<=n. No isotropy of X is assumed.

Round each coordinate of Y to its nearest sign. A wrong coordinate
costs at least one in squared error, so the mean fraction of errors is
at most r/s+sp/n=e. Bitwise conditional entropy, followed by concavity
of h, gives H(X|Y,S)<=n h(e). Ties may use either sign. Therefore

```math
H(X)=H(S)+I(X;Y\mid S)+H(X\mid Y,S)
```

is bounded by the right side of (1). This proves the claim. This is the
Gaussian-observation mechanism reconstructed from El Alaoui--Montanari,
used here to count an actual set of Boolean states, not to replace it
by independent posterior coordinates.

If k=k_n=o(n), r=r_n=o(k), and s=ceil(sqrt(rk)), (1) implies
log|F_r|=o(n). Indeed s/k->0, r/s->0, and log binomial(k,s)<=k log2.
All estimates concern finite sets before limits. An independent proof
using Euclidean nets and Boolean Hamming balls is preserved in the
localization researcher's balanced-mode artifact.

## 2. An exact sign construction realizing the useful response

Let H_k be a real sign Hadamard matrix and write n=kp initially.
Identify x with a k by p array. Let

```math
P_a=(h_a h_a^T/k)\otimes I_p.
```

Choose q column types a_j deterministically with counts q_a differing
from q/k by less than one. Independently for each column draw p fair
signs g_j and put C_j=h_(a_j) tensor g_j. Thus C is a FULL sign matrix.
For arbitrary old word x define B_x=max_y|x^T C y|=sum_j|x^T C_j|.
Two exact estimates hold:

```math
\mathbb E B_x\le(q+k)\sqrt{nR(x)/k},                         \tag{2}
```

```math
\mathbb E\exp\{t(B_x-\mathbb EB_x)\}
 \le\exp\{t^2(q+k)n/2\}.                                  \tag{3}
```

For (2), the field on mode a is a sum of independent signs with
squared coefficient norm k||P_a x||^2, so its mean absolute value is
at most sqrt(k)||P_a x||. Sum with q_a<=q/k+1.
For (3), flipping scalar sign g_(j,t) changes B_x by at most twice
the corresponding field coefficient. The bounded-difference MGF is
therefore exp(t^2 V_x/2), with

```math
V_x=k\sum_aq_a\|P_ax\|_2^2\le(q+k)n.
```

No random-mode averaging, independently paid channels, or cancellation
between child and bridge energies is used. The absolute new-spin
maximization is performed EXACTLY before taking concentration bounds.
If k divides q then V_x=qn for every x, even though the JOINT response
process is not the independent-edge process.

For arbitrary n take p=floor(n/k), n0=kp and ell=n-n0. Fill the ell
leftover entries of each column with independent fair signs. Then
(2) becomes

```math
\mathbb EB_x\le(q+k)\sqrt{n_0R(x_{[n_0]})/k}+q\sqrt\ell,
```

and (3) remains valid with n in place of n0. The entropy bound gains
ell log2. Thus all-order limits are valid for k->infinity and k=o(n),
using dyadic k if desired. Neither all-order Hadamards nor a change of
the old signing is required.

## 3. Direct near-level entropy criterion

Let A_n be actual hollow full sign matrices with bounded normalized
caps, and put

```math
E_n(\eta)=\{x:Q(A_n)-|H_{A_n}(x)|\le\eta n^{3/2}\}.
```

For the balanced construction above, allow k=k_n with k_n=o(n).
Suppose deterministic functions s(eta), mu(eta) satisfy, at every
sufficiently small fixed positive eta,

```math
\limsup_n n^{-1}\log|E_n(\eta)|\le s(\eta),\qquad
\limsup_n\max_{x\in E_n(\eta)}\sqrt{R(x_{[n_0]})/k_n}
 \le\mu(\eta).
```

Assume mu(eta)->0 and
K_s=limsup_(eta downarrow0) s(eta)/eta < infinity. Then, for EVERY
tau>K_s/2 and every sufficiently small fixed epsilon>0, there are
full-sign parents P_n preserving A_n as an exact principal block,
with q=floor(epsilon n), such that

```math
\limsup_n\frac{Q(P_n)-Q(A_n)}{n^{3/2}}
 \le\tau\epsilon+\epsilon^{3/2}.                         \tag{4}
```

**Proof and limit order.** At a fixed level, (3) and a union bound over
the actual code give, with probability tending to one after a strict
arbitrarily small slack,

```math
\max_{x\in E_n(\eta)}B_x/n^{3/2}
 \le\epsilon\mu(\eta)+\sqrt{2\epsilon s(\eta)}+o_n(1).
```

Choose K'>K_s, then a small outer eta0 and a geometric mesh of ratio
R0>1. On the shell E(eta) minus E(eta/R0), subtract the available
old-energy deficit eta/R0. The supremum of
sqrt(2 epsilon K' eta)-eta/R0 is K' R0 epsilon/2. The first term is
at most epsilon mu0, where mu0 can be made arbitrarily small by
shrinking eta0. Stop the FINITE hierarchy at eta_J of order epsilon^2;
the innermost entropy contribution is O(epsilon^(3/2)), absorbed in
the strict tau margin after epsilon is small enough.

Outside the outer code, (2)--(3) and a union over all old words give
a bridge bound epsilon+sqrt(2epsilon log2)+o(1), which is below eta0
for sufficiently small fixed epsilon. Its old deficit pays the bridge.
Fill the new principal block by a sign matrix of cap at most q^(3/2),
which exists by the elementary random-sign union bound. The triangle
inequality pays its entire cap. This respects the exact child-flip
identity. First fix K',R0,eta0 and epsilon and its finite hierarchy;
then send n to infinity. Choose the strict margins so the excess over
K_s epsilon/2 is less than (tau-K_s/2)epsilon. This proves (4).

## 4. What was genuinely combined; what remains open

(1) turns spectral concentration into a rigorous count using Gaussian
information. (2) constructs exact signs whose mean ABSOLUTE responses
exploit that concentration. (3) uses balanced mode allocation to retain
full aggregate variance control without paying modes independently.
Together they can verify the two hypotheses of Section 3 from a single
explicit mode-energy estimate. Neither covariance alone nor a Gaussian
channel alone supplies this construction.

For example, if the normalized effective mode count on E_n(eta) is
at most theta(eta)+o_n(1), (1) gives an explicit entropy upper bound
of order sqrt(theta) log(1/theta), while (2) gives mu<=sqrt(theta).
No such estimate is asserted for actual optimizing sequences. The
current recursive-Hadamard stress tests specifically challenge it.

If A_n is a liminf-realizing minimizing sequence with normalized limit
c_*, (4) implies K_s>=3c_* whenever its hypotheses hold: otherwise
choose K_s/2<tau<3c_*/2 and then small epsilon to contradict the
definition of liminf at the enlarged orders. This is a conditional
landscape restriction, not an improved value of c_* or convergence.

Status distinction: the finite entropy inequality and construction are
unconditional; the optimizer-facing implication requires the explicitly
displayed entropy/response bounds. Growing n is not the same operation
as growing a polynomial degree, and no Bohnenblust--Hille conclusion is
silently being imported here.
