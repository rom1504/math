# Twisted chiral doubling: joint selection and exact convergence obligations

Campaign: 2026-09-18, 17:09:11--20:09:11 UTC. This is the independent
uniform-theorem track. Convention throughout:

\[
H_A(x)=\sum_{i<j}a_{ij}x_ix_j,\qquad Q(A)=\max_x|H_A(x)|,
\qquad M_n=\min_A Q(A),\qquad f(n)=M_n/n^{3/2}.
\]

The older `chiral_scale_preserving_lift.md` and
`dependent_lift_analytic_audit.md` use twice this norm. Their compressed-lift
obstructions do **not** automatically apply to the twisted double studied here.
The two-multiplier power-saving convergence implication was already proved in
`fresh_limit_algebra_2026_09_05.md`; section 5 below sharpens its accounting and
records explicit countermodels to weaker claims.

## 1. Two exact forms of the objective

Let B be switching/permutation-equivalent to A, let d have sign entries, put
C=B+diag(d), and set

\[
D=\begin{pmatrix}A&C\\C&-A\end{pmatrix}.
\]

For arbitrary Boolean x,y,

\[
H_D(x,y)=H_A(x)-H_A(y)+x^TCy.
\]

Flipping all coordinates of x leaves its internal energy unchanged and
negates the bridge. Therefore

\[
\boxed{Q(D)=\max_{x,y}\bigl(|H_A(x)-H_A(y)|+|x^TCy|\bigr).} \tag{1}
\]

Equivalently write y=t*x, I={i:t_i=1}, J=I^c, z=x. Then

\[
H_A(x)-H_A(y)=2z_I^TA_{IJ}z_J,
\]

\[
x^TCy=2H_{B[I]}(z_I)-2H_{B[J]}(z_J)+d(I)-d(J).
\]

This proves the proposed cut/internal-energy identity exactly, with the two
terms evaluated at the **same** partition and spin configuration.

Let D_0 denote the same matrix with d=0. The diagonal matching has n edges,
so for every choice of d,

\[
\boxed{|Q(D)-Q(D_0)|\le n.} \tag{2}
\]

Thus diagonal optimization is important at finite order but cannot remove a
fixed positive leading-scale defect.

For the complex representation, put u=(x+y)/2, v=(x-y)/2, w=u+iv. Every
w_i lies in {1,-1,i,-i}, and

\[
\boxed{H_D(x,y)=\operatorname{Re}\bigl(w^T(C-iA)w\bigr).} \tag{3}
\]

The real linear operator D represents the antilinear map
z -> (A+iC) conjugate(z). Consequently

\[
\|D\|_{op}=\|A+iC\|_{op}.
\]

This identity does not replace the four-phase Boolean constraint by an
arbitrary Euclidean vector.

## 2. A genuine quantitative joint-selection lemma

The following sufficient condition involves only two **marginal** profiles of
the seed. It neither assumes a favorable twist nor defines a new name for the
desired minimum.

For 0<=k<=n and nonnegative integers a,b, define

\[
K_k(a)=\#\{(x,y):x,y\in\{\pm1\}^n,
\ d_H(x,y)=k,\ |H_A(x)-H_A(y)|=a\},
\]

\[
J_k(b)=\#\{(u,v):u,v\in\{\pm1\}^n,
\ d_H(u,v)=k,\ |u^TAv|=b\}.
\]

Write J_k(>r)=sum_{b>r}J_k(b). For L>=2Q(A), set

\[
\mathcal E_A(L)=
\sum_{k=0}^n\sum_a
\frac{K_k(a)J_k(>L-a)}{2^n\binom nk}. \tag{4}
\]

**Selection lemma.** If E_A(L)<8, a signed permutation g exists such that,
with B=g^TAg,

\[
\boxed{Q\begin{pmatrix}A&B\\B&-A\end{pmatrix}\le L.} \tag{5}
\]

Every diagonal sign completion then has cap at most L+n. A uniform random
signed permutation succeeds with probability at least 1-E_A(L)/8.

**Proof.** The signed permutation group acts transitively on the ordered
Boolean pairs at a fixed Hamming distance. The orbit at distance k has
2^n binom(n,k) elements. For fixed (x,y), a uniform g therefore makes
(gx,gy) uniform on that orbit, and

\[
\Pr_g\{|H_A(x)-H_A(y)|+|x^Tg^TAg y|>L\}
=\frac{J_k(>L-|H_A(x)-H_A(y)|)}{2^n\binom nk}.
\]

Summing proves that (4) is the exact expected number of violating ordered
pairs, without any independence assumption between different pairs. No pair
with y=+/-x violates when L>=2Q(A). On all other pairs, independent global
reversals of x,y and their interchange form an eight-element symmetry orbit
on which the tested expression is constant. Thus the number of violating
orbits is an integer, with expectation E_A(L)/8. If that expectation is below
one, some g has no violation; Markov gives the probability claim. QED.

The finite group permits conditional-expectation derandomization, although
this statement makes no polynomial-time complexity claim. Permutations enter
materially: random switching alone is transitive only after the entire
coordinatewise product x*y is fixed, rather than just its Hamming weight.

## 3. Exact tests and their limitation

`computations/twisted_chiral_uniform_2026_09_18.py` computes (4) in rational
arithmetic from exact integer profiles. Its separate order-four check averages
the actual violation count over all 384 signed permutations, at 13 thresholds,
and exactly matches (4). Output and all seed matrices are in
`computations/results/twisted_chiral_uniform_2026_09_18.json`.

Selected values:

| seed | Q(A) | smallest L certified by E_A(L)<8 | 2 sqrt(2) Q(A) |
|---|---:|---:|---:|
| stored order 5 minimizer | 4 | 16 | 11.314 |
| stored order 6 minimizer | 5 | 22 | 14.142 |
| stored order 10 minimizer | 13 | 54 | 36.770 |
| dependent chiral order 12 witness | 20 | 76 | 56.569 |

These are **sufficient certificate values**, not minima over twists. For the
order-12 witness the exact expected count at L=56 is

\[
\mathcal E_A(56)=72777149/16128>8.
\]

Thus this first-moment certificate does not establish the target there. It
does not imply that all twists fail: bad configurations may occur in large,
strongly dependent clusters, and favorable twists can exist despite a large
average number of violations. A uniform theorem would require a substantially
stronger profile inequality or a different selection argument.

## 4. What a doubling theorem actually supplies

Suppose an admissible selection rule is available at every step, and

\[
Q(A_{j+1})\le 2^{3/2}Q(A_j)+C n_j,
\qquad n_{j+1}=2n_j.
\]

Then direct summation gives, for every depth r,

\[
\frac{Q(A_r)}{(2^r n)^{3/2}}
\le\frac{Q(A_0)}{n^{3/2}}
+\frac{C(1+\sqrt2)}{2\sqrt n}. \tag{6}
\]

No compression identity is used here. Conversely, if the selection theorem
only holds for exact global minimizers, it gives a recurrence for M_n by
reselecting an actual minimizer at each size; one may not claim it literally
iterates the particular constructed children.

More generally let normalized one-step defects be bounded by epsilon(n), let
e(n)=sup_{m>=n}epsilon(m), and suppose

\[
S(n):=\sum_{j=0}^{\infty}e(2^j n)<\infty,
\qquad S(n)\longrightarrow0. \tag{7}
\]

Every increasing multiplier word whose letters are at least two has total
normalized defect at most S(n). In particular a power saving suffices. A bare
epsilon(n)->0 assertion does not imply (7).

## 5. Coverage: sufficient supplements and explicit countermodels

### Two multiplicatively independent integer multipliers

If, for k=2 and k=3 and all sufficiently large n,

\[
f(kn)\le f(n)+\epsilon(n), \tag{8}
\]

with (7), then f(n) converges. Indeed every descendant n 2^a 3^b has
normalized cap at most f(n)+S(n). The semigroup {2^a 3^b:a,b>=0} has relative
gaps tending to one: a finite set of residues b log 3 modulo log 2 forms an
arbitrarily fine net, and the remaining sufficiently large logarithm is
covered by a nonnegative multiple of log 2. Principal deletion gives
M_m<=M_N for m<=N. Choosing the least descendant N>=m therefore yields

\[
\limsup_{m\to\infty} f(m)\le f(n)+S(n).
\]

Taking n along a liminf subsequence proves the claim. The same proof uses any
two multiplicatively independent integer multipliers, with the corresponding
minimum-generator Dini tail.

### One doubling rule plus a selected proportional-thinning rule

Instead of a second multiplier, it is enough that the constructed dyadic
descendants A_N admit, uniformly for N/2<=m<=N, a principal set S of size m
with

\[
Q(A_N[S])\le (m/N)^{3/2}Q(A_N)+\eta(N)N^{3/2},
\qquad \eta(N)\longrightarrow0. \tag{9}
\]

Only **one** thinning step is used: for each target m choose the first dyadic
descendant N>=m. The normalized thinning error is at most 2^{3/2}eta(N), so
(6), or (7), again implies limsup f<=liminf f. No summability of eta is needed.
Condition (9) is needed only on the chosen descendants, not on all signings
and not on a uniform average of principal subsets. Earlier counterexamples to
average proportional contraction must not be overstated as refuting this
selected-subset statement.

### Why doubling alone is insufficient, even with an O(1) error

Set alpha=3/2 and, for t>=1,

\[
F(t)=t^\alpha\{c+\rho\sin(2\pi\log_2 t)\},
\]

where c>rho>0 and

\[
\rho\sqrt{\alpha^2+(2\pi/\log2)^2}<\alpha c.
\]

Then F is strictly increasing, F(2t)=2^alpha F(t), and F(n)/n^alpha does not
converge. Thus even exact doubling plus monotonicity leaves multiplicative
phase oscillation. One may take c=.46 and rho=.005, keeping this countermodel
inside the currently reported asymptotic interval. Rounding F(n) upward to an
integer with the parity of binom(n,2) changes it by less than two. For all
sufficiently large n the rounded sequence is increasing, has increments
O(sqrt(n))<=n, and satisfies the doubling upper inequality with additive error
two. Hence parity and the elementary one-vertex extension bound do not repair
the coverage gap. This is a countermodel to a proposed implication, **not** a
claim that the actual M_n behaves this way.

### Why two multipliers with only o(n^{3/2}) error are insufficient

For large t let

\[
G(t)=t^\alpha\{c+\rho\sin(\log\log(t+e))\},\qquad 0<\rho<c.
\]

This is eventually increasing and its normalized values do not converge.
For each fixed k, the mean-value theorem gives

\[
G(kn)/(kn)^\alpha-G(n)/n^\alpha=O(1/\log n).
\]

Consequently it satisfies normalized o(1) upper defects simultaneously for
k=2 and k=3. The corresponding dyadic Dini tail diverges. This countermodel
separates the two independent obligations: multiplicatively dense coverage
and summable accumulated error.

## 6. Complexification: exact useful theorem, wrong input norm

The primary paper by
[Muñoz, Sarantopoulos and Tonge, *Complexifications of real Banach spaces,
polynomials and multilinear maps* (1999)](https://www.researchgate.net/profile/Andrew-Tonge-2/publication/228582051_Complexifications_of_real_Banach_spaces_polynomials_and_multilinear_maps/links/0c9605227ec723f09c000000/Complexifications-of-real-Banach-spaces-polynomials-and-multilinear-maps.pdf)
was checked directly, especially Propositions 18 and 20. Proposition 20 does
give norm-preserving extension of a real quadratic polynomial in the
Lindenstrauss--Tzafriri complexification norm. It is not norm preservation on
the ordinary complex cube. For a four-phase vector w=u+iv with both disjoint
supports nonempty, the relevant squared norm is

\[
\|w\|_{LT}^2=
\sup_\theta\{\|u\cos\theta-v\sin\theta\|_\infty^2
+\|u\sin\theta+v\cos\theta\|_\infty^2\}=2.
\]

The ordinary complex cube/Taylor norm of w is one. Thus applying the
norm-preserving theorem to the actual parent states pays a factor two in the
quadratic value. This is an input-domain loss, not a removable normalization
convention.

An elementary refinement for the **untwisted** choice B=A is

\[
Q(D_0)\le 2(1+\sqrt2)W(A)\le 2(1+\sqrt2)Q(A). \tag{10}
\]

Here W is the half-width of the seed energy interval. To verify it, at a
disjoint pair u,v set d=|H_A(u)-H_A(v)| and t=|u^TAv|. Comparing the two cube
points au+bv and bu-av cancels the average of H_A(u),H_A(v), so

\[
W(A)\ge
\begin{cases}
t,&d\le t,\\
(d^2+t^2)/(2d),&d\ge t.
\end{cases}
\]

The parent value is 2(d+t). Optimizing the displayed ratio gives (10), with
t/d=sqrt(2)-1 in its nontrivial branch. This refinement still has a fixed
leading loss and is recorded only to prevent misuse of complexification as a
2 sqrt(2) transfer.

A separate natural chiral identity is also unfavorable. If J is a signed
complex structure with J^2=-I and AJ=-JA, the symmetric companion C=JA gives

\[
H_{\left[\begin{smallmatrix}A&JA\\JA&-A\end{smallmatrix}\right]}(x,y)
=H_A(x-Jy),
\qquad Q\left(\begin{smallmatrix}A&JA\\JA&-A\end{smallmatrix}\right)=4Q(A).
\tag{11}
\]

The upper bound follows because (x-Jy)/2 lies in the real cube; equality uses
x=z, y=Jz for a seed extremizer z. This weighted companion has zeros at the
matching positions and is not automatically in the stipulated orbit family.
Filling its missing matching entries changes the norm by only O(n), so this
particular algebraic completion retains leading factor four, not 2 sqrt(2).

### Bilinear slack does not supply an exact finite isometry

The independent computation track exhaustively checked the order-six global
minimizer. It is centered, has A^2=5I, Q(A)=5 and beta(A)=12, hence
beta(A)<2 sqrt(2) Q(A) by a strict margin. Nevertheless its 384 distinct orbit
bridges have core-cap histogram

\[
\#\{B:Q(D_0)=18\}=104,\qquad
\#\{B:Q(D_0)=22\}=280.
\]

Even independent arbitrary shifts in every partition have minimum possible
maximal half-width 16, still larger than 10 sqrt(2). See
`computations/results/twisted_chiral_2026_09_18_independent_r6_core.json`.
Thus centering, exact conference spectrum, and strict beta slack do not imply
the zero-error joint inequality. This finite statement does not refute an
O(n)-defect asymptotic theorem.

## 7. A variance-flat matching twist, and its increment obstruction

For this section only, A has independent standard Gaussian upper-triangle
entries and zero diagonal. Take a signed permutation G with
G^T=-G and G^2=-I, and let B=G^T A G. Such G are signed perfect matchings.
The parent core process is

\[
X_G(x,y)=H_A(x)-H_A(y)+(Gx)^T A(Gy).
\]

For arbitrary vectors p,q,u,v, direct coefficient multiplication gives

\[
\operatorname{Cov}(H_A(p),u^TAv)
=(p\cdot u)(p\cdot v)-\sum_i p_i^2u_iv_i.
\]

For Boolean x,y the diagonal terms cancel in the difference, and skew
symmetry gives

\[
\operatorname{Cov}(H_A(x)-H_A(y),(Gx)^TA(Gy))=0.
\]

The two separate variances are n^2-(x dot y)^2 and
n^2+(x dot y)^2-2n. Consequently

\[
\boxed{\operatorname{Var}X_G(x,y)=2n(n-1)\quad\text{for every }(x,y).}
\tag{12}
\]

This equals the variance of the independent Gaussian signing core on 2n
vertices with its n matching edges omitted. Equal variances do not give a
Gaussian maximum comparison. To see the exact remaining issue, put

\[
a=x\cdot x',\ b=x\cdot y',\ c=y\cdot x',\ d=y\cdot y',
\]

\[
E=x\cdot Gx',\ F=x\cdot Gy',\ K=y\cdot Gx',\ H=y\cdot Gy',
\qquad r=\sum_i x_i y_i x'_i y'_i.
\]

The covariance is

\[
\operatorname{Cov}(X_G(x,y),X_G(x',y'))
=\tfrac12\{(a+d)^2-(b-c)^2\}+(E-H)(F+K)-2r. \tag{13}
\]

Subtracting the independent-core covariance gives

\[
\boxed{\Delta=n-r-\tfrac12(b-c)^2+(E-H)(F+K).} \tag{14}
\]

The sign is genuinely indefinite even on balanced pairs x dot y=0 and
x' dot y'=0. Use G=diag(J_2,J_2), J_2=[[0,-1],[1,0]], and repeat the following
four-coordinate patterns r times, giving n=4r:

| sign | x | y | x' | y' | covariance difference |
|---|---|---|---|---|---:|
| negative | ++++ | ++-- | +++- | -+-- | -(3/2)n^2+n |
| positive | ++++ | ++-- | +-++ | --+- | n^2+n |

Neither second state belongs to the first state's global/chiral symmetry
orbit. The disagreement is therefore not just a redundant antipodal index,
and its magnitude is of order n^2, not a negligible O(n) covariance error.
This refutes a direct covariance-order/Slepian proof for this natural
reference process; it does not refute every possible Gaussian comparison or
the desired deterministic seed-selection theorem.

`computations/twisted_chiral_gaussian_matching_2026_09_18.py` independently
checks all 65,536 covariance entries at n=4 by integer coefficient dot
products. The explicit replicated witnesses are retained in its JSON output.

### Actual sign-matrix subfamily check

An independent direct integer pair-energy enumeration of all signed perfect
matchings modulo G->-G gives:

| child | matching G count | distinct B | minimum core cap |
|---|---:|---:|---:|
| order 4 minimizer | 6 | 5 | 12 |
| order 6 minimizer | 60 | 36 | 18 |
| order 8, class 0 | 840 | 753 | 32 |
| order 8, class 1 | 840 | 808 | 32 |

At order six only one of the 36 bridges attains core cap 18; all other 35
have cap 22. The computation track also exhausts every diagonal matching in
this subfamily: both order-eight classes have minimum completed cap 32,
whereas their unrestricted signed-permutation families attain 30. The
variance-flat subclass is therefore not even finite-optimal there.

All stored best completed witnesses at child orders 3--10 were independently
checked by the pair-energy formula in
`computations/twisted_chiral_matching_core_audit_2026_09_18.py`; its output is
`computations/results/twisted_chiral_matching_core_audit_2026_09_18.json`.
Odd-order partial matchings are explicitly marked as not skew-symmetric;
the exact variance identity (12) is not asserted for them.

## Linear cube-pullback certificates have an asymptotic factor-four floor

This is a proof-class obstruction, not an obstruction to the twisted family.
Let D_0 be any order-2n twisted core, so its n matching edges are zero and
all its other off-diagonal entries have modulus one. Suppose, modulo arbitrary
diagonal matrices,

\[
D_0=\sum_r\lambda_r R_r^T A_r R_r+E,
\tag{15}
\]

where every A_r is an order-n hollow full signing, every real linear map
R_r:ell_infinity^{2n}->ell_infinity^n has norm at most one, and E is hollow
and symmetric. The A_r may be signed permutations of the selected seed and
all choices may depend on that seed. Then

\[
\boxed{\sum_r|\lambda_r|
\ge 4-\frac{4\sqrt{6n}\,Q(E)}{n(n-1)}.}
\tag{16}
\]

In particular, a remainder Q(E)=o(n^{3/2}) still forces total triangle
coefficient at least 4-o(1), rather than 2 sqrt(2).

**Proof.** Write |F|_1=sum_{ij}|F_ij| for entrywise mass, including the
diagonal. Each row of R_r has absolute sum at most one. Hence

\[
|R_r^T A_rR_r|_1
\le\sum_{i\ne j}|(A_r)_{ij}|
 \left(\sum_u|(R_r)_{iu}|\right)
 \left(\sum_v|(R_r)_{jv}|\right)
\le n(n-1).
\]

Meanwhile |D_0|_1=4n(n-1). Ignoring the diagonal in (15) and taking
entrywise masses yields

\[
4n(n-1)\le n(n-1)\sum_r|\lambda_r|+|E|_1.
\]

For any order-N hollow symmetric E, let S_i=sum_j E_ij epsilon_j with
independent uniform signs. Its variance is sigma_i^2=sum_j E_ij^2 and its
fourth moment is at most 3 sigma_i^4. Interpolation of first, second, and
fourth absolute moments gives E|S_i|>=sigma_i/sqrt(3). Therefore

\[
\beta(E)=\max_y\sum_i|(Ey)_i|
\ge\frac1{\sqrt3}\sum_i\|E_{i,*}\|_2
\ge\frac{|E|_1}{\sqrt{3N}}.
\]

Polarization and multilinear maximization on the real cube give
beta(E)<=4Q(E), so |E|_1<=4 sqrt(3N)Q(E). Substitute N=2n to obtain (16).

The argument allows arbitrary cancellation at the matrix-identity level;
it only says that subsequently bounding every pulled-back seed energy by
Q(A_r) and adding absolute coefficients cannot certify the desired factor.
It does not address nonlinear/state-dependent cube maps, bounds that exploit
joint cancellation of the seed energies, or a successful direct mixed-energy
inequality. Diagonal constants in the pullbacks must also be accounted for
in any actual certificate; discarding them in the mass calculation only
weakens this necessary condition. Filling the matching changes the leading
factor-four statement only by O(1/n).

## A positive uniform theorem for all bipolar two-macrotype models

This section was developed on resumption, 2026-09-19. It is a genuine
selection theorem on a specified class, not a claim for near-minimizers.
The director proposed the equal-block model; the argument below proves its
exact formula and extends it to arbitrary block proportions.

Consider the homogeneous two-macrotype kernel

\[
T=\begin{pmatrix}1&t\\t&-1\end{pmatrix},\qquad 0\le t\le1,
\]

on macroblocks of masses p,q, where p+q=1. Swapping macroblocks and negating
the kernel preserves the relevant norms, so assume p>=q. Choose a switching
sign s balanced **within each macroblock**, and use B=SAS, with no
permutation. Then the exact continuous aggregate-model cap of its twisted
core, in units N^2, is

\[
\boxed{F_{p,q,t}=\max\{\beta_{p,q,t},\ p^2+\tfrac12q^2+tpq\},\qquad
\beta_{p,q,t}=p^2+tpq+q|q-tp|.}
\tag{17}
\]

The corresponding seed cap is

\[
Q_{p,q,t}=
\begin{cases}
\tfrac12p^2(1+t^2),&tp\le q,\\
\tfrac12(p^2-q^2)+tpq,&tp\ge q.
\end{cases}
\tag{18}
\]

Consequently F_{p,q,t}<=4Q_{p,q,t} for every permitted p,q,t. The constant
four is sharp at p=q=1/2,t=0. At equal proportions, (17) simplifies to

\[
F_{1/2,1/2,t}=\max\{\tfrac12,\tfrac38+\tfrac14t\}.
\tag{19}
\]

In particular, throughout 0<=t<=1/2 this chosen switching attains the
twist-invariant bilinear lower bound exactly in the homogeneous model.
At t=sqrt(2)-1, it removes the untwisted model's factor 2(1+sqrt(2)) loss.

### Short analytic proof by four diamond faces

For a spin x, let a_i be its ordinary signed aggregate on macroblock i and
b_i its aggregate weighted by s. Balanced switching gives the diamond
constraints |a_1|+|b_1|<=p and |a_2|+|b_2|<=q. For y write c_i,d_i
similarly. Twice the parent energy is

\[
a_1^2-a_2^2-c_1^2+c_2^2
+2t(a_1a_2-c_1c_2)+2b^TTd.
\]

For fixed magnitudes a_i,c_i, their signs can make both displayed mixed
seed terms nonnegative, independently of the bridge signs. The bridge
maximum is attained with |b|=(p-a_1,q-a_2) and
|d|=(p-c_1,q-c_2). For nonnegative budgets B_1,B_2,D_1,D_2 its exact value is

\[
\beta_T(B,D)=B_1D_1+B_2D_2+tB_1D_2+tB_2D_1
-2\min\{B_1D_1,B_2D_2,tB_1D_2,tB_2D_1\}.
\]

Indeed the four edge signs of T have negative product, so an extremal
choice leaves exactly the lightest edge unsatisfied. Equivalently this
bridge norm is a maximum of four bilinear functions of the budgets.
The resulting full objective is separately convex in a_1 and c_2:
their square coefficients are positive, and the bridge is a maximum of
affine functions of either variable. It therefore suffices to examine
(a_1,c_2)=(p,q),(p,0),(0,q),(0,0). Put a=a_2 and c=c_1.

For the first three faces, twice the parent energy is respectively at most

\[
\begin{array}{ll}
(p,q):&p^2+q^2+2tpq-a^2-c^2+2tac,\\
(p,0):&p^2+2q^2+2tpq-a^2-c^2+2tac-2qa-2tqc,\\
(0,q):&2p^2+q^2+2tpq-a^2-c^2+2tac-2pc-2tpa.
\end{array}
\]

All three are maximized at a=c=0, because
-a^2-c^2+2tac<=0. On the final face, p>=q and t<=1 imply

\[
\beta_T=p(p-c)+tpq+(q-a)|q-t(p-c)|.
\]

Thus its objective is

\[
-a^2-c^2+2p(p-c)+2tpq+2(q-a)|q-t(p-c)|.
\]

It decreases with a. After setting a=0, its one-sided derivatives in c
are -2c-2p plus or minus 2tq, both nonpositive. Its maximum is therefore
also at a=c=0 and equals 2 beta_{p,q,t}. Of the other three endpoint
values, 2p^2+q^2+2tpq is largest. Every retained endpoint is feasible in
the diamonds, proving (17), including attainment.

To obtain (18), maximize (u^2-v^2)/2+tuv on
|u|<=p,|v|<=q. The positive maximum has |u|=p and
|v|=min(tp,q). It dominates the negative maximum when p>=q.
Finally, writing r=q/p, if t<=r then beta=p^2+q^2<=4Q and

\[
4Q-(p^2+q^2/2+tpq)
=p^2(1-r^2/2-tr+2t^2)>0.
\]

The bracket is at least 1/2-t+2t^2>0, since r<=1. If t>=r then
beta=2Q and the same remaining difference is
p^2-(5/2)q^2+3tpq>=p^2+q^2/2>0.

### Scope and finite full-signing realizations

The formulas are exact for the continuous weighted macrotype model. For a
finite block matrix, diagonal removal and balancing-rounding cost O(N).
Replacing its t-valued cross block by full signs with spectral error
||E|| gives seed-cap error at most N||E||/2 and core-cap error at most
2N||E||. Hence the positive factor-four transfer becomes

\[
Q(D)\le4Q(A)+4N\|E\|+O(N).
\tag{20}
\]

Full-sign cross-block approximations with ||E||=O(sqrt(N)) therefore give
an actual full-signing family with asymptotic ratio at most four, sharp in
the balanced t=0 limit. This is an O(N^{3/2}) error, not a subleading error
at the near-minimizer scale; here the seed cap itself has order N^2.
Thus this theorem neither supplies the requested near-minimizer transfer
nor closes the global convergence problem. It does give a nontrivial,
explicit switch selection that repairs the entire two-macrotype source
of the untwisted 2(1+sqrt(2)) bound.

## Universal diamond-resource theorem and a full-signing clone corollary

The two-macrotype proof suggests a more general positive theorem. Unlike
the preceding four-face proof, the following identity works for an arbitrary
number of types and has no spectral assumption.

**Diamond theorem.** Let A be any real hollow symmetric matrix. For any
real vectors a,b,c,d with

\[
|a_i|+|b_i|\le1,\qquad |c_i|+|d_i|\le1\quad\text{for every }i,
\]

one has

\[
\boxed{|H_A(a)-H_A(c)+b^TAd|\le4Q(A).}
\tag{21}
\]

More generally the same assertion holds for any symmetric A if Q(A) is
interpreted as the supremum of |x^TAx/2| on the entire real cube. For hollow
A this equals the Boolean-cube maximum by coordinatewise multilinearity.
The coefficient four is sharp, as shown by disconnected positive and
negative equal-size clique blocks in the large-block limit.

For hollow A there is the stronger unconditional invariant estimate

\[
\boxed{\Delta(A)\le2W(A)+\tfrac12\beta(A)\le4W(A).}
\tag{21a}
\]

To see this, reduce to diamond vertices by coordinatewise affinity and use
the complementary supports I of a+h,a-h and J of v in the proof below.
The six-form identity is at most U+V+U_I+V_J, where U=max H_A and
V=-min H_A. Polarization gives
max_(I partition J)(U_I+V_J)=beta(A)/2: write x=u+v and y=u-v for signs u,v
on disjoint complementary supports, so x^TAy=2(H_A(u)-H_A(v)). Applying
the same argument to -A controls the negative side. This proves (21a).

**Proof by nonlinear clipping.** Switch coordinates so b_i>=0; this changes
neither the hypotheses nor Q(A). Set, coordinatewise,

\[
h_i=(b_i+d_i-1)_+,\qquad
u=b+d-h,\quad v=b-h,\quad w=d-h.
\]

Then 0<=h_i<=b_i and |w_i|<=|d_i|. In particular each of
u,v,a+h,a-h,c+w,c-w belongs to [-1,1]^n. The exact algebraic identity

\[
b^TAd=H_A(u)-H_A(v)-H_A(w)+H_A(h)
\]

and parallelogram identities give

\[
\begin{split}
H_A(a)-H_A(c)+b^TAd
={}&H_A(u)-H_A(v)\\
&+\tfrac12H_A(a+h)+\tfrac12H_A(a-h)\\
&-\tfrac12H_A(c+w)-\tfrac12H_A(c-w).
\end{split}
\tag{22}
\]

The absolute coefficients sum to four. This proves (21).

The clipping h is essential: these are nonlinear, state-dependent cube
maps, not a fixed linear pullback decomposition of the original doubled
matrix. The construction is consistent with the earlier factor-four
linear-certificate floor but is not an instance of its restricted proof
class. An exploratory finite LP first exposed (22); the proof above is
exact and independent of that computation.

### A selectable full-signing seed class with a uniform factor-four lift

Given any full signing A of order n, form a full signing of order N=2n
by duplicating every vertex into a twin pair:

\[
\widetilde A=A\otimes J_2+
\operatorname{diag}_{i=1}^n
\begin{pmatrix}0&e_i\\e_i&0\end{pmatrix},\qquad e_i\in\{-1,1\}.
\]

Choose S to have opposite signs on each twin pair and use the explicit
allowed twist B=S tilde(A) S. Let D be its order-2N completed chiral double,
with any matching signs. Then

\[
Q(D)\le16Q(A)+6n
\le4Q(\widetilde A)+10n
=4Q(\widetilde A)+5N.
\tag{23}
\]

To verify this, first omit the n twin edges and the 2n parent matching
edges. If x is a spin on the 2n twin coordinates, its average a_i and
signed average b_i on each pair satisfy |a_i|+|b_i|=1. For y use c,d.
The ideal parent energy is exactly
4[H_A(a)-H_A(c)+b^TAd], bounded by 16Q(A) using (21).
Restoring the child twin-edge matrix has operator norm one and perturbs
the parent block matrix by operator norm at most two, so its cap cost is
at most 4n. The parent matching costs at most 2n. Finally the unfilled
child has cap exactly 4Q(A), and filling its n twin edges changes the cap
by at most n. These estimates prove both inequalities in (23).

This is a uniform theorem for an explicit, selectable full-signing seed
class, with a genuinely O(N) defect and an explicit switching. It does
not prove coefficient four for every seed. Nor does it supply the desired
2 sqrt(2) coefficient: exact twin cloning multiplies the normalized child
cap by sqrt(2). In fact the preserved lower bound
L=0.4333221116640807 and upper bound C<0.493609 rigorously exclude this
clone class from near-minimality: Q(tilde(A))>=4M_n-n implies its normalized
liminf is at least sqrt(2)L>0.6128, strictly above C. An o(N^(3/2))
quadratic-norm perturbation cannot close that leading gap. The lower and
upper inputs are audited in
`artifacts/decisive_audit_fresh_full_lower_chain_2026_09_07.md` and
`artifacts/flatify_independent_2026_09_07_ternary_upper_proof.md`.
Its value is a general positive joint-energy
mechanism, independent of the seed's spectrum or its number of macrotypes.

## Uniform all-seed factor four at the coarse dense scale

The diamond theorem also yields a uniform statement for every full signing,
but its error is much too large for the near-minimizer problem. The precise
error scale matters.

**Coarse uniform selection theorem.** For every order-n full signing A and
every epsilon>0, there is a diagonal switching S such that its completed
twisted double, with B=SAS and any matching signs, satisfies

\[
\boxed{Q(D)\le4Q(A)+16\epsilon n^2+4kn+n,
\qquad k\le4^{\lceil\epsilon^{-2}\rceil}.}
\tag{24}
\]

In particular, for n>1 choose
epsilon^2=2 log(4)/log(n). Then k<=4 sqrt(n), giving

\[
Q(D)\le4Q(A)
+16\sqrt{\frac{2\log4}{\log n}}\,n^2
+16n^{3/2}+n.
\tag{25}
\]

Thus the leading factor four is uniformly achievable up to
O(n^2/sqrt(log n)); the bipolar clique obstruction shows that the leading
factor cannot be reduced below four on all dense, quadratic-cap seeds.
This does **not** prove 4Q(A)+O(n), nor any useful transfer at the
n^(3/2) cap scale.

### Elementary approximation and balancing proof

Define the rectangular cut norm by

\[
\|E\|_\square=\max_{U,V\subseteq[n]}
\left|\sum_{i\in U,j\in V}E_{ij}\right|.
\]

There exists a vertex partition into at most
4^{ceil(epsilon^{-2})} cells whose symmetric block-average matrix P obeys
||A-P||_square<=epsilon n^2. For completeness, start with one cell. If
the residual has a violating rectangle U times V, refine every cell by
membership in U and V, increasing the cell count by at most four. In the
normalized Frobenius inner product n^{-2} sum X_ij Y_ij, the old residual
is orthogonal to the old block-constant space. The new space contains
the rectangle indicator, so its projection captures squared norm at
least epsilon^2. These orthogonal energy increments can occur fewer than
epsilon^{-2} times, since ||A||_F^2/n^2<=1. This proves the stated bound
without an external regularity theorem.

Put E=A-P. Decomposing two sign vectors into their positive and negative
coordinate sets gives beta(E)<=4||E||_square. By independent coordinate
rounding this bilinear bound also holds on the real cube. Hence

\[
\sup_{z\in[-1,1]^n}|z^TEz/2|\le2\epsilon n^2,
\qquad Q_{\rm cube}(P)\le Q(A)+2\epsilon n^2.
\]

Delete one vertex from each odd partition cell, at most k vertices in all,
and balance S inside every remaining cell. Because P is constant on each
cell pair, the aggregate ordinary and S-weighted spin vectors obey the
diamonds of (21). Consequently the uncompleted parent on the undeleted
vertices has cap at most 4Q_cube(P). Restoring the deleted vertices costs
at most 4kn, as its parent has at most 2k restored vertices and every
entry of P has modulus at most one. This estimate also safely accounts
for any diagonal/cross-matching entries of P; opposite parent diagonals
have zero total contribution to its quadratic energy.

The ideal parent matrices constructed from A and P differ by
[[E,SES],[SES,-E]]. On Boolean pairs its energy difference is bounded by
2Q_cube(E)+beta(E)<=8epsilon n^2. Finally, completing the actual hollow-A
parent matching costs at most n. Combining these inequalities proves
(24). The switching is chosen using A's block approximation and is an
allowed native switching; no clone enlargement is performed in this
all-seed theorem.

This result identifies exactly where the positive diamond mechanism
loses the desired scale: approximating an arbitrary n-vertex signing by
large twin cells currently costs n^2/sqrt(log n), not o(n^(3/2)). It is
a uniform positive theorem with a quantified coarse remainder, not a
renaming of the full joint optimization condition.

## A strict low-cap improvement of the diamond constant

This section uses an already archived input, not a new discrepancy theorem.
The normalization-correct translation of Bollobas--Scott is proved in
`artifacts/one_sided_energy_product.md` and
`artifacts/macroscopic_closure_block_dichotomy.md`, Section3. The primary
source is [Discrepancy in graphs and hypergraphs, Theorem1](https://people.maths.ox.ac.uk/~scott/Papers/disc.pdf). In the present
half-energy convention, for a full signing E on m vertices with
Q(E)=O(m^(3/2)), it gives

\[
\ell(E):=\min\{\max H_E,-\min H_E\}
\ge(1-o(1))\frac{m^3}{12800Q(E)}.
\tag{26}
\]

Write W(A)=(max H_A-min H_A)/2, q=Q(A), and let Delta(A) denote the
left-hand supremum in (21). For every fixed K and every full signing with
q<=K n^(3/2), the diamond theorem improves uniformly to

\[
\boxed{\Delta(A)\le4W(A)
 -(1-o_K(1))\frac{n^3}{51200q}.}
\tag{27}
\]

In particular Delta(A)<= [4-1/(51200K^2)+o_K(1)]q. This is a strict but
very small improvement of four, not the desired 2 sqrt(2).

**Proof.** Since A is hollow, the diamond objective is affine in each
coordinate pair (a_i,b_i), and in each (c_i,d_i), with all other variables
fixed. A maximizer therefore exists at diamond vertices. After switching
b nonnegative, b_i is zero or one and h_i is an indicator. In the clipping
identity (22), v=b-h and a+h,a-h have complementary supports J and I.
The latter two vectors have the same support I. Write
U=max H_A and V=-min H_A, and similarly U_I,V_I for principal blocks.
The identity gives the sharper bound

\[
\Delta(A)\le U+V+U_I+V_J
\le2(U+V)-U_J-V_I
\le4W(A)-\ell(A[I])-\ell(A[J]).
\]

The middle inequality follows from the exact one-sided superadditivity
U>=U_I+U_J and V>=V_I+V_J: choose the two internal maximizers and then
their relative global sign so the cross term has the favorable sign.
The same reasoning applies to -A and hence controls the absolute value.

For both |I|,|J|>=n/4, the archived bound (26), principal monotonicity,
and q<=K n^(3/2) give, uniformly,

\[
\ell(A[I])+\ell(A[J])
\ge(1-o_K(1))\frac{|I|^3+|J|^3}{12800q}
\ge(1-o_K(1))\frac{n^3}{51200q}.
\]

If one block has size less than n/4, apply (26) only to the other,
whose cubed size exceeds 27n^3/64>n^3/4. This proves the same bound.

A finite version is also immediate. Put m0=floor(n/4)>=2 and
rho=2q/[m0(m0-1)]. If (1-rho^2)/4>=1/m0, the exact archived product
inequalities imply

\[
\Delta(A)\le4W(A)-(1-\rho^2)\frac{n^3}{51200q}.
\tag{28}
\]

The proof uses the exact lower bound
ell(A[S])>=(1-rho^2)|S|^3/(12800q) for |S|>=m0. No claim is made when
the displayed finite hypothesis fails.

For comparison, direct bilinear polarization has two disjoint complementary
supports and yields the slightly stronger invariant consequence

\[
\beta(A)\le4W(A)
 -(1-o_K(1))\frac{n^3}{25600q}.
\tag{29}
\]

Both consequences are genuine low-cap joint-energy improvements. They do
not supply a native twist of A: turning A itself into the diamond resource
still requires actual twins or an approximation whose error is controlled
at the relevant scale. The coarse approximation (25) has a much larger
error than the saving in (27).

### A precise limitation of the bilinear-slack condition

Even beta(A)<2 sqrt(2)Q(A), together with an exactly centered energy
interval, is not enough by itself to give Delta(A)<=2 sqrt(2)Q(A) in the
weighted class. In the equal-proportion bipolar model with t=3/4,

\[
\frac{\beta}{Q}=\frac{64}{25}=2.56,
\qquad
\frac{\Delta}{Q}=\frac{72}{25}=2.88>2\sqrt2.
\]

These exact values follow from (17)--(19). This is a quadratic-scale
weighted/clone-resource example; it is not a near-minimizer counterexample
and not a lower bound for optimization over every twist.

### Independent checks

`computations/twisted_chiral_diamond_2026_09_19.py` checks (22) on 4,400
rational-vector cases using integer arithmetic, checks (21) against exact
Boolean caps on 1,400 hollow cases, and checks 24 completed full-signing
twin-clone witnesses by exhaustive parent-spin enumeration. Every check
passes. Results are in
`computations/results/twisted_chiral_diamond_2026_09_19.json`.
The adversarial track independently audited the six-vector identity and
the clone-corollary constants; its audit is recorded in
`artifacts/twisted_chiral_symmetry_followup_2026_09_19.md`, Section8.

## A sharp bilinear comparison for the diamond norm

There is a substantially stronger universal estimate in terms of beta:

\[
\boxed{\Delta(A)\le\frac43\beta(A)
\quad\text{for every real hollow symmetric }A.}
\tag{30}
\]

This theorem is proved by an exact rational certificate below, not by a
floating-point optimization. Its constant is sharp over hollow real
forms, including asymptotic full-signing realizations at quadratic scale.
No claim of literature novelty is made.

### Exact eight-type bilinear certificate

By hollow-coordinate affinity, maximize the diamond objective at vertices.
Switch coordinates by sigma_i=a_i+b_i in {+-1}. Then a_i+b_i=1 and the
four entries at coordinate i have one of the eight types

\[
(a_i,b_i,c_i,d_i)=(e,1-e,fr,(1-f)r),\qquad
(e,f,r)\in\{0,1\}^2\times\{-1,1\}.
\]

Order these types lexicographically, with -1 before +1. On this eight-type
space define

\[
K=aa^T-cc^T+bd^T+db^T.
\]

The following table gives positive integers w and two eight-coordinate
sign vectors X,Y. The exact integer identity is

\[
\sum_r w_r(X_rY_r^T+Y_rX_r^T)=24K,
\qquad\sum_r w_r=32.
\tag{31}
\]

| w | X | Y |
|---:|---|---|
| 2 | -------+ | +-+-+--- |
| 2 | ------+- | +--++--- |
| 4 | -----+++ | +-+++-++ |
| 4 | ----+-+- | +---+--+ |
| 2 | ----++-+ | +-+-+-++ |
| 2 | ----+++- | +--++-++ |
| 1 | ---+-+-+ | +++-+-+- |
| 2 | ---++--+ | --+-+-+- |
| 4 | ---+++-+ | +-+-+++- |
| 4 | --+----- | +--+---- |
| 1 | --+--++- | ++-++--+ |
| 2 | -+------ | +-++---- |
| 1 | -+++---- | -+++---- |
| 1 | -+++++++ | -+++++++ |

Lift each X_r,Y_r to n coordinates by looking up the coordinate's type.
Arbitrary type multiplicities are allowed, including zero. Symmetry of A
and (31) yield exactly

\[
H_A(a)-H_A(c)+b^TAd
=\sum_r\frac{w_r}{24}X_r^TAY_r
\le\frac{32}{24}\beta(A).
\]

Apply the same argument to -A for the absolute value. This proves (30).
In finite-dimensional norm language, (31) is a projective bilinear
decomposition of the symmetric eight-type coefficient kernel; the native
norm on A is its real cube-to-cube bilinear norm. It does not identify
that bilinear norm with the Boolean quadratic norm.

The verifier
`computations/twisted_chiral_diamond_bilinear_2026_09_19.py` uses only integer
matrix arithmetic, contains no optimizer, and checks all 64 coefficient
identities as well as 1,900 lifted random identities. The adversarial track
independently checked the same 64 identities by plain scalar loops and
audited the type reduction and normalization. The discovery LP is not
needed for the proof.

### Sharpness and a nonvacuous sufficient condition

Consider the symmetric six-type integer matrix

\[
M=\begin{pmatrix}
-3&-1&-1&-1&-2&2\\
-1&3&2&2&-1&1\\
-1&2&-2&1&0&0\\
-1&2&1&-2&0&0\\
-2&-1&0&0&2&1\\
2&1&0&0&1&2
\end{pmatrix}.
\]

Exhaustive integer enumeration of its 64-by-64 bilinear spin pairs gives
beta(M)=18. Its Boolean quadratic extrema are +/-9, and beta/2 also bounds
the entire real-cube quadratic norm, so that cube norm is exactly 9.
The diamond vectors

\[
a=(0,0,0,0,1,1),\quad b=(1,1,1,1,0,0),
\]
\[
c=(0,0,-1,1,0,0),\quad d=(-1,1,0,0,-1,1)
\]

give energy 24. To obtain hollow examples, replace M by
M tensor J_r and remove its diagonal. The replicated diamond value is
24r^2-10r, while beta differs from 18r^2 by at most 14r. Thus the ratio
tends to 4/3. Scaling M by 1/3 and replacing its bounded block values by
full signs with those means gives full-signing sharpness as well: an
elementary union bound over bilinear sign pairs gives residual beta
O(N^(3/2)), negligible compared with the N^2 leading values. These remain
quadratic-cap examples, not near-minimizers.

The exact, independently verifiable seed condition

\[
\beta(A)\le\frac{3\sqrt2}{2}Q(A)
\]

therefore suffices for Delta(A)<=2 sqrt(2)Q(A). It is a genuine smaller
invariant condition, not the full paired landscape. For example beta=2Q
gives the strictly better diamond factor 8/3, and the sharpness family
above has exactly that ratio asymptotically.

### Scaling caveat for the actual cloned double

For the full-signing two-clone child C of order 2n constructed in (23),
the balanced parent of order 4n now satisfies

\[
Q(\text{parent})\le\frac{16}{3}\beta(A)+6n.
\tag{32}
\]

In parent-size normalized units this is
(2/3) beta(A)/n^(3/2)+o(1). Since beta(A)>=2Q(A), this certificate alone
cannot preserve the normalized cap of the original seed A. In particular,
the bilinear-slack condition solves the diamond norm comparison, not the
missing native twist selection or the convergence problem.

The coarse all-seed theorem (24) also improves in its beta formulation:

\[
\min_S Q(D)\le\tfrac43\beta(A)
+\tfrac{40}{3}\epsilon n^2+4kn+\tfrac{10}{3}n,
\qquad k\le4^{\lceil\epsilon^{-2}\rceil}.
\tag{33}
\]

Indeed beta(P)<=beta(A)+4epsilon n^2, and the double approximation costs
8epsilon n^2 as before. One must also account for the diagonal of the
block average P: the diamond vertex proof assumes a hollow matrix.
Writing P0=P-diag(P), the triangle bounds give

```math
 \Delta(P)\le\Delta(P0)+\sum_i|P_{ii}|
 \le\tfrac43\beta(P)+\tfrac73\sum_i|P_{ii}|
 \le\tfrac43\beta(P)+\tfrac73 n.
```

Here the diamond norm of a diagonal matrix is its entrywise diagonal
mass, since each scalar diamond contribution has absolute value at most
one. This is the extra 7n/3 correction in (33); a version with only the
last matching cost n silently applied the hollow theorem to P. Its
optimized remainder is still
O(n^2/sqrt(log n)), not O(n). Sharpness of (30) does not assert sharpness
of the optimized native constant in (33).

### Exact finite resource/native comparisons

The diamond norm itself has the exact formula

\[
\Delta(A)=\max_{b\in\{0,\pm1\}^n,\,J\subseteq[n]}
\left\{U(A[\operatorname{supp}(b)^c])+V(A[J])
+\sum_{i\notin J}|(Ab)_i|\right\}.
\tag{34}
\]

For fixed b, the a-vector occupies its complementary support and can
independently maximize its positive internal energy. For fixed support J
of c, optimize its negative internal energy, while each d_i outside J
independently takes the sign of (Ab)_i. This proves (34). Precomputing the
one-sided principal profiles from 3^n partial spins leaves a 6^n integer
enumeration; it does not search the native twist family.

`computations/twisted_chiral_diamond_profile_2026_09_19.py` performs this
enumeration and reconstructs an attaining four-vector witness for every
record. Among the results:

| seed | Q | beta | Delta | native completed family minimum |
|---|---:|---:|---:|---:|
| order6 minimizer | 5 | 12 | 12 | 18 |
| order7 class0 | 9 | 18 | 21 | 21 |
| order7 class1 | 9 | 18 | 20 | 21 |
| order7 class2 | 9 | 18 | 18 | 25 |
| order8 class0 | 10 | 24 | 26 | 30 |
| order8 class1 | 10 | 24 | 24 | 30 |
| stored order9 minimizer | 12 | 28 | 30 | 33 |
| stored order10 minimizer | 13 | 40 | 40 | 44 |
| order10 conference | 15 | 30 | 30 | 40 |

The last column is the separately audited computation track's full family
result, not a conclusion from (34). In particular, among the three order7
classes the best diamond resource has the worst native family minimum.
The independently exhausted order6 **core** family has minimum18, whereas
(4/3)beta=16. Thus the proposed zero-error purification
min_switch Q(D0)<=4beta/3 is already false at order6, even if permutations
are allowed. This finite example does not refute an O(n) remainder.

There is also a leading-order obstruction to lossless native purification.
For every A, the operator estimate on the diamond domain gives

\[
\Delta(A)\le n\|A\|.
\]

Indeed the two diamond constraints imply
||a||_2^2+||b||_2^2<=n and ||c||_2^2+||d||_2^2<=n. Bound the two internal
quadratics by half the operator norm times their squared lengths, and the
bridge by ||A||||b||||d||; arithmetic-geometric mean finishes the estimate.
On conference seeds this is n sqrt(n-1). Therefore a hypothetical native
bound F(A)<=Delta(A)+O(n) would force parent normalized caps at most
1/(2 sqrt(2))+o(1), contradicting the preserved universal lower bound
L>0.43332. Resource-to-native conversion must lose at leading order on
some conference sequence. This argument does not rule out a larger
conversion constant or a selected restricted-seed theorem.

## 15. Native 4/3-beta audit: two necessary nontrivial selection effects

The resource theorem (30) does **not** presently prove

```math
 F(A):=\min_{g,d}Q(D(A,g^TAg+\operatorname{diag}d))
       \le {4\over3}\beta(A)+O(n).                         \tag{36}
```

The exact order-six example already rules out the zero-error version,
but does not settle (36). Neither of the following results is a
counterexample to (36); both rule out concrete proposed selection shortcuts.

### Identity twisting fails at the leading competitive scale

Let H_n be Sylvester of order n=4^k and A_n=H_n-diag(H_n). Bent signs give
beta(H_n)=n^(3/2), and diagonal removal changes beta by at most n. The
identity-twist zero-matching core differs from H_2 tensor H_n only on
O(n) entries of bounded magnitude. The primary Schmidt theorem and the
same-spin embedding audited in
`fresh_schmidt_odd_walsh_regularization_2026_09_05.md` imply

```math
 Q(D_0(A_n,A_n))=(\sqrt2-o(1))n^{3/2},
 \qquad \beta(A_n)=(1+o(1))n^{3/2}.                       \tag{37}
```

Here the same-spin consequence is important: a statement only about
bilinear Walsh excess would not prove (37). Any diagonal matching changes
the cap by at most n. Hence (36), if true, must choose genuinely nontrivial
twists even on a flat Sylvester family. No lower bound for all twists is
deduced from this identity-twist failure.

### A bilinear-ground-pair switching rule can lose a factor two

For m>=2 set

```math
 A_m=\begin{pmatrix}J_m-I_m&J_m\\J_m&-J_m+I_m\end{pmatrix}.
```

Then, exactly,

```math
 \beta(A_m)=2m^2,\qquad Q(D_0(A_m,A_m))=4m^2-2m.          \tag{38}
```

Moreover, every pair attaining |x^T A_m y|=beta(A_m) has coordinate
product x*y constant. Thus the rule S=diag(x*y), even with free choice
among all bilinear maximizing pairs, supplies only B=A_m and is trapped
at `(2-1/m) beta(A_m)`. Diagonal matchings save at most 2m, so this rule
does not achieve (36) even at the easier m^2 scale.

To prove the assertions, write u and v for the two block sums of y.
The exact optimized bilinear value for this y is

```math
 \sum_{i\le m}|u+v-y_i|+\sum_{i>m}|u-v+y_i|.              \tag{39}
```

When |u|>|v|, this is `(2m-1)|u|+v sign(u)<=2m^2-2`;
when |v|>|u| it is `(2m-1)|v|-u sign(v)<=2m^2-2`.
When |u|=|v|>0 it is `(2m-1)|u|+m`, reaching 2m^2 only
at |u|=|v|=m. At u=v=0 it is 2m<2m^2. Thus the only
maximizing y are the four macroblock-constant corners. All their fields
are nonzero, and their uniquely aligned x equals either y or -y.
Allowing a negative bilinear maximum simply negates x.

For the second equality in (38), x=(1,1) and y=(1,-1), with blocks of
length m, attain 4m^2-2m. The eigenvalues on the block-constant subspace
are +-sqrt(2m^2-2m+1), and all remaining eigenvalues are +-1. Therefore

```math
 Q(D_0)\le2m\sqrt{2(2m^2-2m+1)}<4m^2-2m+1.
```

The cap is an integer, proving equality. In particular the elementary
universal untwisted bound `Q(D_0)<=2 beta(A)` is asymptotically sharp.
That upper follows from `|H(x)-H(y)|<=2W<=beta` and
`|x^TAy|<=beta`.

This is **not** an obstruction to other twists: the previously proved
two-macrotype theorem explicitly supplies a balanced switch with leading
cap `(5/2)m^2`, versus beta=2m^2, for this t=1 family.

The exact program
`computations/twisted_chiral_native_selection_2026_09_19.py` checks all
root-gauged seeds through order six and all 2^(2m) spins of this family
for 2<=m<=8. On an optimal order-six seed the bilinear-ground-product
rule gives core22 for every eligible product, whereas unrestricted
switching gives18 (beta=12). Thus the failure is visible even before
the explicit quadratic-scale counterfamily. The JSON retains matrices,
switches, enumeration counts, and the integer spectral certificate.
An independent audit by the adversarial track checked the block-sum
classification, the attaining pair, and the strict integer spectral bound;
the squared gap in the last bound is exactly (2m-1)^2.

### Scope of the remaining native claim

The coarse regularity result above already gives
`F(A)<=(4/3)beta(A)+o(n^2)` uniformly. Consequently a fixed finite
block model cannot disprove (36) by a leading n^2 gap. Conversely,
the global parent lower bound L>0.433322 would disprove (36) if one
could produce an infinite hollow full-sign symmetric family with

```math
 \limsup {\beta(A_n)\over n^{3/2}}
 <{3\over4}\,2\sqrt2\,(0.4333221116640807)
 =0.9192150107\ldots.
```

This is only a sufficient falsifier criterion, not an asserted family.
Finite low-beta seeds do not supply it: bilinear norm is not
multiplicative under tensoring, and unrestricted permutations may mix
tensor fibers. The archived primary-source Gale--Berlekamp bounds and
Hadamard regularization audits do not provide such a family. No
native 4/3-beta selection theorem or asymptotic counterexample has
been obtained here, and the bound is not promoted as a justified
next construction without an additional effective selection mechanism.

## Current conclusion

The exact objective, a computable joint-selection lemma, and the convergence
accounting are proved. No uniform low-cap/minimizer selection theorem has been
proved. In particular a successful O(n)-defect twisted doubling estimate would
be a meaningful step but would still need an order-coverage supplement such as
(8) for a second multiplier or the selected thinning statement (9).
