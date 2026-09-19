# Independent adversarial audit of twisted chiral doubles

Campaign: 2026-09-18, 17:09:11--20:09:11 UTC. This is an independent
audit in the original normalization
\(H_A(x)=\sum_{i<j}a_{ij}x_ix_j\), \(Q(A)=\max_x|H_A(x)|\).
Statements about an explicitly checked matrix, a twist family, and the
minimum over all signings are kept separate.

## 1. Exact identity and its scope

Let \(A,B\) be any symmetric hollow real matrices and put
\[
D=\begin{pmatrix}A&B+\operatorname{diag}d\\
B+\operatorname{diag}d&-A\end{pmatrix}.
\]
No switch/permutation hypothesis is needed for the following identity.
For \(x,y\in\{\pm1\}^n\), set \(z=x\),
\(I=\{i:y_i=x_i\}\), \(J=I^c\). With
\(C_A(I,J;z)=\sum_{i\in I,j\in J}a_{ij}z_iz_j\), direct expansion gives
\[
H_D(x,y)=2C_A(I,J;z)+2H_{B[I]}(z_I)-2H_{B[J]}(z_J)
                 +d(I)-d(J).
\]
Negating every \(z_i\) in \(I\) negates the first summand and fixes
the second. Consequently
\[
\boxed{Q(D)=\max_{I,z}\left(2|C_A(I,I^c;z)|+
 \left|2H_{B[I]}(z_I)-2H_{B[I^c]}(z_{I^c})+d(I)-d(I^c)\right|\right).}
\]
The maximizing partition and spin are shared by both terms. Separate
upper bounds on their individual maxima need not capture their interaction.

The signed permutation \(J_0=\left(\begin{smallmatrix}0&-I\\I&0\end{smallmatrix}\right)\)
satisfies \(J_0^2=-I\), \(J_0^TDJ_0=-D\). Thus every such double is
chiral and its positive and negative extrema agree in magnitude.
This centering statement by itself gives no target upper bound.

## 2. Exact weighted-shell rank obstruction

Let \(W=(w_{ij})\) be hollow symmetric with \(|w_{ij}|\le1\), and
choose \(\sigma\in\{\pm1\}\), \(x_*\in\{\pm1\}^n\) such that
\(\sigma H_W(x_*)=Q(W)\). Write \(x_*^T\) for the spin obtained by
flipping the coordinates in \(T\). Then
\[
0\le Q(W)-\sigma H_W(x_*^T)
=2\sigma\sum_{i\in T,j\notin T}w_{ij}x_{*,i}x_{*,j}
\le2|T|(n-|T|).
\]
For \(n\ge3\), the signed nearmaximum shell of width \(4(n-2)\)
contains the unflipped spin and every single and double flip. A width
\(4n\) works without small-order qualifications. The absolute shell
\(\{|H_W|\ge Q(W)-\eta\}\) also contains those spins whenever
\(\eta\ge4n\).

Let \(\phi(x)=(x_ix_j)_{i<j}\). For every edge \(ij\),
\[
\phi(x_*)-\phi(x_*^{\{i\}})-\phi(x_*^{\{j\}})
+\phi(x_*^{\{i,j\}})=4x_{*,i}x_{*,j}e_{ij}.
\]
The shell feature vectors therefore span the entire
\(\binom n2\)-dimensional edge space. In particular, for every fixed
\(\delta>0\), an exact linear perturbation preserving all energies
on a shell of width \(\delta n^{3/2}\) must vanish for all sufficiently
large \(n\); the sufficient threshold \(n\ge16/\delta^2\) follows
from the convenient width \(4n\).

Allowing a common energy shift does not evade the result: the same
mixed difference annihilates a constant. Hence a quadratic edge
perturbation constant on this shell has every coefficient zero and
its constant is zero.

### Quantitative statement, and why it is weaker asymptotically

If \(V=(v_{ij})\) obeys
\[
|H_V(x)-c|\le\varepsilon
\quad\text{on the unflipped, single-flip, and double-flip spins},
\]
the four-term difference gives exactly
\[
\boxed{|v_{ij}|\le\varepsilon\quad(i<j).}
\]
Thus rounding an edge with distance at least \(\rho\) from
\(\{\pm1\}\) requires \(\varepsilon\ge\rho\). This is an exact
finite-tolerance obstruction, **not** an obstruction to a tolerated
\(o(n^{3/2})\) error.

Nor does it bound how many order-one coefficients can change once
\(\varepsilon\) is order one. After switching \(x_*\) to all plus,
any hollow symmetric \(V\) whose row sums vanish satisfies
\[
H_V(x_*)=H_V(x_*^{\{i\}})=0,\qquad
H_V(x_*^{\{i,j\}})=4v_{ij}.
\]
Dense examples exist: for even \(m\), choose balanced sign vectors
\(u,v\in\{\pm1\}^m\) and take
\(V=\left(\begin{smallmatrix}0&uv^T\\vu^T&0\end{smallmatrix}\right)\).
This changes \(m^2\) edge coordinates by one while keeping all those
local energies within four. Controlling the *entire* wider shell might
give stronger inequalities, but the local rank identity alone does not.

## 3. Escape from a protected shell

Even when a narrower shell has a nonzero annihilator, preserving it
does not control the new cap. An explicit bounded-coefficient example
has order four. Let \(w_{12}=w_{13}=w_{14}=1/10\), all other
\(w_{ij}=0\). Then \(Q(W)=3/10\). Let
\(v_{23}=1,v_{24}=-1\), all other \(v_{ij}=0\).
Both positive and negative old extremizers have all three leaf spins
equal, so \(H_V=0\) there. For every absolute shell width
\(0\le\eta<1/5\), these are all the protected states. Nevertheless
at \(x=(1,1,1,-1)\),
\[
H_W(x)=1/10,\qquad H_V(x)=2,\qquad H_{W+V}(x)=21/10>3/10.
\]
Both \(W\) and \(W+V\) have coefficients in \([-1,1]\).

For a shell \(S_\eta=\{x:|H_W(x)|\ge Q(W)-\eta\}\), the sufficient
conditions
\[
\sup_{S_\eta}|H_V|\le\varepsilon,\qquad
\sup_{S_\eta^c}|H_V|\le\eta+\varepsilon
\]
ensure \(Q(W+V)\le Q(W)+\varepsilon\). The second condition, or some
equally strong one-sided margin condition, must not be omitted.
For example \(n\|V\|_{\rm op}/2\le\eta+\varepsilon\) certifies
the needed global perturbation bound, but may be too costly in the
proposed application.

### Full rounding with an \(O(n)\)-protected shell can still escape

There is a stronger asymptotic example. Let \(n=3m\) and partition the
vertices into three equal blocks. Let \(A\) be the full signing with
negative within-block edges and positive cross-block edges. Every row
sum is \(m+1\). Put
\[
\rho=\frac{m+1}{3m-1},\qquad W=\rho(J-I),\qquad V=A-W.
\]
Then every coefficient of \(W\) is genuinely fractional, \(W+V\)
rounds **all** edges to signs, and \(V\) has row sums zero. At a spin
obtained by flipping \(T\) from all plus,
\[
H_V(x^T)=4\sum_{i<j\in T}v_{ij},\qquad
|H_V(x^T)|\le2(1+\rho)|T|^2.
\]

Fix \(\delta>0\). For sufficiently large \(n\), the absolute
nearmaximum shell of \(W\) with width \(\delta n^{3/2}\) consists
only of positive-energy states close to one of the two constant
spins, since the most negative \(W\)-energy has magnitude at most
\(\rho n/2\), whereas \(Q(W)=\Theta(n^2)\). If \(k\le n/2\)
is the distance to the nearer constant spin, then
\[
2\rho k(n-k)\le\delta n^{3/2},\qquad
k\le\frac{\delta\sqrt n}{\rho}.
\]
Because \(\rho\ge1/3\), the **entire** protected shell satisfies
\[
\boxed{\sup_{\{|H_W|\ge Q(W)-\delta n^{3/2}\}}|H_V|
\le24\delta^2 n.}
\]
Yet
\[
Q(W)=\frac{3m^2+3m}{2},\qquad
Q(A)=\frac{5m^2-3m}{2}\quad(m\ge3),
\]
so the cap increases by \(m^2-3m\), of order \(n^2\).
The negative extremizer has block spins \((+,+,-)\). To verify the
exact cap, let \(M_1,M_2,M_3\) be the three block magnetizations;
\[
H_A=\frac{(M_1+M_2+M_3)^2-2\sum M_i^2+3m}{2}.
\]
The quadratic part is concave in each individual \(M_i\), so its
minimum occurs at a cube vertex; evaluating the eight vertices gives
the stated minimum. Its maximum is at most
\((\sum M_i^2+3m)/2\le(3m^2+3m)/2\), attained at all plus.

Thus the exact shell-rank obstruction cannot be promoted to a claim
that approximate whole-shell protection with \(O(n)\) error prevents
extensive rounding. It also cannot substitute for an escape bound.
This example has a high-cap fractional seed \(Q(W)=\Theta(n^2)\);
it makes no assertion about approximate shell geometry under an
additional \(Q(W)=O(n^{3/2})\) hypothesis.

## 4. Scope of the prior chiral exclusions

`artifacts/chiral_scale_preserving_lift.md` uses
\(Q_{\rm old}(A)=\max|x^TAx|=2Q(A)\). Its order-12 seed has
\(Q_{\rm old}=40\), hence original \(Q=20\). Its order-48 equality
target \(320\) means original \(Q=160\), and its rejected candidate
separator \(440\) means original energy \(220\).

The completely excluded class is the explicitly parameterized
complementary-support inherited-chiral Clifford **four-lift** family
for that fixed order-12 seed: 1008 masks, 912 rejected by a uniform-fibre
profile, and 96 remaining masks in four microcoordinate-permutation
orbits excluded with exact subset-sum certificates. The unrestricted
compressed four-lift problem remains unresolved in that artifact.
An initial symmetry-restricted 1507-constraint candidate was rejected;
failure of the subsequent time-bounded solve is not an UNSAT proof.

The order-two chiral example excludes a theorem promising exact
scale-preserving four-lifts for **every** chiral seed. It does not
exclude an additive-defect statement, favorable minimizer selection,
or a theorem restricted to growing orders.

The present twist \(B=P^TD_sAD_sP\) generally does not have the
old block-compression constraints, and it is an order-two multiplier.
None of the old finite-family UNSAT conclusions transfers to all
twisted doubles without a proved embedding into that excluded family.

## 5. A twist-invariant obstruction: the Boolean bilinear norm

Define
\[
\beta(A)=\max_{x,y\in\{\pm1\}^n}|x^TAy|
        =\max_y\sum_i|(Ay)_i|.
\]
For symmetric \(C\), swapping \(x\) and \(y\) negates
\(H_A(x)-H_A(y)\) but fixes \(x^TCy\). Therefore the stronger identity
\[
Q\!\begin{pmatrix}A&C\\C&-A\end{pmatrix}
=\max_{x,y}\bigl(|H_A(x)-H_A(y)|+|x^TCy|\bigr)
\ge\beta(C)
\]
holds. Signed permutation congruence preserves \(\beta\), so every
allowed twisted double satisfies
\[
\boxed{Q(D)\ge\beta(B+\operatorname{diag}d)\ge\beta(A)-n.}
\]
This lower bound is insensitive to the choice of twist.

There is a still stronger invariant
\[
\gamma(A)=\min_{d\in\{\pm1\}^n}\beta(A+\operatorname{diag}d),
\qquad \boxed{Q(D)\ge\gamma(A)}.
\]
Indeed, under signed permutation congruence a diagonal sign matrix
is simply permuted. These statements do not claim that either lower
bound is achievable by a twisted double.

### Deterministic full-signing counterexample to an arbitrary-seed theorem

Let \(m\) be a power of two and \(R_m\) the Sylvester Hadamard matrix,
so \(R_mR_m^T=mI\). Define the actual hollow full signing of order \(n=2m\)
\[
A_m=\begin{pmatrix}J_m-I_m&R_m\\R_m^T&-J_m+I_m\end{pmatrix}.
\]
For \(u,v\in\{\pm1\}^m\),
\[
H_{A_m}(u,v)=\tfrac12\bigl((\sum u_i)^2-(\sum v_i)^2\bigr)+u^TR_mv.
\]
Consequently
\[
Q(A_m)\le m^2/2+m^{3/2}.
\]
The explicit bilinear pair
\(x=(\mathbf1,\mathbf1)\), \(y=(\mathbf1,-\mathbf1)\) has
\[
x^TA_my=2m(m-1),
\]
because the two crossblock terms cancel exactly. Hence **every**
switch/permutation twist and **every** diagonal filling obeys
\[
\boxed{Q(D_m)\ge2m^2-4m.}
\]
In particular,
\[
\liminf_{m\to\infty}\frac{\min_{P,s,d}Q(D_m)}{Q(A_m)}\ge4>2\sqrt2.
\]
The gap from \(2\sqrt2 Q(A_m)\) is of order \(m^2\), not merely
\(m^{3/2}\). Thus a uniform bound
\(Q(D)\le2\sqrt2 Q(A)+O(n^{3/2-\delta})\), or even one with
an \(o(n^2)\) defect, is false when quantified over **all** full
signings.

This is a high-cap seed family: \(Q(A_m)=\Theta(m^2)\). It does not
disprove a theorem restricted to minimizing or uniformly
\(O(n^{3/2})\)-cap seeds. The lower bound \(\beta(A)\le4Q(A)\)
from polarization shows that the ratio-four bilinear obstruction is
asymptotically maximal: write \(u=(x+y)/2\), \(v=(x-y)/2\), use
\(x^TAy=2(H_A(u)-H_A(v))\), and use multilinearity to bound each
\(|H_A(u)|,|H_A(v)|\le Q(A)\) on \([-1,1]^n\).

At \(m=64\), the entirely explicit order-128 seed has certified
\(Q(A)\le2560\) and every order-256 twisted double has
\(Q(D)\ge7936\), whereas \(2\sqrt2\cdot2560<7241\).
No exhaustive order-128 cap computation or global optimum claim is
needed for this certificate.

### The obstruction can persist at bounded normalized cap

The following strengthening was proposed by the campaign director and
independently reconstructed here. It shows that merely adding the
assumption \(Q(A)=O(N^{3/2})\), with an arbitrary fixed implied constant,
does not rescue the arbitrary-seed theorem.

Let \(N\) run through powers of two, take a symmetric Sylvester
Hadamard \(H_N\), and start with the full signing
\(A_0=H_N-\operatorname{diag}(H_N)\). Fix \(K>0\), put
\(m=\lfloor K N^{3/4}\rfloor\), and for sufficiently large \(N\)
choose two disjoint vertex sets of size \(m\). Replace their internal
edges by a positive clique and a negative clique, leaving all other
edges unchanged. The resulting \(A\) remains an actual full signing.

Write \(R\) for the block diagonal matrix formed by the two old
principal blocks of \(A_0\), and \(C\) for their clique replacements.
Since \(\|A_0\|_{\rm op}\le\sqrt N+1\),
\[
Q(A_0)\le\tfrac12N^{3/2}+\tfrac12N,\qquad
Q(R)\le m(\sqrt N+1).
\]
The two new blocks have equal size, so their diagonal constants cancel
and \(Q(C)=m^2/2\) when \(m\) is even, with the upper bound
\(Q(C)\le m^2/2\) valid for every \(m\). Therefore
\[
Q(A)\le\tfrac12N^{3/2}+\tfrac12N+\tfrac12m^2+m(\sqrt N+1)
=\left(\tfrac12+\tfrac12K^2+o(1)\right)N^{3/2}.
\]
The Boolean bilinear norm of a principal submatrix never exceeds
that of the whole matrix: extend its two Boolean test vectors by zero,
then maximize the resulting bilinear form over the full cube. On the
two patched blocks the same explicit pair as above cancels the cross
terms, regardless of what those cross edges are. Hence
\[
\beta(A)\ge2m(m-1),\qquad
\min_{P,s,d}Q(D)\ge2m(m-1)-N
=(2K^2-o(1))N^{3/2}.
\]
Whenever \(K^2>1+\sqrt2\), the lower bound exceeds
\(2\sqrt2 Q(A)+c_KN^{3/2}\) for some \(c_K>0\) and all sufficiently
large \(N\). Taking \(K=2\) gives
\[
Q(A)\le(5/2+o(1))N^{3/2},\qquad
\min_{P,s,d}Q(D)\ge(8-o(1))N^{3/2},
\]
whose ratio is at least \(16/5-o(1)>2\sqrt2\).

This **does not** settle the near-minimizer regime. In particular it
does not construct counterexamples with normalized cap near the
project's approximately \(0.43\)--\(0.50\) range. A theorem restricted
to exact minimizers or a sufficiently stringent low-cap class remains
unrefuted by this construction.

### Sharper padded obstruction using the unmodified complement

The constant can be improved while retaining the same construction.
Restrict to \(N=4^k\). For the standard symmetric Sylvester matrix,
\[
z=(1,1,1,-1)^{\otimes k}\in\{\pm1\}^N,\qquad H_Nz=\sqrt N\,z.
\]
Also \(\operatorname{tr}H_N=0\), so
\(Q(A_0)=N^{3/2}/2\) exactly. The previous upper bound sharpens to
\[
Q(A)\le\tfrac12N^{3/2}+\tfrac12m^2+m\sqrt N+m.
\tag{P1}
\]

The bilinear norm is superadditive over disjoint principal blocks.
To see this for a partition \(T\sqcup U\), take maximizing pairs
on the two blocks with nonnegative bilinear values and combine them.
Simultaneously negating both vectors on \(T\) keeps both internal
values fixed and negates the sum of cross terms. One of the two
choices makes that sum nonnegative. Thus
\[
\beta(A)\ge\beta(A[T])+\beta(A[U]).
\]

Here \(|T|=2m\) consists of the patched blocks, and
\(\beta(A[T])\ge2m(m-1)\). On \(U=T^c\), let \(w\) be \(z\)
with the coordinates in \(T\) replaced by zero. Then
\[
\begin{aligned}
w^TH_Nw
&=N^{3/2}-2\sqrt N|T|+z_T^TH_N[T]z_T\\
&\ge N^{3/2}-6m\sqrt N.
\end{aligned}
\]
The diagonal correction on \(U\) has absolute value at most \(2m\),
because the full trace is zero. Hence
\[
\boxed{\beta(A)\ge N^{3/2}+2m^2-6m\sqrt N-4m.}
\tag{P2}
\]
Combining (P1)--(P2) with \(Q(D)\ge\beta(A)-N\) gives the entirely
explicit finite lower bound
\[
\begin{aligned}
Q(D)-2\sqrt2 Q(A)\ge{}&
(1-\sqrt2)N^{3/2}+(2-\sqrt2)m^2\\
&-(6+2\sqrt2)m\sqrt N-(4+2\sqrt2)m-N.
\end{aligned}
\tag{P3}
\]
For \(m\sim K N^{3/4}\), its leading coefficient is positive exactly
when \(K^2>1/\sqrt2\). Thus no theorem of the proposed form can hold
on **all** signings with \(Q(A)\le C N^{3/2}\) for any fixed
\[
\boxed{C>\frac{1+1/\sqrt2}{2}=0.853553\ldots.}
\]
Indeed, choose \(K^2\) strictly between \(1/\sqrt2\) and \(2C-1\).
This is stronger than the first padded counterexample but still does
not reach the near-minimizer regime.

For completeness, the cap estimate is asymptotically exact. If \(m\)
is even, retain the bent spin on \(U\), use all plus on the positive
clique and a balanced spin on the negative clique. The two-clique
internal energy is \(m^2/2\), its mutual cross term has magnitude at
most \(m\sqrt N\), and a common sign flip on the patched coordinates
makes their interaction with \(U\) nonnegative. This gives
\[
Q(A)\ge\tfrac12N^{3/2}+\tfrac12m^2-4m\sqrt N-m.
\]
Choosing even \(m\sim K N^{3/4}\) therefore yields
\(Q(A)=(1/2+K^2/2+o(1))N^{3/2}\). Similarly the triangle inequality
\(\beta(A)\le\beta(A_0)+\beta(C)+\beta(R)\) gives
\(\beta(A)\le N^{3/2}+2m^2+2m\sqrt N+N\), so
\(\beta(A)=(1+2K^2+o(1))N^{3/2}\).

## 6. Independent exact finite checks

`computations/twisted_chiral_adversary_2026_09_18.py` checks 48 arbitrary
\(A,B,d\) identities at orders \(1\le n\le6\), checks the edge mixed
differences at \(2\le n\le14\), and exhausts the escape example.
The results and every generated input are retained in
`computations/results/twisted_chiral_adversary_2026_09_18.json`.

`computations/twisted_chiral_bilinear_audit_2026_09_18.py` independently
exhausts all projective Boolean spins of stored witnesses to compute
\(Q\) and \(\beta\). It also exhausts every diagonal sign vector to
compute \(\gamma\). This is a certificate for the stated finite
bilinear optimization, not for the minimum over all signing matrices.
For local fields \(a_i=y_i(Ay)_i\), it uses the exact integer identity
\[
|(Ay)_i+d_iy_i|=|a_i|+d_i\operatorname{sgn}(a_i)
                      +\mathbf1_{a_i=0}.
\]
The sum is evaluated using exactly representable integers. The same
script reconstructs and verifies the full-signing Hadamard obstruction.

For the stored `exact_m3.json` through `exact_m10.json` representatives:

| Order | Exact witness \(Q\) | \(\beta\) | \(\gamma\) |
| --- | --- | --- | --- |
| 3 | 3 | 6 | 5 |
| 4 | 4 | 8 | 8 |
| 5 | 4 | 8 | 11 |
| 6 | 5 | 12 | 14 |
| 7 | 9 | 18 | 17 |
| 8 | 10 | 24 | 20 |
| 9 | 12 | 28 | 27 |
| 10 | 13 | 40 | 30 |

All archived orbit representatives through order eight were also
checked: the other two order-seven classes have \(\gamma=19\),
and the other order-eight class has \(\gamma=22\). No complete
order-nine or order-ten orbit classification is claimed here.

The order-ten minimizing witness has \(\beta/Q=40/13>2\sqrt2\),
so the unfilled bilinear norm can exceed the target even on an actual
small minimizing example. But \(\gamma=30<2\sqrt2\cdot13\): diagonal
filling removes that finite obstruction. An additive \(O(n)\) theorem
is not contradicted. The older order-twelve chiral seed has
\((Q,\beta,\gamma)=(20,52,44)\).

Additional archived witness checks (without asserting any newly
certified global optimum) give:

| Order | Exact witness \(Q\) | \(\beta\) | \(\gamma\) |
| --- | --- | --- | --- |
| 11 | 17 | 38 | 35 |
| 12 | 18 | 40 | 36 |
| 13 | 20 | 40 | 49 |
| 14 | 21 | 46 | 52 |
| 15 | 27 | 58 | 61 |
| 16 | 30 | 68 | 64 |

The retained order-sixteen \(Q=32\) witness instead gives
\((\beta,\gamma)=(72,68)\). Source paths, exact spins attaining each
\(\beta\), minimizing diagonal vectors, and the number of diagonals
checked are retained in the bilinear audit JSON.

The discarded eight-type graphon optimization is retained as
`computations/twisted_chiral_graphon_adversary_2026_09_18.py`.
Its initial exact Boolean-type evaluation exposed the ratio-four
obstruction; the invariant above supersedes numerical minimax search.
It is not used as a certificate of a graphon optimum.

`computations/twisted_chiral_witness_recheck_2026_09_18.py` reconstructs
the submitted switched/permuted bridges and parents, then directly
enumerates \(x_0=+1\) and every \(y\), using
\(H_A(x)-H_A(y)+x^T(B+\operatorname{diag}d)y\). This is independent of
the cut-profile search implementation. The first nine submitted
witnesses, with seed orders seven through ten, passed reconstruction
and complete integer-energy histogram checks. Their reported caps
were respectively \(21,21,25,30,30,37,44,37,44\). These are witness
caps, not conclusions about the best twist or all-signing optimum.

### Independent profile-width certificate for the third order-seven class

For fixed \(B\), write the unfilled parent energy as
\(E_B(x,y)=H_A(x)-H_A(y)+x^TBy\). At a fixed relative spin
\(t_i=x_iy_i\), the diagonal filling adds only the scalar
\(d\cdot t\). Thus every filling satisfies
\[
Q(D)\ge W(A,B):=\max_t\frac{\max_{x\odot y=t}E_B(x,y)
                                    -\min_{x\odot y=t}E_B(x,y)}2.
\]
This lower bound survives even if each profile is allowed an
independently chosen *real* centering scalar, rather than a scalar
of the constrained form \(d\cdot t\).

The C++ search uses \(x_0=y_0=+1\), which meets every absolute-energy
orbit: global negation and the chiral map \((x,y)\mapsto(y,-x)\)
reduce every state to that subset. Its signed automorphism and
antiautomorphism quotient preserves the fixed-seed twist family;
an antiautomorphism is accompanied by exchange of the two parent
halves. Those facts justify the stated reduction, but were not relied
upon for the following independent full-orbit enumeration.

`computations/twisted_chiral_width_recheck_2026_09_18.py` separately
enumerates all 420 rooted permutation copies of archived order-seven
class 2 and all 64 switches of each. Direct integer matrix products,
without Gray recurrence or the automorphism quotient, reproduce
\[
\boxed{W(A,B)=24\quad\text{for all 26880 distinct allowed }B.}
\]
Every full order-fourteen signing has odd energies. Therefore every
twisted double of this seed has \(Q(D)\ge25\). The independently
reconstructed cap-25 witness attains this lower bound, certifying
the optimum **within this fixed-seed twisted family**. Neither the
global order-fourteen optimum nor an asymptotic exclusion follows.
The full enumeration digest and explicit high/low profile witnesses
are retained in
`computations/results/twisted_chiral_width_recheck_2026_09_18.json`.

## 7. Broader order-sixteen chiral feasibility relaxation

At the director's request, a bounded CP-SAT feasibility test asks
whether **any** chiral order-sixteen full signing has cap at most 28,
without requiring \(B\) to be a switched/permuted copy of \(A\).
Thus infeasibility would exclude every order-sixteen twisted double,
while feasibility need not produce a member of the twist family.

The gauge reduction is complete. On the \(i\)-th chiral coordinate
pair, a quarter rotation
\(R=\left(\begin{smallmatrix}0&-1\\1&0\end{smallmatrix}\right)\)
changes the internal matching sign \(d_i\) to \(-d_i\). Rotating
one endpoint of a cross-pair block replaces
\[
\begin{pmatrix}a&b\\b&-a\end{pmatrix}
\quad\text{by}\quad
\begin{pmatrix}b&-a\\-a&-b\end{pmatrix};
\]
rotating both negates the whole block. Full signs and chirality are
preserved. Independent pair rotations can therefore set every
\(d_i=+1\). Subsequent simultaneous sign switches on both members
of each pair preserve \(d\) and can set every \(a_{0i}=+1\).
There remain 21 free \(A\)-edges and 28 free \(B\)-edges, hence
49 Boolean variables.

Permuting the seven nonroot pairs can sort \((b_{01},\ldots,b_{07})\).
Exchanging the two parent halves sends \(A\mapsto-A\), fixes \(B,d\),
and, after restoring \(a_{0i}=+1\), complements all seven root
\(B\)-signs. Consequently at most three sorted root signs need be
negative. These are legitimate symmetries of the *relaxed* chiral
class; no assertion is made that pairwise quarter rotations preserve
the restricted relation \(B\sim A\).

Every absolute-energy orbit meets \(x_0=y_0=+1\), giving 16384
representatives. Their two-sided cap constraints reduce to 16320
distinct integer linear rows after exact duplicate merging. The
complete coefficient matrix and bounds are stored as a compressed
NPZ beside the run JSON. The model uses eight CP-SAT workers, seed
20260918, and a 600-second initial limit. The gauge-normalized
cap-30 witness is used only as a search hint, not as feasible input.

`computations/twisted_chiral_chiral16_model_audit_2026_09_18.py`
independently checks all 65536 parent spins for eight random gauge
transformations and verifies exact equality of their full energy
histograms. It also checks 24 assignments against the generated
linear systems at targets 28, 30, and 48, comparing feasibility with
direct parent quadratic-form evaluation. All checks pass.

The initial solve ended `UNKNOWN` after 600.088 solver seconds
(61724 conflicts, 92624 branches in the reported primary search
statistics). It returned neither a feasible cap-28 witness nor an
infeasibility result. The broader chiral question therefore remains
unresolved; the cap-30 witnesses do not establish its optimum.
The exact model, inputs, search seed, and response statistics are
preserved. No further generic CP-SAT continuation is inferred from
this bounded failure.

For evidence discipline, even a future `INFEASIBLE` CP-SAT response
would be reported as solver-certified without a standalone proof
object, and not as a lower bound over nonchiral order-sixteen signings.

## Resumed follow-up, 2026-09-19

The continuation is recorded in
`artifacts/twisted_chiral_symmetry_followup_2026_09_19.md`.
It proves an arbitrary-baseline low-cap padding theorem with explicit
subleading errors, independently audits the factor-four linear
cube-pullback certificate obstruction, and checks the completed
order-ten minimizing-class and fixed-child-family results. The separate
question of a stronger asymptotic lower bound from global chirality
remains unresolved there.
