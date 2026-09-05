# Fresh limit attack: algebraic scaling (2026-09-05)

This document starts from the problem statement only. No steering/state/ledger/archive material was read before the candidate below was frozen.

For a symmetric hollow sign matrix, write
\[
 P_A(x)=\sum_{i<j}a_{ij}x_ix_j,\qquad
 Q(A)=\max_{x\in\{-1,1\}^n}|P_A(x)|,
 \qquad M_n=\min_AQ(A).
\]

## Frozen candidate and exact convergence criterion

**Candidate amplification theorem.** There exist absolute constants
\(C<\infty\) and \(\delta>0\) such that, for each \(k\in\{2,3\}\)
and all sufficiently large integers \(n\),
\[
 M_{kn}\le k^{3/2}M_n+C n^{3/2-\delta}.
 \tag{A}
\]
It is enough that the theorem hold for one minimizer at each order; a universal theorem for all signed matrices is not asserted.

The reason to isolate two incommensurable integer multipliers is that a single dyadic construction is insufficient: it controls individual dyadic orbits but does not bridge all multiplicative phases. A proof of (A), unlike an exact identity at a fixed finite order, would imply convergence. A full proof of that implication, including the summable accumulated defect and the use of deletion, is given below.

### First algebraic construction and a warning about its constants

The most immediate candidate at multiplier 2 is
\[
 L(A)=\begin{pmatrix}A&A+I\\A+I&-A\end{pmatrix}.
\]
This is again a symmetric hollow sign matrix. For \(x,y\in\{-1,1\}^n\), put
\(u=(x+y)/2\), \(v=(x-y)/2\). Their supports are disjoint and cover all coordinates. Direct expansion gives
\[
 P_{L(A)}(x,y)
 =2\{P_A(u)-P_A(v)+u^TAv\}+x^Ty.
 \tag{B}
\]
The factor 2 in (B) is essential. Writing \(z=u+iv\), its first term is
\(2\sqrt2\,\Re(e^{-i\pi/4}P_A(z))\). The general degree-two complexification bound \(|P_A(z)|\le2Q(A)\) therefore gives only \(4\sqrt2Q(A)+n\), not the needed \(2\sqrt2Q(A)+O(n)\).

Moreover, a universal amplification estimate with leading factor \(2\sqrt2\) is false without a near-optimality hypothesis. For the hollow matrix \(A=K_s\oplus(-K_s)\), with zero cross entries, choosing \(u\) on the first component and \(v\) on the second gives lifted energy \(2s(s-1)\), whereas \(Q(A)=s^2/2+O(s)\). Filling the zero cross block with a sign matrix having bilinear norm \(O(s^{3/2})\) preserves this asymptotic ratio 4. Thus a universal all-sign-matrix version also fails.

The intended next step is to test the lift specifically on exact small-order minimizers, then decide whether its failure is merely a finite defect or a structural obstruction. No claim that (A) is true is being made at this checkpoint.

For completeness, (A) implies convergence as follows. Put \(f(n)=M_n/n^{3/2}\). Along any word in multipliers 2 and 3, the sum of normalized defects is at most
\[
 Dn^{-\delta},\qquad D=\frac{C2^{-3/2}}{1-2^{-\delta}}.
\]
The semigroup \(2^a3^b\), \(a,b\ge0\), has consecutive multiplicative gaps tending to 1. Indeed, for any \(\epsilon>0\), finitely many residues \(b\log3\pmod{\log2}\) form an \(\epsilon\)-net, and sufficiently large target logarithms can be approximated from above using these finitely many \(b\)'s and a nonnegative \(a\). Principal deletion gives \(M_m\le M_N\) for \(m\le N\): extend a maximizing configuration of the principal submatrix with independent fair spins and take expectations. Thus, for each sufficiently large fixed \(n\),
\[
 \limsup_{m\to\infty}f(m)\le f(n)+Dn^{-\delta}.
\]
Taking \(n\) through a liminf sequence proves convergence. This criterion overlaps existing dense-good-order arguments and is not the main result below.

## 1. A genuine convergent tensor-regularized surrogate

The subsequent construction does **not** assume that the original sequence converges, nor that a tensor lift preserves every seed.

### 1.1 Fixed outer matrices

Let \(H_4=J_4-2I_4\), so \(H_4^2=4I\) and \(H_4\mathbf1=2\mathbf1\). We also require one symmetric regular Hadamard matrix of order 144. Here is a direct construction.

Start with any order-12 Hadamard matrix \(F\). One explicit choice is the Paley matrix over \(\mathbb F_{11}\): take \(C_{ij}=\chi(i-j)\), with \(\chi(0)=0\), and
\[
 F=\begin{pmatrix}1&\mathbf1^T\\-\mathbf1&I+C\end{pmatrix}.
\]
The identities \(C^T=-C\), \(C\mathbf1=0\), and \(CC^T=11I-J\) directly give \(FF^T=12I\).

Index a 144-by-144 matrix by pairs and set
\[
 K_{(i,j),(p,q)}=F_{iq}F_{pj},\qquad
 D_{(i,j),(i,j)}=F_{ij},\qquad H_{144}=DKD.
\]
Then \(K^T=K\), \(K^2=144I\), and \(K\operatorname{vec}(F)=12\operatorname{vec}(F)\); hence \(H_{144}\) is symmetric, has sign entries, and has every row sum 12.

For \(a,b\ge0\), let
\[
 s=4^a144^b,\qquad H_s=H_4^{\otimes a}\otimes H_{144}^{\otimes b}.
\]
These outer matrices satisfy \(H_s^2=sI\) and \(H_s\mathbf1=\sqrt s\mathbf1\). Their orders form a multiplicatively dense semigroup at infinity, because \(\log144/\log4\) is irrational. The proof is the same finite-residue argument as above.

### 1.2 Definition and basic properties

For any real symmetric hollow matrix \(A\), including matrices with some zero off-diagonal entries, define
\[
 \mathcal R(A)=\sup_{a,b\ge0}
 \frac{Q(H_{4^a144^b}\otimes A)}{(4^a144^b)^{3/2}}.
 \tag{R}
\]

The numerator uses the same-spin **absolute quadratic** objective; it is not a bilinear norm. No diagonal-fiber sign completion is included in (R).

The normalized quantities in (R) are coordinatewise nondecreasing in \((a,b)\). To see this, tensor any witness with an all-one vector in an additional outer factor. The energy then multiplies by that factor's order to the power \(3/2\). Also,
\[
 Q(A)\le\mathcal R(A)\le\frac n2\|A\|_{\mathrm{op}}<\infty.
 \tag{R1}
\]
Thus the supremum is also the cofinal limit as both indices tend to infinity. It follows by shifting the indices that, for every outer order \(s\),
\[
 \mathcal R(H_s\otimes A)=s^{3/2}\mathcal R(A).
 \tag{R2}
\]
The functional is subadditive and positively homogeneous, as a supremum of the corresponding normalized seminorms. It is monotone under principal deletion: tensor the smaller principal matrix with each \(H_s\), then use the elementary principal-deletion inequality for \(Q\).

### 1.3 Exact completion payment

For a hollow sign matrix \(A\) of order \(n\), the core \(H_s\otimes A\) has zero off-diagonal entries exactly when the two inner coordinates agree. Complete these entries using
\[
 E_s=(H_s-\operatorname{diag}H_s)\otimes I_n,
 \qquad A'=H_s\otimes A+E_s.
\]
This \(A'\) is a symmetric hollow sign matrix of order \(sn\). Since
\(\|E_s\|_{\mathrm{op}}\le\sqrt s+1\), (R1) and (R2) give
\[
 \frac{\mathcal R(A')}{(sn)^{3/2}}
 \le\frac{\mathcal R(A)}{n^{3/2}}
       +\frac1{2\sqrt n}+\frac1{2\sqrt{sn}}.
 \tag{R3}
\]
This is an upper bound with the completion paid absolutely. It does not assume cancellation between the core and filling.

### 1.4 The stabilized minimum converges

Set
\[
 r_n=\min_{A\text{ hollow symmetric sign, order }n}
             \frac{\mathcal R(A)}{n^{3/2}}.
\]
The minimum is attained because the set of such \(A\)'s is finite. For each fixed \(n\), combine (R3), the dense outer-order mesh, and principal deletion to obtain
\[
 \limsup_{m\to\infty}r_m\le r_n+\frac1{2\sqrt n}.
 \tag{R4}
\]
Here one chooses an outer order \(s\) with \(sn\ge m\) and \(sn/m\to1\); the last term of (R3) tends to zero. Taking \(n\) through a liminf sequence proves the unconditional theorem
\[
 \boxed{\lim_{n\to\infty}r_n=c_{\mathcal R}\text{ exists}.}
 \tag{R5}
\]

Since \(Q\le\mathcal R\), this constant only gives
\[
 \limsup_n\frac{M_n}{n^{3/2}}\le c_{\mathcal R}.
 \tag{R6}
\]
It does **not** identify the original liminf. At each **fixed seed order** \(n\), the set of sign matrices is finite and the normalized tensor caps are coordinatewise nondecreasing, so the finite minimum actually does commute with the directed supremum:
\[
 r_n=\sup_{a,b\ge0}\ \min_{A\text{ of order }n}
          \frac{Q(H_{4^a144^b}\otimes A)}{(4^a144^b n)^{3/2}}.
\]
For example, to prove this interchange, choose for each of the finitely many seeds an outer index within \(\epsilon\) of its own supremum, then take one coordinatewise common upper index. The unproved comparison concerns **unbounded growing seed orders** and unrestricted parent minimizers versus tensor-restricted families, not this valid fixed-order identity. Even convergence of the original objective to a smaller constant would not contradict the surrogate theorem.

## 2. The absolute PSD-majorant relaxation is exactly tensor-stable, but collapses to 1/2

For a real symmetric matrix \(A\), define
\[
 \mathcal T(A)=\frac12\min_{B\succeq A,\ B\succeq-A}
                       \max_{x\in\{-1,1\}^n}x^TBx.
 \tag{T}
\]
The minimum exists: feasible \(B\)'s are positive semidefinite, their objective bounds \(\operatorname{tr}B\), and a finite objective sublevel is compact. The spectral choice \(B=\|A\|_{\mathrm{op}}I\) is feasible.

### 2.1 Tensor stability and comparison

For every symmetric regular Hadamard matrix \(H_s\),
\[
 \mathcal T(H_s\otimes A)=s^{3/2}\mathcal T(A).
 \tag{T1}
\]
For the upper bound, if \(B\succeq\pm A\), then
\(\sqrt s I_s\otimes B\succeq\pm H_s\otimes A\), as seen by decomposing \(H_s/\sqrt s\) into its positive and negative spectral projections. The maximum quadratic value of this majorant is \(s^{3/2}\max_xx^TBx\).

For the lower bound, let \(D\succeq\pm H_s\otimes A\), and put \(V=\mathbf1_s\otimes I_n\). Then
\[
 B=s^{-3/2}V^TDV\succeq\pm A,
\]
because \(\mathbf1^TH_s\mathbf1=s^{3/2}\). Also \(Vx\) is Boolean whenever \(x\) is Boolean, so
\[
 \max_z z^TDz\ge s^{3/2}\max_xx^TBx\ge2s^{3/2}\mathcal T(A).
\]

Clearly \(Q(A)\le\mathcal T(A)\). Applying this after tensoring and using (T1) gives
\[
 \boxed{Q(A)\le\mathcal R(A)\le\mathcal T(A).}
 \tag{T2}
\]

### 2.2 Universal half-scale floor for \(\mathcal T\)

**Theorem.** Every symmetric hollow sign matrix of order \(n\ge2\) satisfies
\[
 \boxed{\mathcal T(A)\ge\frac12n\sqrt{n-1}.}
 \tag{T3}
\]

Proof. Let \(B\succeq\pm A\) and \(C=\max_xx^TBx\). Write \(a_i=Ae_i\). Each vector \(a_i+e_i\) and \(a_i-e_i\) is Boolean, because \(A\) is hollow with all off-diagonal entries of magnitude one. Averaging their quadratic values and summing over \(i\) gives
\[
 nC\ge\operatorname{tr}B(A^2+I).
 \tag{T4}
\]
If \(\lambda_1,\ldots,\lambda_n\) are the eigenvalues of \(A\), then \(v_j^TBv_j\ge|\lambda_j|\) for an eigenbasis \(v_j\). Consequently
\[
 nC\ge S_3+S_1,
 \qquad S_p=\sum_j|\lambda_j|^p.
 \tag{T5}
\]
Here \(S_2=n(n-1)\), \(S_1\le n\sqrt{n-1}\), and Cauchy--Schwarz gives \(S_1S_3\ge S_2^2\). The function \(s+S_2^2/s\) is decreasing for \(0<s\le n\sqrt{n-1}\), since \(n\ge2\). Hence
\[
 S_3+S_1\ge S_1+\frac{S_2^2}{S_1}
 \ge n\sqrt{n-1}+n(n-1)^{3/2}
 =n^2\sqrt{n-1}.
\]
This proves \(C\ge n\sqrt{n-1}\), and minimizing over \(B\) proves (T3).

The bound is attained whenever \(A^2=(n-1)I\), by \(B=\sqrt{n-1}I\). More generally, the dense symmetric Hadamard orders constructed above, with their diagonals removed and then principal deletion, prove
\[
 \boxed{
 \lim_{n\to\infty}\min_A\frac{\mathcal T(A)}{n^{3/2}}=\frac12.
 }
 \tag{T6}
\]
Indeed, if \(m\le s\) and \(s/m\to1\), a principal order-\(m\) submatrix of \(H_s-\operatorname{diag}H_s\) has operator norm at most \(\sqrt s+1\), so its \(\mathcal T\)-value is at most \(m(\sqrt s+1)/2\).

In particular, \(c_{\mathcal R}\le1/2\), but (T3) gives no lower bound on \(\mathcal R\), since the inequality in (T2) goes in the opposite direction. Proving a half-scale floor for the true \(\mathcal R\) remains open in this attack.

### 2.3 What this rules out

Any argument that replaces the genuine tensor-regularized objective by the majorant functional (T), even while optimizing the majorant and the growing seed jointly, has asymptotic minimum exactly \(1/2\). It therefore cannot transport a putative original liminf below \(1/2\) with vanishing loss. This does **not** rule out an exact Boolean tensor argument strictly stronger than (T), nor a completion exploiting the seed.

## 3. Rigorous lower bound for the true regularization

Write \(\mathcal B(A)=\max_{f,g\in\{-1,1\}^n}|f^TAg|\). The matrix \(H_4=J-2I\) has orthogonal Boolean eigenvectors
\(u_+=\mathbf1\) and \(u_-=(1,1,-1,-1)\), with eigenvalues \(+2\) and \(-2\). For arbitrary Boolean \(f,g\), set
\[
 z_i=\frac12\{(f_i+g_i)u_+ +(f_i-g_i)u_-\}\in\{-1,1\}^4.
\]
A direct calculation gives
\[
 \frac{z_i^TH_4z_j}{8}=\frac{f_ig_j+g_if_j}{2},
 \qquad
 Q(H_4\otimes A)\ge4\mathcal B(A).
\]
Therefore
\[
 \mathcal R(A)\ge\frac12\mathcal B(A).
 \tag{R7}
\]
For a hollow sign matrix, a uniform random \(g\) gives
\[
 \mathcal B(A)\ge\mathbb E\|Ag\|_1
 =n\,\mathbb E|\varepsilon_1+\cdots+\varepsilon_{n-1}|.
\]
The exact binomial expression, or the elementary central limit asymptotic, yields
\[
 \frac1{\sqrt{2\pi}}\le c_{\mathcal R}\le\frac12.
 \tag{R8}
\]
This is an interval for the surrogate, not a new lower bound on the original sequence.

## 4. Exact seed audits and their scope

The independent script `computations/fresh_limit_algebra_lift.py` enumerates all first-row-normalized order-six signs. It finds \(M_6=5\), with 12 normalized minimizers, and the completed Sylvester lift from (B) has cap 18 for all of them.

After that candidate was frozen, the named archive note `artifacts/cross_order_selectable_hadamard_replica_no_go.md` was inspected. It supplied the order-five representative
\[
 A_5=\begin{pmatrix}
0&-1&1&-1&1\\
-1&0&-1&1&1\\
1&-1&0&1&-1\\
-1&1&1&0&-1\\
1&1&-1&-1&0
\end{pmatrix},\qquad A_5^2=5I-J.
\]
For this seed, \(|A_5|=\sqrt5(I-J/5)\). Thus, for **every** symmetric Hadamard outer matrix of order \(s\),
\[
 |H_s\otimes A_5|=\sqrt s I_s\otimes|A_5|,
\]
and applying the absolute-value majorant to both sides of the quadratic form proves the improved uniform core bound
\[
 Q(H_s\otimes A_5)\le\frac{12}{\sqrt5}s^{3/2}.
 \tag{S1}
\]
Every Boolean vector in five coordinates has odd coordinate sum, so
\(\max_xx^T|A_5|x=\sqrt5(5-1/5)=24/\sqrt5\). In fact the lower-bound proof (T4)--(T5) is exact here and gives
\[
 \mathcal T(A_5)=\frac{12}{\sqrt5},\qquad
 \frac{\mathcal R(A_5)}{5^{3/2}}\le\frac{12}{25}=0.48.
 \tag{S2}
\]
This does not improve the original asymptotic upper bound: (R3) still pays \(1/(2\sqrt5)\) for completion. A finite seed below \(1/2\) is not enough; the sum of its regularized ratio and its completion payment must be below \(1/2\), or a growing-seed gap / sharper completion theorem is needed.

The independent Gray-code enumerator `computations/fresh_limit_algebra_core.cpp` checks all configurations, modulo global reversal, for the regular outer \(H_4=J_4-2I_4\). It gives

| Seed | \(Q(A)\) | \(Q(H_4\otimes A)\) | preserving target \(8Q(A)\) |
|---|---:|---:|---:|
| \(A_5\) above | 4 | 36 | 32 |
| order-six conference minimizer, normalized code 220 | 5 | 48 | 40 |

Both objective signs were exhaustively checked: the energy ranges are \([-36,36]\) and \([-48,48]\). The regular outer improves the order-five cap from the archived standard Sylvester value 40 to 36, but does not preserve the seed coefficient. This is direct evidence that selecting the outer matrix matters, without resolving growing-seed regularization.

For exploratory lower witnesses only, `computations/fresh_limit_algebra_anneal.cpp` uses the fixed RNG seed 912734 and simulated annealing. On the regular \(H_4\)-powers it finds:

| Seed | Outer order | Certified witness energy | Divided by \(s^{3/2}\) |
|---|---:|---:|---:|
| \(A_5\) | 16 | 324 | 5.0625 |
| \(A_5\) | 64 | 2612 | 5.1015625 |
| \(A_6\) | 16 | 404 | 6.3125 |
| \(A_6\) | 64 | 3220 | 6.2890625 |

These are not exact caps. In the last row the weaker normalized value is simply an optimization failure: the order-16 witness can always be amplified to normalized value 6.3125 at order 64. No empirical convergence of these searches is asserted. All displayed exact and heuristic spin witnesses are saved in `computations/results/fresh_limit_algebra_witnesses.json`; their energies were recomputed independently in JavaScript from the seed and outer definitions before saving.

## 5. Full-sign seeds remove the persistent completion payment

This strengthening was suggested by the root agent after the hollow-seed theorem above was written. It changes the finite-seed criterion substantially.

Extend the definition (R) to **all** real symmetric matrices, using
\(Q(B)=\tfrac12\max_x|x^TBx|\), including the diagonal contribution. The comparisons, subadditivity, exact tensor stability, and spectral bound remain valid. In fact this extension is a norm: (R7) remains valid for diagonal terms as well, so \(\mathcal R(B)\ge\mathcal B(B)/2\), and \(\mathcal B\) is a norm. For a diagonal matrix,
\[
 \mathcal R(\operatorname{diag}d)=\frac12\sum_i|d_i|.
 \tag{F1}
\]
The lower bound follows from (R7); the upper bound follows by splitting into one-coordinate matrices and using the spectral bound for each one.

Let \(B\) now be an order-\(n\) symmetric **full sign** matrix: diagonal entries are also \(\pm1\). Then \(H_s\otimes B\) has no missing entries. Removing only its parent diagonal gives a hollow sign matrix
\[
 A_s=H_s\otimes B-\operatorname{diag}(H_s\otimes B)
\]
with
\[
 \frac{\mathcal R(A_s)}{(sn)^{3/2}}
 \le\frac{\mathcal R(B)}{n^{3/2}}+\frac1{2\sqrt{sn}}.
 \tag{F2}
\]
Thus the entire completion payment now tends to zero for **each fixed full seed**.

Define
\[
 c_* = \inf_{n\ge1}\ \min_{B\in\{-1,1\}^{n\times n},\ B=B^T}
                   \frac{\mathcal R(B)}{n^{3/2}}.
 \tag{F3}
\]
Dense outer orders and hollow principal deletion applied to (F2) prove \(\limsup_mr_m\le c_*\). Conversely, for each hollow seed \(A\), add any sign diagonal \(D\). By (F1),
\[
 c_*\le\frac{\mathcal R(A+D)}{n^{3/2}}
 \le\frac{\mathcal R(A)}{n^{3/2}}+\frac1{2\sqrt n}.
\]
Minimize over \(A\) and let \(n\to\infty\). Consequently the earlier limit has the stronger finite-seed characterization
\[
 \boxed{c_{\mathcal R}=c_*.
 }
 \tag{F4}
\]
The normalized minima over full sign seeds also converge to this number, since their difference from the hollow minima is at most \(1/(2\sqrt n)\). In particular,
\[
 \boxed{\limsup_m M_m/m^{3/2}
       \le \mathcal R(B)/n^{3/2}\quad
       \text{for every fixed full symmetric sign seed }B.}
 \tag{F5}
\]

So a single certified full seed with \(\mathcal R(B)<n^{3/2}/2\) **would** improve the original asymptotic upper bound. This does not contradict the earlier warning about the hollow order-five seed: its missing diagonal fibers are a positive fraction of each large core, whereas a full seed has only the final parent diagonal to remove.

For clarity, ordinary quadratic principal monotonicity fails on general symmetric matrices with diagonals (e.g. \(\operatorname{diag}(1,-1)\)). The proof above avoids that issue by deleting the parent diagonal first. If desired, \(\mathcal R\) itself is principal-monotone on the full domain: tensor a principal-deletion expectation inequality with \(H_s\); its trace correction is at most \(s\sum_{i\notin I}|B_{ii}|/2\), which vanishes after division by \(s^{3/2}\) and taking a cofinal outer limit.

### 5.1 The full-sign PSD relaxation has an exact half floor

If \(B\) is a full symmetric sign matrix and \(D\succeq\pm B\), set \(C=\max_xx^TDx\). Every column \(Be_i\) is Boolean, so
\[
 nC\ge\operatorname{tr}(DB^2)
     \ge\sum_j|\lambda_j(B)|^3
     \ge n^{5/2}.
\]
The last step is the power-mean inequality with \(\sum_j\lambda_j(B)^2=n^2\). Therefore
\[
 \boxed{\mathcal T(B)\ge\tfrac12n^{3/2}
          \quad\text{for every full symmetric sign matrix }B.}
 \tag{F6}
\]
Equality is attained by symmetric Hadamard matrices. Hence replacing \(\mathcal R\) by \(\mathcal T\) in the finite-seed characterization forces the value exactly \(1/2\).

### 5.2 A concrete structural realization problem

On the uniform four-point space, \(U=H_4/2=2\mathbb E-I\). Its tensor power acts by \((-1)^{|S|}\) on the Efron--Stein component supported on coordinate set \(S\). For a Boolean vector map \(f:\{1,2,3,4\}^d\to\{-1,1\}^n\), put
\[
 G=\mathbb E(ff^T),\qquad
 K=\sum_{|S|\text{ even}}\mathbb E(f_Sf_S^T)
   -\sum_{|S|\text{ odd}}\mathbb E(f_Sf_S^T).
\]
Then \(G\) belongs to the classical cut-correlation polytope and \(-G\preceq K\preceq G\). The normalized same-spin tensor objective is \(\tfrac12|\operatorname{tr}(AK)|\). Convex duality gives the PSD-majorant target body
\[
 \mathcal T(A)=\frac12\max\{\operatorname{tr}(AK):
     G\in\operatorname{conv}\{xx^T:x\in\{-1,1\}^n\},\ -G\preceq K\preceq G\}.
 \tag{F7}
\]
For an explicit duality derivation, write the primal as minimizing \(C/2\) subject to \(C\ge x^TDx\) for every Boolean \(x\), and \(D\succeq\pm A\). Its dual variables are weights \(p_x\ge0\) with \(\sum_xp_x=1\), and PSD matrices \(P,Q\) with \(P+Q=\sum_xp_xxx^T\). Set \(G=P+Q\), \(K=P-Q\). The objective becomes \(\operatorname{tr}(AK)/2\), and the PSD constraints become \(-G\preceq K\preceq G\). Strict primal feasibility gives strong duality.

Thus proving \(\mathcal R_4=\mathcal T\) would require an actual parity-covariance realization theorem, or equality of support functions after taking the closure/convex hull of the attainable covariance pairs. It is not justified merely by the PSD inequalities. Conversely, a separating quadratic functional would give a genuine tensor upper certificate stronger than \(\mathcal T\).

The simplest full-sign candidate is
\(B=H_2=\begin{psmallmatrix}1&1\\1&-1\end{psmallmatrix}\):
\(Q(B)=1\) and \(\mathcal T(B)=\sqrt2\). A uniform strict upper bound \(\mathcal R(B)<\sqrt2\), if proved, would improve the original upper coefficient. No such bound is established here.

For this full seed, exact enumeration gives \(Q(H_4\otimes H_2)=10\); fixed-seed heuristic searches give lower witnesses 80 at outer order 16 and 676 at outer order 64. These are normalized seed values 1.25, 1.25, and 1.3203125, still below the relaxation ceiling \(\sqrt2\). Only the first cap is exact. The full spin witnesses and independent energy checks are saved in `computations/results/fresh_limit_algebra_h2_witnesses.json`.

The optional numerical screening script `computations/fresh_limit_algebra_full_seeds.py` checks all 32 or 64 sign diagonals added to the exact order-five or order-six seed. It finds numerical minima \(\min_D\mathcal T(A_5+D)\approx6.6114484\) and \(\min_D\mathcal T(A_6+D)\approx8.0138769\), above the respective half floors \(5.5901699\) and \(7.3484692\). These are solver estimates, not certified bounds or a no-go theorem for the true regularization.

## 6. Remaining precise questions

1. Does the genuine regularization satisfy \(\mathcal R(A)\ge\tfrac12n^{3/2}-o(n^{3/2})\) uniformly over hollow sign matrices? A proof would show \(c_{\mathcal R}=1/2\), but would not settle the original problem.
2. Is \(\mathcal R(A)=\mathcal T(A)\) for all real symmetric hollow \(A\), or at least asymptotically at the minimizing sign matrices? This would be a nontrivial Boolean-Hadamard realization theorem; it is not an ordinary same-spin SDP or bilinear completely-bounded-norm identity.
3. Can original minimizers along an order sequence attaining the original liminf be selected with \(\mathcal R(A_n)-Q(A_n)=o(n^{3/2})\)? Together with (R5)--(R6), this **would** prove convergence of the original sequence. No such selection theorem is proved here.

The deliverable is therefore an unconditional convergent tensor-stable surrogate and an exact half-scale collapse theorem for its natural PSD relaxation, not a proof of convergence or nonconvergence of \(M_n/n^{3/2}\).

## Reproduction

Run `.venv/bin/python -B computations/fresh_limit_algebra_verify.py` to check the explicit order-144 outer and every saved exact/heuristic witness by integer arithmetic. The C++ sources compile into the workspace's `tmp` directory; binaries are not part of the deliverable. No commits were made.
