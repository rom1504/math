# Phase 2: true tensor regularization and the full H2 seed

This continues `fresh_limit_algebra_2026_09_05.md`. The original convergence problem remains open in this work. Here the target is a genuine upper certificate below the PSD-majorant relaxation, or a structural theorem explaining why that relaxation is unavoidable.

## 1. The same-spin regularization equals a regularized bilinear norm

Write \(\mathcal B(A)=\max_{x,y\in\{-1,1\}^n}|x^TAy|\). The elementary inequalities
\[
 Q(A)\le\tfrac12\mathcal B(A),\qquad
 \mathcal R(A)\ge\tfrac12\mathcal B(A)
\]
have different left sides. The second is the exact one-step \(H_4=J-2I\) construction from the first note. Applying it to every outer tensor and using exact \(\mathcal R\)-stability proves
\[
 \boxed{\mathcal R(A)=\frac12\sup_{s=4^a144^b}
       \frac{\mathcal B(H_s\otimes A)}{s^{3/2}}.}
 \tag{1}
\]
This identity was also independently derived by the variational agent. It does not replace the original objective with an unregularized bilinear norm.

## 2. Explicit congruence with odd-dimensional Walsh matrices

Let \(A=H_4\otimes H_2\), with \(H_2=\begin{psmallmatrix}1&1\\1&-1\end{psmallmatrix}\), and let \(W_8(i,j)=(-1)^{\operatorname{popcount}(i\mathbin\&j)}\). With
\[
 p=(0,3,5,6,7,4,2,1),\qquad
 d=(1,-1,-1,-1,-1,-1,-1,1),
\]
direct integer checking gives
\[
 W_8(i,j)=-d_id_j A(p_i,p_j).
 \tag{2}
\]
Tensor induction therefore shows that \(H_4^{\otimes a}\otimes H_2\) is a signed-permutation/global-sign congruence of \(W_{2^{2a+1}}\). The finite certificate is saved in `computations/results/fresh_limit_algebra_h2_r4_congruence.json`.

This only identifies the \(H_4\)-only branch. The order-144 outer allowed in \(\mathcal R\) is an additional family, so an upper certificate restricted to Walsh powers would not bound the full regularization.

## 3. The known self-dual-bent recursion can stall permanently

For a symmetric Walsh matrix \(W_N\), let \(f\) be Boolean and set \(g=\operatorname{sign}(W_Nf)\), choosing \(+1\) at zero. The recursion
\[
 F=(f,g,g,-f)
\]
under \(W_4\otimes W_N\) satisfies
\[
 F^T(W_4\otimes W_N)F=8\|W_Nf\|_1
 \ge8|f^TW_Nf|.
 \tag{3}
\]
This is exactly the monotonicity mechanism in Theorem 3.4 of Carlet--Danielsen--Parker--Sol\'e, *Self-Dual Bent Functions* (2010); the primary manuscript is author-uploaded at [ResearchGate](https://www.researchgate.net/publication/210230164_Self-dual_bent_functions), DOI [10.1504/IJICOT.2010.032864](https://doi.org/10.1504/IJICOT.2010.032864). The paper does not claim that the limiting normalized value is one.

There is a precise obstruction to extracting strict progress from (3). If \(\operatorname{sign}(W_Nf)=f\) and no coordinate of \(W_Nf\) is zero, then \(g=f\), and
\[
 (W_4\otimes W_N)(f,f,f,-f)
   =2(W_Nf,W_Nf,W_Nf,-W_Nf).
\]
Hence the new Boolean vector is again its own sign best response, has exactly the same normalized Rayleigh quotient, and the same statement repeats indefinitely.

For an explicit non-spectral fixed point, take
\[
 f=(1,-1,-1,1,-1,1,1,1).
\]
Then
\[
 W_8f=(2,-2,-2,2,-2,2,2,6),\qquad
 f^TW_8f=20,
\]
so all aligned fields are strictly positive, but the normalized Rayleigh quotient is
\(20/(8\sqrt8)=5/(4\sqrt2)<1\). The recursion stalls forever on this starting vector. Thus monotonicity plus the orthogonality of the outer matrices alone cannot prove a half-floor or \(\mathcal R(H_2)=\sqrt2\).

The variational agent has additionally certified the exact next finite plateau \(Q(H_4^{\otimes2}\otimes H_2)=80\), so the first two normalized regularized seed values are both 1.25. This stronger finite statement is proved in that agent's profile certificate, not inferred from heuristic failure.

## 4. Ordinary degree-two SDP is already exactly spectral on the H2 family

There is also an exact limitation on a natural upper-certificate class. Let \(H_s\) be a symmetric Hadamard matrix with constant diagonal \(\delta\in\{-1,1\}\), and assume both spectral signs occur. Put
\[
 U=H_s/\sqrt s,\qquad V=H_2/\sqrt2,\qquad
 \tau=\delta/\sqrt s.
\]
The matrix
\[
 G=\frac{I\otimes I+U\otimes V-\tau(U\otimes I+I\otimes V)}{1-\tau^2}
 \tag{4}
\]
is positive semidefinite, has every diagonal entry equal to one, and is supported in the positive eigenspace of \(U\otimes V\). Indeed, writing \(P_\pm\) and \(Q_\pm\) for the spectral projections of \(U\) and \(V\), respectively,
\[
 G=\frac{2}{1+\tau}P_+\otimes Q_+
       +\frac{2}{1-\tau}P_-\otimes Q_-.
\]
The diagonal of \(P_\pm\) is \((1\pm\tau)/2\), and the diagonals of \(Q_++Q_-\) sum to one, which proves the diagonal claim. Consequently the ordinary same-spin correlation SDP attains
\[
 \frac12\operatorname{tr}[(H_s\otimes H_2)G]
   =\frac12(2s)^{3/2}
\]
exactly. A certificate below that value must therefore use genuinely higher-degree Boolean constraints, not merely a better diagonal majorant or the ordinary degree-two elliptope.

## 5. A concrete higher-moment obstruction target

By (1) and the identity \(\|H_2v\|_1=2\|v\|_\infty\),
\[
 \mathcal R(H_2)=\sup_s\sup_{f,g\text{ Boolean}}
          \mathbb E\max(|U_sf|,|U_sg|).
 \tag{5}
\]
The upper bound \(\sqrt2\) follows from \(\mathbb E(U_sf)^2=\mathbb E(U_sg)^2=1\). Near equality requires both:

- the two transformed profiles have almost disjoint support;
- their squared magnitudes sum to almost the constant 2.

This has an exact nonquadratic defect formulation. Let
\[
 a=\max(|u|,|v|),\qquad b=\min(|u|,|v|),\qquad
 \Delta(u,v)=\mathbb E\{b^2+(a-\sqrt2)^2\}.
\]
Since \(\mathbb E(a^2+b^2)=2\),
\[
 \boxed{\mathbb E a=\sqrt2-\frac{\Delta(u,v)}{2\sqrt2}.}
 \tag{6}
\]
Thus a uniform positive lower bound on this explicit defect, for Boolean input pairs and every allowed outer, would be exactly a strict upper certificate for the full seed. The elementary lattice obstruction rules out exact zero at finite square outer orders, but its normalized defect can tend to zero. No dimension-independent lower bound is proved here.

There is also a useful phase-cube formulation. Put \(h=(f+ig)/\sqrt2\), so every input coordinate belongs to \(e^{i\pi/4}\{1,i,-1,-i\}\). Let
\[
 d_s^2=\min_{h\in e^{i\pi/4}\{1,i,-1,-i\}^s}
        \min_{w\in\{1,i,-1,-i\}^s}
            \mathbb E|U_sh-w|^2.
\]
Choosing the nearest coordinate phase gives
\[
 d_s^2=2-\sqrt2\sup_{f,g}\mathbb E\max(|U_sf|,|U_sg|),
 \qquad
 \boxed{\mathcal R(H_2)=\sqrt2-\frac1{\sqrt2}\inf_s d_s^2.}
 \tag{7}
\]
This is a concrete uniform phase-separation problem for two finite phase cubes under the fixed signed orthogonal outer transforms. A positive separation must hold for the order-144 factors as well as the Walsh branch.

For Walsh outers, the second condition can be expressed through autocorrelations. If \(u=W_Nf/\sqrt N\), then
\[
 \mathbb E u^4=\sum_t C_f(t)^2,\qquad
 \mathbb E(u^2v^2)=\sum_t C_f(t)C_g(t),
 \quad C_f(t)=\mathbb E_xf(x)f(x+t).
\]
These identities suggest a quartic/Boolean-convolution certificate. However, a positive fourth-moment defect alone does not imply a gap in (6), because a small exceptional set can have large transformed magnitudes; tail control or a bounded/truncated certificate is necessary. The nonquadratic defect in (6) avoids this hidden assumption. No uniform positive defect has been established.

## 6. Operator-theory scope audit

The norm in (1) is an endpoint vector-valued Fourier norm:
\(\sup_s\|A U_s\|_{L_\infty(\ell_\infty^n)\to L_1(\ell_1^n)}\), with uniform probability normalization externally. It is not the usual \(L_2\)-to-\(L_2\) Fourier-type norm. Hinrichs, *Hilbert space factorization and Fourier type of operators*, Studia Math. 145 (2001), shows nonuniform equivalence of finite Fourier/Walsh and Hilbert-factorization gradations; its counterexamples use other operator domains and do not settle this endpoint norm. See the [publisher's paper](https://www.impan.pl/shop/en/publication/transaction/download/product/89984). No Hilbert-factorization identity is imported from that result.

## 7. Quadratic quaternary phase constructions cannot improve the H2 baseline

Here is a rigorous obstruction to one apparently promising construction family. It is **not** an upper bound on unrestricted \(\mathcal R\).

Let \(r\) be even, \(N=2^r\), and consider the quaternary quadratic phase
\[
 z(x)=i^{q(x)},\qquad
 q(x)=c+\sum_{j=1}^r\alpha_jx_j
               +2\sum_{j<k}\beta_{jk}x_jx_k\pmod4,
 \tag{8}
\]
where \(c,\alpha_j\in\mathbb Z/4\mathbb Z\) and \(\beta_{jk}\in\mathbb F_2\). The even coefficient on every cross term is part of the hypothesis; allowing an odd cross coefficient is a genuinely larger, non-Clifford family.

Put \(h=e^{i\pi/4}z\), and \(U=W_N/\sqrt N\). Then
\[
 \boxed{\mathbb E\max_{w\in\{1,i,-1,-i\}}
               \Re(\overline w\,Uh)\le\frac1{\sqrt2}.}
 \tag{9}
\]
Equivalently, the input pairs \(h=(f+ig)/\sqrt2\) obtained from (8) have H2 seed objective at most 1 in (5). They cannot even improve the unamplified H2 baseline, whereas arbitrary inputs already attain 1.25 at outer order four.

Proof. Associate to \(q\) the symmetric binary matrix \(L\) with diagonal \(\alpha_j\bmod2\) and off-diagonal entries \(\beta_{jk}\). The derivative \(z(x+t)/z(x)\) is a constant fourth root of unity times the binary character \((-1)^{x^TLt}\). Hence the Fourier magnitude of \(z\) is supported on an affine subspace of size \(2^{\operatorname{rank}L}\), and every nonzero normalized Fourier value has magnitude \(2^{k/2}\), where \(k=r-\operatorname{rank}L\).

For completeness, elimination of the binary quadratic form gives one-dimensional blocks with odd diagonal, two-dimensional alternating blocks, and radical coordinates. Their nonzero unnormalized Gauss sums are respectively fourth-root multiples of \(1+i\), of 2, and of 2. Thus every nonzero Fourier phase is a fourth-root multiple of \(e^{i\pi(r-k)/4}\). This also follows by performing the one-variable sum \(1\pm i\) and the nondegenerate alternating two-variable sum explicitly at each elimination step.

If \(k\) is even, then \(r-k\) is even; multiplying by \(e^{i\pi/4}\) leaves every nonzero output halfway between adjacent fourth roots. The left side of (9) is then \(2^{-k/2}/\sqrt2\le1/\sqrt2\). If \(k\) is odd, the output phases are fourth roots, but their support fraction is \(2^{-k}\), so the left side is \(2^{-k/2}\le1/\sqrt2\). This proves (9).

The same restricted-family bound holds for \(H_4^{\otimes a}\) after its quadratic sign gauge and binary linear reindexing. It does not analyze the order-144 outer, and it does not constrain nonquadratic quaternary phase inputs. Its use is to rule out importing an ordinary quadratic generalized-bent construction as a phase-conversion proof.

## 8. Literature correction: the H2 phase defect does tend to zero

The preceding finite-orbit and restricted-quadratic obstructions remain
valid, but they do **not** describe the global optimum. Kai-Uwe Schmidt,
*Asymptotically optimal Boolean functions*, Theorem 1 and Section 2,
prove that
\[
 \mu_r:=\min_{f\in\{\pm1\}^{2^r}}
       \frac{\|W_{2^r}f\|_\infty}{2^{r/2}}\longrightarrow1
\]
as the number of variables tends to infinity, including odd numbers of
variables. This statement and its normalization were read directly in
the [primary manuscript](https://math.uni-paderborn.de/fileadmin-eim/mathematik/AG-Diskrete_Mathematik/Publications-schmidt/pw.pdf).
The literature agent identified the missing implication; the following
mapping was independently checked by the algebra agent.

Put `N=2^r`, and select a minimizing Boolean `f`. Parseval gives
`||W_N f||_2^2=N^2`, hence
\[
 \|W_Nf\|_1\ge
 \frac{\|W_Nf\|_2^2}{\|W_Nf\|_\infty}
 \ge\frac{N^{3/2}}{\mu_r}.
\]
Choose `g=sign(W_N f)` and use the exact lift (3). Its same-spin
normalized Rayleigh quotient at order `4N` is at least `1/mu_r`, while
orthogonality bounds every such quotient by one. Therefore
\[
 \lim_{r\to\infty}
 \frac{\max_{x\in\{\pm1\}^{2^r}}|x^TW_{2^r}x|}{2^{3r/2}}=1.
\]
For the odd subsequence use odd `r` in the construction, which remains
odd after the two-variable lift. The congruence (2) now proves
\[
 \boxed{\mathcal R_4(H_2)=\mathcal R(H_2)=\sqrt2.}
\]
Indeed the lower bound already occurs on the `H4` branch, and the
operator norm gives the matching upper bound on every admitted outer.

This chooses globally appropriate starting functions at increasing
orders. It does not assert improvement along the single permanently
stagnating orbit in Section 3, and it does not contradict the quadratic
phase restriction in Section 7. It resolves this seed's regularized
value, not the original all-order minimax convergence problem or the
general identity `R=T`.
